"""
quality_assessment.py - Step 4.1 quality assessment of one parcel's
Sentinel-2 NDVI / NDMI daily time series.

This module reads the time-series CSV produced in Steps 2/3 (read-only),
preserves every original observation and column, and appends quality-assurance
columns plus a human-readable Markdown report. It never removes, replaces,
smooths or interpolates observations, and it does not touch
src/sentinel_stats.py, src/plots.py, the source CSV, the PNG figures or the
satellite-processing methodology.

Statistical rules (all documented in the report and configurable here):

    * 1.5 x IQR fences  -> distribution outliers,
    * MAD robust z-score (0.6745 * (x - median) / MAD, |z| > 3.5)
      -> distribution outliers (Iglewicz & Hoaglin),
    * robust rules are skipped safely when IQR == 0 or MAD == 0,
    * coverage completeness is a definitional flag (clear_fraction < 1.0 for
      in-polygon masking) and is kept distinct from statistical outliers.

Flags mark *potential anomalies requiring review*; they are never treated as
confirmed invalid observations.

Run the offline self-check (no network, no files written):
    python -m src.quality_assessment --check

Generate the QA dataset and report:
    python -m src.quality_assessment
"""

from __future__ import annotations

import argparse
import csv
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

try:  # works for "python -m src.quality_assessment"
    from .plots import parse_timestamp
except ImportError:  # fallback for "python src/quality_assessment.py"
    from plots import parse_timestamp  # type: ignore

# --------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CSV_PATH = PROJECT_ROOT / "data" / "results" / "time_series_p001_2023_2024.csv"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "data" / "results"

# Robust-outlier rules (Iglewicz & Hoaglin). Documented in the report.
IQR_MULTIPLIER = 1.5
MAD_SCALE = 0.6745
MAD_Z_THRESHOLD = 3.5

# Tolerance used only to detect the "complete coverage" definitional case.
COMPLETE_COVERAGE_TOL = 1e-9

# Columns appended to the original CSV (original columns are kept untouched).
QA_FIELDS = [
    "is_valid",
    "days_since_prev_valid",
    "ndvi_change",
    "ndmi_change",
    "ndvi_change_per_day",
    "ndmi_change_per_day",
    "coverage_status",
    "valid_count_outlier",
    "clear_fraction_outlier",
    "ndvi_change_outlier",
    "ndmi_change_outlier",
    "ndvi_rate_outlier",
    "ndmi_rate_outlier",
    "low_quality",
    "qa_flags",
]


# --------------------------------------------------------------------------
# Small helpers
# --------------------------------------------------------------------------

def _num(value) -> float | None:
    """Parse a numeric CSV cell; blank / unparsable cells become None."""
    if value is None:
        return None
    text = str(value).strip()
    if text == "":
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _fmt(value: float | None, digits: int = 6) -> str:
    """Format a number for the QA CSV, or '' when it is missing."""
    if value is None:
        return ""
    return f"{value:.{digits}f}"


def load_rows(csv_path: str | Path) -> tuple[list[str], list[dict]]:
    """Read a time-series CSV (read-only). Returns (fieldnames, rows)."""
    csv_path = Path(csv_path)
    if not csv_path.exists():
        raise FileNotFoundError(f"Input CSV not found: {csv_path}")
    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)
    return fieldnames, rows


# --------------------------------------------------------------------------
# Robust statistics (pure, offline-testable)
# --------------------------------------------------------------------------

def metric_stats(values: list[float]) -> dict:
    """Quartiles, MAD and IQR fences for one metric.

    Returns a dict with q1/median/q3/iqr/mad/lower_fence/upper_fence and the
    booleans iqr_usable / mad_usable. When IQR == 0 (a (near-)constant metric)
    the IQR rule is unusable; when MAD == 0 the robust z-score is unusable.
    """
    arr = np.asarray(values, dtype=float)
    q1 = float(np.percentile(arr, 25))
    q3 = float(np.percentile(arr, 75))
    median = float(np.median(arr))
    iqr = q3 - q1
    mad = float(np.median(np.abs(arr - median)))
    return {
        "n": int(arr.size),
        "min": float(arr.min()),
        "max": float(arr.max()),
        "mean": float(arr.mean()),
        "q1": q1,
        "median": median,
        "q3": q3,
        "iqr": iqr,
        "mad": mad,
        "lower_fence": q1 - IQR_MULTIPLIER * iqr,
        "upper_fence": q3 + IQR_MULTIPLIER * iqr,
        "iqr_usable": iqr > 0.0,
        "mad_usable": mad > 0.0,
    }


def flag_metric(values: list[float], stats: dict | None = None) -> list[str | None]:
    """Return 'low'/'high'/None per value using IQR fences OR MAD robust z.

    A value is flagged when it falls outside a 1.5 x IQR fence or when its MAD
    robust z-score exceeds the threshold. IQR takes precedence on the (rare)
    conflicting case. Degenerate rules are skipped rather than guessed.
    """
    if stats is None:
        stats = metric_stats(values)
    mad = stats["mad"]
    median = stats["median"]
    flags: list[str | None] = []
    for value in values:
        flag: str | None = None
        if stats["iqr_usable"]:
            if value < stats["lower_fence"]:
                flag = "low"
            elif value > stats["upper_fence"]:
                flag = "high"
        if flag is None and stats["mad_usable"]:
            z = MAD_SCALE * (value - median) / mad
            if abs(z) > MAD_Z_THRESHOLD:
                flag = "high" if z > 0 else "low"
        flags.append(flag)
    return flags


# --------------------------------------------------------------------------
# Analysis (pure, offline-testable)
# --------------------------------------------------------------------------

def analyze(rows: list[dict]) -> tuple[list[dict], dict]:
    """Append QA columns to every row and return (qa_rows, analysis).

    The output has exactly the same number of rows as the input; nothing is
    dropped or modified. Change columns are computed only between consecutive
    *valid* observations (valid_count > 0), using the elapsed days between
    them so that multi-day gaps are never reported as daily changes.
    """
    n_total = len(rows)

    valid_positions = [
        i for i, row in enumerate(rows) if (_num(row.get("valid_count")) or 0) > 0
    ]
    valid_rows = [rows[i] for i in valid_positions]

    # ---- per-valid-observation change series -----------------------------
    days_since_prev: list[int | None] = [None]
    d_ndvi: list[float | None] = [None]
    d_ndmi: list[float | None] = [None]
    r_ndvi: list[float | None] = [None]
    r_ndmi: list[float | None] = [None]
    for prev, cur in zip(valid_rows, valid_rows[1:]):
        prev_date = parse_timestamp(prev["interval_from"])
        cur_date = parse_timestamp(cur["interval_from"])
        dt = (cur_date - prev_date).days
        days_since_prev.append(dt)

        prev_ndvi, cur_ndvi = _num(prev.get("ndvi_mean")), _num(cur.get("ndvi_mean"))
        prev_ndmi, cur_ndmi = _num(prev.get("ndmi_mean")), _num(cur.get("ndmi_mean"))

        ndvi_delta = cur_ndvi - prev_ndvi if None not in (cur_ndvi, prev_ndvi) else None
        ndmi_delta = cur_ndmi - prev_ndmi if None not in (cur_ndmi, prev_ndmi) else None
        d_ndvi.append(ndvi_delta)
        d_ndmi.append(ndmi_delta)
        r_ndvi.append(ndvi_delta / dt if ndvi_delta is not None and dt > 0 else None)
        r_ndmi.append(ndmi_delta / dt if ndmi_delta is not None and dt > 0 else None)

    # ---- robust statistics over valid observations ----------------------
    metric_series = {
        "valid_count": [float(_num(r.get("valid_count")) or 0) for r in valid_rows],
        "valid_coverage": [float(_num(r.get("valid_coverage")) or 0) for r in valid_rows],
        "clear_fraction": [float(_num(r.get("clear_fraction")) or 0) for r in valid_rows],
        "ndvi_change": [x for x in d_ndvi if x is not None],
        "ndmi_change": [x for x in d_ndmi if x is not None],
        "ndvi_rate": [x for x in r_ndvi if x is not None],
        "ndmi_rate": [x for x in r_ndmi if x is not None],
    }
    stats = {name: metric_stats(values) for name, values in metric_series.items() if values}

    def _flags_positions(name: str, series: list[float | None]) -> list[str | None]:
        """Flag a series with blanks by scoring only its present values."""
        present = [x for x in series if x is not None]
        if not present or name not in stats:
            return [None] * len(series)
        present_flags = flag_metric(present, stats[name])
        out: list[str | None] = []
        it = iter(present_flags)
        for x in series:
            out.append(next(it) if x is not None else None)
        return out

    count_flags = flag_metric(metric_series["valid_count"], stats["valid_count"])
    cf_flags = flag_metric(metric_series["clear_fraction"], stats["clear_fraction"])
    ndvi_ch_flags = _flags_positions("ndvi_change", d_ndvi)
    ndmi_ch_flags = _flags_positions("ndmi_change", d_ndmi)
    ndvi_rt_flags = _flags_positions("ndvi_rate", r_ndvi)
    ndmi_rt_flags = _flags_positions("ndmi_rate", r_ndmi)

    # ---- assemble QA rows (every original row preserved) ----------------
    qa_by_position: dict[int, dict] = {}

    for row in rows:
        valid = (_num(row.get("valid_count")) or 0) > 0
        if valid:
            continue
        qa_by_position[id(row)] = {
            "is_valid": "no",
            "days_since_prev_valid": "",
            "ndvi_change": "",
            "ndmi_change": "",
            "ndvi_change_per_day": "",
            "ndmi_change_per_day": "",
            "coverage_status": "masked",
            "valid_count_outlier": "",
            "clear_fraction_outlier": "",
            "ndvi_change_outlier": "",
            "ndmi_change_outlier": "",
            "ndvi_rate_outlier": "",
            "ndmi_rate_outlier": "",
            "low_quality": "yes",
            "qa_flags": "fully masked (valid_count=0)",
        }

    coverage_counts = {"complete": 0, "partial": 0, "masked": 0}
    flagged: list[dict] = []

    for k, row in enumerate(valid_rows):
        cf = _num(row.get("clear_fraction"))
        complete = cf is not None and abs(cf - 1.0) <= COMPLETE_COVERAGE_TOL
        status = "complete" if complete else "partial"
        coverage_counts[status] += 1

        reasons: list[str] = []
        if not complete:
            reasons.append(
                f"reduced in-polygon clear coverage (clear_fraction={cf:.3f})"
            )

        flag_fields = {
            "valid_count_outlier": count_flags[k],
            "clear_fraction_outlier": cf_flags[k],
            "ndvi_change_outlier": ndvi_ch_flags[k],
            "ndmi_change_outlier": ndmi_ch_flags[k],
            "ndvi_rate_outlier": ndvi_rt_flags[k],
            "ndmi_rate_outlier": ndmi_rt_flags[k],
        }
        for metric, flag in flag_fields.items():
            if flag:
                reasons.append(f"statistical outlier: {metric[:-8]} ({flag})")

        qa_by_position[id(row)] = {
            "is_valid": "yes",
            "days_since_prev_valid": (
                "" if days_since_prev[k] is None else str(days_since_prev[k])
            ),
            "ndvi_change": _fmt(d_ndvi[k]),
            "ndmi_change": _fmt(d_ndmi[k]),
            "ndvi_change_per_day": _fmt(r_ndvi[k]),
            "ndmi_change_per_day": _fmt(r_ndmi[k]),
            "coverage_status": status,
            "valid_count_outlier": "yes" if flag_fields["valid_count_outlier"] else "no",
            "clear_fraction_outlier": "yes" if flag_fields["clear_fraction_outlier"] else "no",
            "ndvi_change_outlier": "yes" if flag_fields["ndvi_change_outlier"] else "no",
            "ndmi_change_outlier": "yes" if flag_fields["ndmi_change_outlier"] else "no",
            "ndvi_rate_outlier": "yes" if flag_fields["ndvi_rate_outlier"] else "no",
            "ndmi_rate_outlier": "yes" if flag_fields["ndmi_rate_outlier"] else "no",
            "low_quality": "no" if complete else "yes",
            "qa_flags": "; ".join(reasons) if reasons else "",
        }

        if reasons:
            flagged.append(
                {
                    "date": row.get("interval_from", "")[:10],
                    "ndvi": _num(row.get("ndvi_mean")),
                    "ndmi": _num(row.get("ndmi_mean")),
                    "valid_count": int(_num(row.get("valid_count")) or 0),
                    "clear_fraction": cf,
                    "coverage_status": status,
                    "days_since_prev_valid": days_since_prev[k],
                    "ndvi_change": d_ndvi[k],
                    "ndmi_change": d_ndmi[k],
                    "ndvi_change_per_day": r_ndvi[k],
                    "ndmi_change_per_day": r_ndmi[k],
                    "reasons": reasons,
                }
            )

    qa_rows = [dict(row, **qa_by_position[id(row)]) for row in rows]
    coverage_counts["masked"] = n_total - len(valid_rows)
    analysis = {
        "n_total": n_total,
        "n_valid": len(valid_rows),
        "n_masked": n_total - len(valid_rows),
        "stats": stats,
        "coverage_counts": coverage_counts,
        "flagged": flagged,
        "valid_rows": valid_rows,
        "days_since_prev": days_since_prev,
        "d_ndvi": d_ndvi,
        "d_ndmi": d_ndmi,
        "r_ndvi": r_ndvi,
        "r_ndmi": r_ndmi,
    }
    return qa_rows, analysis


# --------------------------------------------------------------------------
# Output writers
# --------------------------------------------------------------------------

def write_qa_csv(
    qa_rows: list[dict],
    original_fields: list[str],
    output_path: str | Path,
) -> Path:
    """Write the QA dataset: original columns verbatim + appended QA columns."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(original_fields) + QA_FIELDS
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, restval="")
        writer.writeheader()
        for row in qa_rows:
            writer.writerow(row)
    return output_path


def _stats_row(name: str, values: dict, unit: str = "") -> str:
    return (
        f"| {name}{unit} | {values['n']} | {values['min']:.4f} | {values['q1']:.4f} | "
        f"{values['median']:.4f} | {values['q3']:.4f} | {values['max']:.4f} | "
        f"{values['iqr']:.4f} | {values['mad']:.4f} | "
        f"{'yes' if values['iqr_usable'] else 'no'} | "
        f"{'yes' if values['mad_usable'] else 'no'} |"
    )


def build_report(analysis: dict, csv_path: str | Path) -> str:
    """Render the reproducible Markdown quality report."""
    stats = analysis["stats"]
    cc = analysis["coverage_counts"]
    geometry_values = sorted(
        {round(_num(r.get("geometry_coverage")) or -1, 12) for r in analysis["valid_rows"]}
    )

    lines: list[str] = []
    lines.append("# Step 4.1 - Sentinel-2 time-series quality assessment (parcel P001)")
    lines.append("")
    lines.append(f"- **Input (read-only):** `{csv_path}`")
    lines.append("- **Reproduce:** `python -m src.quality_assessment`")
    lines.append(
        f"- **Observations:** {analysis['n_total']} total; "
        f"**{analysis['n_valid']} valid** (`valid_count > 0`); "
        f"{analysis['n_masked']} fully masked."
    )
    lines.append(
        "- **No observation is removed, replaced, smoothed or interpolated**; the "
        "original CSV, methodology, masking rules and PNG figures are unchanged."
    )
    lines.append("")

    lines.append("## 1. Coverage metrics: what each one measures")
    lines.append("")
    lines.append(
        "- `valid_count` - number of 20 m pixels inside the polygon that passed "
        "`dataMask == 1 AND SCL in {4, 5}`."
    )
    lines.append(
        "- `geometry_coverage = geometry_pixel_count / sample_count` - the "
        "polygon's share of the request bounding box. This is a **pure geometry "
        "ratio** and is **constant** for a fixed parcel: here "
        f"{geometry_values[0]:.6f} on every observation."
    )
    lines.append(
        "- `valid_coverage = valid_count / sample_count` - valid pixels as a "
        "fraction of the **bounding box** (cover + cloud + non-agricultural SCL)."
    )
    lines.append(
        "- `clear_fraction = valid_count / geometry_pixel_count` - the "
        "**cloud-free, SCL-valid fraction inside the polygon**. For P001 the "
        "bounding box is 23 x 9 = 207 pixels while the polygon contains 96, so "
        "`valid_coverage` and `clear_fraction` differ by the constant geometric "
        "factor 96/207 = 0.4638. They measure **different aspects of coverage**: "
        "`valid_coverage` mixes geometry with data availability, whereas "
        "`clear_fraction` isolates in-polygon masking."
    )
    lines.append("")

    lines.append("## 2. Statistical rules (all computed from the data)")
    lines.append("")
    lines.append(
        f"- **IQR fences:** a value is an outlier when it falls below "
        f"`Q1 - {IQR_MULTIPLIER} * IQR` or above `Q3 + {IQR_MULTIPLIER} * IQR`."
    )
    lines.append(
        f"- **MAD robust z-score:** `z = {MAD_SCALE} * (x - median) / MAD`; "
        f"`|z| > {MAD_Z_THRESHOLD}` is an outlier (Iglewicz & Hoaglin)."
    )
    lines.append(
        "- A value is flagged if **either** rule fires; IQR takes precedence if "
        "the two ever conflict."
    )
    lines.append(
        "- **Degenerate rules are skipped safely:** when `IQR == 0` the IQR rule "
        "is unusable and when `MAD == 0` the z-score is unusable. This is "
        "expected for `clear_fraction`, which equals 1.0 for most valid "
        "observations, so it receives **no** robust outlier flags."
    )
    lines.append(
        "- **Coverage completeness is a definitional flag, not a statistical "
        "one:** `clear_fraction < 1.0` means some in-polygon pixels were masked "
        "(cloud / SCL / no-data). It is reported separately so that *low-quality* "
        "observations are distinguished from *statistical outliers*."
    )
    lines.append(
        "- **Changes** are computed only between consecutive valid observations. "
        "`*_change` is the difference across the elapsed interval and "
        "`*_change_per_day` divides it by the elapsed days, so a gap of several "
        "days is never described as a daily change."
    )
    lines.append("")
    lines.append("Computed dispersion thresholds per metric:")
    lines.append("")
    lines.append("| Metric | n | min | Q1 | median | Q3 | max | IQR | MAD | IQR usable | MAD usable |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|---|")
    for name, label in (
        ("valid_count", "valid_count (px)"),
        ("valid_coverage", "valid_coverage"),
        ("clear_fraction", "clear_fraction"),
        ("ndvi_change", "dNDVI (interval)"),
        ("ndmi_change", "dNDMI (interval)"),
        ("ndvi_rate", "dNDVI per day"),
        ("ndmi_rate", "dNDMI per day"),
    ):
        if name in stats:
            lines.append(_stats_row(label, stats[name]))
    lines.append("")

    lines.append("## 3. Coverage summary")
    lines.append("")
    lines.append(
        f"- Complete, fully clear valid observations (`clear_fraction = 1.0`): "
        f"**{cc['complete']}**."
    )
    lines.append(
        f"- Valid but partially masked observations (`0 < clear_fraction < 1`): "
        f"**{cc['partial']}** (low-quality candidates)."
    )
    lines.append(f"- Fully masked observations (`valid_count = 0`): **{cc['masked']}**.")
    lines.append("")

    lines.append("## 4. Flagged observations (potential anomalies - review, not rejection)")
    lines.append("")
    if analysis["flagged"]:
        lines.append(
            "| Date | NDVI | NDMI | valid_count | clear_fraction | coverage | "
            "gap (d) | dNDVI | dNDMI | dNDVI/day | dNDMI/day | Reasons |"
        )
        lines.append(
            "|---|---|---|---|---|---|---|---|---|---|---|---|"
        )
        for item in analysis["flagged"]:
            lines.append(
                "| {date} | {ndvi:.3f} | {ndmi:.3f} | {vc} | {cf:.3f} | {cov} | "
                "{gap} | {dn} | {dm} | {rn} | {rm} | {why} |".format(
                    date=item["date"],
                    ndvi=item["ndvi"],
                    ndmi=item["ndmi"],
                    vc=item["valid_count"],
                    cf=item["clear_fraction"],
                    cov=item["coverage_status"],
                    gap="" if item["days_since_prev_valid"] is None
                    else item["days_since_prev_valid"],
                    dn="-" if item["ndvi_change"] is None else f"{item['ndvi_change']:+.3f}",
                    dm="-" if item["ndmi_change"] is None else f"{item['ndmi_change']:+.3f}",
                    rn="-" if item["ndvi_change_per_day"] is None
                    else f"{item['ndvi_change_per_day']:+.3f}",
                    rm="-" if item["ndmi_change_per_day"] is None
                    else f"{item['ndmi_change_per_day']:+.3f}",
                    why="; ".join(item["reasons"]),
                )
            )
    else:
        lines.append("No observations were flagged.")
    lines.append("")

    lines.append("## 5. Interpretation")
    lines.append("")
    lines.append(
        "- Statistical outliers (IQR / MAD) are values far from the parcel's own "
        "distribution of that metric; they are flagged for review only."
    )
    lines.append(
        "- Low-quality observations are those with reduced in-polygon clear "
        "coverage. Because most valid days have `clear_fraction = 1.0`, the "
        "robust rules are degenerate there and these days are identified by the "
        "definitional coverage flag instead."
    )
    lines.append(
        "- Outliers concentrated on low-`clear_fraction` days should be treated as "
        "less reliable, not as confirmed land-cover change; no automatic removal "
        "or arbitrary thresholding is applied."
    )
    lines.append("")
    return "\n".join(lines)


# --------------------------------------------------------------------------
# Offline self-check (no network, no files written)
# --------------------------------------------------------------------------

def _self_check() -> None:
    # Robust stats and degenerate handling.
    constant = metric_stats([1.0, 1.0, 1.0, 1.0, 0.5])
    assert constant["iqr"] == 0.0 and not constant["iqr_usable"], constant
    assert constant["mad"] == 0.0 and not constant["mad_usable"], constant
    assert flag_metric([1.0, 1.0, 1.0, 1.0, 0.5], constant) == [None] * 5

    spread = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 100.0]
    flags = flag_metric(spread)
    assert flags[-1] == "high" and all(f is None for f in flags[:-1]), flags

    # Change computation + elapsed-days handling.
    def _row(date, ndvi, ndmi, vc, cf):
        return {
            "interval_from": f"{date}T00:00:00Z",
            "interval_to": "",
            "ndvi_mean": ndvi,
            "ndvi_stdev": "0.05",
            "ndmi_mean": ndmi,
            "ndmi_stdev": "0.05",
            "sample_count": "207",
            "no_data_count": "0",
            "valid_count": vc,
            "geometry_pixel_count": "96",
            "geometry_coverage": "0.463768115942029",
            "valid_coverage": str(int(vc) / 207),
            "clear_fraction": cf,
        }

    rows = [
        _row("2024-01-17", "0.6446", "0.1286", "96", "1.0"),
        _row("2024-01-18", "", "", "0", "0.0"),  # masked, must be preserved
        _row("2024-01-22", "0.2576", "0.6321", "61", "0.6354"),
    ]
    qa_rows, analysis = analyze(rows)

    # Requirement: every original observation preserved.
    assert len(qa_rows) == len(rows) == 3
    assert qa_rows[1]["coverage_status"] == "masked"
    assert qa_rows[1]["qa_flags"] == "fully masked (valid_count=0)"
    assert qa_rows[0]["days_since_prev_valid"] == ""  # first valid -> no change
    assert qa_rows[2]["days_since_prev_valid"] == "5"
    assert qa_rows[2]["ndvi_change"] == "-0.387000"
    assert qa_rows[2]["ndmi_change"] == "0.503500"
    assert qa_rows[2]["ndvi_change_per_day"] == "-0.077400"
    assert qa_rows[2]["ndmi_change_per_day"] == "0.100700"
    assert qa_rows[2]["low_quality"] == "yes"  # partial coverage
    assert "reduced in-polygon clear coverage" in qa_rows[2]["qa_flags"]

    # Original columns are all still present on every QA row.
    original_fields = list(rows[0].keys())
    for field in original_fields:
        assert field in qa_rows[0], field

    # Report renders and documents the key concepts.
    report = build_report(analysis, "sample.csv")
    assert "geometry_coverage" in report
    assert "clear_fraction" in report
    assert "IQR" in report and "MAD" in report

    print("quality_assessment.py self-check passed (offline, no files written).")


# --------------------------------------------------------------------------
# Command-line interface
# --------------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m src.quality_assessment",
        description="Step 4.1 quality assessment for a parcel time-series CSV.",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Run the offline self-check and exit (no files written).",
    )
    parser.add_argument(
        "--csv",
        default=str(DEFAULT_CSV_PATH),
        help=f"Input time-series CSV (default: {DEFAULT_CSV_PATH}).",
    )
    parser.add_argument(
        "--outdir",
        default=str(DEFAULT_OUTPUT_DIR),
        help=f"Directory for the QA outputs (default: {DEFAULT_OUTPUT_DIR}).",
    )
    parser.add_argument(
        "--suffix",
        default="_qa",
        help="Suffix appended to the input stem for the QA CSV (default: _qa).",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    if args.check:
        _self_check()
        return 0

    csv_path = Path(args.csv)
    outdir = Path(args.outdir)
    fieldnames, rows = load_rows(csv_path)
    qa_rows, analysis = analyze(rows)

    qa_csv_path = outdir / f"{csv_path.stem}{args.suffix}.csv"
    report_path = outdir / "quality_report_p001_2023_2024.md"

    write_qa_csv(qa_rows, fieldnames, qa_csv_path)
    report_path.write_text(build_report(analysis, csv_path), encoding="utf-8")

    flagged = len(analysis["flagged"])
    cc = analysis["coverage_counts"]
    print(
        f"Quality assessment: {analysis['n_valid']} valid / "
        f"{analysis['n_total']} total observation(s); "
        f"{flagged} flagged for review "
        f"({cc['partial']} partial-coverage, {cc['masked']} masked)."
    )
    print(f"Wrote QA dataset: {qa_csv_path}")
    print(f"Wrote report:     {report_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
