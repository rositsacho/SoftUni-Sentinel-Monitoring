"""
crop_signature.py - Step 4.2 retrospective crop-signature analysis for P001.

Compares monthly and seasonal NDVI / NDMI profiles of the observed Sentinel-2
time series against *reference* crop labels supplied by the user. This module
is descriptive only: it performs NO crop classification and produces NO
satellite-derived predictions.

Three things are kept strictly separate:
    1. OBSERVATIONS - the real NDVI/NDMI values from the Step 2/3/4.1 workflow,
    2. REFERENCE LABELS - the user-provided ground-truth crop per harvest year
       (auditable config in data/reference_crops.csv, not hard-coded here),
    3. AGRONOMIC EXPECTATIONS - qualitative published phenology, not data.

Reference harvest years are interpreted as the year of harvest; crop-season
windows are configurable analytical assumptions (see data/reference_crops.csv),
NOT known sowing or harvest dates.

Statistics are computed on real observations only - no interpolation, no
fabrication. Two quality tiers are reported everywhere:
    * "all_valid"       - every observation with valid_count > 0,
    * "complete"        - valid observations with clear_fraction == 1.0
                          (coverage_status == "complete" in the Step 4.1 QA).

No existing file is modified by this module.

Run the offline self-check (no network, no files written):
    python -m src.crop_signature --check

Generate the crop-signature CSV, plots and report:
    python -m src.crop_signature
"""

from __future__ import annotations

import argparse
import csv
import io
import sys
from datetime import date, datetime
from pathlib import Path

import numpy as np

import matplotlib

matplotlib.use("Agg")  # headless, reproducible

import matplotlib.pyplot as plt  # noqa: E402

try:  # works for "python -m src.crop_signature"
    from .plots import parse_timestamp
    from .quality_assessment import load_rows
except ImportError:  # fallback for "python src/crop_signature.py"
    from plots import parse_timestamp  # type: ignore
    from quality_assessment import load_rows  # type: ignore

# --------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_QA_CSV_PATH = PROJECT_ROOT / "data" / "results" / "time_series_p001_2023_2024_qa.csv"
DEFAULT_REFERENCE_CSV = PROJECT_ROOT / "data" / "reference_crops.csv"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "data" / "results"
DEFAULT_PARCEL_ID = "P001"
DEFAULT_DPI = 150

TIERS = ("all_valid", "complete")
TIER_LABELS = {
    "all_valid": "All valid observations",
    "complete": "Complete-coverage observations",
}

# Small, verifiable bibliography. Keys are referenced from reference_crops.csv.
REFERENCE_LIBRARY = {
    "rouse1974": (
        "Rouse, J.W., Haas, R.H., Schell, J.A., Deering, D.W. (1974). "
        "Monitoring vegetation systems in the Great Plains with ERTS. "
        "NASA SP-351, 309-317."
    ),
    "tucker1979": (
        "Tucker, C.J. (1979). Red and photographic infrared linear combinations "
        "for monitoring vegetation. Remote Sensing of Environment, 8(2), 127-150."
    ),
    "gao1996": (
        "Gao, B.-C. (1996). NDWI - A normalized difference water index for "
        "remote sensing of vegetation from space. Remote Sensing of Environment, "
        "58(3), 257-266."
    ),
    "zadoks1974": (
        "Zadoks, J.C., Chang, T.T., Konzak, C.F. (1974). A decimal code for the "
        "growth stages of cereals. Weed Research, 14(6), 415-421."
    ),
    "schneiter1981": (
        "Schneiter, A.A., Miller, J.F. (1981). Description of sunflower growth "
        "stages. Crop Science, 21(6), 901-903."
    ),
    "beckerreshef2010": (
        "Becker-Reshef, I., Vermote, E., Lindeman, M., Justice, C. (2010). A "
        "generalized regression-based model for forecasting winter wheat yields "
        "in Kansas and Ukraine using MODIS data. Remote Sensing of Environment, "
        "114(6), 1312-1323."
    ),
}

# Tidy output columns (one row per record_type x tier).
OUTPUT_COLUMNS = [
    "record_type",
    "parcel_id",
    "tier",
    "crop",
    "harvest_year",
    "period",
    "season_start",
    "season_end",
    "status",
    "n_obs",
    "months_with_obs",
    "months_in_window",
    "observed_start",
    "observed_end",
    "ndvi_mean",
    "ndvi_median",
    "ndvi_min",
    "ndvi_max",
    "ndvi_peak",
    "ndvi_peak_date",
    "ndmi_mean",
    "ndmi_median",
    "ndmi_min",
    "ndmi_max",
    "ndmi_peak",
    "ndmi_peak_date",
    "notes",
]


# --------------------------------------------------------------------------
# Reference-crop configuration
# --------------------------------------------------------------------------

def parse_reference_crops(text: str) -> list[dict]:
    """Parse the reference-crop CSV text into a list of dicts (dates as date)."""
    crops: list[dict] = []
    reader = csv.DictReader(io.StringIO(text))
    required = {"parcel_id", "harvest_year", "crop", "season_start", "season_end"}
    missing = required - set(reader.fieldnames or [])
    if missing:
        raise ValueError(f"reference CSV missing column(s): {sorted(missing)}")
    for row in reader:
        row = {k: (v or "").strip() for k, v in row.items()}
        if not row["parcel_id"]:
            continue
        crops.append(
            {
                "parcel_id": row["parcel_id"],
                "harvest_year": int(row["harvest_year"]),
                "crop": row["crop"],
                "season_start": date.fromisoformat(row["season_start"]),
                "season_end": date.fromisoformat(row["season_end"]),
                "window_rule": row.get("window_rule", ""),
                "expected_phenology": row.get("expected_phenology", ""),
                "reference_keys": [
                    k.strip()
                    for k in row.get("reference_keys", "").split(";")
                    if k.strip()
                ],
            }
        )
    return crops


def load_reference_crops(path: str | Path) -> list[dict]:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Reference-crop config not found: {path}")
    return parse_reference_crops(path.read_text(encoding="utf-8"))


# --------------------------------------------------------------------------
# Observation selection (pure, offline-testable)
# --------------------------------------------------------------------------

def _num(value) -> float | None:
    if value is None:
        return None
    text = str(value).strip()
    if text == "":
        return None
    try:
        return float(text)
    except ValueError:
        return None


def tier_matches(row: dict, tier: str) -> bool:
    """True when a QA-CSV row belongs to the requested quality tier."""
    if str(row.get("is_valid", "")).strip().lower() not in ("yes", "true", "1"):
        return False
    if tier == "all_valid":
        return True
    if tier == "complete":
        return str(row.get("coverage_status", "")).strip() == "complete"
    raise ValueError(f"unknown tier: {tier}")


def _observation(row: dict) -> dict | None:
    """(date, ndvi, ndmi) for a row, or None when either index is missing."""
    ndvi = _num(row.get("ndvi_mean"))
    ndmi = _num(row.get("ndmi_mean"))
    stamp = (row.get("interval_from") or "").strip()
    if ndvi is None or ndmi is None or not stamp:
        return None
    return {"date": parse_timestamp(stamp).date(), "ndvi": ndvi, "ndmi": ndmi}


def select_observations(rows: list[dict], tier: str) -> list[dict]:
    """Real observations of one tier, chronologically sorted. No fabrication."""
    selected = [
        obs
        for row in rows
        if tier_matches(row, tier) and (obs := _observation(row)) is not None
    ]
    selected.sort(key=lambda o: o["date"])
    return selected


# --------------------------------------------------------------------------
# Aggregation (pure, offline-testable)
# --------------------------------------------------------------------------

def _stats(values: list[float]) -> dict:
    arr = np.asarray(values, dtype=float)
    return {
        "n": int(arr.size),
        "mean": float(arr.mean()),
        "median": float(np.median(arr)),
        "min": float(arr.min()),
        "max": float(arr.max()),
    }


def monthly_profile(observations: list[dict]) -> list[dict]:
    """Per-calendar-month summary of real observations (no resampling)."""
    buckets: dict[str, list[dict]] = {}
    for obs in observations:
        buckets.setdefault(obs["date"].strftime("%Y-%m"), []).append(obs)
    profile: list[dict] = []
    for month in sorted(buckets):
        obss = buckets[month]
        ndvi = _stats([o["ndvi"] for o in obss])
        ndmi = _stats([o["ndmi"] for o in obss])
        profile.append({"month": month, "n_obs": len(obss), "ndvi": ndvi, "ndmi": ndmi})
    return profile


def _iter_months(start: date, end: date) -> list[str]:
    months: list[str] = []
    y, m = start.year, start.month
    while (y, m) <= (end.year, end.month):
        months.append(f"{y:04d}-{m:02d}")
        m += 1
        if m > 12:
            m, y = 1, y + 1
    return months


def seasonal_summary(observations: list[dict], window: dict) -> dict:
    """Summary of one crop-season window from real observations only.

    status is "no observations" when nothing falls inside the window,
    "partially observed" when some window months have no usable observation
    (or the window extends beyond the observed span), and "observed" otherwise.
    """
    start, end = window["season_start"], window["season_end"]
    in_window = [o for o in observations if start <= o["date"] <= end]
    window_months = _iter_months(start, end)
    observed_months = sorted({o["date"].strftime("%Y-%m") for o in in_window})

    if not in_window:
        return {
            "crop": window["crop"],
            "harvest_year": window["harvest_year"],
            "season_start": start,
            "season_end": end,
            "status": "no observations",
            "n_obs": 0,
            "months_with_obs": 0,
            "months_in_window": len(window_months),
            "observed_start": None,
            "observed_end": None,
            "ndvi": None,
            "ndmi": None,
            "ndvi_peak": None,
            "ndvi_peak_date": None,
            "ndmi_peak": None,
            "ndmi_peak_date": None,
        }

    ndvi = _stats([o["ndvi"] for o in in_window])
    ndmi = _stats([o["ndmi"] for o in in_window])
    ndvi_peak_obs = max(in_window, key=lambda o: o["ndvi"])
    ndmi_peak_obs = max(in_window, key=lambda o: o["ndmi"])

    complete_months = len(observed_months) == len(window_months)
    spans_window = in_window[0]["date"] <= start and in_window[-1]["date"] >= end
    if observed_months and complete_months and spans_window:
        status = "observed"
    else:
        status = "partially observed"

    return {
        "crop": window["crop"],
        "harvest_year": window["harvest_year"],
        "season_start": start,
        "season_end": end,
        "status": status,
        "n_obs": len(in_window),
        "months_with_obs": len(observed_months),
        "months_in_window": len(window_months),
        "observed_start": in_window[0]["date"],
        "observed_end": in_window[-1]["date"],
        "ndvi": ndvi,
        "ndmi": ndmi,
        "ndvi_peak": ndvi_peak_obs["ndvi"],
        "ndvi_peak_date": ndvi_peak_obs["date"],
        "ndmi_peak": ndmi_peak_obs["ndmi"],
        "ndmi_peak_date": ndmi_peak_obs["date"],
    }


# --------------------------------------------------------------------------
# Consistency text (descriptive only - never a classification)
# --------------------------------------------------------------------------

def consistency_statement(summary: dict) -> str:
    """Qualitative, hedged comparison of observed behaviour vs expectation."""
    if summary["status"] == "no observations":
        return (
            "No usable observations fall inside this window in the existing "
            "dataset, so the seasonal signature cannot be assessed. It is NOT "
            "inferred from other years."
        )
    ndvi = summary["ndvi"]
    ndmi = summary["ndmi"]
    text = (
        f"Observed NDVI rose to a peak of {summary['ndvi_peak']:.2f} on "
        f"{summary['ndvi_peak_date'].isoformat()} and ranged "
        f"{ndvi['min']:.2f}-{ndvi['max']:.2f} across {summary['n_obs']} day(s); "
        f"NDMI peaked at {summary['ndmi_peak']:.2f} on "
        f"{summary['ndmi_peak_date'].isoformat()}. This window is "
        f"**{summary['status']}** ({summary['months_with_obs']}/"
        f"{summary['months_in_window']} window months have data). The observed "
        "seasonal behaviour may be *consistent or inconsistent* with the "
        "reference label, but it is not confirmation: the reference label is "
        "user-provided ground truth, and cloud gaps plus incomplete seasonal "
        "coverage add substantial uncertainty."
    )
    return text


# --------------------------------------------------------------------------
# Plots
# --------------------------------------------------------------------------

def _profiles_by_tier(rows: list[dict]) -> dict:
    return {tier: monthly_profile(select_observations(rows, tier)) for tier in TIERS}


def plot_monthly_profile(profiles: dict, index: str, out_path: str | Path, dpi: int) -> Path:
    """Monthly NDVI or NDMI profile, one line per quality tier."""
    months = sorted({m["month"] for p in profiles.values() for m in p})
    xs = list(range(len(months)))
    fig, ax = plt.subplots(figsize=(12, 5))
    styles = {"all_valid": ("-", "o"), "complete": ("--", "s")}
    for tier in TIERS:
        lookup = {m["month"]: m[index]["mean"] for m in profiles[tier]}
        ys = [lookup.get(mmonth) for mmonth in months]
        line, marker = styles[tier]
        ax.plot(
            xs,
            ys,
            linestyle=line,
            marker=marker,
            markersize=4,
            linewidth=1.2,
            label=TIER_LABELS[tier],
        )
    ax.set_xticks(xs)
    ax.set_xticklabels(months, rotation=45, ha="right")
    ax.set_xlabel("Month")
    ax.set_ylabel(f"{index.upper()} monthly mean (unitless)")
    ax.set_title(
        "Parcel P001 - monthly mean "
        f"{index.upper()} of real observations (2023-09 to 2024-08)"
    )
    ax.set_ylim(-1.0, 1.0)
    ax.grid(True, linestyle=":", alpha=0.5)
    ax.legend(loc="best")
    fig.tight_layout()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    return out_path


def plot_seasonal(seasonal_by_tier: dict, out_path: str | Path, dpi: int) -> Path:
    """Seasonal mean NDVI (top) and NDMI (bottom) per reference crop window."""
    windows = [
        {"crop": s["crop"], "harvest_year": s["harvest_year"], "status": s["status"]}
        for s in seasonal_by_tier["all_valid"]
    ]
    labels = [f"{w['crop']}\n{w['harvest_year']}" for w in windows]
    xs = np.arange(len(windows))
    width = 0.38

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    for ax, index, title in ((ax1, "ndvi", "NDVI"), (ax2, "ndmi", "NDMI")):
        for offset, tier in zip((-width / 2, width / 2), TIERS):
            lookup = {
                (s["crop"], s["harvest_year"]): s
                for s in seasonal_by_tier[tier]
            }
            heights, positions = [], []
            for i, w in enumerate(windows):
                summary = lookup[(w["crop"], w["harvest_year"])]
                if summary[index] is None:
                    continue
                heights.append(summary[index]["mean"])
                positions.append(i + offset)
            ax.bar(positions, heights, width=width, label=TIER_LABELS[tier])
        ax.set_ylabel(f"{title} seasonal mean (unitless)")
        ax.set_ylim(-1.0, 1.0)
        ax.grid(True, axis="y", linestyle=":", alpha=0.5)
        # Mark windows with no observations.
        for i, w in enumerate(windows):
            if w["status"] == "no observations":
                ax.text(i, 0.0, "no\nobservations", ha="center", va="center", fontsize=8)
        ax.legend(loc="best")

    ax2.set_xticks(xs)
    ax2.set_xticklabels(labels)
    ax1.set_title("Parcel P001 - seasonal mean indices per reference crop window")
    fig.tight_layout()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    return out_path


# --------------------------------------------------------------------------
# CSV + report writers
# --------------------------------------------------------------------------

def _fmt(value, digits: int = 4, iso_date: bool = False) -> str:
    if value is None:
        return ""
    if iso_date:
        return value.isoformat()
    return f"{value:.{digits}f}"


def build_output_rows(rows: list[dict], crops: list[dict], parcel_id: str) -> list[dict]:
    """Tidy monthly + seasonal rows for both tiers (CSV records)."""
    records: list[dict] = []
    profiles = _profiles_by_tier(rows)
    for tier in TIERS:
        for m in profiles[tier]:
            records.append(
                {
                    "record_type": "monthly",
                    "parcel_id": parcel_id,
                    "tier": tier,
                    "crop": "",
                    "harvest_year": "",
                    "period": m["month"],
                    "season_start": "",
                    "season_end": "",
                    "status": "",
                    "n_obs": m["n_obs"],
                    "months_with_obs": "",
                    "months_in_window": "",
                    "observed_start": "",
                    "observed_end": "",
                    "ndvi_mean": _fmt(m["ndvi"]["mean"]),
                    "ndvi_median": _fmt(m["ndvi"]["median"]),
                    "ndvi_min": _fmt(m["ndvi"]["min"]),
                    "ndvi_max": _fmt(m["ndvi"]["max"]),
                    "ndvi_peak": "",
                    "ndvi_peak_date": "",
                    "ndmi_mean": _fmt(m["ndmi"]["mean"]),
                    "ndmi_median": _fmt(m["ndmi"]["median"]),
                    "ndmi_min": _fmt(m["ndmi"]["min"]),
                    "ndmi_max": _fmt(m["ndmi"]["max"]),
                    "ndmi_peak": "",
                    "ndmi_peak_date": "",
                    "notes": "",
                }
            )

    all_valid = {tier: select_observations(rows, tier) for tier in TIERS}
    for window in crops:
        for tier in TIERS:
            s = seasonal_summary(all_valid[tier], window)
            records.append(
                {
                    "record_type": "seasonal",
                    "parcel_id": parcel_id,
                    "tier": tier,
                    "crop": s["crop"],
                    "harvest_year": s["harvest_year"],
                    "period": f"{s['season_start'].isoformat()}..{s['season_end'].isoformat()}",
                    "season_start": s["season_start"].isoformat(),
                    "season_end": s["season_end"].isoformat(),
                    "status": s["status"],
                    "n_obs": s["n_obs"],
                    "months_with_obs": s["months_with_obs"],
                    "months_in_window": s["months_in_window"],
                    "observed_start": _fmt(s["observed_start"], iso_date=True),
                    "observed_end": _fmt(s["observed_end"], iso_date=True),
                    "ndvi_mean": "" if s["ndvi"] is None else _fmt(s["ndvi"]["mean"]),
                    "ndvi_median": "" if s["ndvi"] is None else _fmt(s["ndvi"]["median"]),
                    "ndvi_min": "" if s["ndvi"] is None else _fmt(s["ndvi"]["min"]),
                    "ndvi_max": "" if s["ndvi"] is None else _fmt(s["ndvi"]["max"]),
                    "ndvi_peak": _fmt(s["ndvi_peak"]),
                    "ndvi_peak_date": _fmt(s["ndvi_peak_date"], iso_date=True),
                    "ndmi_mean": "" if s["ndmi"] is None else _fmt(s["ndmi"]["mean"]),
                    "ndmi_median": "" if s["ndmi"] is None else _fmt(s["ndmi"]["median"]),
                    "ndmi_min": "" if s["ndmi"] is None else _fmt(s["ndmi"]["min"]),
                    "ndmi_max": "" if s["ndmi"] is None else _fmt(s["ndmi"]["max"]),
                    "ndmi_peak": _fmt(s["ndmi_peak"]),
                    "ndmi_peak_date": _fmt(s["ndmi_peak_date"], iso_date=True),
                    "notes": window.get("window_rule", ""),
                }
            )
    return records


def write_output_csv(records: list[dict], output_path: str | Path) -> Path:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        for record in records:
            writer.writerow(record)
    return output_path


def _coverage_audit_table(rows: list[dict]) -> list[str]:
    months: dict[str, dict] = {}
    for row in rows:
        key = (row.get("interval_from") or "")[:7]
        bucket = months.setdefault(
            key, {"total": 0, "valid": 0, "complete": 0, "partial": 0, "masked": 0}
        )
        bucket["total"] += 1
        if str(row.get("is_valid")) == "yes":
            bucket["valid"] += 1
            if row.get("coverage_status") == "complete":
                bucket["complete"] += 1
            else:
                bucket["partial"] += 1
        else:
            bucket["masked"] += 1
    lines = [
        "| Month | rows | valid | complete | partial | masked |",
        "|---|---|---|---|---|---|",
    ]
    for month in sorted(months):
        b = months[month]
        lines.append(
            f"| {month} | {b['total']} | {b['valid']} | {b['complete']} | "
            f"{b['partial']} | {b['masked']} |"
        )
    return lines


def build_report(
    rows: list[dict],
    crops: list[dict],
    parcel_id: str,
    qa_csv_path: str | Path,
    reference_csv_path: str | Path,
) -> str:
    profiles = _profiles_by_tier(rows)
    valid_all = select_observations(rows, "all_valid")
    valid_complete = select_observations(rows, "complete")
    seasonal = {tier: [seasonal_summary(select_observations(rows, tier), w) for w in crops]
                for tier in TIERS}

    span_start = valid_all[0]["date"].isoformat() if valid_all else "n/a"
    span_end = valid_all[-1]["date"].isoformat() if valid_all else "n/a"

    lines: list[str] = []
    lines.append("# Step 4.2 - Retrospective crop-signature analysis (parcel P001)")
    lines.append("")
    lines.append(f"- **Observed time series (read-only):** `{qa_csv_path}`")
    lines.append(f"- **Reference-crop config (read-only):** `{reference_csv_path}`")
    lines.append("- **Reproduce:** `python -m src.crop_signature`")
    lines.append(
        "- **Descriptive only:** no crop classification and no satellite-derived "
        "prediction is performed."
    )
    lines.append("")
    lines.append("## 0. What is observed vs labelled vs expected")
    lines.append("")
    lines.append(
        "- **Observations** are the real Sentinel-2 NDVI/NDMI values from Steps "
        "2-4.1, computed on real pixels only; nothing is interpolated or "
        "fabricated."
    )
    lines.append(
        "- **Reference labels** (winter wheat 2023, sunflower 2024, winter wheat "
        "2025) are **user-provided ground truth**, supplied in the auditable "
        "config `data/reference_crops.csv`; they are not derived from the "
        "imagery."
    )
    lines.append(
        "- **Agronomic expectations** are qualitative published phenology "
        "(Section 5); they are not data."
    )
    lines.append("")

    lines.append("## 1. Observation coverage audit")
    lines.append("")
    lines.append(
        f"- Real observations span **{span_start} to {span_end}** "
        f"({len(valid_all)} valid days; {len(valid_complete)} complete-coverage "
        "days)."
    )
    lines.append(
        "- The dataset therefore covers essentially only the **2024 growing "
        "season**. The 2023 winter-wheat and 2025 winter-wheat windows lie "
        "outside this span and have **no usable observations**."
    )
    lines.append("")
    lines.extend(_coverage_audit_table(rows))
    lines.append("")

    lines.append("## 2. Reference crop-season windows")
    lines.append("")
    lines.append(
        "| Crop | Harvest year | Window | Rule (assumption) | Status | n (all) | "
        "n (complete) | Months with data |"
    )
    lines.append("|---|---|---|---|---|---|---|---|")
    for i, window in enumerate(crops):
        s_all = seasonal["all_valid"][i]
        s_cmp = seasonal["complete"][i]
        lines.append(
            f"| {window['crop']} | {window['harvest_year']} | "
            f"{window['season_start'].isoformat()}..{window['season_end'].isoformat()} | "
            f"{window['window_rule']} | {s_all['status']} | {s_all['n_obs']} | "
            f"{s_cmp['n_obs']} | {s_all['months_with_obs']}/{s_all['months_in_window']} |"
        )
    lines.append("")

    lines.append("## 3. Monthly profiles (real observations only)")
    lines.append("")
    for tier in TIERS:
        lines.append(f"**Tier: {TIER_LABELS[tier]}**")
        lines.append("")
        lines.append("| Month | n | NDVI mean | NDVI median | NDVI min | NDVI max | "
                     "NDMI mean | NDMI median | NDMI min | NDMI max |")
        lines.append("|---|---|---|---|---|---|---|---|---|---|")
        for m in profiles[tier]:
            lines.append(
                f"| {m['month']} | {m['n_obs']} | {m['ndvi']['mean']:.3f} | "
                f"{m['ndvi']['median']:.3f} | {m['ndvi']['min']:.3f} | "
                f"{m['ndvi']['max']:.3f} | {m['ndmi']['mean']:.3f} | "
                f"{m['ndmi']['median']:.3f} | {m['ndmi']['min']:.3f} | "
                f"{m['ndmi']['max']:.3f} |"
            )
        lines.append("")

    lines.append("## 4. Seasonal summaries and consistency (hedged)")
    lines.append("")
    for i, window in enumerate(crops):
        s = seasonal["all_valid"][i]
        lines.append(f"### {window['crop']} {window['harvest_year']}")
        lines.append("")
        lines.append(f"- Window: {window['season_start'].isoformat()} to "
                     f"{window['season_end'].isoformat()} ({window['window_rule']}).")
        lines.append(f"- Status: **{s['status']}**.")
        if s["n_obs"]:
            lines.append(
                f"- Observed days: {s['n_obs']} ({s['observed_start'].isoformat()} "
                f"to {s['observed_end'].isoformat()}); months with data "
                f"{s['months_with_obs']}/{s['months_in_window']}."
            )
            lines.append(
                f"- NDVI mean/min/max {s['ndvi']['mean']:.3f}/"
                f"{s['ndvi']['min']:.3f}/{s['ndvi']['max']:.3f}; "
                f"NDMI mean/min/max {s['ndmi']['mean']:.3f}/"
                f"{s['ndmi']['min']:.3f}/{s['ndmi']['max']:.3f}."
            )
        lines.append(f"- Expected phenology ({window['reference_keys']}): "
                     f"{window['expected_phenology']}.")
        lines.append(f"- Assessment: {consistency_statement(s)}")
        lines.append("")

    lines.append("## 5. Limitations")
    lines.append("")
    lines.append(
        "- Only the 2024 season has real observations; **no cross-year crop "
        "comparison is possible** without new data (not fetched here)."
    )
    lines.append(
        "- Crop-season windows are **analytical assumptions**, not known sowing/"
        "harvest dates."
    )
    lines.append(
        "- Frequent cloud cover leaves multi-day and multi-week gaps; "
        "fully-masked days are absent, so the seasonal curves are sampled "
        "unevenly."
    )
    lines.append(
        "- The analysis reports whether observed behaviour is *consistent* with "
        "a reference label and never claims satellite confirmation."
    )
    lines.append("")

    lines.append("## 6. References")
    lines.append("")
    used_keys: list[str] = []
    for window in crops:
        used_keys.extend(window["reference_keys"])
    for key in dict.fromkeys(used_keys):
        citation = REFERENCE_LIBRARY.get(key)
        if citation:
            lines.append(f"- {citation}")
    lines.append("")
    return "\n".join(lines)


# --------------------------------------------------------------------------
# High-level workflow
# --------------------------------------------------------------------------

def run_crop_signature(
    qa_csv_path: str | Path = DEFAULT_QA_CSV_PATH,
    reference_csv_path: str | Path = DEFAULT_REFERENCE_CSV,
    output_dir: str | Path = DEFAULT_OUTPUT_DIR,
    parcel_id: str = DEFAULT_PARCEL_ID,
    dpi: int = DEFAULT_DPI,
) -> dict:
    qa_csv_path = Path(qa_csv_path)
    reference_csv_path = Path(reference_csv_path)
    output_dir = Path(output_dir)

    _fieldnames, rows = load_rows(qa_csv_path)
    crops = [c for c in load_reference_crops(reference_csv_path) if c["parcel_id"] == parcel_id]
    if not crops:
        raise ValueError(f"No reference crops for parcel {parcel_id}.")

    profiles = _profiles_by_tier(rows)
    seasonal = {
        tier: [seasonal_summary(select_observations(rows, tier), w) for w in crops]
        for tier in TIERS
    }

    parcel = parcel_id.lower()
    outputs = {
        "ndvi_monthly": plot_monthly_profile(
            profiles, "ndvi", output_dir / f"ndvi_monthly_profile_{parcel}_2023_2024.png", dpi
        ),
        "ndmi_monthly": plot_monthly_profile(
            profiles, "ndmi", output_dir / f"ndmi_monthly_profile_{parcel}_2023_2024.png", dpi
        ),
        "seasonal": plot_seasonal(
            seasonal, output_dir / f"crop_signature_seasonal_{parcel}.png", dpi
        ),
        "csv": write_output_csv(
            build_output_rows(rows, crops, parcel_id),
            output_dir / f"crop_signature_{parcel}_2023_2024.csv",
        ),
        "report": output_dir / f"crop_signature_report_{parcel}.md",
    }
    outputs["report"].write_text(
        build_report(rows, crops, parcel_id, qa_csv_path, reference_csv_path),
        encoding="utf-8",
    )
    return outputs


# --------------------------------------------------------------------------
# Offline self-check (no network, no files written)
# --------------------------------------------------------------------------

def _self_check() -> None:
    # Reference-crop config parsing.
    text = (
        "parcel_id,harvest_year,crop,season_start,season_end,window_rule,"
        "expected_phenology,reference_keys\n"
        "P001,2024,Sunflower,2024-04-01,2024-09-30,Apr-Sep,peak in summer,"
        "schneiter1981;gao1996\n"
    )
    crops = parse_reference_crops(text)
    assert len(crops) == 1
    assert crops[0]["crop"] == "Sunflower"
    assert crops[0]["season_start"] == date(2024, 4, 1)
    assert crops[0]["reference_keys"] == ["schneiter1981", "gao1996"]

    # Tier filtering respects the QA flags.
    valid_complete = {
        "is_valid": "yes", "coverage_status": "complete", "valid_count": "96",
        "interval_from": "2024-06-15T00:00:00Z", "ndvi_mean": "0.80", "ndmi_mean": "0.38",
    }
    valid_partial = dict(valid_complete, coverage_status="partial")
    masked = {
        "is_valid": "no", "coverage_status": "masked", "valid_count": "0",
        "interval_from": "2024-06-16T00:00:00Z", "ndvi_mean": "", "ndmi_mean": "",
    }
    rows = [valid_complete, valid_partial, masked]
    assert len(select_observations(rows, "all_valid")) == 2
    assert len(select_observations(rows, "complete")) == 1

    # Monthly aggregation uses only real observations.
    profile = monthly_profile(select_observations(rows, "all_valid"))
    assert len(profile) == 1 and profile[0]["month"] == "2024-06"
    assert profile[0]["n_obs"] == 2
    assert abs(profile[0]["ndvi"]["mean"] - 0.80) < 1e-9

    # Seasonal windows: no observations / partially observed, peak detection.
    window_none = {
        "crop": "Winter wheat", "harvest_year": 2023,
        "season_start": date(2022, 10, 1), "season_end": date(2023, 7, 15),
    }
    empty = seasonal_summary(select_observations(rows, "all_valid"), window_none)
    assert empty["status"] == "no observations" and empty["n_obs"] == 0
    assert "cannot be assessed" in consistency_statement(empty)

    window_partial = {
        "crop": "Sunflower", "harvest_year": 2024,
        "season_start": date(2024, 4, 1), "season_end": date(2024, 9, 30),
    }
    partial = seasonal_summary(select_observations(rows, "all_valid"), window_partial)
    assert partial["status"] == "partially observed", partial["status"]
    assert partial["n_obs"] == 2
    assert partial["ndvi_peak"] == 0.80

    # Fully-covered window is reported as observed.
    window_full = {
        "crop": "Test", "harvest_year": 2024,
        "season_start": date(2024, 6, 1), "season_end": date(2024, 6, 15),
    }
    filled = seasonal_summary(
        select_observations(rows + [dict(valid_complete, interval_from="2024-06-01T00:00:00Z")], "all_valid"),
        window_full,
    )
    assert filled["status"] == "observed", filled["status"]

    # Output rows: one monthly set per tier + one seasonal per crop per tier.
    records = build_output_rows(rows, crops, "P001")
    assert any(r["record_type"] == "monthly" for r in records)
    assert sum(1 for r in records if r["record_type"] == "seasonal") == 2  # 1 crop x 2 tiers
    assert all(set(r.keys()) == set(OUTPUT_COLUMNS) for r in records)

    # Report mentions the crucial separation and limitation concepts.
    report = build_report(rows, crops, "P001", "qa.csv", "ref.csv")
    assert "user-provided ground truth" in report
    assert "no usable observations" in report
    assert "consistency" in report.lower()

    print("crop_signature.py self-check passed (offline, no files written).")


# --------------------------------------------------------------------------
# Command-line interface
# --------------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m src.crop_signature",
        description="Step 4.2 retrospective crop-signature analysis (descriptive).",
    )
    parser.add_argument("--check", action="store_true",
                        help="Run the offline self-check and exit (no files written).")
    parser.add_argument("--csv", default=str(DEFAULT_QA_CSV_PATH),
                        help=f"QA time-series CSV (default: {DEFAULT_QA_CSV_PATH}).")
    parser.add_argument("--reference", default=str(DEFAULT_REFERENCE_CSV),
                        help=f"Reference-crop config (default: {DEFAULT_REFERENCE_CSV}).")
    parser.add_argument("--outdir", default=str(DEFAULT_OUTPUT_DIR),
                        help=f"Output directory (default: {DEFAULT_OUTPUT_DIR}).")
    parser.add_argument("--parcel", default=DEFAULT_PARCEL_ID,
                        help=f"Parcel id (default: {DEFAULT_PARCEL_ID}).")
    parser.add_argument("--dpi", type=int, default=DEFAULT_DPI,
                        help=f"PNG resolution (default: {DEFAULT_DPI}).")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    if args.check:
        _self_check()
        return 0

    outputs = run_crop_signature(
        qa_csv_path=args.csv,
        reference_csv_path=args.reference,
        output_dir=args.outdir,
        parcel_id=args.parcel,
        dpi=args.dpi,
    )
    for name, path in outputs.items():
        print(f"Wrote {name}: {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
