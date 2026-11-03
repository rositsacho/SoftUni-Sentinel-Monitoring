"""
plots.py - Sentinel-2 seasonal NDVI / NDMI time-series plots for one parcel.

Reads the daily per-parcel CSV produced by src.sentinel_stats (Steps 2/3) and
renders one PNG per index (NDVI, NDMI). Only valid observations are plotted:

    * valid_count > 0, and
    * the index value for that row is present.

Observations are drawn at their true acquisition dates. There is deliberately
no interpolation, resampling or fabrication of missing values - fully masked
and scene-less days are simply absent from the series.

The source CSV is opened read-only and is never modified, and this module
does not touch the satellite-processing methodology (evalscript, SCL masking,
coverage calculations).

Run the offline self-check (logic only, no files written, no network):
    python -m src.plots --check

Generate the NDVI and NDMI PNGs:
    python -m src.plots
"""

from __future__ import annotations

import argparse
import csv
import sys
from datetime import datetime
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")  # headless, reproducible: never opens a GUI window

import matplotlib.pyplot as plt  # noqa: E402  (must follow matplotlib.use)
from matplotlib.dates import AutoDateLocator, ConciseDateFormatter  # noqa: E402

# --------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CSV_PATH = PROJECT_ROOT / "data" / "results" / "time_series_p001_2023_2024.csv"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "data" / "results"
DEFAULT_PARCEL_ID = "P001"
DEFAULT_DPI = 150

# One "figure" definition per index: value column, label text and line colour.
INDEX_SPECS = {
    "ndvi": {"value_column": "ndvi_mean", "ylabel": "NDVI (unitless)", "color": "#2e7d32"},
    "ndmi": {"value_column": "ndmi_mean", "ylabel": "NDMI (unitless)", "color": "#1565c0"},
}


# --------------------------------------------------------------------------
# Parsing helpers
# --------------------------------------------------------------------------

def parse_timestamp(value: str) -> datetime:
    """Parse an ISO-8601 UTC timestamp such as '2023-09-02T00:00:00Z'."""
    text = (value or "").strip()
    if not text:
        raise ValueError("Empty timestamp.")
    # datetime.fromisoformat handles 'Z' from Python 3.11; normalise anyway.
    if text.endswith("Z"):
        text = text[:-1] + "+00:00"
    return datetime.fromisoformat(text)


def _to_float(value) -> float | None:
    """Parse a numeric CSV cell; blank / unparsable cells become None."""
    if value is None:
        return None
    text = str(value).strip()
    if text == "":
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    return number


# --------------------------------------------------------------------------
# Data selection (pure, offline-testable)
# --------------------------------------------------------------------------

def collect_series(
    rows: Iterable[dict],
    value_column: str,
) -> tuple[list[datetime], list[float], dict]:
    """Select the plottable points for one index from parsed CSV rows.

    Returns (dates, values, stats) where the series is sorted chronologically
    and never interpolated. Rows are excluded, with a reason, when:

        * valid_count is missing          -> stats["excluded_missing_valid_count"]
        * valid_count <= 0 (fully masked) -> stats["excluded_masked"]
        * the index value is blank        -> stats["excluded_missing_value"]
    """
    dates: list[datetime] = []
    values: list[float] = []
    stats = {
        "rows_total": 0,
        "plotted": 0,
        "excluded_masked": 0,
        "excluded_missing_valid_count": 0,
        "excluded_missing_value": 0,
    }

    for row in rows:
        stats["rows_total"] += 1

        valid_count = _to_float(row.get("valid_count"))
        if valid_count is None:
            stats["excluded_missing_valid_count"] += 1
            continue
        if valid_count <= 0:
            stats["excluded_masked"] += 1
            continue

        value = _to_float(row.get(value_column))
        stamp = (row.get("interval_from") or "").strip()
        if value is None or not stamp:
            stats["excluded_missing_value"] += 1
            continue

        dates.append(parse_timestamp(stamp))
        values.append(value)

    ordered = sorted(zip(dates, values), key=lambda pair: pair[0])
    stats["plotted"] = len(ordered)
    return [d for d, _ in ordered], [v for _, v in ordered], stats


def load_valid_series(
    csv_path: str | Path,
    value_column: str,
) -> tuple[list[datetime], list[float], dict]:
    """Read a time-series CSV (read-only) and select one index's series."""
    csv_path = Path(csv_path)
    if not csv_path.exists():
        raise FileNotFoundError(f"Input CSV not found: {csv_path}")
    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    return collect_series(rows, value_column)


# --------------------------------------------------------------------------
# Plot rendering
# --------------------------------------------------------------------------

def _season_span(dates: list[datetime]) -> str:
    """Human-readable 'YYYY-MM-DD to YYYY-MM-DD' span of the plotted points."""
    if not dates:
        return "no observations"
    return f"{dates[0].date().isoformat()} to {dates[-1].date().isoformat()}"


def render_plot(
    dates: list[datetime],
    values: list[float],
    *,
    title: str,
    ylabel: str,
    color: str,
    output_path: str | Path,
    dpi: int = DEFAULT_DPI,
) -> Path:
    """Render a single clean time-series line+marker plot to a PNG file."""
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.plot(
        dates,
        values,
        marker="o",
        markersize=4,
        linewidth=1.2,
        color=color,
        label="Valid observations (valid_count > 0)",
    )

    locator = AutoDateLocator(minticks=4, maxticks=12)
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(ConciseDateFormatter(locator))

    ax.set_xlabel("Date (UTC)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_ylim(-1.0, 1.0)
    ax.grid(True, linestyle=":", alpha=0.5)
    ax.legend(loc="best")
    fig.tight_layout()

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)
    return output_path


def build_time_series_plots(
    csv_path: str | Path = DEFAULT_CSV_PATH,
    output_dir: str | Path = DEFAULT_OUTPUT_DIR,
    parcel_id: str = DEFAULT_PARCEL_ID,
    dpi: int = DEFAULT_DPI,
) -> dict:
    """Build the separate NDVI and NDMI PNGs and return a summary dict."""
    csv_path = Path(csv_path)
    output_dir = Path(output_dir)
    parcel = parcel_id.lower()
    summary: dict = {"csv_path": csv_path, "outputs": {}, "stats": {}}

    for index, spec in INDEX_SPECS.items():
        dates, values, stats = load_valid_series(csv_path, spec["value_column"])
        summary["stats"][index] = stats

        if not dates:
            print(
                f"Warning: no valid {index.upper()} observations in {csv_path}; "
                "skipping its plot."
            )
            continue

        title = (
            f"Parcel {parcel_id} - {index.upper()} daily observations "
            f"({_season_span(dates)})"
        )
        out_path = output_dir / f"{index}_{parcel}_2023_2024.png"
        render_plot(
            dates,
            values,
            title=title,
            ylabel=spec["ylabel"],
            color=spec["color"],
            output_path=out_path,
            dpi=dpi,
        )
        summary["outputs"][index] = out_path
        print(f"{index.upper()}: plotted {stats['plotted']} observation(s) -> {out_path}")

    return summary


# --------------------------------------------------------------------------
# Offline self-check (no network, no files written)
# --------------------------------------------------------------------------

def _self_check() -> None:
    # Timestamp parsing: 'Z' and explicit offset must both work.
    assert parse_timestamp("2023-09-02T00:00:00Z") == parse_timestamp(
        "2023-09-02T00:00:00+00:00"
    )
    assert parse_timestamp("2024-08-29T00:00:00Z").year == 2024

    # Numeric parsing: blank cells are missing, not zero.
    assert _to_float("") is None
    assert _to_float("  ") is None
    assert _to_float("0.71") == 0.71
    assert _to_float("0") == 0.0

    sample_rows = [
        {  # valid -> plotted
            "interval_from": "2023-09-02T00:00:00Z",
            "valid_count": "96",
            "ndvi_mean": "0.63",
            "ndmi_mean": "0.23",
        },
        {  # out of order on purpose; valid -> plotted, must be sorted
            "interval_from": "2023-09-05T00:00:00Z",
            "valid_count": "96",
            "ndvi_mean": "0.70",
            "ndmi_mean": "0.29",
        },
        {  # fully masked -> excluded as masked
            "interval_from": "2023-09-08T00:00:00Z",
            "valid_count": "0",
            "ndvi_mean": "",
            "ndmi_mean": "",
        },
        {  # blank valid_count -> excluded (unknown coverage)
            "interval_from": "2023-09-10T00:00:00Z",
            "valid_count": "",
            "ndvi_mean": "0.5",
            "ndmi_mean": "0.1",
        },
        {  # valid coverage but missing index value -> excluded as missing
            "interval_from": "2023-09-12T00:00:00Z",
            "valid_count": "96",
            "ndvi_mean": "",
            "ndmi_mean": "0.15",
        },
    ]

    dates, values, stats = collect_series(sample_rows, "ndvi_mean")
    assert stats["rows_total"] == 5
    assert stats["plotted"] == 2, stats
    assert stats["excluded_masked"] == 1, stats
    assert stats["excluded_missing_valid_count"] == 1, stats
    assert stats["excluded_missing_value"] == 1, stats
    # Points keep their true dates (no interpolation) and are chronological.
    assert [d.date().isoformat() for d in dates] == ["2023-09-02", "2023-09-05"]
    assert values == [0.63, 0.70]

    # NDMI uses the same rows but a different value column: the last row has
    # a valid NDMI even though its NDVI was blank, so it is plotted here.
    _, ndmi_values, ndmi_stats = collect_series(sample_rows, "ndmi_mean")
    assert ndmi_stats["plotted"] == 3, ndmi_stats
    assert ndmi_stats["excluded_missing_value"] == 0, ndmi_stats
    assert ndmi_values == [0.23, 0.29, 0.15]

    # Span helper.
    assert _season_span(dates) == "2023-09-02 to 2023-09-05"
    assert _season_span([]) == "no observations"

    print("plots.py self-check passed (offline, no files written).")


# --------------------------------------------------------------------------
# Command-line interface
# --------------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m src.plots",
        description="Render NDVI/NDMI time-series PNGs from a parcel CSV.",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Run the offline self-check and exit (no plots written).",
    )
    parser.add_argument(
        "--csv",
        default=str(DEFAULT_CSV_PATH),
        help=f"Input time-series CSV (default: {DEFAULT_CSV_PATH}).",
    )
    parser.add_argument(
        "--outdir",
        default=str(DEFAULT_OUTPUT_DIR),
        help=f"Directory for the PNGs (default: {DEFAULT_OUTPUT_DIR}).",
    )
    parser.add_argument(
        "--parcel",
        default=DEFAULT_PARCEL_ID,
        help=f"Parcel id (used in titles/file names; default: {DEFAULT_PARCEL_ID}).",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=DEFAULT_DPI,
        help=f"PNG resolution in dots per inch (default: {DEFAULT_DPI}).",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    if args.check:
        _self_check()
        return 0

    build_time_series_plots(
        csv_path=args.csv,
        output_dir=args.outdir,
        parcel_id=args.parcel,
        dpi=args.dpi,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
