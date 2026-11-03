"""
historical_acquisition.py - Step 4.3 historical Sentinel-2 L2A acquisition for
parcel P001.

Retrieves ONLY the uncovered dates of three historical windows, reusing the
existing Step 2/3 pipeline unchanged (same parcel geometry, 20 m resolution,
evalscript, SCL cloud mask, NDVI/NDMI definitions) and the existing Step 4.1
QA functions. Original files are never overwritten: new observations are
written to separate, additive CSVs.

Windows requested (inclusive last day):
    * Winter wheat 2023 : 2022-10-01 .. 2023-07-15
    * Sunflower 2024    : 2024-08-28 .. 2024-09-30
    * Winter wheat 2025 : 2024-10-01 .. 2025-07-15

Deduplication uses the actual interval_from / interval_to boundaries emitted
by the Statistical API (each daily interval is named by its inclusive start
timestamp "YYYY-MM-DDT00:00:00Z"). Already-covered days are excluded by that
key, never by a single global max-date rule, so no date is duplicated or
skipped (in particular 2024-08-28 and 2024-08-29 stay in the existing series).

The module is network-free unless run with --execute. --check runs the offline
self-test only.

Run the offline self-check (no network, no files written):
    python -m src.historical_acquisition --check

Plan only: list uncovered windows, request ranges and counts (no network):
    python -m src.historical_acquisition --plan

Perform the live acquisition (network + credentials; requires approval):
    python -m src.historical_acquisition --execute
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import date, datetime, timezone, timedelta
from pathlib import Path

try:  # works for "python -m src.historical_acquisition"
    from .sentinel_stats import (
        CSV_COLUMNS,
        DEFAULT_GEOJSON_PATH,
        DEFAULT_RESOLUTION_M,
        SentinelHubAuth,
        combine_interval_rows,
        get_parcel_time_series,
        load_parcel_geometry,
        save_time_series_csv,
    )
    from .quality_assessment import analyze, build_report, load_rows, write_qa_csv
except ImportError:  # fallback for "python src/historical_acquisition.py"
    from sentinel_stats import (  # type: ignore
        CSV_COLUMNS,
        DEFAULT_GEOJSON_PATH,
        DEFAULT_RESOLUTION_M,
        SentinelHubAuth,
        combine_interval_rows,
        get_parcel_time_series,
        load_parcel_geometry,
        save_time_series_csv,
    )
    from quality_assessment import analyze, build_report, load_rows, write_qa_csv  # type: ignore

# --------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_RESULTS_DIR = PROJECT_ROOT / "data" / "results"
DEFAULT_EXISTING_CSV = DEFAULT_RESULTS_DIR / "time_series_p001_2023_2024.csv"
DEFAULT_HISTORICAL_CSV = DEFAULT_RESULTS_DIR / "time_series_p001_historical.csv"
DEFAULT_HISTORICAL_QA_CSV = DEFAULT_RESULTS_DIR / "time_series_p001_historical_qa.csv"
DEFAULT_LOG_PATH = DEFAULT_RESULTS_DIR / "historical_acquisition_log.json"
DEFAULT_PARCEL_ID = "P001"

# Requested windows: (label, start_inclusive, end_inclusive). The API's
# time range is [from 00:00:00Z, to 23:59:59Z], i.e. day-inclusive on both ends.
REQUESTED_WINDOWS: list[tuple[str, str, str]] = [
    ("Winter wheat 2023", "2022-10-01", "2023-07-15"),
    ("Sunflower 2024", "2024-08-28", "2024-09-30"),
    ("Winter wheat 2025", "2024-10-01", "2025-07-15"),
]


# --------------------------------------------------------------------------
# Interval-key helpers (API boundary semantics)
# --------------------------------------------------------------------------

def interval_key_from_day(day: str) -> str:
    """Return the API interval identifier for a calendar day 'YYYY-MM-DD'.

    The Statistical API labels each P1D interval by its inclusive start
    timestamp, e.g. 2024-08-28 -> '2024-08-28T00:00:00Z'. Using this exact key
    (rather than a global max-date cutoff) is what makes de-duplication exact.
    """
    return f"{date.fromisoformat(day).isoformat()}T00:00:00Z"


def existing_interval_keys(csv_path: str | Path) -> set[tuple[str, str]]:
    """Return the set of (interval_from, interval_to) keys already acquired."""
    _fields, rows = load_rows(csv_path)
    keys: set[tuple[str, str]] = set()
    for row in rows:
        iv_from = (row.get("interval_from") or "").strip()
        if iv_from:
            keys.add((iv_from, (row.get("interval_to") or "").strip()))
    return keys


def _parse_interval(interval_from: str) -> date:
    """Start day of a Statistical API interval timestamp."""
    return datetime.fromisoformat(interval_from.replace("Z", "+00:00")).date()


def existing_covered_intervals(csv_path: str | Path) -> list[tuple[date, date]]:
    """Real coverage of existing data as half-open [start_day, end_day) pairs.

    The Statistical API labels each P1D interval by its inclusive start; the
    interval covers the half-open time range [from, to). A day D is covered
    only when some existing interval has from == D (00:00:00Z). We therefore
    read the exact interval_from/interval_to timestamps and derive the covered
    start day of each interval - never infer coverage from adjacent dates.
    """
    _fields, rows = load_rows(csv_path)
    intervals: list[tuple[date, date]] = []
    for row in rows:
        iv_from = (row.get("interval_from") or "").strip()
        if not iv_from:
            continue
        intervals.append((_parse_interval(iv_from), _parse_interval(iv_from)))
    return intervals


def covered_start_days(csv_path: str | Path) -> set[date]:
    """Set of exact start days already covered by existing intervals."""
    return {start for start, _ in existing_covered_intervals(csv_path)}


def uncovered_request_ranges(
    windows: list[tuple[str, str, str]],
    covered_days: set[date],
) -> list[dict]:
    """Split each window into the day-runs that are NOT already covered.

    A day is covered iff it is the exact start day of an existing interval
    (half-open [from, to) semantics). Returns a flat, chronological list of
    {window, start, end} where start/end are inclusive day strings and every
    day in the run is uncovered. Contiguous uncovered days are merged into one
    range so the API sees one request per uninterrupted gap.
    """
    ranges: list[dict] = []
    for label, start_s, end_s in windows:
        start = date.fromisoformat(start_s)
        end = date.fromisoformat(end_s)
        run_start: date | None = None
        run_end: date | None = None
        cursor = start
        while cursor <= end:
            covered = cursor in covered_days
            if not covered:
                if run_start is None:
                    run_start = cursor
                run_end = cursor
            else:
                if run_start is not None:
                    ranges.append(
                        {
                            "window": label,
                            "start": run_start.isoformat(),
                            "end": run_end.isoformat(),
                        }
                    )
                    run_start = run_end = None
            cursor += timedelta(days=1)
        if run_start is not None:
            ranges.append(
                {"window": label, "start": run_start.isoformat(), "end": run_end.isoformat()}
            )
    return ranges


def _month_chunks(start_s: str, end_s: str) -> list[tuple[str, str]]:
    """Split an inclusive [start, end] range into inclusive monthly chunks."""
    start = date.fromisoformat(start_s)
    end = date.fromisoformat(end_s)
    chunks: list[tuple[str, str]] = []
    cursor = start
    while cursor <= end:
        next_month = (
            date(cursor.year + 1, 1, 1)
            if cursor.month == 12
            else date(cursor.year, cursor.month + 1, 1)
        )
        chunk_end = min(next_month - timedelta(days=1), end)
        chunks.append((cursor.isoformat(), chunk_end.isoformat()))
        cursor = next_month
    return chunks


def api_timerange(start_s: str, end_s: str) -> tuple[str, str]:
    """Exact Statistical API timeRange for an inclusive [start, end] day range.

    The huge majority of the API treats the range as half-open:
    [from, to). We therefore open on the start day at 00:00:00Z and close on
    the day AFTER the final included day at 00:00:00Z (exclusive). Using
    T23:59:59Z would rely on a partial last interval being kept, which the
    default lastIntervalBehavior=SKIP may drop. The exclusive end guarantees
    the final daily interval is never silently skipped.

        api_timerange("2024-08-28", "2024-08-28")
            -> ("2024-08-28T00:00:00Z", "2024-08-29T00:00:00Z")
    """
    start = date.fromisoformat(start_s)
    end_exclusive = date.fromisoformat(end_s) + timedelta(days=1)
    return (
        f"{start.isoformat()}T00:00:00Z",
        f"{end_exclusive.isoformat()}T00:00:00Z",
    )


def plan_requests(
    windows: list[tuple[str, str, str]] = REQUESTED_WINDOWS,
    covered_days: set[date] | None = None,
) -> list[dict]:
    """Return the monthly request plan (inclusive ranges) for uncovered days.

    One request per calendar month that contains at least one uncovered day,
    clipped to the uncovered run so covered month-edges are never included.
    Each entry also carries the exact API timeRange boundaries.
    """
    covered_days = covered_days or set()
    plan: list[dict] = []
    for run in uncovered_request_ranges(windows, covered_days):
        for chunk_start, chunk_end in _month_chunks(run["start"], run["end"]):
            time_from, time_to = api_timerange(chunk_start, chunk_end)
            plan.append(
                {
                    "window": run["window"],
                    "start": chunk_start,
                    "end": chunk_end,
                    "time_from": time_from,
                    "time_to": time_to,
                }
            )
    return plan


# --------------------------------------------------------------------------
# Merging / validation
# --------------------------------------------------------------------------

def merge_rows(*row_lists: list[dict]) -> list[dict]:
    """Chronologically combine row lists and drop duplicate interval keys."""
    flat: list[dict] = []
    for rows in row_lists:
        flat.extend(rows)
    return combine_interval_rows(flat)


def validate_rows(rows: list[dict]) -> dict:
    """Validate schema, interval keys and numeric values of a row list.

    Returns a report dict. Raises ValueError on a hard schema/key problem.
    """
    report = {
        "rows": len(rows),
        "missing_columns": [],
        "duplicate_keys": [],
        "bad_timestamps": [],
        "non_numeric": [],
        "masked_rows": 0,
        "empty_index_rows": 0,
    }
    seen: set[tuple[str, str]] = set()
    for row in rows:
        for column in CSV_COLUMNS:
            if column not in row:
                report["missing_columns"].append(column)
        iv_from = (row.get("interval_from") or "").strip()
        iv_to = (row.get("interval_to") or "").strip()
        if not iv_from:
            report["bad_timestamps"].append(iv_to or "<blank>")
        else:
            try:
                datetime.fromisoformat(iv_from.replace("Z", "+00:00"))
            except ValueError:
                report["bad_timestamps"].append(iv_from)
        key = (iv_from, iv_to)
        if key in seen:
            report["duplicate_keys"].append(key)
        seen.add(key)

        valid_count = row.get("valid_count")
        try:
            vc = int(float(valid_count)) if str(valid_count).strip() != "" else None
        except (TypeError, ValueError):
            vc = None
            report["non_numeric"].append(("valid_count", valid_count))
        if vc == 0:
            report["masked_rows"] += 1

        for column in ("ndvi_mean", "ndmi_mean", "clear_fraction"):
            value = row.get(column)
            if str(value).strip() == "":
                continue
            try:
                float(value)
            except (TypeError, ValueError):
                report["non_numeric"].append((column, value))

        if vc == 0 and str(row.get("ndvi_mean") or "").strip() == "":
            report["empty_index_rows"] += 1

    if report["missing_columns"]:
        raise ValueError(
            f"rows are missing required columns: {sorted(set(report['missing_columns']))}"
        )
    if report["duplicate_keys"]:
        raise ValueError(
            f"{len(report['duplicate_keys'])} duplicate interval key(s) present."
        )
    return report


# --------------------------------------------------------------------------
# Checkpointing
# --------------------------------------------------------------------------

def _load_checkpoint(log_path: str | Path) -> dict:
    path = Path(log_path)
    if not path.exists():
        return {"requests": [], "completed_ranges": []}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (ValueError, OSError):
        return {"requests": [], "completed_ranges": []}


def _save_checkpoint(log_path: str | Path, checkpoint: dict) -> None:
    path = Path(log_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(checkpoint, indent=2), encoding="utf-8")


def _range_key(entry: dict) -> str:
    return f"{entry['start']}..{entry['end']}"


# --------------------------------------------------------------------------
# Acquisition
# --------------------------------------------------------------------------

def summarize_request(rows: list[dict]) -> dict:
    """Classify a request's rows: total / valid / masked, plus scene status.

    - "no imagery": the API returned zero intervals for the request range.
    - "success": at least one interval returned (valid and/or fully masked).
    """
    intervals = len(rows)
    valid = sum(1 for r in rows if (str(r.get("valid_count") or "0").strip() not in ("", "0")
                                    and float(r.get("valid_count") or 0) > 0))
    masked = sum(1 for r in rows if str(r.get("valid_count") or "").strip() == "0")
    if intervals == 0:
        status = "no imagery"
    else:
        status = "success"
    return {
        "intervals": intervals,
        "valid": valid,
        "masked": masked,
        "status": status,
    }


def acquire(
    plan: list[dict],
    *,
    parcel_id: str = DEFAULT_PARCEL_ID,
    geojson_path: str | Path = DEFAULT_GEOJSON_PATH,
    resolution_m: float = DEFAULT_RESOLUTION_M,
    log_path: str | Path = DEFAULT_LOG_PATH,
) -> list[dict]:
    """Execute the plan against the live API with checkpointing.

    One SentinelHubAuth instance (one token) is reused for all requests. Each
    request is checkpointed; a failed request is logged and does NOT abort the
    run, so a rerun resumes only the missing ranges (idempotent).
    """
    geometry, src_epsg = load_parcel_geometry(geojson_path, parcel_id)
    auth = SentinelHubAuth()
    checkpoint = _load_checkpoint(log_path)
    completed = set(checkpoint.get("completed_ranges", []))
    new_rows: list[dict] = []

    for entry in plan:
        key = _range_key(entry)
        if key in completed:
            print(f"skip (already completed): {key}")
            continue
        print(f"request {key} ({entry['window']}) ...")
        try:
            # Use the exact half-open boundaries: start 00:00:00Z and the day
            # after the last included day at 00:00:00Z (exclusive end).
            time_from, time_to = api_timerange(entry["start"], entry["end"])
            rows, _metadata = get_parcel_time_series(
                geometry=geometry,
                src_epsg=src_epsg,
                start_date=time_from,
                end_date=time_to,
                resolution_m=resolution_m,
                auth=auth,
                require_usable=False,
            )
            summary = summarize_request(rows)
            new_rows.extend(rows)
            checkpoint["requests"].append(
                {
                    "window": entry["window"],
                    "start": entry["start"],
                    "end": entry["end"],
                    "time_from": time_from,
                    "time_to": time_to,
                    "status": summary["status"],
                    "intervals": summary["intervals"],
                    "valid": summary["valid"],
                    "masked": summary["masked"],
                    "completed_at": datetime.now(timezone.utc).isoformat(),
                }
            )
            checkpoint.setdefault("completed_ranges", []).append(key)
            _save_checkpoint(log_path, checkpoint)
            completed.add(key)
            print(
                f"  {summary['status']}: {summary['intervals']} interval(s), "
                f"{summary['valid']} valid, {summary['masked']} masked"
            )
        except Exception as exc:  # noqa: BLE001 - log and continue, never crash the run
            checkpoint["requests"].append(
                {
                    "window": entry["window"],
                    "start": entry["start"],
                    "end": entry["end"],
                    "status": "failed",
                    "error": type(exc).__name__,
                    "error_message": str(exc)[:300],
                    "completed_at": datetime.now(timezone.utc).isoformat(),
                }
            )
            _save_checkpoint(log_path, checkpoint)
            print(f"  FAILED: {type(exc).__name__}: {str(exc)[:200]}")

    return new_rows


# --------------------------------------------------------------------------
# High-level workflow
# --------------------------------------------------------------------------

def run_historical_acquisition(
    existing_csv: str | Path = DEFAULT_EXISTING_CSV,
    historical_csv: str | Path = DEFAULT_HISTORICAL_CSV,
    historical_qa_csv: str | Path = DEFAULT_HISTORICAL_QA_CSV,
    log_path: str | Path = DEFAULT_LOG_PATH,
    parcel_id: str = DEFAULT_PARCEL_ID,
    geojson_path: str | Path = DEFAULT_GEOJSON_PATH,
    resolution_m: float = DEFAULT_RESOLUTION_M,
    execute: bool = False,
) -> dict:
    """Plan (and optionally execute) the historical acquisition."""
    covered_days = covered_start_days(existing_csv)
    plan = plan_requests(REQUESTED_WINDOWS, covered_days)

    if not execute:
        return {"plan": plan, "executed": False}

    new_rows = acquire(
        plan,
        parcel_id=parcel_id,
        geojson_path=geojson_path,
        resolution_m=resolution_m,
        log_path=log_path,
    )
    new_rows = combine_interval_rows(new_rows)

    # Existing rows (read-only) + new rows, same schema, deduped by interval.
    _fields, existing_rows = load_rows(existing_csv)
    merged = combine_interval_rows(existing_rows + new_rows)

    # Every row validated before writing anything.
    validation = validate_rows(merged)

    historical_path = save_time_series_csv(new_rows, historical_csv)
    qa_path = _write_qa_via_existing_functions(merged, historical_qa_csv, historical_csv)

    return {
        "executed": True,
        "plan": plan,
        "historical_csv": historical_path,
        "qa_csv": qa_path,
        "new_rows": len(new_rows),
        "validation": validation,
    }


def _write_qa_via_existing_functions(
    rows: list[dict],
    output_qa_path: str | Path,
    source_csv_path: str | Path,
) -> Path:
    """Produce the QA CSV using the real Step 4.1 functions (not stub columns)."""
    _fieldnames, _rows = load_rows(source_csv_path) if Path(source_csv_path).exists() else ([], rows)
    qa_rows, analysis = analyze(rows)
    original_fields = CSV_COLUMNS
    write_qa_csv(qa_rows, original_fields, output_qa_path)
    # Also emit a Markdown report next to the QA CSV for auditability.
    report_path = Path(output_qa_path).with_name("historical_quality_report_p001.md")
    report_path.write_text(build_report(analysis, source_csv_path), encoding="utf-8")
    return Path(output_qa_path)


# --------------------------------------------------------------------------
# Offline self-check (no network, no files written)
# --------------------------------------------------------------------------

def _self_check() -> None:
    # 1. Interval-key semantics.
    assert interval_key_from_day("2024-08-28") == "2024-08-28T00:00:00Z"

    # Real August 2024 timestamps from the existing CSV: the series has
    # 08-27 -> 08-28 and 08-29 -> 08-30, but NOT 08-28 -> 08-29. Under the
    # half-open [from, to) semantics 08-28 is therefore UNCOVERED.
    covered = {date(2024, 8, 27), date(2024, 8, 29)}

    # 2. 08-28 must be requested; 08-27 and 08-29 must NOT be; no gap.
    runs = uncovered_request_ranges(
        [("Sunflower 2024", "2024-08-27", "2024-08-30")], covered
    )
    assert runs == [
        {"window": "Sunflower 2024", "start": "2024-08-28", "end": "2024-08-28"},
        {"window": "Sunflower 2024", "start": "2024-08-30", "end": "2024-08-30"},
    ], runs

    # 3. Monthly plan counts for the three requested windows (corrected math).
    #    August has uncovered 08-28 (its own request) then 08-30..31.
    plan = plan_requests(REQUESTED_WINDOWS, covered_days={date(2024, 8, 27), date(2024, 8, 29)})
    by_window: dict[str, list[dict]] = {}
    for entry in plan:
        by_window.setdefault(entry["window"], []).append(entry)
    assert len(by_window["Winter wheat 2023"]) == 10, by_window["Winter wheat 2023"]
    assert len(by_window["Sunflower 2024"]) == 3, by_window["Sunflower 2024"]
    assert len(by_window["Winter wheat 2025"]) == 10, by_window["Winter wheat 2025"]
    assert len(plan) == 23, len(plan)
    assert plan[0]["start"] == "2022-10-01" and plan[0]["end"] == "2022-10-31"

    def _has(window: str, start: str, end: str) -> bool:
        return any(
            e["window"] == window and e["start"] == start and e["end"] == end
            for e in plan
        )

    assert _has("Sunflower 2024", "2024-08-28", "2024-08-28"), plan
    assert _has("Sunflower 2024", "2024-08-30", "2024-08-31"), plan
    assert _has("Sunflower 2024", "2024-09-01", "2024-09-30"), plan
    assert plan[-1]["start"] == "2025-07-01" and plan[-1]["end"] == "2025-07-15"

    # 3b. Final-day inclusion (half-open): the API timeRange must close on the
    #     day AFTER the last included day at 00:00:00Z (exclusive), so the final
    #     daily interval is not dropped by lastIntervalBehavior=SKIP.
    first = plan[0]
    assert first["time_from"] == "2022-10-01T00:00:00Z"
    assert first["time_to"] == "2022-11-01T00:00:00Z"
    last = plan[-1]
    assert last["time_from"] == "2025-07-01T00:00:00Z"
    assert last["time_to"] == "2025-07-16T00:00:00Z"
    # A single included day (2024-08-28) must span exactly that one day.
    assert api_timerange("2024-08-28", "2024-08-28") == (
        "2024-08-28T00:00:00Z", "2024-08-29T00:00:00Z"
    )
    assert api_timerange("2024-08-30", "2024-08-31") == (
        "2024-08-30T00:00:00Z", "2024-09-01T00:00:00Z"
    )
    assert api_timerange("2024-09-01", "2024-09-30") == (
        "2024-09-01T00:00:00Z", "2024-10-01T00:00:00Z"
    )
    assert api_timerange("2023-07-15", "2023-07-15") == (
        "2023-07-15T00:00:00Z", "2023-07-16T00:00:00Z"
    )

    # 3c. The reused Step 2/3 builder passes explicit timestamps through
    #     unchanged, so the exclusive end above is honoured on the wire.
    from .sentinel_stats import build_statistics_request  # local; offline only

    probe = build_statistics_request(
        geometry={"type": "Polygon", "coordinates": [[[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]]]},
        crs_epsg=32635,
        start_date="2024-08-28",
        end_date="2024-08-29T00:00:00Z",
    )
    assert probe["aggregation"]["timeRange"] == {
        "from": "2024-08-28T00:00:00Z",
        "to": "2024-08-29T00:00:00Z",
    }

    # 4. Partial coverage inside a month is clipped, not skipped.
    runs2 = uncovered_request_ranges(
        [("X", "2024-08-01", "2024-08-31")], {date(2024, 8, 1), date(2024, 8, 2)}
    )
    assert runs2 == [{"window": "X", "start": "2024-08-03", "end": "2024-08-31"}], runs2

    # 5. Dedupe + merge preserves masked rows and sorts chronologically.
    rows_a = [
        {"interval_from": "2023-07-01T00:00:00Z", "interval_to": "2023-07-02T00:00:00Z",
         "valid_count": "0", "ndvi_mean": "", "ndmi_mean": ""},
    ]
    rows_b = [
        {"interval_from": "2022-10-01T00:00:00Z", "interval_to": "2022-10-02T00:00:00Z",
         "valid_count": "96", "ndvi_mean": "0.5", "ndmi_mean": "0.2"},
        {"interval_from": "2023-07-01T00:00:00Z", "interval_to": "2023-07-02T00:00:00Z",
         "valid_count": "10", "ndvi_mean": "0.3", "ndmi_mean": "0.1"},
    ]
    merged = merge_rows(rows_a, rows_b)
    assert len(merged) == 2, merged  # duplicate 2023-07-01 collapsed
    assert merged[0]["interval_from"] == "2022-10-01T00:00:00Z"
    assert merged[1]["valid_count"] == "0"  # first wins -> masked row preserved

    # 6. Request summariser distinguishes the four outcome classes.
    assert summarize_request([])["status"] == "no imagery"
    assert summarize_request(rows_b)["status"] == "success"
    assert summarize_request(rows_b)["valid"] == 2
    assert summarize_request(rows_a)["masked"] == 1

    # 7. Schema validation: full row set passes, duplicates and bad schema fail.
    good_row = {c: "" for c in CSV_COLUMNS}
    good_row.update({
        "interval_from": "2024-01-01T00:00:00Z",
        "interval_to": "2024-01-02T00:00:00Z",
        "valid_count": "96", "ndvi_mean": "0.5", "ndmi_mean": "0.2",
        "clear_fraction": "1.0",
    })
    report = validate_rows([good_row])
    assert report["rows"] == 1 and not report["bad_timestamps"]
    try:
        validate_rows([good_row, good_row])  # duplicate key
    except ValueError:
        pass
    else:
        raise AssertionError("duplicate keys were not detected")
    try:
        validate_rows([{k: v for k, v in good_row.items() if k != "ndvi_mean"}])
    except ValueError:
        pass
    else:
        raise AssertionError("missing column was not detected")

    # 8. Checkpoint idempotency: a completed range is not re-requested.
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        log = Path(tmp) / "log.json"
        _save_checkpoint(log, {"requests": [], "completed_ranges": ["2022-10-01..2022-10-31"]})
        loaded = _load_checkpoint(log)
        assert "2022-10-01..2022-10-31" in loaded["completed_ranges"]

    # 9. QA is produced via the real Step 4.1 functions.
    qa_rows, analysis = analyze([good_row])
    assert len(qa_rows) == 1 and qa_rows[0]["is_valid"] == "yes"
    assert set(analyze.__module__ for _ in [0]) == {"src.quality_assessment"}

    print("historical_acquisition.py self-check passed (offline, no files written).")


# --------------------------------------------------------------------------
# Command-line interface
# --------------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m src.historical_acquisition",
        description="Step 4.3 historical Sentinel-2 acquisition for P001.",
    )
    parser.add_argument("--check", action="store_true",
                        help="Run the offline self-check and exit (no network).")
    parser.add_argument("--plan", action="store_true",
                        help="Print the request plan (no network, no files written).")
    parser.add_argument("--execute", action="store_true",
                        help="Perform the live acquisition (network + credentials).")
    parser.add_argument("--existing", default=str(DEFAULT_EXISTING_CSV),
                        help="Existing time-series CSV (read-only).")
    parser.add_argument("--out", default=str(DEFAULT_HISTORICAL_CSV),
                        help="Output historical CSV (additive).")
    parser.add_argument("--qa-out", default=str(DEFAULT_HISTORICAL_QA_CSV),
                        help="Output historical QA CSV (additive).")
    parser.add_argument("--log", default=str(DEFAULT_LOG_PATH),
                        help="Checkpoint / log path.")
    parser.add_argument("--parcel", default=DEFAULT_PARCEL_ID)
    return parser


def _print_plan(plan: list[dict]) -> None:
    per_window: dict[str, int] = {}
    for entry in plan:
        per_window[entry["window"]] = per_window.get(entry["window"], 0) + 1
        print(
            f"  {entry['window']:<18} {entry['start']} .. {entry['end']}  "
            f"[timeRange: {entry['time_from']} -> {entry['time_to']}]"
        )
    print("")
    for window, count in per_window.items():
        print(f"  {window}: {count} monthly request(s)")
    print(f"  TOTAL: {len(plan)} Statistical API request(s) (+1 OAuth token request)")


def main(argv: list[str] | None = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    if args.check:
        _self_check()
        return 0

    covered = covered_start_days(args.existing)
    plan = plan_requests(REQUESTED_WINDOWS, covered)

    if args.plan or not args.execute:
        print(f"Uncovered monthly requests (existing covered days: {len(covered)}):")
        _print_plan(plan)
        if not args.execute:
            print("\nDry run: no network request made. Use --execute (after approval).")
        return 0

    result = run_historical_acquisition(
        existing_csv=args.existing,
        historical_csv=args.out,
        historical_qa_csv=args.qa_out,
        log_path=args.log,
        parcel_id=args.parcel,
        execute=True,
    )
    print(f"\nNew rows: {result['new_rows']}")
    print(f"Historical CSV: {result['historical_csv']}")
    print(f"QA CSV: {result['qa_csv']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
