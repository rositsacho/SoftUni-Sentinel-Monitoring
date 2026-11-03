"""
sentinel_stats.py - per-parcel Sentinel-2 L2A index statistics via the
Copernicus Data Space Ecosystem (CDSE) Sentinel Hub Statistical API.

Endpoint:
    https://sh.dataspace.copernicus.eu/statistics/v1

Workflow (one function call per parcel):
    1. load the parcel geometry (EPSG:4326) from a GeoJSON file,
    2. reproject it to the local UTM zone so that the 20 m resolution is
       expressed in metres (the Statistical API expects resx/resy in the
       units of the request geometry),
    3. build a Statistical API request that runs an Evalscript V3 which
       computes NDVI and NDMI per pixel (units: REFLECTANCE, mosaicking:
       SIMPLE) and exports an explicit dataMask,
    4. POST the request with an OAuth2 Bearer token (src.sentinelhub_auth),
    5. parse the JSON response into one row per daily (P1D) interval with
       mean / standard deviation for NDVI and NDMI, sampleCount, noDataCount
       and the valid-pixel coverage,
    6. write the rows to a CSV under data/results/.

The module never prints or stores credentials or access tokens.

Run the offline self-test (no network):
    python -m src.sentinel_stats

Run the real P001 test for June 2024 (network + credentials):
    python -m src.sentinel_stats --execute

Run full-season monitoring for P001 (network + credentials), one request
per calendar month, combined chronologically into a single CSV:
    python -m src.sentinel_stats --monthly \
        --start 2023-09-01 --end 2024-09-01 \
        --output data/results/time_series_p001_2023_2024.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from datetime import date, timedelta
from pathlib import Path

import requests
from pyproj import Transformer

try:  # works for "python -m src.sentinel_stats"
    from .sentinelhub_auth import SentinelHubAuth
except ImportError:  # fallback for "python src/sentinel_stats.py"
    from sentinelhub_auth import SentinelHubAuth  # type: ignore

# --------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------

STATISTICS_URL = "https://sh.dataspace.copernicus.eu/statistics/v1"
S2_L2A_COLLECTION = "sentinel-2-l2a"

# EPSG:4326 is what GeoJSON always stores.
WGS84_EPSG = 4326

# Mosaicking order used inside dataFilter (which source tile wins when tiles
# overlap). "leastCC" prefers the least cloudy tile for the day.
DATA_FILTER_MOSAICKING_ORDER = "leastCC"

# SCL (Scene Classification Layer) classes considered valid agricultural
# surfaces: 4 = vegetation, 5 = not-vegetated / bare soil. Everything else
# (water, cloud, shadow, snow, no-data, ...) is masked out. See the module
# docstring and README for the rationale.
VALID_SCL_CLASSES: tuple[int, ...] = (4, 5)

# Default test parameters (Step 2 first test).
DEFAULT_PARCEL_ID = "P001"
DEFAULT_START = "2024-06-01"
DEFAULT_END = "2024-06-30"
DEFAULT_RESOLUTION_M = 20.0

# Defaults for Step 3 full-season monitoring (P001). The end date is
# exclusive, so the last processed day is 2024-08-31.
DEFAULT_SEASON_START = "2023-09-01"
DEFAULT_SEASON_END = "2024-09-01"

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_GEOJSON_PATH = PROJECT_ROOT / "data" / "results" / "all_parcels.geojson"
DEFAULT_RESULTS_DIR = PROJECT_ROOT / "data" / "results"

# Columns written to the output CSV, in order.
#
# Two coverage measures are reported on purpose, because they answer
# different questions (they coincide only when no pixel inside the polygon
# is masked):
#   * geometry_coverage = geometry_pixel_count / sample_count
#       -> the polygon's share of the request bounding box (pure geometry,
#          independent of clouds / data availability),
#   * valid_coverage    = valid_count / sample_count
#       -> the raw valid-pixel coverage over the bounding box (this is the
#          "valid-pixel coverage" required by Step 2),
#   * clear_fraction    = valid_count / geometry_pixel_count
#       -> cloud-free, SCL-valid fraction INSIDE the polygon.
# See README "Methodology notes" for the boundary-pixel limitation.
CSV_COLUMNS = [
    "interval_from",
    "interval_to",
    "ndvi_mean",
    "ndvi_stdev",
    "ndmi_mean",
    "ndmi_stdev",
    "sample_count",
    "no_data_count",
    "valid_count",
    "geometry_pixel_count",
    "geometry_coverage",
    "valid_coverage",
    "clear_fraction",
]


# --------------------------------------------------------------------------
# Errors
# --------------------------------------------------------------------------

class SentinelStatsError(RuntimeError):
    """Base class for errors raised by this module."""


class SentinelApiError(SentinelStatsError):
    """Raised when the Statistical API returns an error response."""


class EmptyObservationsError(SentinelStatsError):
    """Raised when the API response contains no usable observations."""


class ParcelNotFoundError(SentinelStatsError):
    """Raised when the requested parcel id is not present in the GeoJSON."""


# --------------------------------------------------------------------------
# Evalscript V3
# --------------------------------------------------------------------------

def build_evalscript(valid_scl_classes: tuple[int, ...] = VALID_SCL_CLASSES) -> str:
    """Return the Evalscript V3 string used for the NDVI/NDMI statistics.

    - REFLECTANCE units for the reflectance bands (B04/B08/B11). SCL is a
      classification layer and only supports DN, so units are given per band
      (SCL and dataMask stay in DN),
    - mosaicking SIMPLE (one mosaicked scene per aggregation interval),
    - explicit 'dataMask' output combining the source dataMask with the
      selected SCL classes.
    """
    scl_checks = " || ".join(f"sample.SCL == {cls}" for cls in valid_scl_classes)
    return f"""//VERSION=3
function setup() {{
  return {{
    input: [{{
      bands: ["B04", "B08", "B11", "SCL", "dataMask"],
      units: ["REFLECTANCE", "REFLECTANCE", "REFLECTANCE", "DN", "DN"]
    }}],
    mosaicking: "SIMPLE",
    output: [
      {{ id: "ndvi", bands: 1, sampleType: "FLOAT32" }},
      {{ id: "ndmi", bands: 1, sampleType: "FLOAT32" }},
      {{ id: "dataMask", bands: 1 }}
    ]
  }};
}}

function evaluatePixel(sample) {{
  // Valid only where the source dataMask is set AND SCL is an agricultural
  // surface (default: vegetation 4 / bare soil 5).
  var sclValid = ({scl_checks});
  var valid = (sample.dataMask == 1 && sclValid) ? 1 : 0;

  var ndvi = 0;
  var ndmi = 0;
  if (valid == 1) {{
    var denomNdvi = sample.B08 + sample.B04;
    var denomNdmi = sample.B08 + sample.B11;
    if (denomNdvi != 0) {{ ndvi = (sample.B08 - sample.B04) / denomNdvi; }}
    if (denomNdmi != 0) {{ ndmi = (sample.B08 - sample.B11) / denomNdmi; }}
  }}

  return {{
    ndvi: [ndvi],
    ndmi: [ndmi],
    dataMask: [valid]
  }};
}}
"""


# --------------------------------------------------------------------------
# Geometry helpers (reprojection without GeoPandas)
# --------------------------------------------------------------------------

def _iter_positions(coords):
    """Yield every [lon, lat] position from a nested GeoJSON coordinate list."""
    first = coords[0]
    if isinstance(first, (int, float)):
        yield coords
        return
    for item in coords:
        yield from _iter_positions(item)


def utm_epsg_for_lon_lat(lon: float, lat: float) -> int:
    """Return the EPSG code of the UTM zone covering a lon/lat point."""
    zone = int((lon + 180.0) // 6.0) + 1
    hemisphere = 32600 if lat >= 0 else 32700
    return hemisphere + zone


def geometry_centroid_lon_lat(geometry: dict) -> tuple[float, float]:
    """Approximate (lon, lat) centroid of a GeoJSON geometry (bounds mid)."""
    positions = list(_iter_positions(geometry["coordinates"]))
    lons = [pos[0] for pos in positions]
    lats = [pos[1] for pos in positions]
    return (min(lons) + max(lons)) / 2.0, (min(lats) + max(lats)) / 2.0


def _transform_coords(coords, transform):
    first = coords[0]
    if isinstance(first, (int, float)):
        x, y = transform(coords[0], coords[1])
        return [x, y]
    return [_transform_coords(item, transform) for item in coords]


def project_geometry(geometry: dict, src_epsg: int, dst_epsg: int) -> dict:
    """Return a copy of a GeoJSON geometry reprojected to dst_epsg."""
    transformer = Transformer.from_crs(src_epsg, dst_epsg, always_xy=True)
    return {
        "type": geometry["type"],
        "coordinates": _transform_coords(geometry["coordinates"], transformer.transform),
    }


# --------------------------------------------------------------------------
# Request builder
# --------------------------------------------------------------------------

def _to_utc_iso(date_str: str, end_of_day: bool = False) -> str:
    """Turn 'YYYY-MM-DD' into an ISO-8601 UTC timestamp."""
    date_str = date_str.strip()
    if "T" in date_str:
        return date_str if date_str.endswith("Z") else f"{date_str}Z"
    suffix = "T23:59:59Z" if end_of_day else "T00:00:00Z"
    return f"{date_str}{suffix}"


def build_statistics_request(
    geometry: dict,
    crs_epsg: int,
    start_date: str,
    end_date: str,
    resolution_m: float = DEFAULT_RESOLUTION_M,
    collection: str = S2_L2A_COLLECTION,
    evalscript: str | None = None,
    mosaicking_order: str = DATA_FILTER_MOSAICKING_ORDER,
) -> dict:
    """Build the JSON body for a Statistical API POST request."""
    if evalscript is None:
        evalscript = build_evalscript()

    return {
        "input": {
            "bounds": {
                "geometry": geometry,
                "properties": {
                    "crs": f"http://www.opengis.net/def/crs/EPSG/0/{crs_epsg}"
                },
            },
            "data": [
                {
                    "type": collection,
                    "dataFilter": {"mosaickingOrder": mosaicking_order},
                }
            ],
        },
        "aggregation": {
            "timeRange": {
                "from": _to_utc_iso(start_date, end_of_day=False),
                "to": _to_utc_iso(end_date, end_of_day=True),
            },
            "aggregationInterval": {"of": "P1D"},
            "evalscript": evalscript,
            "resx": resolution_m,
            "resy": resolution_m,
        },
        "calculations": {"default": {}},
    }


# --------------------------------------------------------------------------
# Response parsing
# --------------------------------------------------------------------------

def _extract_band_stats(output: dict, band_key: str = "B0") -> dict | None:
    """Return the stats dict for one output/band, or None when absent."""
    if not isinstance(output, dict):
        return None
    bands = output.get("bands")
    if not isinstance(bands, dict):
        return None
    band = bands.get(band_key)
    if not isinstance(band, dict):
        return None
    stats = band.get("stats")
    return stats if isinstance(stats, dict) else None


def _safe_float(value) -> float | None:
    """Parse a numeric stat. NaN (all pixels masked) is treated as missing."""
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if math.isnan(number):
        return None
    return number


def _ratio(numerator, denominator) -> float | None:
    """Return numerator/denominator, or None when it cannot be computed."""
    if numerator is None or denominator is None:
        return None
    if denominator <= 0:
        return None
    return numerator / denominator


def parse_statistics_response(payload: dict) -> list[dict]:
    """Convert a Statistical API response into one dict per daily interval.

    Response shape (confirmed against the CDSE examples and the live API):
        {"status": "OK",
         "geometryPixelCount": 96,          # top level, geometry only
         "data": [{"interval": {"from": ..., "to": ...},
                   "outputs": {"ndvi": {"bands": {"B0": {"stats": {...}}}},
                               "ndmi": {...}}}]}

    `geometryPixelCount` (pixels intersecting the request geometry) is
    reported once at the top level and is copied onto every row so that
    geometry coverage and cloud-free coverage stay separable.
    """
    if not isinstance(payload, dict):
        raise SentinelApiError("Unexpected (non-JSON-object) API response.")

    status = payload.get("status")
    if status not in (None, "OK"):
        raise SentinelApiError(f"Statistical API returned status '{status}'.")

    entries = payload.get("data")
    if not isinstance(entries, list):
        raise SentinelApiError("API response is missing the 'data' list.")

    geometry_pixel_count = payload.get("geometryPixelCount")
    if geometry_pixel_count is not None:
        try:
            geometry_pixel_count = int(geometry_pixel_count)
        except (TypeError, ValueError):
            geometry_pixel_count = None

    rows: list[dict] = []
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        interval = entry.get("interval") or {}
        outputs = entry.get("outputs") or {}

        ndvi_stats = _extract_band_stats(outputs.get("ndvi"))
        ndmi_stats = _extract_band_stats(outputs.get("ndmi"))

        # Counts are identical for every output (same mask). Prefer NDVI,
        # fall back to NDMI if NDVI is missing for this interval.
        counts_source = ndvi_stats or ndmi_stats or {}
        sample_count = counts_source.get("sampleCount")
        no_data_count = counts_source.get("noDataCount")

        sample_count = int(sample_count) if sample_count is not None else None
        no_data_count = int(no_data_count) if no_data_count is not None else None

        valid_count = None
        if sample_count is not None and no_data_count is not None:
            valid_count = max(sample_count - no_data_count, 0)

        rows.append(
            {
                "interval_from": interval.get("from"),
                "interval_to": interval.get("to"),
                "ndvi_mean": _safe_float((ndvi_stats or {}).get("mean")),
                "ndvi_stdev": _safe_float((ndvi_stats or {}).get("stDev")),
                "ndmi_mean": _safe_float((ndmi_stats or {}).get("mean")),
                "ndmi_stdev": _safe_float((ndmi_stats or {}).get("stDev")),
                "sample_count": sample_count,
                "no_data_count": no_data_count,
                "valid_count": valid_count,
                "geometry_pixel_count": geometry_pixel_count,
                # Geometry vs bounding box (pure geometry, no clouds).
                "geometry_coverage": _ratio(geometry_pixel_count, sample_count),
                # Raw valid-pixel coverage over the bounding box.
                "valid_coverage": _ratio(valid_count, sample_count),
                # Cloud-free, SCL-valid fraction inside the polygon.
                "clear_fraction": _ratio(valid_count, geometry_pixel_count),
            }
        )
    return rows


# --------------------------------------------------------------------------
# API call
# --------------------------------------------------------------------------

def _safe_api_error_text(response: requests.Response) -> str:
    try:
        payload = response.json()
    except ValueError:
        return (response.text or "").strip()[:300]
    if isinstance(payload, dict):
        for key in ("message", "reason", "error_description", "error"):
            if payload.get(key):
                return str(payload[key])[:300]
        return str(payload)[:300]
    return str(payload)[:300]


def fetch_statistics(
    request_payload: dict,
    auth: SentinelHubAuth,
    timeout: float = 60.0,
) -> dict:
    """POST a Statistical API request and return the parsed JSON response."""
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json",
        **auth.auth_headers(),
    }
    try:
        response = requests.post(
            STATISTICS_URL,
            headers=headers,
            json=request_payload,
            timeout=timeout,
        )
    except requests.RequestException as exc:
        raise SentinelApiError(f"Could not reach the Statistical API: {exc}") from exc

    if response.status_code != 200:
        raise SentinelApiError(
            f"Statistical API error (HTTP {response.status_code}): "
            f"{_safe_api_error_text(response)}"
        )
    try:
        return response.json()
    except ValueError as exc:
        raise SentinelApiError("Statistical API returned a non-JSON response.") from exc


# --------------------------------------------------------------------------
# Seasonal monitoring helpers (Step 3)
# --------------------------------------------------------------------------

def _parse_iso_date(value: str) -> date:
    """Parse a 'YYYY-MM-DD' string, raising a clear error when invalid."""
    try:
        return date.fromisoformat(value.strip())
    except (AttributeError, ValueError) as exc:
        raise ValueError(
            f"Invalid ISO date '{value}' (expected YYYY-MM-DD)."
        ) from exc


def month_ranges(start_date: str, end_date: str) -> list[tuple[str, str]]:
    """Split the half-open range [start_date, end_date) into month chunks.

    Returns a list of (chunk_start, chunk_end) day strings where chunk_end is
    the last day *included* in the chunk (the overall end_date is exclusive).
    The chunks are contiguous and non-overlapping.
    """
    start = _parse_iso_date(start_date)
    end = _parse_iso_date(end_date)
    if end <= start:
        raise ValueError(
            f"end date ({end_date}) must be after start date ({start_date})."
        )

    chunks: list[tuple[str, str]] = []
    cursor = start
    while cursor < end:
        next_month = (
            date(cursor.year + 1, 1, 1)
            if cursor.month == 12
            else date(cursor.year, cursor.month + 1, 1)
        )
        chunk_end_exclusive = min(next_month, end)
        chunks.append(
            (
                cursor.isoformat(),
                (chunk_end_exclusive - timedelta(days=1)).isoformat(),
            )
        )
        cursor = next_month
    return chunks


def combine_interval_rows(rows: list[dict]) -> list[dict]:
    """Sort rows chronologically and drop duplicate intervals.

    Rows are keyed by (interval_from, interval_to); the first occurrence of a
    duplicate interval wins. Missing and fully masked rows (valid_count == 0)
    are preserved unchanged.
    """
    ordered = sorted(
        rows,
        key=lambda row: (str(row.get("interval_from")), str(row.get("interval_to"))),
    )
    combined: list[dict] = []
    seen: set[tuple] = set()
    for row in ordered:
        key = (row.get("interval_from"), row.get("interval_to"))
        if key in seen:
            continue
        seen.add(key)
        combined.append(row)
    return combined


# --------------------------------------------------------------------------
# High-level workflow
# --------------------------------------------------------------------------

def load_parcel_geometry(
    geojson_path: str | Path,
    parcel_id: str,
) -> tuple[dict, int]:
    """Return (geometry, source_epsg) for one parcel id from a GeoJSON file."""
    geojson_path = Path(geojson_path)
    if not geojson_path.exists():
        raise FileNotFoundError(f"GeoJSON not found: {geojson_path}")

    with geojson_path.open("r", encoding="utf-8") as handle:
        collection = json.load(handle)

    for feature in collection.get("features", []):
        properties = feature.get("properties", {}) or {}
        if str(properties.get("parcel_id")) == str(parcel_id):
            geometry = feature.get("geometry")
            if not geometry:
                raise ParcelNotFoundError(
                    f"Parcel {parcel_id} has no geometry in {geojson_path}."
                )
            return geometry, WGS84_EPSG

    raise ParcelNotFoundError(
        f"Parcel {parcel_id} not found in {geojson_path}."
    )


def get_parcel_time_series(
    geometry: dict,
    src_epsg: int,
    start_date: str,
    end_date: str,
    resolution_m: float = DEFAULT_RESOLUTION_M,
    auth: SentinelHubAuth | None = None,
    env_path: str | Path | None = None,
    require_usable: bool = True,
) -> tuple[list[dict], dict]:
    """Fetch and parse the daily index statistics for one parcel.

    Returns (rows, metadata). By default raises EmptyObservationsError when no
    interval contains a usable observation. With require_usable=False the rows
    are returned as-is (possibly empty, or containing only fully masked rows),
    which lets the seasonal runner preserve masked observations instead of
    aborting on an overcast month.
    """
    lon, lat = geometry_centroid_lon_lat(geometry)
    dst_epsg = utm_epsg_for_lon_lat(lon, lat)
    projected = project_geometry(geometry, src_epsg, dst_epsg)

    payload = build_statistics_request(
        geometry=projected,
        crs_epsg=dst_epsg,
        start_date=start_date,
        end_date=end_date,
        resolution_m=resolution_m,
    )

    auth = auth or SentinelHubAuth(env_path=env_path)
    raw = fetch_statistics(payload, auth)
    rows = parse_statistics_response(raw)

    usable = [row for row in rows if (row.get("valid_count") or 0) > 0]
    if require_usable and not usable:
        raise EmptyObservationsError(
            "No valid Sentinel-2 observations for the requested area/period "
            f"({start_date}..{end_date}). Every pixel was masked as no-data or "
            "as a non-agricultural SCL class."
        )

    metadata = {
        "crs_epsg": dst_epsg,
        "resolution_m": resolution_m,
        "intervals_returned": len(rows),
        "intervals_with_data": len(usable),
    }
    return rows, metadata


def _csv_cell(value):
    """Render a cell: missing values become an empty CSV field (not 'nan')."""
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return ""
    return value


def save_time_series_csv(rows: list[dict], output_path: str | Path) -> Path:
    """Write parsed rows to CSV with the canonical column order.

    Missing statistics (e.g. NDVI/NDMI of a fully masked interval) are
    written as empty cells rather than the literal string 'nan'.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {column: _csv_cell(row.get(column)) for column in CSV_COLUMNS}
            )
    return output_path


def run_first_test(
    parcel_id: str = DEFAULT_PARCEL_ID,
    geojson_path: str | Path = DEFAULT_GEOJSON_PATH,
    start_date: str = DEFAULT_START,
    end_date: str = DEFAULT_END,
    output_path: str | Path | None = None,
    resolution_m: float = DEFAULT_RESOLUTION_M,
    env_path: str | Path | None = None,
) -> Path:
    """End-to-end first test: P001, June 2024, 20 m, daily intervals."""
    if output_path is None:
        output_path = DEFAULT_RESULTS_DIR / f"time_series_{parcel_id.lower()}_202406.csv"

    geometry, src_epsg = load_parcel_geometry(geojson_path, parcel_id)
    rows, metadata = get_parcel_time_series(
        geometry=geometry,
        src_epsg=src_epsg,
        start_date=start_date,
        end_date=end_date,
        resolution_m=resolution_m,
        env_path=env_path,
    )
    saved = save_time_series_csv(rows, output_path)

    print(
        f"Parcel {parcel_id}: {metadata['intervals_with_data']} interval(s) with "
        f"data out of {metadata['intervals_returned']} "
        f"(UTM EPSG:{metadata['crs_epsg']}, {metadata['resolution_m']:.0f} m)."
    )
    print(f"Wrote {len(rows)} row(s) to {saved}")
    return saved


def run_seasonal_monitoring(
    parcel_id: str = DEFAULT_PARCEL_ID,
    start_date: str = DEFAULT_SEASON_START,
    end_date: str = DEFAULT_SEASON_END,
    output_path: str | Path | None = None,
    resolution_m: float = DEFAULT_RESOLUTION_M,
    geojson_path: str | Path = DEFAULT_GEOJSON_PATH,
    env_path: str | Path | None = None,
) -> Path:
    """Full-season monitoring: one request per calendar month, one CSV.

    The half-open period [start_date, end_date) is split into monthly chunks
    (Step 3, requirement 5). Daily P1D intervals are preserved (no monthly
    averaging), combined chronologically and de-duplicated. Months without any
    scene, and fully masked observations, are preserved rather than dropped.
    A single SentinelHubAuth instance is reused so the OAuth2 token is fetched
    once and reused for every monthly request.
    """
    private = parcel_id.lower()
    if output_path is None:
        first_year = _parse_iso_date(start_date).year
        last_year = (_parse_iso_date(end_date) - timedelta(days=1)).year
        output_path = (
            DEFAULT_RESULTS_DIR / f"time_series_{private}_{first_year}_{last_year}.csv"
        )

    chunks = month_ranges(start_date, end_date)
    geometry, src_epsg = load_parcel_geometry(geojson_path, parcel_id)
    auth = SentinelHubAuth(env_path=env_path)

    all_rows: list[dict] = []
    print(
        f"Parcel {parcel_id}: monitoring {start_date}..{end_date} "
        f"(exclusive) in {len(chunks)} monthly request(s)."
    )
    for chunk_start, chunk_end in chunks:
        rows, metadata = get_parcel_time_series(
            geometry=geometry,
            src_epsg=src_epsg,
            start_date=chunk_start,
            end_date=chunk_end,
            resolution_m=resolution_m,
            auth=auth,
            require_usable=False,
        )
        all_rows.extend(rows)
        print(
            f"  {chunk_start}..{chunk_end}: {metadata['intervals_with_data']} with "
            f"data / {metadata['intervals_returned']} interval(s)"
        )

    combined = combine_interval_rows(all_rows)
    duplicates_removed = len(all_rows) - len(combined)
    usable_total = sum(1 for row in combined if (row.get("valid_count") or 0) > 0)
    masked_total = sum(
        1
        for row in combined
        if row.get("valid_count") is not None and (row.get("valid_count") or 0) == 0
    )

    saved = save_time_series_csv(combined, output_path)
    print(
        f"Parcel {parcel_id}: {usable_total} observation(s) with data, "
        f"{masked_total} fully masked, {duplicates_removed} duplicate(s) removed, "
        f"{len(combined)} row(s) total ({len(chunks)} monthly request(s))."
    )
    print(f"Wrote {saved}")
    return saved


# --------------------------------------------------------------------------
# Offline self-test (no network): python -m src.sentinel_stats
# --------------------------------------------------------------------------

def _self_test() -> None:
    evalscript = build_evalscript()
    for fragment in (
        "//VERSION=3",
        'mosaicking: "SIMPLE"',
        'units: ["REFLECTANCE", "REFLECTANCE", "REFLECTANCE", "DN", "DN"]',
        'id: "dataMask"',
        "sample.SCL == 4 || sample.SCL == 5",
        "(sample.B08 - sample.B04) / denomNdvi",
        "(sample.B08 - sample.B11) / denomNdmi",
    ):
        assert fragment in evalscript, f"evalscript missing: {fragment}"

    # A small synthetic square near Stara Zagora (EPSG:4326).
    geometry = {
        "type": "Polygon",
        "coordinates": [
            [
                [25.30, 42.40],
                [25.31, 42.40],
                [25.31, 42.41],
                [25.30, 42.41],
                [25.30, 42.40],
            ]
        ],
    }
    lon, lat = geometry_centroid_lon_lat(geometry)
    dst_epsg = utm_epsg_for_lon_lat(lon, lat)
    assert dst_epsg == 32635, f"unexpected UTM zone: {dst_epsg}"

    projected = project_geometry(geometry, WGS84_EPSG, dst_epsg)
    xs = [pos[0] for pos in _iter_positions(projected["coordinates"])]
    assert min(xs) > 100_000, "coordinates were not reprojected to metres"

    payload = build_statistics_request(projected, dst_epsg, DEFAULT_START, DEFAULT_END)
    assert payload["aggregation"]["resx"] == DEFAULT_RESOLUTION_M
    assert payload["aggregation"]["aggregationInterval"]["of"] == "P1D"
    assert payload["aggregation"]["timeRange"]["from"] == "2024-06-01T00:00:00Z"
    assert payload["aggregation"]["timeRange"]["to"] == "2024-06-30T23:59:59Z"
    assert payload["input"]["data"][0]["type"] == S2_L2A_COLLECTION

    # Synthetic API response resembling the confirmed CDSE structure,
    # including the top-level geometryPixelCount.
    synthetic = {
        "status": "OK",
        "geometryPixelCount": 80,
        "data": [
            {
                "interval": {"from": "2024-06-05T00:00:00Z", "to": "2024-06-06T00:00:00Z"},
                "outputs": {
                    "ndvi": {"bands": {"B0": {"stats": {
                        "min": 0.1, "max": 0.9, "mean": 0.55, "stDev": 0.12,
                        "sampleCount": 100, "noDataCount": 20}}}},
                    "ndmi": {"bands": {"B0": {"stats": {
                        "min": -0.1, "max": 0.5, "mean": 0.22, "stDev": 0.07,
                        "sampleCount": 100, "noDataCount": 20}}}},
                },
            },
            {
                # A fully masked interval (e.g. complete cloud cover).
                "interval": {"from": "2024-06-12T00:00:00Z", "to": "2024-06-13T00:00:00Z"},
                "outputs": {
                    "ndvi": {"bands": {"B0": {"stats": {
                        "mean": float("nan"), "stDev": float("nan"),
                        "sampleCount": 100, "noDataCount": 100}}}},
                },
            },
        ],
    }
    rows = parse_statistics_response(synthetic)
    assert len(rows) == 2
    assert rows[0]["ndvi_mean"] == 0.55
    assert rows[0]["valid_count"] == 80
    assert rows[0]["geometry_pixel_count"] == 80
    # Geometry coverage and cloud-free coverage are reported separately.
    assert abs(rows[0]["geometry_coverage"] - 0.8) < 1e-9
    assert abs(rows[0]["valid_coverage"] - 0.8) < 1e-9
    assert abs(rows[0]["clear_fraction"] - 1.0) < 1e-9
    # Fully masked interval: zero valid pixels, no usable statistics.
    assert rows[1]["sample_count"] == 100
    assert rows[1]["no_data_count"] == 100
    assert rows[1]["valid_count"] == 0
    assert rows[1]["geometry_pixel_count"] == 80
    assert rows[1]["valid_coverage"] == 0.0
    assert rows[1]["clear_fraction"] == 0.0
    assert rows[1]["ndvi_mean"] is None, "NaN must be normalised to None"

    # Empty / error responses must raise.
    try:
        parse_statistics_response({"status": "FAILED", "data": []})
    except SentinelApiError:
        pass
    else:
        raise AssertionError("expected SentinelApiError for a failed status")

    # CSV writer round-trip: header, empty cells instead of 'nan'.
    scratch = PROJECT_ROOT / "data" / "results" / "_selftest_stats.csv"
    save_time_series_csv(rows, scratch)
    content = scratch.read_text(encoding="utf-8")
    scratch.unlink()
    header = content.splitlines()[0].split(",")
    assert header == CSV_COLUMNS, f"unexpected CSV header: {header}"
    assert "nan" not in content.lower(), "'nan' leaked into the CSV"
    masked_line = [
        line for line in content.splitlines() if line.startswith("2024-06-12")
    ]
    assert masked_line, "masked interval missing from CSV"
    masked_cells = masked_line[0].split(",")
    # NDVI/NDMI mean+stDev (columns 3-6) must be empty for the masked day,
    # while the counts and coverage columns are still populated.
    assert masked_cells[2:6] == ["", "", "", ""], (
        f"masked row not blanked: {masked_line[0]}"
    )

    # --- Step 3: monthly chunking -----------------------------------------
    chunks = month_ranges("2023-09-01", "2024-09-01")
    assert len(chunks) == 12, f"expected 12 monthly chunks, got {len(chunks)}"
    assert chunks[0] == ("2023-09-01", "2023-09-30"), chunks[0]
    assert chunks[-1] == ("2024-08-01", "2024-08-31"), chunks[-1]
    # Contiguous, non-overlapping, and covering the whole half-open range.
    for (_, prev_end), (next_start, _) in zip(chunks, chunks[1:]):
        expected_next = (
            _parse_iso_date(prev_end) + timedelta(days=1)
        ).isoformat()
        assert next_start == expected_next, (prev_end, next_start)
    covered = sum(
        (_parse_iso_date(end) - _parse_iso_date(start)).days + 1
        for start, end in chunks
    )
    assert covered == (_parse_iso_date("2024-09-01") - _parse_iso_date("2023-09-01")).days
    # Leap year: 2024-02 must span 29 days.
    assert ("2024-02-01", "2024-02-29") in chunks

    # A single-day range still yields exactly one chunk.
    assert month_ranges("2024-01-01", "2024-01-02") == [("2024-01-01", "2024-01-01")]
    for bad in (("2024-01-01", "2024-01-01"), ("2024-02-01", "2024-01-01")):
        try:
            month_ranges(*bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f"expected ValueError for range {bad}")

    # --- Step 3: chronological combine + de-duplication -------------------
    duplicate = dict(rows[0])  # same interval as a previous month (overlap)
    out_of_order = [rows[1], rows[0], duplicate]
    combined = combine_interval_rows(out_of_order)
    assert len(combined) == 2, "duplicate interval was not removed"
    assert [r["interval_from"] for r in combined] == [
        "2024-06-05T00:00:00Z",
        "2024-06-12T00:00:00Z",
    ], "rows were not sorted chronologically"
    # The fully masked row survives de-duplication and stays masked.
    assert combined[1]["valid_count"] == 0
    assert combined[1]["ndvi_mean"] is None

    # --- Step 3: CLI parser ----------------------------------------------
    parser = build_arg_parser()
    args = parser.parse_args(
        [
            "--monthly",
            "--start", "2023-09-01",
            "--end", "2024-09-01",
            "--output", "data/results/time_series_p001_2023_2024.csv",
        ]
    )
    assert args.monthly is True
    assert args.start == "2023-09-01"
    assert args.end == "2024-09-01"
    assert args.parcel == DEFAULT_PARCEL_ID
    assert args.resolution == DEFAULT_RESOLUTION_M
    assert args.execute is False
    # `--execute` remains the backward-compatible single-window mode.
    args_exec = parser.parse_args(["--execute"])
    assert args_exec.execute is True
    assert args_exec.monthly is False

    print("Statistical API module self-test passed (offline, no network).")


# --------------------------------------------------------------------------
# Command-line interface
# --------------------------------------------------------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    """Build the command-line parser (dates configurable, req 3)."""
    parser = argparse.ArgumentParser(
        prog="python -m src.sentinel_stats",
        description=(
            "Sentinel-2 L2A NDVI/NDMI statistics for CAP parcels via the "
            "CDSE Sentinel Hub Statistical API."
        ),
    )
    parser.add_argument(
        "--execute",
        action="store_true",
        help="Backward-compatible single-window test (defaults P001, June 2024).",
    )
    parser.add_argument(
        "--monthly",
        action="store_true",
        help="Full-season monitoring: one request per calendar month, one CSV.",
    )
    parser.add_argument(
        "--parcel",
        default=DEFAULT_PARCEL_ID,
        help=f"Parcel id to process (default: {DEFAULT_PARCEL_ID}).",
    )
    parser.add_argument(
        "--start",
        default=None,
        help="Start date YYYY-MM-DD, inclusive "
        f"(monthly default: {DEFAULT_SEASON_START}; execute default: {DEFAULT_START}).",
    )
    parser.add_argument(
        "--end",
        default=None,
        help="End date YYYY-MM-DD, exclusive "
        f"(monthly default: {DEFAULT_SEASON_END}; execute default: {DEFAULT_END}).",
    )
    parser.add_argument(
        "--resolution",
        type=float,
        default=DEFAULT_RESOLUTION_M,
        help=f"Pixel resolution in metres (default: {DEFAULT_RESOLUTION_M:.0f}).",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output CSV path (defaults per mode).",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Entry point. No args runs the offline self-test (no network)."""
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    if args.monthly and args.execute:
        parser.error("--monthly and --execute are mutually exclusive.")

    if args.monthly:
        run_seasonal_monitoring(
            parcel_id=args.parcel,
            start_date=args.start or DEFAULT_SEASON_START,
            end_date=args.end or DEFAULT_SEASON_END,
            output_path=args.output,
            resolution_m=args.resolution,
        )
        return 0

    if args.execute:
        run_first_test(
            parcel_id=args.parcel,
            start_date=args.start or DEFAULT_START,
            end_date=args.end or DEFAULT_END,
            output_path=args.output,
            resolution_m=args.resolution,
        )
        return 0

    _self_test()
    return 0


if __name__ == "__main__":
    sys.exit(main())
