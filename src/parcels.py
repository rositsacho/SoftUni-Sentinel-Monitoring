"""
parcels.py - read agricultural parcel boundaries from KML files.

Pipeline:
    KML file(s) -> GeoDataFrame in EPSG:4326 (ALL layers/folders merged)
                -> cleaned geometries (every skipped feature is reported)
                -> unique parcel IDs (across files when using the directory API)
                -> area in hectares measured in EPSG:3035 (equal-area)
                -> export to GeoJSON and CSV

Run the synthetic self-test with:
    python -m src.parcels
"""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import pandas as pd
from shapely import make_valid
from shapely.geometry import GeometryCollection, LineString, Polygon, box
from shapely.ops import unary_union

# pyproj.Geod is only used by the self-test at the bottom of this file.
from pyproj import Geod

# KML always stores coordinates as WGS84 lat/lon (EPSG:4326) by specification.
KML_CRS = "EPSG:4326"

# Equal-area projection for Europe, used ONLY for measuring areas.
# Bulgaria lies fully inside its domain; the EU itself uses EPSG:3035
# for agricultural statistics, so our hectares are comparable with CAP data.
AREA_CRS = "EPSG:3035"


# --------------------------------------------------------------------------
# Geometry helpers
# --------------------------------------------------------------------------

def _polygon_parts(geom) -> list[Polygon]:
    """Recursively collect every non-empty polygon inside a geometry.

    make_valid() can turn one broken polygon into a GeometryCollection that
    mixes polygons with lines/points. We keep every polygon component and
    ignore only the zero-area debris (points and lines).
    """
    if geom is None or geom.is_empty:
        return []
    if geom.geom_type == "Polygon":
        return [geom]
    if geom.geom_type == "MultiPolygon":
        return [part for part in geom.geoms if not part.is_empty]
    if geom.geom_type == "GeometryCollection":
        parts: list[Polygon] = []
        for part in geom.geoms:
            parts.extend(_polygon_parts(part))
        return parts
    return []  # Point / LineString debris carries no area


def _clean_geometry_with_reason(geom):
    """Return (cleaned_geometry_or_None, skip_reason_or_None)."""
    if geom is None:
        return None, "missing geometry"
    if geom.is_empty:
        return None, "empty geometry"
    if not geom.is_valid:
        geom = make_valid(geom)  # repairs self-intersections etc. (Shapely 2.x API)
    parts = _polygon_parts(geom)
    if not parts:
        return None, f"no polygon inside geometry (type was {geom.geom_type})"
    merged = unary_union(parts)  # combine into one valid Polygon/MultiPolygon
    if merged.is_empty:
        return None, "geometry became empty after repair"
    return merged, None


def clean_geometry(geom):
    """Public wrapper: returns a valid Polygon/MultiPolygon, or None."""
    cleaned, _reason = _clean_geometry_with_reason(geom)
    return cleaned


# --------------------------------------------------------------------------
# Reading KML
# --------------------------------------------------------------------------

def _read_all_kml_layers(kml_path: Path) -> gpd.GeoDataFrame:
    """Read EVERY layer of a KML file and merge them into one table.

    GDAL may expose one KML as a single flattened layer, or as one layer
    per (nested) Folder. We enumerate layers first so no Placemark is lost,
    then fall back to a plain single-layer read if enumeration fails.
    """
    frames: list[gpd.GeoDataFrame] = []

    try:
        layers = gpd.list_layers(kml_path)
        for _, layer_row in layers.iterrows():
            if pd.isna(layer_row["geometry_type"]):
                continue  # attribute-only layer, nothing to map
            layer_name = layer_row["name"]
            try:
                frame = gpd.read_file(kml_path, layer=layer_name)
            except Exception as exc:
                print(
                    f"WARNING: {kml_path.name}: layer '{layer_name}' could "
                    f"not be read and was skipped ({exc})"
                )
                continue
            if not isinstance(frame, gpd.GeoDataFrame) or frame.empty:
                continue
            if frame.crs is None:
                frame = frame.set_crs(KML_CRS)
            else:
                frame = frame.to_crs(KML_CRS)
            frames.append(frame)
    except Exception as exc:
        print(
            f"WARNING: {kml_path.name}: layer enumeration failed ({exc}); "
            f"falling back to a single-layer read"
        )

    if not frames:
        frame = gpd.read_file(kml_path, driver="KML")
        if frame.crs is None:
            frame = frame.set_crs(KML_CRS)
        else:
            frame = frame.to_crs(KML_CRS)
        frames = [frame]

    merged = pd.concat(frames, ignore_index=True)
    return gpd.GeoDataFrame(merged, geometry="geometry", crs=KML_CRS)


def _parcel_name_from_row(row) -> str:
    """Original Placemark name, or empty string if the KML has none.

    Different GDAL drivers capitalise the field differently (Name / name).
    """
    for key in ("Name", "name", "NAME"):
        value = row.get(key)
        if value is not None and str(value).strip() and str(value) != "None":
            return str(value).strip()
    return ""


def load_parcels_from_kml(kml_path: str | Path) -> gpd.GeoDataFrame:
    """Read all parcel polygons from ONE KML file (all folders/layers).

    Returns columns: parcel_id, name, source_file, geometry.
    IDs are unique within this file (P001, P002, ...). For IDs that stay
    unique across several files use load_parcels_from_directory().
    Every skipped feature is reported on the console.
    """
    kml_path = Path(kml_path)
    if not kml_path.exists():
        raise FileNotFoundError(f"KML file not found: {kml_path}")

    raw = _read_all_kml_layers(kml_path)

    records = []
    skipped = 0
    for feature_index, (_, row) in enumerate(raw.iterrows(), start=1):
        geom, reason = _clean_geometry_with_reason(row.geometry)
        if geom is None:
            skipped += 1
            print(
                f"WARNING: {kml_path.name}: feature {feature_index} "
                f"skipped - {reason}"
            )
            continue
        records.append(
            {
                "name": _parcel_name_from_row(row),
                "source_file": kml_path.name,
                "geometry": geom,
            }
        )

    if not records:
        raise ValueError(
            f"No valid parcel geometries found in {kml_path} "
            f"({skipped} feature(s) skipped)"
        )
    if skipped:
        print(
            f"NOTE: {kml_path.name}: kept {len(records)} of {len(raw)} "
            f"features ({skipped} skipped)"
        )

    parcels = gpd.GeoDataFrame(records, geometry="geometry", crs=KML_CRS)
    parcels.insert(
        0, "parcel_id", [f"P{i:03d}" for i in range(1, len(parcels) + 1)]
    )
    return parcels.reset_index(drop=True)


def load_parcels_from_directory(kml_dir: str | Path) -> gpd.GeoDataFrame:
    """Read EVERY *.kml file in a folder and assign global unique IDs.

    Files are processed in alphabetical order, so the IDs are reproducible:
    the first parcel of the first file becomes P001, and so on.
    Original parcel names and source filenames are preserved in every row.
    """
    kml_dir = Path(kml_dir)
    kml_files = sorted(kml_dir.glob("*.kml"))
    if not kml_files:
        raise FileNotFoundError(f"No .kml files found in {kml_dir}")

    frames = [load_parcels_from_kml(path) for path in kml_files]
    combined = gpd.GeoDataFrame(
        pd.concat(frames, ignore_index=True), geometry="geometry", crs=KML_CRS
    )
    # Re-assign IDs so they stay unique across ALL files.
    combined["parcel_id"] = [f"P{i:03d}" for i in range(1, len(combined) + 1)]
    return combined


# --------------------------------------------------------------------------
# Areas and exports
# --------------------------------------------------------------------------

def calculate_areas_ha(parcels: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Add an area_ha column, measured in the equal-area CRS (EPSG:3035)."""
    projected = parcels.to_crs(AREA_CRS)               # degrees -> meters
    out = parcels.copy()
    out["area_ha"] = projected.geometry.area / 10_000  # m^2 -> hectares
    return out


def export_geojson(parcels: gpd.GeoDataFrame, output_path: str | Path) -> Path:
    """Save full parcel geometries (EPSG:4326) as a GeoJSON file."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    parcels.to_file(output_path, driver="GeoJSON")
    return output_path


def export_csv(parcels: gpd.GeoDataFrame, output_path: str | Path) -> Path:
    """Save a small attribute table (no geometry blob) as CSV."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Centroids are computed in the projected CRS (avoids the
    # "geographic CRS" warning), then converted back to lat/lon.
    centroids = parcels.to_crs(AREA_CRS).geometry.centroid.to_crs(KML_CRS)

    table = parcels.drop(columns="geometry").copy()
    table["centroid_lon"] = centroids.x
    table["centroid_lat"] = centroids.y
    table.to_csv(output_path, index=False)
    return output_path


# --------------------------------------------------------------------------
# Full pipelines (used by the future notebook)
# --------------------------------------------------------------------------

def process_kml_file(
    kml_path: str | Path,
    results_dir: str | Path | None = None,
) -> gpd.GeoDataFrame:
    """Full pipeline for ONE KML file: load -> clean -> area -> export."""
    parcels = calculate_areas_ha(load_parcels_from_kml(kml_path))
    if results_dir is not None:
        stem = Path(kml_path).stem
        export_geojson(parcels, Path(results_dir) / f"{stem}_parcels.geojson")
        export_csv(parcels, Path(results_dir) / f"{stem}_parcels.csv")
    return parcels


def process_kml_directory(
    kml_dir: str | Path,
    results_dir: str | Path | None = None,
) -> gpd.GeoDataFrame:
    """Full pipeline for EVERY KML in a folder (global unique IDs)."""
    parcels = calculate_areas_ha(load_parcels_from_directory(kml_dir))
    if results_dir is not None:
        results_dir = Path(results_dir)
        export_geojson(parcels, results_dir / "all_parcels.geojson")
        export_csv(parcels, results_dir / "all_parcels.csv")
    return parcels


# --------------------------------------------------------------------------
# Synthetic self-test (run with: python -m src.parcels)
# --------------------------------------------------------------------------

if __name__ == "__main__":
    root = Path(__file__).resolve().parents[1]
    scratch = root / "data" / "results"
    scratch.mkdir(parents=True, exist_ok=True)

    # ---- 1) Area test with synthetic rectangles -------------------------
    # Near Stara Zagora (~42.4 deg N):
    #   field_1: 0.01 deg lon (~822 m) x 0.01 deg lat (~1106 m)  ~= 91 ha
    #   field_2: 0.02 deg lon (~1644 m) x 0.01 deg lat (~1106 m) ~= 182 ha
    # The exact assertion below does NOT trust these comments: it compares
    # EPSG:3035 areas against independent geodesic (ellipsoid) areas.
    synthetic = gpd.GeoDataFrame(
        {
            "name": ["test_field_1", "test_field_2"],
            "source_file": ["synthetic", "synthetic"],
            "geometry": [
                box(25.30, 42.40, 25.31, 42.41),
                box(25.35, 42.42, 25.37, 42.43),
            ],
        },
        crs=KML_CRS,
    )
    synthetic.insert(0, "parcel_id", ["T001", "T002"])

    result = calculate_areas_ha(synthetic)
    print(result[["parcel_id", "name", "area_ha"]].to_string(index=False))

    geod = Geod(ellps="WGS84")
    for _, row in result.iterrows():
        reference_m2, _perimeter = geod.geometry_area_perimeter(row.geometry)
        reference_ha = abs(reference_m2) / 10_000
        relative_error = abs(row["area_ha"] - reference_ha) / reference_ha
        assert relative_error < 0.01, f"{row['parcel_id']}: area off by >1%"
    print("Area cross-check passed (EPSG:3035 vs geodesic).")

    # ---- 2) Geometry validation test ------------------------------------
    bowtie = Polygon([(0, 0), (1, 1), (1, 0), (0, 1)])  # self-intersecting
    assert clean_geometry(bowtie) is not None

    # GeometryCollection with a usable polygon AND useless line debris:
    # the polygon must survive, the line must be dropped.
    mixed = GeometryCollection(
        [box(25.0, 42.0, 25.01, 42.01), LineString([(0, 0), (1, 1)])]
    )
    fixed = clean_geometry(mixed)
    assert fixed is not None
    assert fixed.geom_type in ("Polygon", "MultiPolygon")

    assert clean_geometry(None) is None
    assert clean_geometry(Polygon()) is None
    print("Geometry validation checks passed.")

    # ---- 3) End-to-end test: synthetic KML with a nested Folder ---------
    synthetic_kml = scratch / "_synthetic_fields.kml"
    synthetic_kml.write_text(
        """<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
  <Document>
    <name>synthetic_doc</name>
    <Folder>
      <name>fields_folder</name>
      <Placemark>
        <name>alpha</name>
        <Polygon><outerBoundaryIs><LinearRing><coordinates>
          25.30,42.40,0 25.31,42.40,0 25.31,42.41,0 25.30,42.41,0 25.30,42.40,0
        </coordinates></LinearRing></outerBoundaryIs></Polygon>
      </Placemark>
      <Placemark>
        <name>beta</name>
        <Polygon><outerBoundaryIs><LinearRing><coordinates>
          25.35,42.42,0 25.37,42.42,0 25.37,42.43,0 25.35,42.43,0 25.35,42.42,0
        </coordinates></LinearRing></outerBoundaryIs></Polygon>
      </Placemark>
    </Folder>
  </Document>
</kml>
""",
        encoding="utf-8",
    )

    parcels = process_kml_file(synthetic_kml, scratch)
    assert len(parcels) == 2, f"expected 2 parcels, got {len(parcels)}"
    assert sorted(parcels["name"]) == ["alpha", "beta"], "original names lost"
    assert parcels["parcel_id"].is_unique
    assert set(parcels["source_file"]) == {"_synthetic_fields.kml"}
    print("Synthetic KML read passed (2 Placemarks inside a nested Folder).")

    # Keep the example exports, remove only the scratch KML input.
    export_csv(result, scratch / "synthetic_parcels.csv")
    export_geojson(result, scratch / "synthetic_parcels.geojson")
    synthetic_kml.unlink()

    # ---- 4) Real files, if any were dropped into data/kml/ --------------
    kml_dir = root / "data" / "kml"
    if sorted(kml_dir.glob("*.kml")):
        combined = process_kml_directory(kml_dir, scratch)
        assert combined["parcel_id"].is_unique
        print(
            f"Processed {len(combined)} parcels from data/kml/, "
            f"total {combined['area_ha'].sum():.2f} ha"
        )
    else:
        print("No .kml files in data/kml/ yet - drop some in and rerun.")

    print("Self-test passed.")
