# Sentinel-2 Crop Monitoring for CAP - Bulgaria

Retrospective crop verification using Sentinel-2 L2A imagery.


## Progress

- [x] Step 1 - Parcel boundaries from KML and area calculations
- [x] Step 2 - Sentinel-2 L2A retrieval and parcel-level statistics using Sentinel Hub Statistical API
- [x] Step 3 - Vegetation indices: NDVI and NDMI, with quality assessment
- [x] Step 4 - Phenology and time-series analysis, including descriptive statistics, regression, anomaly detection, and temporal changes
- [x] Step 5 - Retrospective analysis of historical Sentinel-2 observations for parcel P001 (2023–2024)

### Scope and limitations

The project demonstrates retrospective Sentinel-2 crop monitoring using Python and mathematical time-series analysis.

The historical analysis includes 141 observation intervals, of which 70 are valid and 58 have complete coverage.

NDYI is not included in the verified analysis. Crop identity cannot be conclusively established from NDVI and NDMI alone. A validated preliminary-versus-historical crop comparison has not yet been demonstrated.


## Run the self-test (Step 1)

```
pip install -r requirements.txt
python -m src.parcels
```

Drop your parcel `.kml` files into `data/kml/` and rerun the command.
Results are written to `data/results/` (not committed to git).

## Methodology notes

### NDVI / NDMI statistics (Step 2)

Daily per-parcel statistics are computed with the CDSE Sentinel Hub
Statistical API over Sentinel-2 L2A. Each parcel geometry is reprojected to
its local UTM zone so that `resx`/`resy` are true metres, and a pixel is
treated as valid only when `dataMask == 1` and `SCL ∈ {4, 5}` (vegetation /
bare soil); clouds, shadow, snow, water and no-data are excluded.

Per interval the CSV reports two independent coverage measures:

- `geometry_coverage = geometry_pixel_count / sample_count` — the parcel's
  share of the request bounding box (pure geometry, cloud-independent);
- `clear_fraction = valid_count / geometry_pixel_count` — the cloud-free,
  SCL-valid fraction **inside** the polygon.

`valid_coverage = valid_count / sample_count` is the raw valid-pixel
coverage over the bounding box. Intervals with no scene over the parcel are
omitted; fully masked intervals keep their counts but leave the NDVI/NDMI
cells empty (no literal `nan`).

**Known limitation — partially intersecting boundary pixels.** Every
20 m pixel that intersects the polygon, however slightly, is counted in
`geometry_pixel_count` and contributes its full-pixel value to the mean.
For small parcels the boundary pixels can be a large share of the total
(e.g. parcel P001: ~96 intersecting pixels vs ~73 pixel-equivalents by area,
so ~24 partial edge pixels are included). This dilutes the reported means
towards neighbouring land cover. The pixel-inclusion rule is intentionally
left unchanged for now; a future step may weight or clip edge pixels.

