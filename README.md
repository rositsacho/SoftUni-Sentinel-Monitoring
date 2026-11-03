# Sentinel-2 Crop Monitoring for CAP - Bulgaria

Retrospective crop verification using Sentinel-2 L2A imagery.

## Progress

- [x] Step 1 - parcel boundaries from KML, areas in hectares
      (synthetic test passed; real KML files processed successfully)
- [ ] Step 2 - Sentinel-2 L2A retrieval
- [ ] Step 3 - vegetation indices (NDVI, NDMI, NDYI)
- [ ] Step 4 - phenology and time series
- [ ] Step 5 - preliminary vs historical comparison

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

