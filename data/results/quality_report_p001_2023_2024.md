# Step 4.1 - Sentinel-2 time-series quality assessment (parcel P001)

- **Input (read-only):** `C:\Users\User\Documents\SoftUni_Sentinel_Monitoring\data\results\time_series_p001_2023_2024.csv`
- **Reproduce:** `python -m src.quality_assessment`
- **Observations:** 141 total; **70 valid** (`valid_count > 0`); 71 fully masked.
- **No observation is removed, replaced, smoothed or interpolated**; the original CSV, methodology, masking rules and PNG figures are unchanged.

## 1. Coverage metrics: what each one measures

- `valid_count` - number of 20 m pixels inside the polygon that passed `dataMask == 1 AND SCL in {4, 5}`.
- `geometry_coverage = geometry_pixel_count / sample_count` - the polygon's share of the request bounding box. This is a **pure geometry ratio** and is **constant** for a fixed parcel: here 0.463768 on every observation.
- `valid_coverage = valid_count / sample_count` - valid pixels as a fraction of the **bounding box** (cover + cloud + non-agricultural SCL).
- `clear_fraction = valid_count / geometry_pixel_count` - the **cloud-free, SCL-valid fraction inside the polygon**. For P001 the bounding box is 23 x 9 = 207 pixels while the polygon contains 96, so `valid_coverage` and `clear_fraction` differ by the constant geometric factor 96/207 = 0.4638. They measure **different aspects of coverage**: `valid_coverage` mixes geometry with data availability, whereas `clear_fraction` isolates in-polygon masking.

## 2. Statistical rules (all computed from the data)

- **IQR fences:** a value is an outlier when it falls below `Q1 - 1.5 * IQR` or above `Q3 + 1.5 * IQR`.
- **MAD robust z-score:** `z = 0.6745 * (x - median) / MAD`; `|z| > 3.5` is an outlier (Iglewicz & Hoaglin).
- A value is flagged if **either** rule fires; IQR takes precedence if the two ever conflict.
- **Degenerate rules are skipped safely:** when `IQR == 0` the IQR rule is unusable and when `MAD == 0` the z-score is unusable. This is expected for `clear_fraction`, which equals 1.0 for most valid observations, so it receives **no** robust outlier flags.
- **Coverage completeness is a definitional flag, not a statistical one:** `clear_fraction < 1.0` means some in-polygon pixels were masked (cloud / SCL / no-data). It is reported separately so that *low-quality* observations are distinguished from *statistical outliers*.
- **Changes** are computed only between consecutive valid observations. `*_change` is the difference across the elapsed interval and `*_change_per_day` divides it by the elapsed days, so a gap of several days is never described as a daily change.

Computed dispersion thresholds per metric:

| Metric | n | min | Q1 | median | Q3 | max | IQR | MAD | IQR usable | MAD usable |
|---|---|---|---|---|---|---|---|---|---|---|
| valid_count (px) | 70 | 8.0000 | 96.0000 | 96.0000 | 96.0000 | 96.0000 | 0.0000 | 0.0000 | no | no |
| valid_coverage | 70 | 0.0386 | 0.4638 | 0.4638 | 0.4638 | 0.4638 | 0.0000 | 0.0000 | no | no |
| clear_fraction | 70 | 0.0833 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0.0000 | 0.0000 | no | no |
| dNDVI (interval) | 69 | -0.3870 | -0.0463 | -0.0077 | 0.0332 | 0.2269 | 0.0795 | 0.0405 | yes | yes |
| dNDMI (interval) | 69 | -0.6366 | -0.0416 | -0.0059 | 0.0328 | 0.5035 | 0.0744 | 0.0378 | yes | yes |
| dNDVI per day | 69 | -0.0774 | -0.0112 | -0.0015 | 0.0088 | 0.0756 | 0.0200 | 0.0103 | yes | yes |
| dNDMI per day | 69 | -0.2122 | -0.0109 | -0.0007 | 0.0114 | 0.1144 | 0.0223 | 0.0103 | yes | yes |

## 3. Coverage summary

- Complete, fully clear valid observations (`clear_fraction = 1.0`): **58**.
- Valid but partially masked observations (`0 < clear_fraction < 1`): **12** (low-quality candidates).
- Fully masked observations (`valid_count = 0`): **71**.

## 4. Flagged observations (potential anomalies - review, not rejection)

| Date | NDVI | NDMI | valid_count | clear_fraction | coverage | gap (d) | dNDVI | dNDMI | dNDVI/day | dNDMI/day | Reasons |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 2023-09-17 | 0.223 | -0.279 | 38 | 0.396 | partial | 3 | +0.022 | +0.017 | +0.007 | +0.006 | reduced in-polygon clear coverage (clear_fraction=0.396) |
| 2023-10-04 | 0.462 | -0.107 | 96 | 1.000 | complete | 10 | +0.180 | +0.155 | +0.018 | +0.016 | statistical outlier: ndvi_change (high); statistical outlier: ndmi_change (high) |
| 2023-10-14 | 0.623 | 0.006 | 96 | 1.000 | complete | 10 | +0.160 | +0.114 | +0.016 | +0.011 | statistical outlier: ndvi_change (high) |
| 2023-10-22 | 0.733 | 0.166 | 45 | 0.469 | partial | 8 | +0.110 | +0.160 | +0.014 | +0.020 | reduced in-polygon clear coverage (clear_fraction=0.469); statistical outlier: ndmi_change (high) |
| 2023-10-24 | 0.595 | 0.068 | 83 | 0.865 | partial | 2 | -0.137 | -0.098 | -0.069 | -0.049 | reduced in-polygon clear coverage (clear_fraction=0.865); statistical outlier: ndvi_rate (low); statistical outlier: ndmi_rate (low) |
| 2024-01-22 | 0.258 | 0.632 | 61 | 0.635 | partial | 5 | -0.387 | +0.503 | -0.077 | +0.101 | reduced in-polygon clear coverage (clear_fraction=0.635); statistical outlier: ndvi_change (low); statistical outlier: ndmi_change (high); statistical outlier: ndvi_rate (low); statistical outlier: ndmi_rate (high) |
| 2024-01-25 | 0.427 | -0.004 | 41 | 0.427 | partial | 3 | +0.169 | -0.637 | +0.056 | -0.212 | reduced in-polygon clear coverage (clear_fraction=0.427); statistical outlier: ndvi_change (high); statistical outlier: ndmi_change (low); statistical outlier: ndvi_rate (high); statistical outlier: ndmi_rate (low) |
| 2024-02-24 | 0.408 | -0.074 | 93 | 0.969 | partial | 15 | +0.001 | +0.029 | +0.000 | +0.002 | reduced in-polygon clear coverage (clear_fraction=0.969) |
| 2024-02-26 | 0.393 | -0.008 | 55 | 0.573 | partial | 2 | -0.016 | +0.065 | -0.008 | +0.033 | reduced in-polygon clear coverage (clear_fraction=0.573) |
| 2024-03-22 | 0.309 | -0.142 | 72 | 0.750 | partial | 2 | +0.014 | +0.033 | +0.007 | +0.016 | reduced in-polygon clear coverage (clear_fraction=0.750) |
| 2024-04-06 | 0.374 | 0.023 | 96 | 1.000 | complete | 2 | +0.121 | +0.229 | +0.060 | +0.114 | statistical outlier: ndmi_change (high); statistical outlier: ndvi_rate (high); statistical outlier: ndmi_rate (high) |
| 2024-04-09 | 0.274 | -0.191 | 96 | 1.000 | complete | 3 | -0.100 | -0.214 | -0.033 | -0.071 | statistical outlier: ndmi_change (low); statistical outlier: ndmi_rate (low) |
| 2024-05-14 | 0.490 | 0.061 | 96 | 1.000 | complete | 8 | +0.206 | +0.219 | +0.026 | +0.027 | statistical outlier: ndvi_change (high); statistical outlier: ndmi_change (high) |
| 2024-05-16 | 0.515 | 0.195 | 96 | 1.000 | complete | 2 | +0.024 | +0.134 | +0.012 | +0.067 | statistical outlier: ndmi_rate (high) |
| 2024-05-21 | 0.678 | 0.241 | 96 | 1.000 | complete | 5 | +0.164 | +0.046 | +0.033 | +0.009 | statistical outlier: ndvi_change (high) |
| 2024-06-10 | 0.630 | 0.232 | 96 | 1.000 | complete | 2 | -0.107 | -0.051 | -0.054 | -0.025 | statistical outlier: ndvi_rate (low) |
| 2024-06-15 | 0.796 | 0.379 | 96 | 1.000 | complete | 5 | +0.167 | +0.147 | +0.033 | +0.029 | statistical outlier: ndvi_change (high); statistical outlier: ndmi_change (high) |
| 2024-07-20 | 0.467 | 0.110 | 69 | 0.719 | partial | 2 | -0.074 | -0.007 | -0.037 | -0.004 | reduced in-polygon clear coverage (clear_fraction=0.719) |
| 2024-07-23 | 0.693 | 0.315 | 96 | 1.000 | complete | 3 | +0.227 | +0.205 | +0.076 | +0.068 | statistical outlier: ndvi_change (high); statistical outlier: ndmi_change (high); statistical outlier: ndvi_rate (high); statistical outlier: ndmi_rate (high) |
| 2024-08-04 | 0.675 | 0.283 | 96 | 1.000 | complete | 2 | +0.045 | +0.091 | +0.022 | +0.046 | statistical outlier: ndmi_rate (high) |
| 2024-08-07 | 0.476 | 0.075 | 8 | 0.083 | partial | 3 | -0.199 | -0.208 | -0.066 | -0.069 | reduced in-polygon clear coverage (clear_fraction=0.083); statistical outlier: ndvi_change (low); statistical outlier: ndmi_change (low); statistical outlier: ndvi_rate (low); statistical outlier: ndmi_rate (low) |
| 2024-08-09 | 0.606 | 0.185 | 70 | 0.729 | partial | 2 | +0.130 | +0.110 | +0.065 | +0.055 | reduced in-polygon clear coverage (clear_fraction=0.729); statistical outlier: ndvi_rate (high); statistical outlier: ndmi_rate (high) |
| 2024-08-14 | 0.414 | -0.036 | 96 | 1.000 | complete | 2 | -0.146 | -0.158 | -0.073 | -0.079 | statistical outlier: ndmi_change (low); statistical outlier: ndvi_rate (low); statistical outlier: ndmi_rate (low) |
| 2024-08-17 | 0.218 | -0.162 | 96 | 1.000 | complete | 3 | -0.196 | -0.126 | -0.065 | -0.042 | statistical outlier: ndvi_change (low); statistical outlier: ndvi_rate (low) |
| 2024-08-27 | 0.113 | -0.139 | 31 | 0.323 | partial | 3 | -0.051 | +0.010 | -0.017 | +0.003 | reduced in-polygon clear coverage (clear_fraction=0.323) |

## 5. Interpretation

- Statistical outliers (IQR / MAD) are values far from the parcel's own distribution of that metric; they are flagged for review only.
- Low-quality observations are those with reduced in-polygon clear coverage. Because most valid days have `clear_fraction = 1.0`, the robust rules are degenerate there and these days are identified by the definitional coverage flag instead.
- Outliers concentrated on low-`clear_fraction` days should be treated as less reliable, not as confirmed land-cover change; no automatic removal or arbitrary thresholding is applied.
