# Step 4.1 - Sentinel-2 time-series quality assessment (parcel P001)

- **Input (read-only):** `C:\Users\User\Documents\SoftUni_Sentinel_Monitoring\data\results\time_series_p001_historical.csv`
- **Reproduce:** `python -m src.quality_assessment`
- **Observations:** 397 total; **196 valid** (`valid_count > 0`); 201 fully masked.
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
| valid_count (px) | 196 | 1.0000 | 96.0000 | 96.0000 | 96.0000 | 96.0000 | 0.0000 | 0.0000 | no | no |
| valid_coverage | 196 | 0.0048 | 0.4638 | 0.4638 | 0.4638 | 0.4638 | 0.0000 | 0.0000 | no | no |
| clear_fraction | 196 | 0.0104 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0.0000 | 0.0000 | no | no |
| dNDVI (interval) | 195 | -0.3870 | -0.0291 | -0.0031 | 0.0321 | 0.2382 | 0.0612 | 0.0305 | yes | yes |
| dNDMI (interval) | 195 | -0.6366 | -0.0355 | -0.0018 | 0.0293 | 0.5035 | 0.0648 | 0.0313 | yes | yes |
| dNDVI per day | 195 | -0.1283 | -0.0087 | -0.0008 | 0.0084 | 0.1716 | 0.0171 | 0.0089 | yes | yes |
| dNDMI per day | 195 | -0.2122 | -0.0105 | -0.0004 | 0.0092 | 0.1144 | 0.0197 | 0.0097 | yes | yes |

## 3. Coverage summary

- Complete, fully clear valid observations (`clear_fraction = 1.0`): **158**.
- Valid but partially masked observations (`0 < clear_fraction < 1`): **38** (low-quality candidates).
- Fully masked observations (`valid_count = 0`): **201**.

## 4. Flagged observations (potential anomalies - review, not rejection)

| Date | NDVI | NDMI | valid_count | clear_fraction | coverage | gap (d) | dNDVI | dNDMI | dNDVI/day | dNDMI/day | Reasons |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 2022-10-04 | 0.071 | -0.121 | 24 | 0.250 | partial | 2 | -0.098 | +0.044 | -0.049 | +0.022 | reduced in-polygon clear coverage (clear_fraction=0.250); statistical outlier: ndvi_rate (low) |
| 2022-10-09 | 0.166 | -0.138 | 88 | 0.917 | partial | 2 | +0.001 | +0.029 | +0.000 | +0.015 | reduced in-polygon clear coverage (clear_fraction=0.917) |
| 2022-10-12 | 0.157 | -0.142 | 78 | 0.812 | partial | 3 | -0.009 | -0.004 | -0.003 | -0.001 | reduced in-polygon clear coverage (clear_fraction=0.812) |
| 2022-10-17 | 0.155 | -0.098 | 50 | 0.521 | partial | 5 | -0.001 | +0.044 | -0.000 | +0.009 | reduced in-polygon clear coverage (clear_fraction=0.521) |
| 2022-11-03 | 0.162 | -0.161 | 63 | 0.656 | partial | 2 | +0.018 | +0.011 | +0.009 | +0.006 | reduced in-polygon clear coverage (clear_fraction=0.656) |
| 2022-11-18 | 0.116 | -0.150 | 57 | 0.594 | partial | 2 | -0.028 | +0.066 | -0.014 | +0.033 | reduced in-polygon clear coverage (clear_fraction=0.594) |
| 2022-12-21 | 0.307 | -0.153 | 96 | 1.000 | complete | 25 | +0.140 | +0.074 | +0.006 | +0.003 | statistical outlier: ndvi_change (high) |
| 2022-12-26 | 0.302 | -0.177 | 63 | 0.656 | partial | 5 | -0.005 | -0.023 | -0.001 | -0.005 | reduced in-polygon clear coverage (clear_fraction=0.656) |
| 2023-01-02 | 0.366 | -0.137 | 96 | 1.000 | complete | 2 | +0.079 | +0.015 | +0.040 | +0.008 | statistical outlier: ndvi_rate (high) |
| 2023-02-01 | 0.680 | 0.154 | 96 | 1.000 | complete | 25 | +0.238 | +0.233 | +0.010 | +0.009 | statistical outlier: ndvi_change (high); statistical outlier: ndmi_change (high) |
| 2023-02-09 | 0.532 | 0.318 | 96 | 1.000 | complete | 8 | -0.147 | +0.164 | -0.018 | +0.020 | statistical outlier: ndvi_change (low); statistical outlier: ndmi_change (high) |
| 2023-02-11 | 0.638 | 0.123 | 96 | 1.000 | complete | 2 | +0.106 | -0.194 | +0.053 | -0.097 | statistical outlier: ndmi_change (low); statistical outlier: ndvi_rate (high); statistical outlier: ndmi_rate (low) |
| 2023-02-14 | 0.530 | 0.043 | 33 | 0.344 | partial | 3 | -0.109 | -0.081 | -0.036 | -0.027 | reduced in-polygon clear coverage (clear_fraction=0.344); statistical outlier: ndvi_rate (low) |
| 2023-02-16 | 0.642 | 0.125 | 96 | 1.000 | complete | 2 | +0.112 | +0.082 | +0.056 | +0.041 | statistical outlier: ndvi_rate (high); statistical outlier: ndmi_rate (high) |
| 2023-03-06 | 0.807 | 0.345 | 96 | 1.000 | complete | 18 | +0.165 | +0.220 | +0.009 | +0.012 | statistical outlier: ndvi_change (high); statistical outlier: ndmi_change (high) |
| 2023-04-25 | 0.799 | 0.519 | 18 | 0.188 | partial | 23 | -0.044 | +0.046 | -0.002 | +0.002 | reduced in-polygon clear coverage (clear_fraction=0.188) |
| 2023-04-30 | 0.819 | 0.607 | 1 | 0.010 | partial | 5 | +0.020 | +0.088 | +0.004 | +0.018 | reduced in-polygon clear coverage (clear_fraction=0.010) |
| 2023-05-02 | 0.887 | 0.562 | 96 | 1.000 | complete | 2 | +0.069 | -0.045 | +0.034 | -0.023 | statistical outlier: ndvi_rate (high) |
| 2023-05-05 | 0.674 | 0.502 | 37 | 0.385 | partial | 3 | -0.214 | -0.060 | -0.071 | -0.020 | reduced in-polygon clear coverage (clear_fraction=0.385); statistical outlier: ndvi_change (low); statistical outlier: ndvi_rate (low) |
| 2023-05-07 | 0.856 | 0.525 | 96 | 1.000 | complete | 2 | +0.183 | +0.023 | +0.091 | +0.012 | statistical outlier: ndvi_change (high); statistical outlier: ndvi_rate (high) |
| 2023-05-22 | 0.701 | 0.356 | 96 | 1.000 | complete | 15 | -0.155 | -0.168 | -0.010 | -0.011 | statistical outlier: ndvi_change (low); statistical outlier: ndmi_change (low) |
| 2023-06-01 | 0.596 | 0.157 | 96 | 1.000 | complete | 10 | -0.104 | -0.199 | -0.010 | -0.020 | statistical outlier: ndmi_change (low) |
| 2023-06-09 | 0.454 | 0.120 | 91 | 0.948 | partial | 3 | -0.028 | +0.025 | -0.009 | +0.008 | reduced in-polygon clear coverage (clear_fraction=0.948) |
| 2023-06-11 | 0.329 | 0.057 | 96 | 1.000 | complete | 2 | -0.126 | -0.063 | -0.063 | -0.032 | statistical outlier: ndvi_change (low); statistical outlier: ndvi_rate (low) |
| 2023-09-02 | 0.153 | -0.278 | 96 | 1.000 | complete | 50 | -0.029 | -0.181 | -0.001 | -0.004 | statistical outlier: ndmi_change (low) |
| 2023-09-17 | 0.223 | -0.279 | 38 | 0.396 | partial | 3 | +0.022 | +0.017 | +0.007 | +0.006 | reduced in-polygon clear coverage (clear_fraction=0.396) |
| 2023-10-04 | 0.462 | -0.107 | 96 | 1.000 | complete | 10 | +0.180 | +0.155 | +0.018 | +0.016 | statistical outlier: ndvi_change (high); statistical outlier: ndmi_change (high) |
| 2023-10-14 | 0.623 | 0.006 | 96 | 1.000 | complete | 10 | +0.160 | +0.114 | +0.016 | +0.011 | statistical outlier: ndvi_change (high) |
| 2023-10-22 | 0.733 | 0.166 | 45 | 0.469 | partial | 8 | +0.110 | +0.160 | +0.014 | +0.020 | reduced in-polygon clear coverage (clear_fraction=0.469); statistical outlier: ndmi_change (high) |
| 2023-10-24 | 0.595 | 0.068 | 83 | 0.865 | partial | 2 | -0.137 | -0.098 | -0.069 | -0.049 | reduced in-polygon clear coverage (clear_fraction=0.865); statistical outlier: ndvi_change (low); statistical outlier: ndvi_rate (low); statistical outlier: ndmi_rate (low) |
| 2024-01-22 | 0.258 | 0.632 | 61 | 0.635 | partial | 5 | -0.387 | +0.503 | -0.077 | +0.101 | reduced in-polygon clear coverage (clear_fraction=0.635); statistical outlier: ndvi_change (low); statistical outlier: ndmi_change (high); statistical outlier: ndvi_rate (low); statistical outlier: ndmi_rate (high) |
| 2024-01-25 | 0.427 | -0.004 | 41 | 0.427 | partial | 3 | +0.169 | -0.637 | +0.056 | -0.212 | reduced in-polygon clear coverage (clear_fraction=0.427); statistical outlier: ndvi_change (high); statistical outlier: ndmi_change (low); statistical outlier: ndvi_rate (high); statistical outlier: ndmi_rate (low) |
| 2024-02-24 | 0.408 | -0.074 | 93 | 0.969 | partial | 15 | +0.001 | +0.029 | +0.000 | +0.002 | reduced in-polygon clear coverage (clear_fraction=0.969) |
| 2024-02-26 | 0.393 | -0.008 | 55 | 0.573 | partial | 2 | -0.016 | +0.065 | -0.008 | +0.033 | reduced in-polygon clear coverage (clear_fraction=0.573) |
| 2024-03-17 | 0.312 | -0.155 | 96 | 1.000 | complete | 20 | -0.080 | -0.147 | -0.004 | -0.007 | statistical outlier: ndmi_change (low) |
| 2024-03-22 | 0.309 | -0.142 | 72 | 0.750 | partial | 2 | +0.014 | +0.033 | +0.007 | +0.016 | reduced in-polygon clear coverage (clear_fraction=0.750) |
| 2024-04-06 | 0.374 | 0.023 | 96 | 1.000 | complete | 2 | +0.121 | +0.229 | +0.060 | +0.114 | statistical outlier: ndmi_change (high); statistical outlier: ndvi_rate (high); statistical outlier: ndmi_rate (high) |
| 2024-04-09 | 0.274 | -0.191 | 96 | 1.000 | complete | 3 | -0.100 | -0.214 | -0.033 | -0.071 | statistical outlier: ndmi_change (low); statistical outlier: ndmi_rate (low) |
| 2024-05-14 | 0.490 | 0.061 | 96 | 1.000 | complete | 8 | +0.206 | +0.219 | +0.026 | +0.027 | statistical outlier: ndvi_change (high); statistical outlier: ndmi_change (high) |
| 2024-05-16 | 0.515 | 0.195 | 96 | 1.000 | complete | 2 | +0.024 | +0.134 | +0.012 | +0.067 | statistical outlier: ndmi_change (high); statistical outlier: ndmi_rate (high) |
| 2024-05-21 | 0.678 | 0.241 | 96 | 1.000 | complete | 5 | +0.164 | +0.046 | +0.033 | +0.009 | statistical outlier: ndvi_change (high) |
| 2024-06-10 | 0.630 | 0.232 | 96 | 1.000 | complete | 2 | -0.107 | -0.051 | -0.054 | -0.025 | statistical outlier: ndvi_rate (low) |
| 2024-06-15 | 0.796 | 0.379 | 96 | 1.000 | complete | 5 | +0.167 | +0.147 | +0.033 | +0.029 | statistical outlier: ndvi_change (high); statistical outlier: ndmi_change (high) |
| 2024-07-20 | 0.467 | 0.110 | 69 | 0.719 | partial | 2 | -0.074 | -0.007 | -0.037 | -0.004 | reduced in-polygon clear coverage (clear_fraction=0.719); statistical outlier: ndvi_rate (low) |
| 2024-07-23 | 0.693 | 0.315 | 96 | 1.000 | complete | 3 | +0.227 | +0.205 | +0.076 | +0.068 | statistical outlier: ndvi_change (high); statistical outlier: ndmi_change (high); statistical outlier: ndvi_rate (high); statistical outlier: ndmi_rate (high) |
| 2024-08-04 | 0.675 | 0.283 | 96 | 1.000 | complete | 2 | +0.045 | +0.091 | +0.022 | +0.046 | statistical outlier: ndmi_rate (high) |
| 2024-08-07 | 0.476 | 0.075 | 8 | 0.083 | partial | 3 | -0.199 | -0.208 | -0.066 | -0.069 | reduced in-polygon clear coverage (clear_fraction=0.083); statistical outlier: ndvi_change (low); statistical outlier: ndmi_change (low); statistical outlier: ndvi_rate (low); statistical outlier: ndmi_rate (low) |
| 2024-08-09 | 0.606 | 0.185 | 70 | 0.729 | partial | 2 | +0.130 | +0.110 | +0.065 | +0.055 | reduced in-polygon clear coverage (clear_fraction=0.729); statistical outlier: ndvi_change (high); statistical outlier: ndvi_rate (high); statistical outlier: ndmi_rate (high) |
| 2024-08-14 | 0.414 | -0.036 | 96 | 1.000 | complete | 2 | -0.146 | -0.158 | -0.073 | -0.079 | statistical outlier: ndvi_change (low); statistical outlier: ndmi_change (low); statistical outlier: ndvi_rate (low); statistical outlier: ndmi_rate (low) |
| 2024-08-17 | 0.218 | -0.162 | 96 | 1.000 | complete | 3 | -0.196 | -0.126 | -0.065 | -0.042 | statistical outlier: ndvi_change (low); statistical outlier: ndvi_rate (low); statistical outlier: ndmi_rate (low) |
| 2024-08-27 | 0.113 | -0.139 | 31 | 0.323 | partial | 3 | -0.051 | +0.010 | -0.017 | +0.003 | reduced in-polygon clear coverage (clear_fraction=0.323) |
| 2024-09-11 | 0.144 | -0.005 | 35 | 0.365 | partial | 3 | -0.017 | +0.113 | -0.006 | +0.038 | reduced in-polygon clear coverage (clear_fraction=0.365) |
| 2024-09-23 | 0.178 | -0.227 | 96 | 1.000 | complete | 12 | +0.034 | -0.222 | +0.003 | -0.019 | statistical outlier: ndmi_change (low) |
| 2024-10-13 | 0.196 | -0.268 | 12 | 0.125 | partial | 5 | +0.007 | +0.006 | +0.001 | +0.001 | reduced in-polygon clear coverage (clear_fraction=0.125) |
| 2024-12-10 | 0.355 | -0.035 | 54 | 0.562 | partial | 15 | +0.078 | +0.198 | +0.005 | +0.013 | reduced in-polygon clear coverage (clear_fraction=0.562); statistical outlier: ndmi_change (high) |
| 2024-12-12 | 0.429 | 0.032 | 14 | 0.146 | partial | 2 | +0.074 | +0.067 | +0.037 | +0.034 | reduced in-polygon clear coverage (clear_fraction=0.146); statistical outlier: ndvi_rate (high) |
| 2024-12-17 | 0.412 | -0.139 | 96 | 1.000 | complete | 5 | -0.017 | -0.171 | -0.003 | -0.034 | statistical outlier: ndmi_change (low) |
| 2025-01-24 | 0.531 | 0.131 | 87 | 0.906 | partial | 5 | -0.094 | +0.090 | -0.019 | +0.018 | reduced in-polygon clear coverage (clear_fraction=0.906) |
| 2025-01-26 | 0.529 | 0.131 | 75 | 0.781 | partial | 2 | -0.002 | +0.000 | -0.001 | +0.000 | reduced in-polygon clear coverage (clear_fraction=0.781) |
| 2025-02-08 | 0.620 | 0.197 | 33 | 0.344 | partial | 10 | -0.010 | +0.067 | -0.001 | +0.007 | reduced in-polygon clear coverage (clear_fraction=0.344) |
| 2025-03-05 | 0.766 | 0.293 | 96 | 1.000 | complete | 18 | +0.153 | +0.131 | +0.008 | +0.007 | statistical outlier: ndvi_change (high); statistical outlier: ndmi_change (high) |
| 2025-03-20 | 0.893 | 0.535 | 96 | 1.000 | complete | 10 | +0.088 | +0.159 | +0.009 | +0.016 | statistical outlier: ndmi_change (high) |
| 2025-04-13 | 0.884 | 0.566 | 74 | 0.771 | partial | 22 | +0.013 | +0.041 | +0.001 | +0.002 | reduced in-polygon clear coverage (clear_fraction=0.771) |
| 2025-04-21 | 0.869 | 0.552 | 59 | 0.615 | partial | 5 | -0.001 | +0.014 | -0.000 | +0.003 | reduced in-polygon clear coverage (clear_fraction=0.615) |
| 2025-05-11 | 0.839 | 0.552 | 7 | 0.073 | partial | 5 | -0.015 | +0.011 | -0.003 | +0.002 | reduced in-polygon clear coverage (clear_fraction=0.073) |
| 2025-05-24 | 0.686 | 0.417 | 8 | 0.083 | partial | 1 | -0.128 | -0.077 | -0.128 | -0.077 | reduced in-polygon clear coverage (clear_fraction=0.083); statistical outlier: ndvi_change (low); statistical outlier: ndvi_rate (low); statistical outlier: ndmi_rate (low) |
| 2025-05-29 | 0.662 | 0.413 | 85 | 0.885 | partial | 5 | -0.024 | -0.004 | -0.005 | -0.001 | reduced in-polygon clear coverage (clear_fraction=0.885) |
| 2025-05-31 | 0.788 | 0.490 | 30 | 0.312 | partial | 2 | +0.127 | +0.077 | +0.063 | +0.038 | reduced in-polygon clear coverage (clear_fraction=0.312); statistical outlier: ndvi_change (high); statistical outlier: ndvi_rate (high) |
| 2025-06-02 | 0.602 | 0.427 | 96 | 1.000 | complete | 2 | -0.187 | -0.063 | -0.093 | -0.031 | statistical outlier: ndvi_change (low); statistical outlier: ndvi_rate (low) |
| 2025-06-03 | 0.773 | 0.407 | 96 | 1.000 | complete | 1 | +0.172 | -0.020 | +0.172 | -0.020 | statistical outlier: ndvi_change (high); statistical outlier: ndvi_rate (high) |
| 2025-06-05 | 0.679 | 0.335 | 96 | 1.000 | complete | 2 | -0.094 | -0.072 | -0.047 | -0.036 | statistical outlier: ndvi_rate (low) |
| 2025-06-10 | 0.490 | 0.243 | 64 | 0.667 | partial | 2 | -0.123 | -0.053 | -0.061 | -0.026 | reduced in-polygon clear coverage (clear_fraction=0.667); statistical outlier: ndvi_change (low); statistical outlier: ndvi_rate (low) |
| 2025-06-12 | 0.426 | 0.156 | 96 | 1.000 | complete | 2 | -0.064 | -0.088 | -0.032 | -0.044 | statistical outlier: ndmi_rate (low) |
| 2025-06-15 | 0.347 | 0.097 | 96 | 1.000 | complete | 2 | -0.082 | -0.039 | -0.041 | -0.019 | statistical outlier: ndvi_rate (low) |

## 5. Interpretation

- Statistical outliers (IQR / MAD) are values far from the parcel's own distribution of that metric; they are flagged for review only.
- Low-quality observations are those with reduced in-polygon clear coverage. Because most valid days have `clear_fraction = 1.0`, the robust rules are degenerate there and these days are identified by the definitional coverage flag instead.
- Outliers concentrated on low-`clear_fraction` days should be treated as less reliable, not as confirmed land-cover change; no automatic removal or arbitrary thresholding is applied.
