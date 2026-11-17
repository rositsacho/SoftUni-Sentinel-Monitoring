"""Mathematical analysis utilities for unevenly spaced satellite observations."""
import numpy as np
import pandas as pd

def normalized_difference(a,b):
    a=np.asarray(a,dtype=float);b=np.asarray(b,dtype=float)
    with np.errstate(divide="ignore",invalid="ignore"):
        return np.where(np.abs(a+b)>1e-12,(a-b)/(a+b),np.nan)


def rates_of_change(times, values):
    x = pd.to_datetime(times, utc=True)
    y = np.asarray(values, dtype=float)
    delta = np.diff(x.astype("int64")) / 86_400_000_000
    if np.any(delta <= 0):
        raise ValueError("Dates must be strictly increasing")
    return np.diff(y) / delta


def rolling_time_mean(df, value="ndvi_mean", window="21D"):
    series=df.set_index("date")[value].sort_index()
    return series.rolling(window,min_periods=2).mean()

def summary(df, value="ndvi_mean"):
    v=pd.to_numeric(df[value],errors="coerce").dropna()
    return {"n":int(v.size),"mean":float(v.mean()) if len(v) else np.nan,
            "std":float(v.std(ddof=1)) if len(v)>1 else np.nan,
            "min":float(v.min()) if len(v) else np.nan,
            "max":float(v.max()) if len(v) else np.nan}

def rmse(observed,predicted):
    a=np.asarray(observed,dtype=float);b=np.asarray(predicted,dtype=float)
    mask=np.isfinite(a)&np.isfinite(b)
    return float(np.sqrt(np.mean((a[mask]-b[mask])**2))) if np.any(mask) else np.nan

def extended_summary(df, value="ndvi_mean"):
    """Calculate extended statistics for NDVI/NDMI observations."""
    if value not in df.columns:
        raise KeyError(f"Missing column: {value}")

    values = pd.to_numeric(df[value], errors="coerce")
    values = values.replace([np.inf, -np.inf], np.nan)
    valid = values.dropna()

    if valid.empty:
        return {
            "count": 0,
            "missing": int(len(values)),
            "mean": np.nan,
            "median": np.nan,
            "std": np.nan,
            "min": np.nan,
            "max": np.nan,
            "q25": np.nan,
            "q75": np.nan,
            "cv_percent": np.nan,
        }

    mean = float(valid.mean())
    std = float(valid.std(ddof=1)) if len(valid) > 1 else np.nan

    return {
        "count": int(len(valid)),
        "missing": int(len(values) - len(valid)),
        "mean": mean,
        "median": float(valid.median()),
        "std": std,
        "min": float(valid.min()),
        "max": float(valid.max()),
        "q25": float(valid.quantile(0.25)),
        "q75": float(valid.quantile(0.75)),
        "cv_percent": (
            float(100 * std / abs(mean))
            if np.isfinite(std) and abs(mean) > 1e-12
            else np.nan
        ),
    }



def linear_trend(times, values):
    """Estimate a linear trend for unevenly spaced observations."""
    dates = pd.to_datetime(times, utc=True)
    y = np.asarray(values, dtype=float)

    if len(dates) != len(y):
        raise ValueError("Dates and values must have equal lengths")

    mask = np.isfinite(y) & ~pd.isna(dates)
    if mask.sum() < 3:
        raise ValueError("At least 3 valid observations are required")

    
    x = (
        (dates[mask] - dates[mask].min())
        / pd.Timedelta(days=1)
    ).to_numpy(dtype=float)

    x = x - x.min()
    y = y[mask]

    if len(np.unique(x)) < 2:
        raise ValueError("At least 2 distinct dates are required")

    slope, intercept = np.polyfit(x, y, 1)
    predicted = slope * x + intercept

    ss_res = np.sum((y - predicted) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)

    r_squared = (
        1 - ss_res / ss_tot
        if ss_tot > 1e-12
        else np.nan
    )

    return {
        "slope_per_day": float(slope),
        "intercept": float(intercept),
        "r_squared": float(r_squared),
        "n_observations": int(len(y)),
    }
