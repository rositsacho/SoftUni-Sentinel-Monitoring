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
