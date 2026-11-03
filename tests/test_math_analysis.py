
import numpy as np
import pandas as pd
import pytest

from src.math_analysis import normalized_difference, rates_of_change, rmse


def test_normalized_difference():
    assert np.isclose(normalized_difference(0.6, 0.2), 0.5)


def test_zero_denominator():
    assert np.isnan(normalized_difference(0.0, 0.0))


def test_daily_rates():
    dates = pd.to_datetime(
        ["2024-06-01", "2024-06-03", "2024-06-08"],
        utc=True,
    )
    values = [0.30, 0.38, 0.43]

    result = rates_of_change(dates, values)

    np.testing.assert_allclose(result, [0.04, 0.01])


def test_invalid_date_order():
    dates = pd.to_datetime(
        ["2024-06-03", "2024-06-01"],
        utc=True,
    )

    with pytest.raises(ValueError):
        rates_of_change(dates, [0.30, 0.38])


def test_rmse_identical_predictions():
    assert np.isclose(rmse([1, 2], [1, 2]), 0.0)
