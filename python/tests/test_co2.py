import pyfprops as fp
import pytest
import math

def test_co2_valid():
    co2 = fp.CarbonDioxide()
    state = co2.p_T(1e6, 280)
    assert state.p == 1e6
    assert state.T == 280
    assert math.isclose(state.rho, 20.199309000812121, abs_tol=1e-14)
