import fprops as fp
import pytest
import math

def test_air_valid():
    air = fp.Air()
    state = air.p_T(101325, 300)
    assert state.p == 101325
    assert state.T == 300
    assert math.isclose(state.rho, 1.1769510785919943, abs_tol=1e-14)
