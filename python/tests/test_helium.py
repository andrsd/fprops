import fprops as fp
import pytest
import math

def test_he_valid():
    n2 = fp.Helium()
    state = n2.p_T(1e6, 280)
    assert state.p == 1e6
    assert state.T == 280
    assert math.isclose(state.rho, 1.7109055009783694, abs_tol=1e-14)
