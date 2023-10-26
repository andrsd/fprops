import fprops as fp
import pytest
import math

def test_n2_valid():
    n2 = fp.Nitrogen()
    state = n2.p_T(1e6, 280)
    assert state.p == 1e6
    assert state.T == 280
    assert math.isclose(state.rho, 12.074993450711256, abs_tol=1e-14)

@pytest.mark.skip(reason="Waiting on range checking")
def test_n2_invalid():
    with pytest.raises(RuntimeError):
        n2 = fp.Nitrogen()
        state = n2.p_T(0.1e6, 60)
