import fprops as fp
import pytest
import math
import re


def test_air_valid():
    air = fp.Air()
    state = air.p_T(101325, 300)
    assert state.p == 101325
    assert state.T == 300
    assert math.isclose(state.rho, 1.1769510785919943, abs_tol=1e-14)


def test_air_repr(capfd):
    air = fp.Air()
    state = air.p_T(101325, 300)
    print(state)
    out, err = capfd.readouterr()
    lines = out.splitlines()
    assert re.match("rho = [0-9.e-]+ kg/m\\^3", lines[0])
    assert re.match("p = [0-9.e-]+ Pa", lines[1])
    assert re.match("T = [0-9.e-]+ K", lines[2])
    assert re.match("e = [0-9.e-]+ J/kg", lines[3])
    assert re.match("v = [0-9.e-]+ m\\^3/kg", lines[4])
    assert re.match("cp = [0-9.e-]+ J/\\(kg-K\\)", lines[5])
    assert re.match("cv = [0-9.e-]+ J/\\(kg-K\\)", lines[6])
    assert re.match("s = [0-9.e-]+ J/\\(kg-K\\)", lines[7])
    assert re.match("h = [0-9.e-]+ J/kg", lines[8])
    assert re.match("c = [0-9.e-]+ m/s", lines[9])
    assert re.match("mu = [0-9.e-]+ Pa-s", lines[10])
    assert re.match("k = [0-9.e-]+ W/\\(m-K\\)", lines[11])
