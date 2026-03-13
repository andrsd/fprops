from CoolProp import AbstractState
from CoolProp.CoolProp import PhaseSI, PropsSI, get_global_param_string
import CoolProp.CoolProp as CoolProp
import math

fluid = "oxygen"
p = 101325
T = 300.
try:
    print("p =", PropsSI("P", "P", p, "T", T, fluid))
    print("T =", PropsSI("T", "P", p, "T", T, fluid))
    print("rho =", PropsSI("D", "P", p, "T", T, fluid))
    print("h =", PropsSI("H", "P", p, "T", T, fluid))
    print("e =", PropsSI("U", "P", p, "T", T, fluid))
    print("s =", PropsSI("S", "P", p, "T", T, fluid))
    print("c =", PropsSI("A", "P", p, "T", T, fluid))
    print("cv =", PropsSI("O", "P", p, "T", T, fluid))
    print("cp =", PropsSI("C", "P", p, "T", T, fluid))

    print("mu =", PropsSI("V", "P", p, "T", T, fluid))
    print("k =", PropsSI("L", "P", p, "T", T, fluid))
except:
    state = None
