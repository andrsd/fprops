#!/usr/bin/env python3

import argparse
import yaml
from CoolProp import AbstractState
from CoolProp.CoolProp import PhaseSI, PropsSI, get_global_param_string
import CoolProp.CoolProp as CoolProp

parser = argparse.ArgumentParser(description='Compute a range of thermodynamical states using CoolProp')
parser.add_argument('cfg_file', help='Config file for generating values (.yml)')
args = parser.parse_args()


def compute_q(q, p, T):
    """
    Compute a quantity (q) for a given state (pressure, Temperature)

    :param q: Quantity name
    :param p: pressure [Pa]
    :param T: Temperature [K]
    :return: The computed quantity or NaN for non-physical state
    """
    try:
        value = PropsSI(q, "P", p, "T", T, cfg['fluid'])
    except ValueError:
        value = float("nan")
    return value


def compute_state(p, T):
    """
    Compute a thermodynamical state for a given point (pressure, Temperature)
    """
    state = {
        "p": p,
        "T": T,
        "rho": compute_q("D", p, T),
        "c_p": compute_q("C", p, T),
        "c_v": compute_q("O", p, T),
        "h": compute_q("H", p, T),
        "s": compute_q("S", p, T),
        "u": compute_q("U", p, T),
        "w": compute_q("A", p, T),
        "mu": compute_q("V", p, T),
        "k": compute_q("L", p, T)
    }
    return state


with open(args.cfg_file) as f:
    inp = yaml.load(f, Loader=yaml.FullLoader)
    cfg = inp['fprops']

    states = []
    for p in cfg['p']:
        for T in cfg['T']:
            s = compute_state(p * 1e6, T)
            states.append(s)

    out = {
        'coolprop_version': get_global_param_string("version"),
        'fluid': cfg['fluid'],
        'data': states
    }
    print(yaml.dump(out))
