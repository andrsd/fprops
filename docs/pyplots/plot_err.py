#!/usr/bin/env python3

import argparse
import yaml
import sys
from CoolProp import AbstractState
from CoolProp.CoolProp import PhaseSI, PropsSI, get_global_param_string
import CoolProp.CoolProp as CoolProp
import pyfprops as fp
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker

parser = argparse.ArgumentParser(description='Compute a range of thermodynamical states using CoolProp')
parser.add_argument('cfg_file', help='Config file for generating values (.yml)')
args = parser.parse_args()

quantities = [
    ('Density', 'rho'),
    ('C_p', 'c_p'),
    ('C_v', 'c_v'),
    ('Enthalpy', 'h'),
    ('Entropy', 's'),
    ('Internal energy', 'u'),
    ('Sound speed', 'w'),
    ('Viscosity', 'mu'),
    ('Thermal conductivity', 'k')
]

def fmt(x):
    s = f"{x:.1f}"
    if s.endswith("0"):
        s = f"{x:.0f}"
    return rf"{s} \%" if plt.rcParams["text.usetex"] else f"{s} %"


def compute_q(q, p, T):
    """
    Compute a quantity (q) for a given state (pressure, Temperature)

    :param q: Quantity name
    :param p: pressure [Pa]
    :param T: Temperature [K]
    :return: The computed quantity or NaN for non-physical state
    """
    try:
        value = PropsSI(q, "P", p, "T", T, cfg['coolprop'])
    except ValueError:
        value = float("nan")
    return value


def compute_state_coolprop(p, T):
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


def compute_state_fprops(props, p, T):
    try:
        state = props.p_T(p, T)
        return {
            "p": p,
            "T": T,
            "rho": state.rho,
            "c_p": state.cp,
            "c_v": state.cv,
            "h": state.h,
            "s": state.s,
            "u": state.u,
            "w": state.w,
            "mu": state.mu,
            "k": state.k
        }
    except:
        return {
            "p": p,
            "T": T,
            "rho": float("nan"),
            "c_p": float("nan"),
            "c_v": float("nan"),
            "h": float("nan"),
            "s": float("nan"),
            "u": float("nan"),
            "w": float("nan"),
            "mu": float("nan"),
            "k": float("nan")
        }


def compute_error_vals():
    vals = {}
    for q in quantities:
        vals[q[1]] = np.zeros(shape=(len(T_range), len(p_range)))

    for p_idx, p in enumerate(p_range):
        for T_idx, T in enumerate(T_range):
            s1 = compute_state_coolprop(p * 1e6, T)
            s2 = compute_state_fprops(props, p * 1e6, T)
            for q in quantities:
                quantity = q[1]
                if not math.isnan(s1[quantity]) and not math.isnan(s2[quantity]):
                    # rel_err = math.fabs(s2[quantity] - s1[v]) / s1[quantity]
                    # err[T_idx, p_idx] = rel_err
                    abs_err = math.fabs(s2[quantity] - s1[quantity])
                    vals[quantity][T_idx, p_idx] = abs_err
                else:
                    vals[quantity][T_idx, p_idx] = float("nan")
    return vals


def create_plot(Quantity, vals):
    fig, ax = plt.subplots(1, 1)
    cp = ax.contourf(X, Y, vals, locator=ticker.LogLocator())
    fig.colorbar(cp, format='%g')
    ax.set_title('{}: Absolute error'.format(Quantity))
    ax.set_xlabel('Pressure [MPa]')
    ax.set_xscale('log', base=10)
    ax.set_ylabel('Temperature [K]')
    ax.set_yscale('log', base=10)


with open(args.cfg_file) as f:
    inp = yaml.load(f, Loader=yaml.FullLoader)
    cfg = inp['fprops']
    props = getattr(fp, cfg['fprops'])()
    p_range = np.array(cfg['p'])
    T_range = np.array(cfg['T'])
    X, Y = np.meshgrid(p_range, T_range)

    vals = compute_error_vals()
    for q in quantities:
        quantity = q[1]
        create_plot(q[0], vals[quantity])
        fname = '{}_err_{}.png'.format(cfg['name'], quantity)
        plt.savefig(fname)
