#!/usr/bin/env python3

import yaml
import sys
import fprops as fp
import math
import numpy as np
import matplotlib.pyplot as plt

def fmt(x):
    s = f"{x:.1f}"
    if s.endswith("0"):
        s = f"{x:.0f}"
    return rf"{s} \%" if plt.rcParams["text.usetex"] else f"{s} %"


with open(sys.argv[1]) as f:
    doc = yaml.load(f, Loader=yaml.FullLoader)

    p_range = set()
    T_range = set()
    for d in doc['data']:
        p = d['p']
        T = d['T']
        p_range.add(p)
        T_range.add(T)

    p_range = sorted(p_range)
    T_range = sorted(T_range)
    err = np.zeros(shape=(len(T_range), len(p_range)))

    n2 = fp.Nitrogen()
    for d in doc['data']:
        p = d['p']
        T = d['T']
        p_idx = p_range.index(p)
        T_idx = T_range.index(T)

        try:
            state = n2.p_T(p, T)
        except:
            state = None

        if state is not None:
            rho = d['rho']
            if not math.isnan(rho):
                rel_err = (state.rho - rho) / rho
                err[T_idx, p_idx] = rel_err
            else:
                err[T_idx, p_idx] = float("nan")
        else:
            err[T_idx, p_idx] = float("nan")

p_range = np.array(p_range)
T_range = np.array(T_range)
X, Y = np.meshgrid(1e-6 * p_range, T_range)

plt.rcParams['contour.negative_linestyle'] = 'solid'
fig, ax = plt.subplots(1, 1)
cp = ax.contourf(X, Y, err)
# cp = ax.contour(X, Y, err)
# cp = ax.contour(X, Y, err, 3, colors='k')
# manual_locations = [
#     (-1, -1.4), (-0.62, -0.7), (-2, 0.5), (1.7, 1.2), (2.0, 1.4), (2.4, 1.7)]
fig.colorbar(cp)
ax.set_title('Density: Relative error')
ax.set_xlabel('Pressure [MPa]')
ax.set_xscale('log', base=10)
ax.set_ylabel('Temperature [K]')
ax.set_yscale('log', base=10)
# ax.clabel(cp, cp.levels, inline=True, fmt=fmt, fontsize=10, manual=manual_locations)
# ax.clabel(cp, cp.levels, inline=True, fontsize=10, manual=manual_locations)
# ax.clabel(cp, cp.levels, inline=True, fontsize=10)

plt.show()
