#!/usr/bin/env python3
"""

Test plotting values.

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


def mindex_table_vernekar(ifile):
    """Conver table to mindex."""
    # values from tables (from present towrds past)
    v25, v20, v0 = np.loadtxt(ifile, skiprows=1, unpack=True)
    nyear = v25.size

    # invert: from past towrds present
    v25 = v25[::-1]
    v20 = v20[::-1]
    v0 = v0[::-1]

    # get absolute values by adding PI values
    for i in range(0, nyear - 1):
        v25[i] += v25[-1]
        v20[i] += v20[-1]
        v0[i] += v0[-1]

    # interpolate to get tropic of cancer
    lin = interp1d([25, 20], np.vstack([v25, v20]), axis=0)
    v23 = lin(23.45)

    # monsoon index gradient
    grad = 2 * v23 - v0
    grad0 = grad[-1]

    return grad - grad0


# read program output
data = pd.read_csv('test_mindex.res', sep='\s+')
time = data.year
mind = data.mindex

# theory data
theory = mindex_table_vernekar('data/vernekar72.dat')

fig, ax = plt.subplots(1, 1, figsize=(8, 4), layout='constrained')

ax.set(xlim=[-500e3, 0], xlabel='Time', ylabel='Monsoon index (Lang d-1)')
ax.grid(True, linestyle='--', color='lightgray', lw=0.75)

ax.axhline(0, linestyle='--', color='k', lw=0.75)

ax.plot(time, mind, lw=1, color='firebrick', label='Monsindex')
ax.plot(time, theory, lw=1, color='goldenrod', label='Table V72')

ax.legend(handlelength=0.9)

plt.savefig('test_mindex.pdf')
