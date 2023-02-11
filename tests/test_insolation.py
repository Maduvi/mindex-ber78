#!/usr/bin/env python3
"""

Test plotting values.

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# read monsindex results 65N
dat0 = pd.read_csv('test_insolation_N65.res', sep='\s+')
tim0 = dat0.year
sas0 = dat0.sastro
was0 = dat0.wastro
sso0 = dat0.ssols
wso0 = dat0.wsols

# read monsindex results 65S
dat1 = pd.read_csv('test_insolation_S65.res', sep='\s+')
tim1 = dat1.year
sas1 = dat1.sastro
was1 = dat1.wastro
sso1 = dat1.ssols
wso1 = dat1.wsols

# read palinsol data
pali = pd.read_csv('data/palinsol.csv')
tim2 = pali.time
sas2 = pali.N65sum
was2 = pali.N65win
sso2 = pali.N65sso
wso2 = pali.N65wso
sas3 = pali.S65sum
was3 = pali.S65win
sso3 = pali.S65sso
wso3 = pali.S65wso

fig, ax = plt.subplots(8, 1, figsize=(210 / 25.4, 580 / 25.4),
                       layout="constrained")

for xx in ax:
    xx.set(xlim=[-500e3, 0])
    xx.grid(True, linestyle='--', color='lightgray', lw=0.75)
    xx.tick_params(axis='x',label1On=False)

ax[0].plot(tim0, sas0, lw=1, clip_on=False, label='monsindex')
ax[0].plot(tim2, sas2, lw=1, clip_on=False, label='palinsol')
ax[0].set_ylabel('65N astrosum (Wm-2)')
ax[0].legend(handlelength=0.9)

ax[1].plot(tim0, was0, lw=1, clip_on=False)
ax[1].plot(tim2, was2, lw=1, clip_on=False)
ax[1].set_ylabel('65N astrowin (Wm-2)')

ax[2].plot(tim0, sso0, lw=1, clip_on=False)
ax[2].plot(tim2, sso2, lw=1, clip_on=False)
ax[2].set_ylabel('65N sumsolst (Wm-2)')

ax[3].plot(tim0, wso0, lw=1, clip_on=False)
ax[3].plot(tim2, wso2, lw=1, clip_on=False)
ax[3].set_ylabel('65N winsolst (Wm-2)')

ax[4].plot(tim1, sas1, lw=1, clip_on=False, label='monsindex')
ax[4].plot(tim2, sas3, lw=1, clip_on=False, label='palinsol')
ax[4].set_ylabel('65S astrosum (Wm-2)')

ax[5].plot(tim1, was1, lw=1, clip_on=False)
ax[5].plot(tim2, was3, lw=1, clip_on=False)
ax[5].set_ylabel('65S astrowin (Wm-2)')

ax[6].plot(tim1, sso1, lw=1, clip_on=False)
ax[6].plot(tim2, wso3, lw=1, clip_on=False) # JJA summer is winter SH
ax[6].set_ylabel('65S sumsolst (Wm-2)')

ax[7].plot(tim1, wso1, lw=1, clip_on=False)
ax[7].plot(tim2, sso3, lw=1, clip_on=False) # JJA summer is winter SH
ax[7].set_ylabel('65S winsolst (Wm-2)')
ax[7].tick_params(axis='x',label1On=True)
ax[7].set_xlabel('Time')

plt.savefig('./test_insolation.pdf')
