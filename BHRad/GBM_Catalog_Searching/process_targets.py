import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


"""
Conditions imposed:

    Test statistic > 0
    Flux_ratio must not be infinite or NaN
    GMB-LAT Time Difference >= 0

Noteable caveats:
    Flux ratio & time are values extracted from maximum likelihood fits performed by Fermi.

"""

targets = pd.read_csv('Query-04-22-2022/BrowseTargets.4859.1650668127', skiprows = 4, sep='|')


# First just drop the null rows
del targets['Unnamed: 23']
del targets['Unnamed: 0']

print(targets.keys()+ '\n')

# Parse out all GRBs which have 0 test statistic
print('Dropping %s targets for TS=0'%(len(targets.loc[targets['like_best_ts'] == 0])))
targets = targets.loc[targets['like_best_ts'] != 0]

targets['flux_ratio'] = targets['like_gbm_flux']/targets['like_lat_flux']
targets['flux_ene_ratio'] = targets['like_gbm_flux_ene']/targets['like_lat_flux_ene']

targets['flux_ratio_error'] = targets['like_gbm_flux_error']/targets['like_lat_flux_error']
targets['flux_ene_ratio_error'] = targets['like_gbm_flux_ene_error']/targets['like_lat_flux_ene_error']

# Drop targets with infinite or NaN flux
print('Dropping %s targets for inf flux ratio'% len(targets.loc[targets['flux_ratio'] == np.inf]) )
targets = targets.loc[ targets['flux_ratio'] != np.inf ]
print('Dropping %s targets for inf flux ratio'% len(targets.loc[ np.isnan(targets['flux_ratio']) ]) )
targets = targets.loc[ ~np.isnan(targets['flux_ratio']) ]

# Copy all targets
all_targets = targets.copy()

# Impose that start of GBM is before start of LAT
targets['GBM-LAT Time Difference'] = targets['like_gbm_t0'] - targets['like_lat_t0']
print('Dropping %s targets for GMB-LAT>=0'%(len(targets.loc[targets['GBM-LAT Time Difference'] < 0])))
targets = targets.loc[targets['GBM-LAT Time Difference'] >= 0]


fig, (ax1,ax2) = plt.subplots(1,2)
#ax1.errorbar(x = targets['tl100    '], y = targets['flux_ratio'], yerr = targets['flux_ratio_error'], fmt = 'o', linewidth = 2, capsize = 6, label = 'Photon Flux Ratio')
#ax2.errorbar(x = targets['tl100    '], y = targets['flux_ene_ratio'], yerr = targets['flux_ene_ratio_error'], fmt = 'o', linewidth = 2, capsize = 6, label = 'Energy Flux Ratio', color = 'orange')
ax1.errorbar(x = targets['GBM-LAT Time Difference'], y = targets['flux_ratio'], yerr = targets['flux_ratio_error'], fmt = 'o', linewidth = 2, capsize = 6, label = 'Photon Flux Ratio')
ax2.errorbar(x = targets['GBM-LAT Time Difference'], y = targets['flux_ene_ratio'], yerr = targets['flux_ene_ratio_error'], fmt = 'o', linewidth = 2, capsize = 6, label = 'Energy Flux Ratio', color = 'orange')

ax1.set_xscale('log')
ax2.set_xscale('log')

ax1.legend()
ax2.legend()
#ax1.set_xlabel('Duration (s)')
ax1.set_xlabel('Time difference of GBM - LAT (s)')
ax1.set_ylabel('Ratio (cm**-2 s**-1)')
ax1.set_title('Photon Flux Ratio vs. Time Difference (Like_GBM_t0 - Like_LAT_t0)')
ax2.set_title('Energy Flux Ratio vs. Time Difference (Like_GBM_t0 - Like_LAT_t0)')
fig.set_figwidth(20)
fig.set_figheight(10)
#ax1.set_ylim(0, .02)
#ax1.set_xlim(0,400)
#ax2.set_xlim(0,400)
#ax2.set_ylim(0, .02)
#ax1.set_xlim(0,100)
#ax2.set_xlim(0,100)
ax1.set_yscale('log')
ax2.set_yscale('log')
fig.savefig('targets.png')

targets[['name        ','time                   ','ra        ','dec      ','GBM-LAT Time Difference','flux_ratio','flux_ratio_error','flux_ene_ratio','flux_ene_ratio_error']].to_csv('parsed_candidates.csv')

####

fig, (ax1, ax2) = plt.subplots(1,2)
ax1.errorbar(x = targets['lle_t90 '], y=targets['flux_ratio'], yerr = targets['flux_ratio_error'], fmt = 'o', linewidth = 2, capsize = 6)
ax1.set_xlabel('Burst duration llet90 (window of 90% flux accumulation)')
ax1.set_ylabel('Flux ratio: hardness in GBM (10 keV - 25 MeV) / LAT (20 MeV - 300 MeV)')
fig.set_figwidth(20)
fig.set_figheight(10)
ax1.set_title('Burst duration vs. Flux ratio (hardness in GBM/LAT) with main targets')

ax2.set_title('Same, with ALL sources')
ax2.errorbar(x = all_targets['lle_t90 '], y=all_targets['flux_ratio'], yerr = all_targets['flux_ratio_error'], fmt = 'o', linewidth = 2, capsize = 6)

fig.savefig('Burst_duration_vs_flux_ratio_hardness.png')


"""
# Let's parse out the targets in window.
def inWindow(flux, error, threshold=50):
    # Target bounds straddle lower bound
    if flux-error<-threshold and flux+error>-threshold:
        return True
    # Target bounds straddle upper bound
    elif flux+error>threshold and flux-error<threshold:
        return True
    # Target fully contained by bounds
    elif flux+error<threshold and flux-error>threshold:
        return True
    # Target has no datapoints within bounds
    else:
        return False

import pdb;pdb.set_trace();

targets = targets[np.where(inWindow(targets['flux_ratio'], targets['flux_ene_ratio']) == True)]

targets.to_csv('targets_in_window.csv')
"""
