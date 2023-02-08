from astropy.table import Table
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Catalogs download https://heasarc.gsfc.nasa.gov/cgi-bin/W3Browse/w3catindex.pl

def generate_mollweise(targets, cmap_var = 'T90', log_cmap = False, calculate_isotropy=True):

    fig = plt.figure(figsize=(12,6))
    ax = fig.add_subplot(111, projection='mollweide')
    plt.grid(True)

    # change the >180s to negative
    modified = targets.copy()
    modified.loc[modified['LII_x']>180, 'LII_x'] -= (360)
    
    if log_cmap is True:
        modified[cmap_var] = np.log(modified[cmap_var])
    sc = ax.scatter(modified['LII_x']*(np.pi/180), targets['BII_x']*(np.pi/180), 50, 
                c=modified[cmap_var], cmap='Greens', edgecolors='lightgray')
    plt.colorbar(sc, label = 'T90')

    ax.set_ylabel('Gal. Longitude')
    ax.set_xlabel('Gal. Latitude')

    caption = 'In Gal. Lat/Long coords, the galactic center is at (0,0) and the plane is on the x-axis.'
    
    if calculate_isotropy:
        center = modified.query('LII_x > -90 and LII_x < 90').shape[0]
        anticenter = modified.query('LII_x > 90 or LII_x < -90').shape[0]
        caption += '\ncenter/anticenter numerical ratio of %.2f'%(center/anticenter)
        northern = modified.query('BII_x > 0 and BII_x < 90').shape[0]
        southern = modified.query('BII_x < 0 and BII_x > -90').shape[0]
        caption += '\nnothern/southern num. ratio of %.2f'%(northern/southern)
    
    plt.figtext(0.5, 0.01, 
                caption, 
                horizontalalignment='center', fontsize=10)
    return fig

def propagate_error(a, b, e_a, e_b):
    return np.abs(a/b) * np.sqrt( (e_a/a)**2 + (e_b/b) )

def main():

    # Load datasets
    GBMTriggerCatalog = Table.read('GBMTriggerCatalog.fits').to_pandas()
    GBMBurstCatalog = Table.read('GBMBurstCatalog.fits').to_pandas()
    LLECatalog = Table.read('LLECatalog.fits').to_pandas()

    # First merge GBM Trigger/Burst catalogs
    GBM = GBMBurstCatalog.merge(
            right = GBMTriggerCatalog, 
            on='TRIGGER_NAME',
            how = 'left'
            )

    # Theta > 70 + signal_time * 10deg/min
    # 10deg/min * 1min/60s = (1/6)
    GBM['condition'] = (GBM['THETA'] > (70 + GBM['T90']*(1/6)))

    # Let's also quickly calculate T90/50
    GBM['T9050'] = GBM['T90']/GBM['T50']
    GBM['T9050_ERROR'] = propagate_error(GBM['T90'], GBM['T50'], GBM['T90_ERROR'], GBM['T50_ERROR'])
    
    # See how many meet/fail condition
    print(GBM.condition.value_counts())

    # Select those that meet condition
    non_LAT_sources = GBM.query('condition == True')

    names = non_LAT_sources['TRIGGER_NAME'].to_list()

    # Pick out sources still cross-listed
    cross_listed = LLECatalog.query('TRIGGER_NAME in @names').TRIGGER_NAME.to_list()

    # Show their boresight angles
    cross_listed = non_LAT_sources.query('TRIGGER_NAME in @cross_listed')[['TRIGGER_NAME', 'THETA', 'T90']]
    print(cross_listed)

    # THETA vs T90
    fig, ax = plt.subplots()
    ax.scatter(non_LAT_sources.THETA, non_LAT_sources.T90, marker='x', color = 'C0')
    ax.scatter(cross_listed.THETA, cross_listed.T90, marker='x', color='darkred', zorder=1, label = 'also listed in LLE')
    ax.set_xlabel('THETA (deg)')
    ax.set_ylabel('T90 (s)')
    ax.set_title('Sources meeting condition THETA > 70 deg - (T90 * 10 deg/min)')
    ax.set_yscale('log')
    #ax.set_xscale('log')
    ax.legend()
    fig.savefig('THETAvsT90.png')

    # Angular distribution
    fig = generate_mollweise(non_LAT_sources)
    fig.suptitle('Sources meeting condition THETA > 70 deg - (T90 * 10 deg/min)')
    fig.savefig('THETA>Isotropy.png')
    
    # Constraining inferred plaw index
    fig, ax = plt.subplots(figsize = (10,7))

    # Pull out sources within bounds
    sources_IN_index_range = non_LAT_sources.query(' (T9050 + T9050_ERROR > 2.4) & (T9050 - T9050_ERROR < 2.6) ')
    selected_names = sources_IN_index_range['TRIGGER_NAME'].to_list()
    ax.scatter(sources_IN_index_range['T9050'], sources_IN_index_range['FLUENCE'], marker='x', s=30, zorder=3, color='darkred', alpha=.7,
               label = r'Source T90/50 $\in$ [2.4, 2.6]')

    # Sources not within bounds
    sources_NOT_in_index_range = non_LAT_sources.query('TRIGGER_NAME not in @selected_names')
    ax.scatter(sources_NOT_in_index_range['T9050'], sources_NOT_in_index_range['FLUENCE'], marker='x', s=30, zorder=3, color='C0', alpha=.3,
               label = r'T90/50 $\notin$')

    ax.set_xscale('log')
    ax.set_yscale('log')
    xl = ax.get_xlim()
    yl = ax.get_ylim()

    ax.errorbar(non_LAT_sources['T9050'], non_LAT_sources['FLUENCE'], xerr = non_LAT_sources['T9050_ERROR'], yerr = non_LAT_sources['FLUENCE_ERROR'],
                    fmt='none', linewidth=1.5, capsize=3, color='black', markersize=3, alpha=.2, zorder=1
                )

    ax.axvspan(2.4, 2.6, facecolor='gold', alpha=.1, zorder=2)
    ax.axvline(2.4, color = 'gold', alpha = .3, zorder = 2)
    ax.axvline(2.6, color = 'gold', alpha = .3, zorder = 2)

    # Add marginal histogram
    divider = make_axes_locatable(ax)
    ax_histy = divider.append_axes("right", 1.2, pad=0.1, sharey=ax)
    ax_histy.yaxis.set_tick_params(labelleft=False)
    ax_histy.hist( [sources_NOT_in_index_range.FLUENCE.to_list(), sources_IN_index_range.FLUENCE.to_list()] , bins = np.logspace(-8 -3, 12),
                  fill=False, histtype='step', stacked=False, orientation='horizontal', color = ['C0', 'darkred'])
    ax_histy.set_yscale('log')
    ax_histy.set_ylim(yl)
    ax_histy.set_xscale('log')
    ax_histy.set_xlabel('Counts')
       
    ax.set_xlim(xl)
    ax.set_ylim(yl)
    ax.set_ylabel('GBM Fluence (erg cm-2)')
    ax.set_xlabel('T90/50 ratio (unitless)')
    ax.legend()

    fig.savefig('FLUENCEvsINDEX.png')

    import pdb;pdb.set_trace()

    
    
    
if __name__ == "__main__":
    main()
