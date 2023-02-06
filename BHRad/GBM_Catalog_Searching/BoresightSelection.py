from astropy.table import Table
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

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

    fig = generate_mollweise(non_LAT_sources)
    fig.suptitle('Sources meeting condition THETA > 70 deg - (T90 * 10 deg/min)')
    fig.savefig('THETA>Isotropy.png')
    
    import pdb;pdb.set_trace()
    
if __name__ == "__main__":
    main()
