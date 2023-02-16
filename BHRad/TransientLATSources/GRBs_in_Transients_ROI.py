import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def generate_mollweise(targets, cmap_var = 'T90', log_cmap = False):

    fig = plt.figure(figsize=(12,6))
    ax = fig.add_subplot(221, projection='mollweide')
    plt.grid(True)

    # change the >180s to negative
    modified = targets.copy()
    modified.loc[modified['LII_x']>180, 'LII_x'] -= (360)
    
    if log_cmap is True:
        modified[cmap_var] = np.log(modified[cmap_var])
    sc = ax.scatter(modified['LII_x']*(np.pi/180), targets['BII_x']*(np.pi/180), 50, 
                c=modified[cmap_var], cmap='Greens', edgecolors='lightgray', zorder=1)
    plt.colorbar(sc, label = 'T90')

    ax.set_ylabel('Gal. Longitude')
    ax.set_xlabel('Gal. Latitude')

    caption = 'In Gal. Lat/Long coords, the galactic center is at (0,0) and the plane is on the x-axis.'
    
    return fig, ax

# Shouldn't need to use this again
def calculate_proximity(grbs, transients):
    # Adds 36 columns to grbs which is kind of funny, but I don't see any more efficient way to do this
    for i in range(transients.shape[0]):

        proximity = []
        #time_diff = []

        for j in range(grbs.shape[0]):

            GRB_RA = grbs.iloc[j]['RA_x']
            GRB_DEC = grbs.iloc[j]['DEC_x'] 
            T_RA = transients.iloc[i]['RA']
            T_DEC = transients.iloc[i]['DEC']

            # Calculate the degrees of proximity to the Transient source for each GRB
            proximity.append( np.sqrt( (T_RA - GRB_RA)**2 + (T_DEC - GRB_DEC)**2 ) )
            #time_diff.append()


        # Place a col for each transient 
        col_name = 'deg_prox_%s'%transients.iloc[i]['SRCNumber']
        grbs[ col_name ] = proximity
    return grbs





####

grbs = pd.read_csv('./../GBM_Catalog_Searching/GBM_BurstTrig_merge.csv')

# Pick sources that have fitted filename & fitted params
transients = pd.read_csv('./TransientSources_fitted_params.csv').dropna()

# You cannot fathom how cursed this is
transients_in_deg = pd.read_csv( './UnassocTransientsCatalogInDeg.csv' )
transients_in_deg.columns = transients_in_deg.columns.str.strip(' ') # Clean columns
transients_in_deg['SRCNumber'] = transients_in_deg['SRCNumber'].str.strip(' ').str.strip('"').astype('int') # Clean SRCNumber
# Clean up, RA DEC
transients_in_deg['RA'] = transients_in_deg['RA (J2000.0)'].str.replace(' ','').str.strip('"').astype('float')
transients_in_deg['DEC'] = transients_in_deg['Dec (J2000.0)'].str.replace(' ','').str.strip('"').astype('float')

# Merge in the RA, DEC in degrees
transients = transients.merge( transients_in_deg , how='left', on = 'SRCNumber')

# proximities = calculate_proximity(grbs, transients)
# proximities.to_csv('GRBs_Transients_proximity.csv')
proximities = pd.read_csv('GRBs_Transients_proximity.csv')

####



def generate_transient_ROI_cut(transient:pd.Series, proximities:pd.DataFrame, ROIcut = 5):


    targets = proximities.query(' deg_prox_%s < @ROIcut'%transient['SRCNumber'])

    # Generate ROI mollweise
    fig, ax = generate_mollweise(targets = targets)

    # Amend mollweise with transient location
    if transient['LII']>180:
        LII = transient['LII'] - 360
    else:
        LII = transient['LII']
    ax.scatter(LII, transient['BII'], color = 'darkred', marker='x', zorder = 3, label = 'Transient location', s=50)
    ax.legend()

    ax = fig.add_subplot(222)
    ax.errorbar(x = targets['deg_prox_%s'%transient['SRCNumber']], y = targets['T9050'], yerr = targets['T9050_ERROR'],
        fmt='o', linewidth=1.5, capsize=3, color='green', markersize=5)
    ax.set_xlabel('Proximity to transient (deg)')
    ax.set_ylabel(r'$\frac{t_{90}}{t_{50}}$')
    ax.set_xlim(0,5)
    ax.axhspan(2.4, 2.6, color='gold', alpha=.5)

    ax = fig.add_subplot(223)
    ax.errorbar(x = targets['deg_prox_%s'%transient['SRCNumber']], y = targets['T90'], yerr = targets['T90_ERROR'],
        fmt='o', markersize=5, color = 'green', capsize=3, linewidth=1.5)
    ax.set_xlabel('Proximity to transient (deg)')
    ax.set_ylabel('T90 (s)')
    ax.set_xlim(0,5)

    ax = fig.add_subplot(224)
    # There's also fluxes 64, 256, 1024 --> decreasing flux
    ax.scatter(x = targets['TIME'], y = targets['FLUX_64'], s=5, color = 'green')
    ax.set_xlabel('Time in MJD')
    ax.set_ylabel('Flux (values are incorrect)')

    # Now plot transient data

    d = pd.read_csv(transient['filename'], header=None)
    # Reformat
    x = []
    y = []
    yerr = []
    # For each row
    for i in range(d.shape[0]):
        # Odds are data, evens are upper/lower error bound
        if i % 2:
            yerr.append(d.iloc[i,1])
        else:
            x.append(d.iloc[i,0])
            y.append(d.iloc[i,1])

    # Make error always positive
    yerr = np.abs(np.array(y) - np.array(yerr))

    # Plot transients data
    ax.errorbar(x, y, yerr = yerr, fmt='o', color = 'darkred')

    # Plot model #TODO
    #ax[i].plot(xaxis, y_model, '-', label = 'd = %.2epc\nT0 = %.2es\nMSE %.2f'%(distance_A(10**coefs[0]), 10**coefs[1], mean_squared_error(y, model(x, *coefs))))

    return fig, ax 

fig, ax = generate_transient_ROI_cut(transients.iloc[0], proximities)




import pdb;pdb.set_trace()