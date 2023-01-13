from threeML import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# lat_catalog = FermiLATSourceCatalog()

# def generate_lat_lightcurve(name = "GRB141222298", time_in = -30, time_out = 30, dt=.1, returnfig = False):
#     """
#     """
#     lat_catalog.query_sources(name) # form 'nFGL J0000.0+0000' 

#     dload = download_LAT_data(name)

#     ts_lat = TimeSeriesBuilder

"""
As I understand it, every GRB event has an associated LLE event created from what LAT notices abt the event. This LLE isn't 
"""

lle_catalog = FermiLLEBurstCatalog()

def generate_lle_lightcurve(name = "GRB141222298", time_in = -30, time_out = 30, dt=.1, returnfig = False):
    """
    Parameters: name, time_in, time_out, dt, 
    returnfig = True to return output as tuple with figure

    Returns a dataframe with binning and counts;
    figure optionally.
    """
    lle_catalog.query_sources(name)

    dload = download_LLE_trigger_data(name)

    ts_lle = TimeSeriesBuilder.from_lat_lle(
        "GRB141222298",
        lle_file=dload["lle"],
        ft2_file=dload["ft2"],
        rsp_file=dload['rsp']
    )

    ts_lle.create_time_bins(time_in, time_out, dt=dt)

    data = ts_lle.bins._create_pandas()
    data['Counts'] = ts_lle.total_counts_per_interval

    if returnfig is True:
        return data, ts_lle.view_lightcurve(time_in, time_out, dt=dt)
    else:
        return data

gbm_catalog = FermiGBMBurstCatalog()

def generate_gbm_lightcurve(name = "GRB141222298", time_in = -30, time_out = 30, dt=.1, gbm_detectors = None, returnfig = False):
    """
    Parameters: name, time_in, time_out, dt, 
    detectors = None for default mask; list for specific detectors
    returnfig = True to return output as tuple with figure

    Returns a dictionary with dataframes with binning and counts for each detector;
    figures optionally.
    """
    gbm_catalog.query_sources(name)

    grb_info = gbm_catalog.get_detector_information()[name]
    background_interval = grb_info["background"]["full"]

    if gbm_detectors is None:
        gbm_detectors = grb_info["detectors"]
    print(gbm_detectors)

    dload = download_GBM_trigger_data('bn%s'%name[3:], detectors=gbm_detectors)

    if returnfig is True:
        figures = []

    time_series = {}
    for det in gbm_detectors:

        # Calculate background
        ts_cspec = TimeSeriesBuilder.from_gbm_cspec_or_ctime(
            det, cspec_or_ctime_file=dload[det]["cspec"], rsp_file=dload[det]["rsp"]
        )
        ts_cspec.set_background_interval(*background_interval.split(","))
        ts_cspec.save_background(f"{det}_bkg.h5", overwrite=True)

        # Build gbm timeseries
        ts_tte = TimeSeriesBuilder.from_gbm_tte(
            det,
            tte_file=dload[det]["tte"],
            rsp_file=dload[det]["rsp"],
            restore_background=f"{det}_bkg.h5",
        )

        # Bin data
        ts_tte.create_time_bins(time_in, time_out, dt=dt)

        # Generate dataframe
        data = ts_tte.bins._create_pandas()
        data['Counts'] = ts_tte.total_counts_per_interval

        time_series[det] = data

        if returnfig is True:
            figures.append(ts_tte.view_lightcurve(time_in, time_out, dt = dt))

    if returnfig is True:
        return time_series, figures
    else:
        return time_series

def plot_dets(ax, data):
    all = []
    for det in list(data.keys()):
        ax.step(data[det]['Start'], data[det]['Counts'], alpha=.5, label = det)
        all.append(data[det]['Counts'])

    mean = np.mean(all, axis=0)
    ax.plot(data[det]['Start'], mean, label = 'Mean', color='black')

    return pd.DataFrame({'Start':data[det]['Start'], 'AvgCounts':mean})



candidates = [
    #'GRB190731943',
    #'GRB180703876',
    #'GRB220408311',
    #'GRB150523396',
    'GRB190531840',
    'GRB181028590',
    'GRB210928084',
    #'GRB140402007',
    #'GRB180210517',
    'GRB160821857',
    #'GRB120915000'
]

BIGFIGS = False


# For LLE
names = []
for name in candidates:
    try:
        lle_catalog.query_sources(name)
        names.append(name)
    except:
        print('Could not retrieve LLE for %s'%name)
        pass
print(names)

if BIGFIGS: bigfig, axes = plt.subplots(2,2, figsize = (20,10))

if BIGFIGS: axlist = axes.flatten()
# import pdb;pdb.set_trace()
for i,name in enumerate(names):
    fig, ax = plt.subplots(figsize = (10,8))

    data = generate_lle_lightcurve(name = name, time_in = -10, time_out = 100, dt=1)

    if BIGFIGS: 
        axlist[i].set_title(name)
        axlist[i].step(data['Start'], data['Counts'], alpha=.5, color = 'darkred')

    ax.set_title(name)
    ax.step(data['Start'], data['Counts'], alpha=.5, color = 'darkred')
        
    fig.savefig('%s.png'%name, dpi=300)

    # data.to_csv('./Data/LAT_%s'%name)

if BIGFIGS: bigfig.savefig('Candidates_LLE.png', dpi=500)

names = []
det_list = []
for name in candidates:
    try:
        gbm_catalog.query_sources(name)

        grb_info = gbm_catalog.get_detector_information()[name]

        if len(grb_info["detectors"]) == 0:
            # There are no dets
            det_list.append(['n%s'%i for i in range(1,10)] + ['b0', 'b1', 'na', 'nb'])
        else:
            det_list.append(None)

        names.append(name)
    except:
        pass

print(names)
#import pdb;pdb.set_trace()

if BIGFIGS: bigfig, axes = plt.subplots(3,2, figsize=(20,12))

if BIGFIGS: axlist = axes.flatten()
for i,name in enumerate(names):

    fig, ax = plt.subplots(figsize=(10,8))

    try:
        data = generate_gbm_lightcurve(name = name, time_in = -10, time_out = 100, gbm_detectors=det_list[i])
        
        if BIGFIGS: 
            axlist[i].set_title(name)
            mean = plot_dets(axlist[i], data)
            axlist[i].legend()

        ax.set_title(name)
        mean = plot_dets(ax, data)
        ax.legend()

        fig.savefig('%s_GBM.png'%name, dpi=300)
    except:
        print('Failed with GBM %s'%name)
        pass

    

    #mean.to_csv('./Data/GBM_%s'%name)
if BIGFIGS: bigfig.savefig('Candidates_.png', dpi=500)