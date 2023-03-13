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

        try:
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
        except:
            print(f'Det {det} of {name} failed, continuing without it.')
            pass

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


"""
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
"""
# candidates = [
#     'GRB121112806',
#     'GRB130624093',
#     'GRB090227772',
#     'GRB130716442'
# ]
candidates = [
    "GRB141205018",
    "GRB110520302",
    "GRB160326062",
    "GRB110818860",
    "GRB171004857",
    "GRB141013803",
    "GRB140115863",
    "GRB130304658",
    "GRB130617564",
    "GRB110207959",
    "GRB110319815",
    "GRB110801335",
    "GRB160827586",
    "GRB120728934",
    "GRB130830921",
    "GRB100819498",
    "GRB170208758",
    "GRB130403866",
    "GRB140817229",
    "GRB170801690",
    "GRB160411062",
    "GRB110625579",
    "GRB140619490",
    "GRB101205309",
    "GRB150423285",
    "GRB150712846",
    "GRB170614505",
    "GRB110702187",
    "GRB180311074",
    "GRB180606730",
    "GRB170718152",
    "GRB090427688",
    "GRB160310291",
    "GRB110509142",
    "GRB171023097",
    "GRB100306199",
    "GRB140701567",
    "GRB100311518",
    "GRB140907429",
    "GRB170621784",
    "GRB120908938",
    "GRB121028280",
    "GRB120412055",
    "GRB151202565",
    "GRB170817529",
    "GRB091018957",
    "GRB100205490",
    "GRB100727238",
    "GRB100204858",
    "GRB150324319",
    "GRB081101167",
    "GRB141102112",
    "GRB100221368",
    "GRB140408553",
    "GRB130626452",
    "GRB120817057",
    "GRB120429003",
    "GRB140106345",
    "GRB100516369",
    "GRB101031625",
    "GRB130908677",
    "GRB110605780",
    "GRB141003788",
    "GRB120618919",
    "GRB121005030",
    "GRB140724533",
    "GRB171025416",
    "GRB120121251",
    "GRB160826938",
    "GRB170326489",
    "GRB160224911",
    "GRB110206202",
    "GRB110726211",
    "GRB160314473",
    "GRB180111815",
    "GRB100616773",
    "GRB170918139",
    "GRB130906222",
    "GRB090411991",
    "GRB160728337",
    "GRB160527080",
    "GRB170916700",
    "GRB131110373",
    "GRB141221338",
    "GRB170307851",
    "GRB130623699",
    "GRB141208038",
    "GRB180410336",
    "GRB150901924",
    "GRB121102064",
    "GRB150605782",
    "GRB100225249",
    "GRB140122597",
    "GRB141213300",
    "GRB170403583",
    "GRB150131951",
    "GRB140224382",
    "GRB100318611",
    "GRB100410356",
    "GRB101213451",
    "GRB160917921",
    "GRB120327418",
    "GRB120312671",
    "GRB180617872",
    "GRB101130074",
    "GRB100406758",
    "GRB100629801",
    "GRB170125022",
    "GRB170101374",
    "GRB140109771",
    "GRB180524920",
    "GRB120109824",
    "GRB100530737",
    "GRB161129300",
    "GRB140516765",
    "GRB100929315",
    "GRB100831651",
    "GRB131108024",
    "GRB100805300",
    "GRB160721806",
    "GRB141005217",
    "GRB110921577",
    "GRB110101202",
    "GRB101008697",
    "GRB160818230",
    "GRB160628136",
    "GRB140626843",
    "GRB180715087",
    "GRB180513815",
    "GRB151022577",
    "GRB120605453",
    "GRB170310883",
    "GRB140501497",
    "GRB130409960",
    "GRB110106893",
    "GRB100714672",
    "GRB160206430",
    "GRB180130744",
    "GRB140311885",
    "GRB161022114",
    "GRB091223511",
    "GRB130215649",
    "GRB160316573",
    "GRB170808065",
    "GRB110722710",
    "GRB150128624",
    "GRB131230529",
    "GRB161212652",
    "GRB170124238",
    "GRB131209963",
    "GRB180416924",
    "GRB141124277",
    "GRB160518039",
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

if BIGFIGS: bigfig.savefig('LLE_Candidates_0310.png', dpi=500)

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
        time_in = -10
        time_out = 100
        data = generate_gbm_lightcurve(name = name, time_in = time_in, time_out = time_out, gbm_detectors=det_list[i])
        
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
if BIGFIGS: bigfig.savefig('GBM_Candidates_0310.png', dpi=500)