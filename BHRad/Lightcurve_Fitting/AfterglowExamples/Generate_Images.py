from threeML import *
import ultranest
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, norm

#GBM Energy Spectrum
Emin = 50*10**(-6) # GeV
Emax = 300*10**(-6) # GeV

gbm_catalog = FermiGBMBurstCatalog()
def query_catalog(GRBNAME = '150902733', time_in = -30, time_out = 60, dt = .1, gbm_detectors = None):
    gbm_catalog.query_sources('GRB%s'%GRBNAME)

    grb_info = gbm_catalog.get_detector_information()['GRB%s'%GRBNAME]
    
    if gbm_detectors is None:
        gbm_detectors = grb_info["detectors"]
    print(gbm_detectors)
    source_interval = grb_info["source"]["fluence"]
    background_interval = grb_info["background"]["full"]
    dload = download_GBM_trigger_data('bn%s'%GRBNAME, detectors=gbm_detectors)

    fluence_plugins = []
    time_series = {}
    figures = []
    for det in gbm_detectors:

        ts_cspec = TimeSeriesBuilder.from_gbm_cspec_or_ctime(
            det, cspec_or_ctime_file=dload[det]["cspec"], rsp_file=dload[det]["rsp"]
        )

        ts_cspec.set_background_interval(*background_interval.split(","))
        ts_cspec.save_background(f"{det}_bkg.h5", overwrite=True)

        ts_tte = TimeSeriesBuilder.from_gbm_tte(
            det,
            tte_file=dload[det]["tte"],
            rsp_file=dload[det]["rsp"],
            restore_background=f"{det}_bkg.h5",
        )

        time_series[det] = ts_tte

        ts_tte.set_active_time_interval(source_interval)

        figures.append(ts_tte.view_lightcurve(time_in, time_out, dt = dt))
    return figures

def get_curves(fig):
    data = fig.get_axes()[0].get_lines()[0].get_xydata()
    x = data[:,0]
    curve = data[:,1]
    
    bckg_data = fig.get_axes()[0].get_lines()[1].get_xydata()
    bckg_data = bckg_data[:,1]
    
    return x, curve, bckg_data

def retrieve_data(figures):
    y = []
    for fig in figures:
        x, curve, bckg_data = get_curves(fig)

        # Get the background
        #bckg_data = fig.get_axes()[0].get_lines()[1].get_xydata()
        #bckg_data = bckg_data[:,1]

        y.append(curve-np.mean(bckg_data))
    return (y[0]+y[1]+y[2]+y[3])/4, x

def show_curves(figures):
    figya, ax = plt.subplots()
    for i,fig in enumerate(figures):
        x,curve, bckg = get_curves(fig)
        ax.plot(x,curve,label="detector %s"%i)
        ax.plot(x[0:-1],bckg, '--',label='bckg %s'%i)
        ax.legend(prop={'size':15})
    return figya

def prior_transform(cube):
    params = cube.copy()

    # K_lightcurve: uniform
    lo = 1e3#1e-26
    hi = 1e6#1e-24
    params[0] = cube[0] * (hi - lo) + lo

    # K_powerlaw: uniform
    lo = 1
    hi = 1e4
    params[1] = cube[1] * (hi - lo) + lo
    
    # t_m_powerlaw: uniform
    lo = 1e-3
    hi = 1e3
    params[2] = cube[2] * (hi - lo) + lo
    
    # t_p_powerlaw: uniform
    lo = 1e-3
    hi = 1e3 # t_p < t_m
    params[3] = cube[3] * (hi - lo) + lo
    
    # delta_powerlaw: uniform
    lo = -20
    hi = 60
    params[4] = cube[4] * (hi - lo) + lo
    
    # index_lightcurve: gaussian on -.52 
    params[5] = norm.ppf(cube[5], -.52, 1)
    
    return params

# Simplifying lightcurve with powerlaw of index -.52
def lightcurve(tau, energy_range, normalization = 1., index = .52):
    return normalization*tau**(-index)

def afterglow(tau, delta, t_m, t_p, normalization):
    # lifetime only, no spectrum.
    # Norm * \frac{e^{(-t -\Delta)/tp } }{1 - e^{(-t - \Delta)/tm} }
    return normalization*np.exp(-(1/t_p)*(tau-delta))/(1 + np.exp(-(1/t_m)*(tau-delta)))

def lightcurve_model(K_lightcurve, K_powerlaw, t_m_powerlaw, t_p_powerlaw, 
                     delta_powerlaw, index_lightcurve,
                     energy_range = np.linspace(Emin, Emax, 1000)):
    
    # Primary/secondary emission interval
    """
    tau = np.linspace(0.5, peak_index - lightcurve_open -.5, peak_index - lightcurve_open)
    """
    tau = np.linspace(peak_index, 0, peak_index, endpoint=False)
    
    curve = lightcurve(tau=tau, energy_range=energy_range, normalization=K_lightcurve, index = index_lightcurve)

    #curve = np.flip(curve)
    
    # Add afterglow element
    afterglow_curve = afterglow(tau=time_domain, delta = delta_powerlaw, t_p = t_p_powerlaw, t_m = t_m_powerlaw, normalization=K_powerlaw)
    # Assemble model
    model = np.zeros(len(data))
    model[0:peak_index]=curve
    #model[lightcurve_open:peak_index] = curve
    model = model + afterglow_curve

    return model

def likelihood_model(params):
    K_lightcurve, K_powerlaw, t_m_powerlaw, t_p_powerlaw, delta_powerlaw, index_lightcurve = params
    
    curve = lightcurve_model(K_lightcurve, K_powerlaw, t_m_powerlaw, t_p_powerlaw, delta_powerlaw, index_lightcurve)
    
    like = -0.5 * (((curve - data)/error_tolerance)**2).sum()
    #like = pearsonr(data,curve)[0]
    return like

def plot_model(x, result = None, log=False, name='bn150902733', vals=None):
    
    fig, ax1 = plt.subplots()#(ax1, ax2) = plt.subplots(1,2)
    fig.set_figwidth(10)
    fig.set_figheight(8)
    
    if vals is None:
        vals = result['posterior']['mean']
    K_lightcurve = vals[0]
    K_powerlaw = vals[1]
    t_m_powerlaw = vals[2]
    t_p_powerlaw = vals[3]
    delta_powerlaw = vals[4]
    index_lightcurve = vals[5]
    
    ax1.scatter(x,data, s=30, label = 'GBM Photon Counts')
    # ax2.scatter(x,data, s=30)
    
    fitted_model = lightcurve_model(K_lightcurve, K_powerlaw, t_m_powerlaw, t_p_powerlaw, delta_powerlaw, index_lightcurve)
    r = pearsonr(data,fitted_model)
    ax1.plot(x, 
            fitted_model,
            label = 'Powerlaw + Afterglow fitted model', color = 'C1')#'r:%.4e  pval:%.4e'%(r[0],r[1]))

    curve = lightcurve(tau=np.linspace(peak_index, 0, peak_index, endpoint=False), 
                          energy_range=np.linspace(Emin, Emax, 1000), 
                          normalization=K_lightcurve)
    #curve = np.flip(curve)
    # ax2.plot(x[0:peak_index], curve, label = 'simplified direct+frag emission')

    afterglow_curve = afterglow(tau=time_domain, delta = delta_powerlaw, t_p = t_p_powerlaw, t_m = t_m_powerlaw, normalization=K_powerlaw)
    # ax2.plot(x, afterglow_curve, '--', label = r'afterglow $\frac{e^{Ax}}{1-e^{Bx}}$', color='orange')

    #plt.ylim(100,200)
    #ax2.ylim(-3500,6500)
    ax1.legend(prop={'size': 13})
    # ax2.legend(prop={'size': 15})
    ax1.set_ylabel('Count rate (/s binned on 100ms)')
    ax1.set_xlabel('Time wrt to trigger time')
    #ax2.set_ylabel('Count rate')
    # ax2.set_xlabel('Time wrt to trigger time')
    fig.suptitle('Composite Bayesian fit of %s'%name, size=15)
    
    # std = result['posterior']['stdev']
    # lc_txt="Lightcurve norm: %.4e +- %.2e & index: %.4f +- %.2f ; \n"%(vals[0], std[0], vals[5], std[5])
    # pl_txt="Powerlaw norm: %.4e +- %.2e , t_m: %.4f +- %.2f , t_p: %.4f +- %.2f , delta: %.4f +- %.2f \n"%( 
    #     vals[1], std[1], vals[2], std[2], vals[3], std[3], vals[4], std[4] )
    # fit_txt="logZ: %.4e +- %.2e (Z = marginal likelihood)"%(result['logz'], result['logzerr'])
    # plt.figtext(0.5, 0.01, lc_txt + pl_txt + fit_txt, wrap=True, horizontalalignment='center', fontsize=13)
    
    if log:
        ax1.set_yscale('log')
        # ax2.set_yscale('log')
    
    return fig


figures = query_catalog(GRBNAME = '150902733', time_in = -30, time_out = 60, dt=.1)

data, time_domain = retrieve_data(figures)
fig = show_curves(figures)
fig.show()

bckg_data=[]
for fig in figures:
    _,_, b = get_curves(fig)
    bckg_data.append(b)
error_tolerance = np.std(bckg_data)

peak_index = np.where(data == max(data))[0][0]
print('Peak index at %s'%peak_index)

param_names = ['K_lightcurve', 'K_powerlaw', 't_m_powerlaw', 't_p_powerlaw', 
               'delta_powerlaw', 'index_lightcurve'] # Normalization & steepness for exp. decay

sampler = ultranest.ReactiveNestedSampler(param_names, likelihood_model, prior_transform)

# result = sampler.run()
# sampler.print_results()

import pdb;pdb.set_trace()

# fig = plot_model(time_domain, result = result, log=False, name = 'bn150902733')
fig = plot_model(time_domain, log=False, name = 'bn150902733', \
    vals = [3.38e3, 7.60e3, 1.27, 2.33, 9.64, 1.49] )

fig.savefig('bn150902733_single.png', dpi=300)

"""
figures = query_catalog(GRBNAME = '140206275', time_in = -20, time_out = 100, dt=.1)

data, time_domain = retrieve_data(figures)
fig = show_curves(figures)
fig.show()

bckg_data=[]
for fig in figures:
    _,_, b = get_curves(fig)
    bckg_data.append(b)
error_tolerance = np.std(bckg_data)

peak_index = np.where(data == max(data))[0][0]
print('Peak index at %s'%peak_index)

param_names = ['K_lightcurve', 'K_powerlaw', 't_m_powerlaw', 't_p_powerlaw', 
               'delta_powerlaw', 'index_lightcurve'] # Normalization & steepness for exp. decay

# sampler = ultranest.ReactiveNestedSampler(param_names, likelihood_model, prior_transform_14)

# result = sampler.run()
# sampler.print_results()

import pdb;pdb.set_trace()

# fig = plot_model(time_domain, result = result, log=False, name = 'bn140206275')
# fig = plot_model(time_domain, log=False, name = 'bn150902733', \
#     vals = [3.38e3, 7.60e3, 1.27, 2.33, 9.64, 1.49] )

fig, ax = plt.subplots()
fig.set_figwidth(10)
fig.set_figheight(8)
ax.scatter(time_domain, data, s=30, label = 'GBM Photon Counts')

ax.set_ylabel('Count rate (/s binned on 100ms)')
ax.set_xlabel('Time wrt to trigger time')
ax.legend(prop={'size': 15})

fig.suptitle('Composite countrate of %s'%'bn140206275', size=15)


import pdb;pdb.set_trace()

fig.savefig('bn140206275_countrate.png', dpi=300)
"""