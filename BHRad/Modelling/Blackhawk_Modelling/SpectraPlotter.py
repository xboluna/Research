import pandas as pd
from matplotlib import pyplot as plt
from numpy import linspace
import numpy as np
import astropy.units as u
from scipy import integrate

try:
    from ebltable.tau_from_model import OptDepth as OD
    ebl = True
    print('ebltable initiated')
except:
    ebl = "Cannot initiate ebltable -- package missing? \n(pip install ebltable)"
    print(ebl)

###
# Recommended use in command line:
#
# from SpectraPlotter import SpectraPlotter
# folder = {name of folder where output blackhawk files are}
#
# sp = SpectraPlotter.from_folder(folder)
# print(sp.particle_names()) <-- This gives all the species you can apply plotting on.
#
# e.g.
# fig = sp.PSD('photon_total')
# fig.savefig("PSD_photon_total.png")
###

file_headers = { 
        'mass' : 'life_evolutions.txt',
        'photon' : ['photon_primary_spectrum.txt', 'photon_secondary_spectrum.txt'],
        'electron' : ['electron_primary_spectrum.txt', 'electron_secondary_spectrum.txt'],
        'neutrino' : ['neutrinos_primary_spectrum.txt', 'neutrino_e_secondary_spectrum.txt', 'neutrino_mu_secondary_spectrum.txt', 'neutrino_tau_secondary_spectrum.txt']}

class SpectraPlotter:

    def __init__(self, dataset):
        # Dataset is a dictionary of particle/dataframe pairs
        self.data = dataset

    @classmethod
    def from_folder(self, folder:str = 'tot_1010/'):
        # Just give the folder where you're keeping all the files -- don't change name/format of any output files.
        data = {}

        data['mass'] = pd.read_csv(folder + file_headers['mass'], sep='         ', skiprows=[0,1])

        files = file_headers['photon']
        data['photon_primary'] = pd.read_csv(folder +  files[0], sep='    ', skiprows=1)
        data['photon_secondary'] = pd.read_csv(folder +  files[1], sep='    ', skiprows=1)
        
        files = file_headers['electron']
        data['electron_primary'] = pd.read_csv(folder +  files[0], sep='    ', skiprows=1)
        data['electron_secondary'] = pd.read_csv(folder +  files[1], sep='    ', skiprows=1)
        
        files = file_headers['neutrino']
        data['neutrino_primary'] = pd.read_csv(folder +  files[0], sep='    ', skiprows=1)
        data['neutrino_e'] = pd.read_csv(folder +  files[1], sep='    ', skiprows=1)
        data['neutrino_mu'] = pd.read_csv(folder +  files[2], sep='    ', skiprows=1)
        data['neutrino_tau'] = pd.read_csv(folder +  files[3], sep='    ', skiprows=1)

        if data['photon_primary'].shape == data['photon_secondary'].shape:
            
            data['photon_total'] = data['photon_primary'].add(data['photon_secondary'])
            # Above also adds time column -- reset it.
            data['photon_total']['time/energy'] = data['photon_primary']['time/energy']

            data['neutrino_total'] = data['neutrino_tau'].add(data['neutrino_mu']).add(data['neutrino_tau']) 
            data['neutrino_total']['time/energy'] = data['neutrino_tau']['time/energy']

        return SpectraPlotter(data)

    def particle_names(self):
        return self.data.keys()

    def apply_telescope_bands(self, axes:list):
        # Helper method to apply telescope bands to any axes, assuming xlim units are GeV
        telescopes = pd.read_csv('telescope_energy_bands.csv')
        colors = ['red', 'orange', 'green', 'brown']

        for ax in axes:
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            for i in range(len(telescopes.index)):
                arr = list(telescopes.loc[i])
                ax.plot( arr[2:4] , np.mean(ylim)/3*(i+1)*np.mean(ylim)/1000 , label = arr[1], color = colors[i] )
            # Reset axes to original dimensions
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

        #optional return
        return axes

    """ TODO
    def apply_ebl(self, ax, model_name = 'franceschini', redshift = 0.2, replace = False):
        if ebl is not True:
            return ebl 
        data_x = ax.lines.get_xdata() * u.GeV.to('TeV')
        data_y = ax.lines.get_ydata()

        tau = OD(model = model_name)
        t = tau.opt_depth_array(z,data_x)

        if replace is True:
            ax.lines[0].set_ydata(new_data_y)
        else:
            ax.plot(data_x, new_data_y, label = 'Attenuated data -- {model_name}')
        return ax
    https://github.com/me-manu/ebltable/blob/master/examples/example.py
    """

    def mass_time_evolution(self):
        ###
        # Evolution of mass over time.
        # Returns fig object
        ###
        accretion = self.data['mass']
        cols = accretion.columns
        accretion[cols[1]] = pd.to_numeric( accretion[cols[1]].str[:11] )

        fig, ax = plt.subplots()
        ax.scatter(accretion['t'], accretion[cols[1]], marker='+', s=5)
        ax.set_title('Mass accreted over time')
        ax.set_xlabel('t [s]')
        ax.set_ylabel('M [g]')
        #fig.savefig('mass(t).png')
        return fig

    def spectra_time_evolution(self, name):
        ###
        # Evolution of given spectra over time.
        # Takes name of particle
        # Returns fig object
        ###

        assert name in self.data.keys(), 'Particle name must be in %s'%file_headers.keys()

        spectrum = self.data[name]

        fig, ax = plt.subplots()
        for i in spectrum.columns[1:]:
            ax.plot(spectrum['time/energy'], spectrum[i], label=i)
        ax.set_title('%s spectrum, time evolved'%name)
        ax.set_xlabel('t [s]')
        ax.set_ylabel('Particle/ unit Energy [GeV^-1 s^-1 cm^-3]')
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.legend()

        return fig

    def initial_ED(self, name):
        ###
        # Gives energy distribution of given particle at time 0.
        # Takes name of particle
        # Returns fig object
        ###

        fig, ax = plt.subplots() 

        ax = self._apply_initial_ED(ax, name)
        
        ax.legend()

        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_title('Initial %s energy distribution'%name)
        ax.set_xlabel('Energies [GeV]')
        ax.set_ylabel('Emission rate d^2N/dtdE [s^-1]')
        #fig.savefig('energy_dist_t0.png')

        return fig
    def _apply_initial_ED(self, ax, name):
        
        spectrum = self.data[name]

        spectrumx = spectrum.columns[1:].astype('float64').to_list()
        spectrumy = spectrum.loc[0].to_list()[1:]

        ax.scatter(spectrumx, spectrumy, s=5, label=name)

        return ax


    def SED(self, name, row = 0):
        ###
        # Gives Spectral Energy Distribution of given particle and row. 
        # Takes name of particle and a row on dataset
        # Returns fig object
        ###

        fig, ax = plt.subplots()

        ax, time = self._apply_SED(ax, name, row=row)
        #ax.set_xlim(min(secondary_energy), max(secondary_energy))
        #ax.set_ylim(min(secondary_SED), max(secondary_SED))

        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_title('%s SED of spectra at t=%ss'%(name,time))
        ax.set_xlabel('Energies [GeV]')
        ax.set_ylabel('Spectral energy density d^2N/dtdE [s^-1] * E^2')
        #fig.savefig('SED_t%s.png'%time)
        
        return fig
    def _apply_SED(self, ax, name, row=0):
        
        spectrum = self.data[name]

        t = row 
        time = spectrum.loc[t].to_list()[0]
        print('row %s, time %s'%(t, time))

        spectrum_energy = spectrum.columns[1:].astype('float64').to_list()
        
        # gets us flux
        spectrum_flux = spectrum.loc[t].to_list()[1:]
        
        # gets us flux*E^2
        spectrum_SED = np.multiply( np.square(spectrum_energy), spectrum_flux)
        ax.scatter(spectrum_energy, spectrum_SED, s=5, label=name)

        return ax, time

    def PSD(self, name):
        ###
        # Gives Power Spectrum Distribution of given particle. 
        # Takes name of particle
        # Returns fig object
        ###
 
        fig, ax = plt.subplots()

        self._apply_PSD(ax, name)         
        
        ax.set_title('Power Spectral Density of %s'%name)
        ax.set_xlabel('Time s')
        ax.set_ylabel('Power W')
        ax.legend()
        #fig.savefig('Power.png')

        return fig
    def _apply_PSD(self, ax, name):
        spectrum = self.data[name]

        E = np.array(spectrum[:0].columns[1:]).astype('float') * u.GeV

        E = E.to(u.J)

        dE = np.diff(E)

        t = spectrum['time/energy'].to_list() * u.s

        dNdE = []
        P = []

        for N in np.array(spectrum):

            #deriv = np.diff(N[1:])/dE
            #dNdE.append(deriv)

            integrand = np.multiply(E,N[1:])#deriv

            P.append(integrate.simps(integrand, x = E))

        #ax.plot(t/u.s,np.ones(len(t)) * np.mean(P), '-.', label = 'Average Power %s %.e'%(name,np.mean(P)))
        ax.plot(t,P, label='Integrated Power %s: maximum %.e'%(name,np.max(P)), color = 'black')

        # Keeping this here for reference, the expected power out of photon spectrum
        #ax.plot(t,np.ones(len(t))*2*10**(-6),'--', label='Expected Power %.e'%(2*10**(-6)))

        power = self._get_total_power(name, x = t, y = P)

        ax.scatter(t,np.ones(len(t))*power[1],label='Average Power %.2e'%power[1].value, s=5)
        ax.text(60, .025, 'Total time-integrated Power %.2e'%power[0]) 
        return 

    def _get_total_power(self, name, x = None, y = None, ax = None): 
        
        if x is None and y is None:
            assert ax is not None, 'Give data through x, y and ax'
            y = ax.lines.get_ydata()
            x = ax.lines.get_xdata()
        
        """
        There are some parts where np.diff(t) = 0 which means that the Riemmann integral cannot be calculated.
        Sketchy workaround:
        locations = np.where(np.diff(t) == 0)
        new_t = np.delete(t, locations)
        new_y = np.delete(y, locations)
        """
        locations = np.where(np.diff(x) == 0)
        if locations != []:    
            x = np.delete(x, locations)
            y = np.delete(y, locations)

        tot_P = integrate.simps(y, x=x)
        average_P = tot_P / (x[-1] - x[0])

        #ax.plot(t, np.zeros(len(t))*average_P, '-.', label = 'Average Power %s %.e'%(name, average_P))
        
        return tot_P, average_P

    def yearly_energy(self, name, burst_density = 3300):

        # Extract average power per burst.
        fig, ax = plt.subplots()
        PSD = self._apply_PSD(ax, name)

        data_x = ax.lines.get_xdata()
        data_y = ax.lines.get_ydata()

        avg_Power = np.mean(P) # Average can't just be mean of array

        # Calculate power over the year, using inverse square law + EBL suppression to find energy incident

        return 'fcn still under construction'


###
    def get_telescope_band_flux_wrt_time(self, energy_band, spectrum_name = 'photon_secondary'):
        spectrum = self.data[spectrum_name]

        # total energy range (in GeV)
        E = np.array(spectrum[:0].columns[1:]).astype('float')
        
        # Match Emin, Emax to elements on E
        Emin,Emax = energy_band

        # Indices to pick up energy band
        index_min = np.where(E > Emin)[0][0]
        index_max = np.where(E < Emax)[-1][-1]
        
        discrete_energy_band = E[index_min:index_max]
        
        # Calculate flux \propto time in the energy band
        flux = []
        for fluxdtdE in np.array(spectrum): 
        # Goes down the array in time, grabbing rows of flux counts over the given energy band
            flux.append(
                    integrate.simps(fluxdtdE[index_min:index_max], x=discrete_energy_band))
        return flux

    def plot_telescope_flux_ratio(self):

        # In GeV
        LAT_energy_band = (0.1, 1000)
        GBM_energy_band = (50 * 10**(-6), 300 * 10**(-6))

        spectra = ['photon_primary', 'photon_secondary']

        # Get the flux spectra for primary, secondary of each LAT and GBM
        LAT_primary = self.get_telescope_band_flux_wrt_time(LAT_energy_band, spectrum_name = 'photon_primary')
        LAT_secondary = self.get_telescope_band_flux_wrt_time(LAT_energy_band, spectrum_name = 'photon_secondary')

        GBM_primary = self.get_telescope_band_flux_wrt_time(GBM_energy_band, spectrum_name = 'photon_primary')
        GBM_secondary = self.get_telescope_band_flux_wrt_time(GBM_energy_band, spectrum_name = 'photon_secondary')

        # Should have same time sampling, so just add primary/secondary fluxes
        LAT_flux = np.add(LAT_primary, LAT_secondary)
        GBM_flux = np.add(GBM_primary, GBM_secondary)

        #On that note, get the time too
        spectrum = self.data['photon_secondary']
        t = spectrum['time/energy'].to_list()

        # Get ratio of GBM to LAT
        flux_ratio = np.divide(GBM_flux,LAT_flux)

        fig, ax0 = plt.subplots()
        fig.set_figheight(10)
        #fig.set_figwidth(20)
        fig.set_figwidth(10)
        ax0.scatter(t,flux_ratio, s=5, label = 'Blackhawk-computed flux ratio')
        
        ax0.set_xlabel('Time (s)')
        ax0.set_ylabel(r'$\gamma$ Flux Ratio $\frac{\phi_{GBM}}{\phi_{LAT}}$')
        
        #ax1.scatter(t, np.divide(GBM_primary, LAT_primary), s=5, color = 'orange')
        #ax1.plot(t, LAT_flux, '--', label = 'Total flux', color = 'black')
        #ax1.scatter(t, LAT_primary, s=5, label = 'Primary (direct) flux')
        #ax1.scatter(t, LAT_secondary, s=5, label = 'Secondary (fragmentary) flux')
        #ax1.set_xlabel('Time (s)')
        #ax1.set_ylabel('Flux Ratio for Direct Spectrum')
        #ax1.set_yscale('log')
        #ax1.legend()

        #ax2.scatter(t, np.divide(GBM_secondary, LAT_secondary), s=5, color = 'green')
        #ax2.plot(t, GBM_flux, '--', label = 'Total flux', color = 'black')
        #ax2.scatter(t, GBM_primary, s=5, label = 'Primary (direct) flux')
        #ax2.scatter(t, GBM_secondary, s=5, label = 'Secondary (fragmentary) flux')
        #ax2.set_xlabel('Time (s)')
        #ax2.set_ylabel('Flux Ratio for Indirect Spectrum')
        #ax2.set_yscale('log')
        #ax2.legend()


        return fig


"""
Look through F-L catalog to see if there are any candidates that match features

there must be no redshift
luminosity distance must be within 10^17 cm

FIND LAT & GBM ENERGY BANDS
Integrate blackhawk flux over the energy band range
query LAT catalog to see if there are sources which fit the RATIO we expect e.g. GBM/LAT photon fluxes

cross-reference with lifetime
"""
