import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt

"""
This is just a quick calc of the proximity limit for an EBH
"""

flux = lambda T: 1e29 * (T/u.TeV)**(1.6) / u.s

d = lambda f, sens = 1e-9: np.sqrt( (1/sens) * u.cm**2 * u.s * f / (4*3.14159) )

E = np.linspace(0.3, 5, 100) * u.TeV

prox = d( flux( E ) )
print(flux(E))
#print(prox)
prox = prox.to('pc')

fig, ax = plt.subplots()
ax.plot( E, prox , '--' , label = 'Maximum distance for photon detection')
ax.fill_between(E.value, prox.value, alpha = 0.5, label = 'Discoverable point sources')

#ax.plot( E, d(flux(E), sens=1e-15).to('pc'), '--', label ='Former limit')

ax.legend()
ax.set_xlabel('Energy (TeV)')

ax.set_ylabel(r'Distance (pc)')

fig.savefig('proximity.png')
