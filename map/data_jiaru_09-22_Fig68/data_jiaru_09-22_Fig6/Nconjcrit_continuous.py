###
# Critical number of conjunctions as a function of eps1 for a given mu
###
import numpy as np
from tqdm import tqdm

neps1 = 1000
mu = 1e-3 #mp/Mstar

rH = (mu/3)**(1/3)
krange = np.random.random(neps1)*3. + 2.
eps1range = krange*rH # initial eps = (a-ap)/ap

g = 2.24

tconjcrit = 1e5 # max number of simulated orbital periods
nconj = int(tconjcrit*3*max(eps1range)/2)
print(f"Max number of conjunctions corresponding to {tconjcrit} orbital periods: {nconj}")

nconjcrit =  np.zeros(neps1)

for k, eps1 in enumerate(tqdm(eps1range)):

    zc = np.zeros(nconj, dtype=complex) # complex eccentricity
    lc = np.zeros(nconj) #phase
    eps = np.zeros(nconj) #eps = (a-ap)/ap

    lc[0] = np.random.random()*2*np.pi
    zc[0] = mu
    eps[0] = eps1
    j = 0 # j is the number of conjunctions before instability

    if g/(eps1**2)*mu > eps1:
        j = 1
    else:

        while((j<nconj-1) and (abs(zc[j]) < 0.5*eps[j])): #instability criterion
            j += 1
            zctemp = zc[j-1] + 1j*g*np.exp(1j*lc[j-1])/(eps[j-1]**2)*np.sign(eps1)*mu
            epstemp = eps[j-1]*np.sqrt(1+4*(np.absolute(zctemp)**2 - np.absolute(zc[j-1])**2)/(3*eps[j-1]**2))
            epstemp = (epstemp+eps[j-1])/2
            zc[j]= zc[j-1] + 1j*g*np.exp(1j*lc[j-1])/(epstemp**2)*np.sign(eps1)*mu
            eps[j] = eps[j-1]*np.sqrt(1+4*(np.absolute(zc[j])**2 - np.absolute(zc[j-1])**2)/(3*eps[j-1]**2))
            lc[j] = lc[j-1] + 2*np.pi/((1+eps[j])**(3/2)-1) % (2*np.pi)

            # print(abs(zc[j])/(0.5*eps[j]))

    nconjcrit[k] = j

np.savetxt(f"Nconjcrit_continuous_{int(-np.log10(mu))}.dat", [eps1range, nconjcrit])
