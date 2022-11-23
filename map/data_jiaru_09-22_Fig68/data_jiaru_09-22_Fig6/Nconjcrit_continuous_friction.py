###
# Critical number of conjunctions as a function of eps1 for a given mu with friction
###
import numpy as np
from tqdm import tqdm

neps1 = 100
mu = 1e-6
taue = 10000 # in orbital periods

rH = (mu/3)**(1/3)
krange = np.random.random(neps1)*3. + 2.
eps1range = krange*rH # initial eps = (a-ap)/ap

g = 2.24

tconjcrit = 1e5
nconj = int(tconjcrit*3*max(eps1range)/2)
print(f"Number of conjunctions corresponding to {tconjcrit} orbital periods: {nconj}")

nconjcrit =  np.zeros(neps1)

for k, eps1 in enumerate(tqdm(eps1range)):

    zc = np.zeros(nconj, dtype=complex)
    lc = np.zeros(nconj)
    eps = np.zeros(nconj)

    lc[0] = np.random.random()*2*np.pi
    zc[0] = mu
    eps[0] = eps1
    j = 0

    if g/(eps1**2)*mu > eps1:
        j = 1
    else:

        while((j<nconj-1) and (abs(zc[j]) < 0.5*eps[j])):
            j += 1
            zctemp = zc[j-1] + 1j*g*np.exp(1j*lc[j-1])/(eps[j-1]**2)*np.sign(eps1)*mu
            epstemp = eps[j-1]*np.sqrt(1+4*(np.absolute(zctemp)**2 - np.absolute(zc[j-1])**2)/(3*eps[j-1]**2))
            epstemp = (epstemp+eps[j-1])/2
            zc[j]= zc[j-1] + 1j*g*np.exp(1j*lc[j-1])/(epstemp**2)*np.sign(eps1)*mu
            eps[j] = eps[j-1]*np.sqrt(1+4*(np.absolute(zc[j])**2 - np.absolute(zc[j-1])**2)/(3*eps[j-1]**2))
            lc[j] = lc[j-1] + 2*np.pi/((1+eps[j])**(3/2)-1) % (2*np.pi)
            Tn = 1/(1-(1+eps[j])**(-1.5))
            zc[j] *= np.exp(-Tn/taue)

            # print(abs(zc[j])/(0.5*eps[j]))

    nconjcrit[k] = j

np.savetxt(f"Nconjcrit_continuous_friction_{int(-np.log10(mu))}_{taue}.dat", [eps1range, nconjcrit])
