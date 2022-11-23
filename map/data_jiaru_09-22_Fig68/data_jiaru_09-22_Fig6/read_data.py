import numpy as np
import pandas as pd

m = 6 # = -log10(mu)

# eps1range, nconjcrit = np.loadtxt(f"Nconjcrit_continuous_{m}.dat") # nofriction

tau = 10000
eps1range, nconjcrit = np.loadtxt(f"Nconjcrit_continuous_friction_{m}_{tau}.dat") # nofriction

mu = 10**(-m)
neps1 = len(eps1range)

### Switch to the N-body parameters K and t_inst

rH = (float(mu)/3)**(1/3) # Hill radius

k = eps1range/rH

tsyn = 2/(3*eps1range) # initial synodic period in units of P1
tinst = nconjcrit*tsyn

tmax = 1e5 # if we want to cap it to 1e5 P1
tinst = np.minimum(tinst, tmax)

### Because of the nature of the map, the timescales can only be a integer number of the synodic period
### To smooth the results, we can assume the instability developped randomly during the previous synodic period, but not sure if it is very rigourous
# tinst -= np.random.random(nlambda)*tsyn

### Representation
import matplotlib.pyplot as plt

plt.scatter(k, tinst, s=2)

plt.xlabel("K")
plt.ylabel("Time to instability [planet orbit]")
plt.yscale("log")
plt.tight_layout()

plt.show()
