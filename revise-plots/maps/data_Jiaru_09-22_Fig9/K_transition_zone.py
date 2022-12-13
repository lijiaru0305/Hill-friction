###
# First stable eps1 and last unstable eps1 as a function of mu
###
import numpy as np
import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(8,5), sharex="all", sharey="all")

murange = np.array([1e-7, 1e-6, 1e-5, 1e-4, 1e-3])

k0jiaru = np.array([2.6, 2.6, 2.5, 2.5, 2.3])
kcjiaru = np.array([3.4, 3.4, 3.4, 3.6, 3.6])
bjiaru = np.array([13.8, 14.6, 10.9, 17.3, 10.4])

coefDQT = 1.6

tau = 1000 # friction time in orbital period

rH = (murange/3)**(1/3)

tsyn = 2/(3*k0jiaru*rH)

ktau = (1/bjiaru * np.log(tau/tsyn) + 1)*k0jiaru
ktau = np.maximum(ktau, k0jiaru)
ktau = np.minimum(ktau, kcjiaru)

axes[0].plot(murange, kcjiaru, "C0o-", label=r"$K_\mathrm{crit}$ (N-body)")
axes[1].plot(murange, ktau, "C0o-", label=r"$\tilde{K}_\tau$ (N-body)")

mu, eps1critmin, eps1critmax = np.loadtxt(f"Interval_phase.dat")

eps1critDQT = coefDQT*(mu**(2/7))
rH = (mu/3)**(1/3)
kDQT = eps1critDQT/rH

for k, m in enumerate(mu):
    # plt.loglog([m, m], [eps1critmin[k], eps1critmax[k]], "C0")
    kmin = eps1critmin[k]/rH[k]
    kmax = eps1critmax[k]/rH[k]
    if k == 0:
        axes[0].plot([m, m], [kmin, kmax], "0.5", label="Map (grey zone)")
        # axes[0].plot([m, m], [2., kmin], "r", label="Map (unstable)")
    else:
        axes[0].plot([m, m], [kmin, kmax], "0.5")
        # axes[0].plot([m, m], [2., kmin], "r")


axes[0].plot(mu, kDQT, 'C2', label="DQT89")

mu, eps1critmin, eps1critmax = np.loadtxt(f"Interval_phase_friction_{tau}.dat")

rH = (mu/3)**(1/3)

for k, m in enumerate(mu):
    # plt.loglog([m, m], [eps1critmin[k], eps1critmax[k]], "C0")
    kmin = eps1critmin[k]/rH[k]
    kmax = eps1critmax[k]/rH[k]
    if k == 0:
        axes[1].plot([m, m], [kmin, kmax], "0.5", label="Map (grey zone)")
        # axes[1].plot([m, m], [2., kmin], "r", label="Map (unstable)")
    else:
        axes[1].plot([m, m], [kmin, kmax], "0.5")
        # axes[1].plot([m, m], [2., kmin], "r")

axes[0].set_xscale("log")
axes[0].set_ylim(ymin=2)
axes[0].legend(fontsize=12)
axes[1].legend(fontsize=12)

axes[0].set_xlabel(r"$\mu$")
axes[0].set_ylabel(r"K")
axes[0].set_title("No friction")
axes[1].set_title(r"$\tau = %s~P_1$" % tau)

axes[1].set_xlabel(r"$\mu$")

# plt.margins(0)
plt.tight_layout()

plt.show()
