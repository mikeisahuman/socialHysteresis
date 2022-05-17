##  Checking fixed point graphical solution & critical point
##  by Mike Phillips, 8/18/2018
##
##  Transcendental equations are plotted for a few situations, and corresponding predicted critical values
##      are recorded for comparison -- given critical value is for appearance of FP near y=0.

import numpy as np
import matplotlib.pyplot as plt

sech = lambda x: 1/np.cosh(x)
tanh = lambda x: np.tanh(x)

x = 0
y = np.linspace(-1,1,100)
dlst = (0, 0.4, 0.8)
for d in dlst:
    g=-d

    kacrit = 2/((1+x)*(sech(d)**2) + 0.5*(1-x)*(sech(g)**2))     #   for kb=ka/2
    ka = np.array( (kacrit/2, kacrit, kacrit*1.5) )
    kb = ka/2

    print("\n" + "critical value for dlt=%1.2f, gam=%1.2f: ka = %1.2f" % (d, g, kacrit) + "\n")

    fig = plt.figure(num=dlst.index(d), figsize=(6,6))       # "fixed point critical pt. check"
    ax = fig.add_subplot(111, xlabel="y", xlim=(-1,1), ylabel="y ; rhs(y)", ylim=(-1,1))
    ax.grid()
    ax.tick_params(axis="both", direction="in", which="both")
    ax.plot(y, y, "k-")

    stys = ("--", "-.", ":")
    for i in range(len(ka)):
        rhs = .5*(1+x)*tanh(d + ka[i]*y) + .5*(1-x)*tanh(g + kb[i]*y)
        ax.plot(y, rhs, linestyle=stys[i],
                label=r"$\kappa_A$" + "/" + r"$\kappa^\ast_A$" + "=%1.2f" % (ka[i]/kacrit) )

    plt.legend(loc="upper left", frameon=True, edgecolor="inherit")
    plt.tight_layout()

plt.show()
plt.close()

