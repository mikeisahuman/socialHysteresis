## Check phase diagram & transition curves: dlt vs. kapa
## by Mike Phillips, 9/18/2019
##
## The phase diagram class "Phases" is used (from socSolver) to generate phase diagrams.
## Predicted transition lines are drawn atop these diagrams for comparison.

import socSolver as soc
from socSolver import np, plt

save = None
colors, stys = ["#0000FF", "#CC6600", "#336633"], ["-","--","-."]

##  Plot font sizes
SMALL_SIZE = 10
MEDIUM_SIZE = 12
LARGE_SIZE = 14

plt.rc('font', size=MEDIUM_SIZE)        # controls default text sizes
plt.rc('axes', titlesize=LARGE_SIZE)    # fontsize of the axes title
plt.rc('axes', labelsize=LARGE_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)   # legend fontsize
plt.rc('figure', titlesize=LARGE_SIZE)  # fontsize of the figure title

##  Setup
modes = ("mult", "diff")
fr = {modes[0]:0.5, modes[1]:2.2}
cosh = np.cosh
def ka_m(d, g, x, f):
    res = 2*( cosh(d)*cosh(g) )**2
    res = res / ((1+x)*(cosh(g))**2 + f*(1-x)*(cosh(d))**2)
    return res
def ka_d(d, g, x, c):
    res = 2*(cosh(g))**2 - c*(1-x)
    res *= (cosh(d))**2
    res = res / ((1+x)*(cosh(g))**2 - (1-x)*(cosh(d))**2)
    return res
func = {modes[0]:ka_m, modes[1]:ka_d}

##  Plots
interp = None
yzi = (-0.4,0.1)
nmax = 250
dmax = 0.8
kamax = 2.0
gam = -0.4
x = 0.4

d = np.linspace(-dmax, dmax, 200)

for m in modes:
    fc = fr[m]
    f = func[m]
    p = soc.Phases(yzi, nmax, respts=12, pars=("dlt","kapa"), dlt=dmax, kapa=kamax, gam=gam, x=x,
               parFrac = [{"kapb":fc}, m])
    pax = p.phasePlot(interp=interp, show=False, save=None, title=False)
    ka = f(d, gam, x, fc)
    pax.plot(ka, d, "k--")
    pax.set_xlim(0, kamax)
    pax.set_ylim(-dmax, dmax)
    if type(save) == str:
        filepath = os.path.join(os.getcwd(), "socSolver_figs", "phases", save)
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        if os.path.isfile(filepath):
            print("ERROR: requested file name '%s' already exists is the directory '%s'."
                  % (save, os.path.dirname(filepath)) )
        else:
            plt.savefig(filepath)
    print("\n" + "plotting mode: '" + m + "' w/ yzi=(%1.1f,%1.1f)...\n" % yzi)
    plt.show()
    plt.close()


del soc, np, plt
