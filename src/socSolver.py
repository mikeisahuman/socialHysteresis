##  Object Class Definition for solving the Social Dynamics ODEs (2 population / 2 behaviors)
##  by Mike Phillips, 7/1/2018
##
##  All features of the problem are well-specified (see "meanValueMain.py"), and we will seek a simple
##      way to manipulate solutions to the ODEs: find, present, record, allow for par. changes/looping.
##  Many built-in attributes (like plots) should be available, as well as the ability to retrieve
##      and compile the various data sets (e.g. for different par.s) to plot them, etc. (separately).


import os
import numpy as np
import sympy as sy
from itertools import combinations
#from scipy import optimize
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib import cm
##from mpl_toolkits.axes_grid1 import make_axes_locatable     # use: matplotlib ver. 2.x (NOT 3.x)
##from mpl_toolkits.mplot3d import axes3d

##  Plot font sizes
SMALL_SIZE = 10
MEDIUM_SIZE = 13
LARGE_SIZE = 14

plt.rc('font', size=MEDIUM_SIZE)        # controls default text sizes
plt.rc('axes', titlesize=LARGE_SIZE)    # fontsize of the axes title
plt.rc('axes', labelsize=LARGE_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)   # legend fontsize
plt.rc('figure', titlesize=LARGE_SIZE)  # fontsize of the figure title

##  Plot ticks / frame
def myTickParams(ax, dr="in"):
    ax.tick_params(axis="both", which="both", direction=dr, left=True, bottom=True, right=True, top=True,
                   labelleft=True, labelbottom=True, labelright=False, labeltop=False)
    return

##  Symbol definitions
time, omega, amp = sy.symbols("t om A", real=True)

#####   Class Definition: critical arguments & parameters with default values   #####
#####   
#####   MODES: "const" = time-indep., "x-osc" = oscillating infection
#####       "dlt-osc" = osc. healthy opinion, "gam-osc" = osc. infected opinion
#####       "kapa-osc" = osc. healthy interactions, "kapb-osc" = osc. interactions
class Soln:
    """
    Solution object for 2-population, 2-behavior model with parasitic behavior modifications.

    Attributes / Args:
        trng (2-tuple): Effective range of time to be used and plotted (real range always begins at t=0).
        yzi (2-tuple):  Co-ordinate pair of initial conditions in y-z plane.
        nmax (int):     Maximum number of time steps used in the solution (i.e. time-space resolution).
        mode (str):     Solution mode: "const" for time-independent parameters, or one of
                            ("x","dlt","gam","kapa","kapb") followed by "-" and one of
                            ("osc","shift","step","smooth") to have the chosen single dynamic parameter
                            following the chosen time-dependent pattern -- 
                            oscillatory (cos), quick shift (tanh), near-step function, smooth shift.
                            In the dynamic case, the given parameter value (see below) is used as a range limit.
        omg_frac (float): [dynamic mode "osc" only] Controls the angular frequency of oscillations, given as
                            a fraction of the system frequency (interaction rate) for population A.
        offset (float):   [dynamic modes only] Controls the vertical offset of the function,
                            e.g. to allow for oscillations about a non-zero value.
        phase (float):    [dynamic modes only] Controls the phase of the cosine function ("osc" mode),
                            or re-scales the amplitude of the shift or step ("shift","step","smooth" modes).
        na, nb (float, positive):  Frequencies (interaction rates) for populations A and B.
        dlt, gam (float):          Internal opinions for populations A and B.
        kapa, kapb (float):        Indirect interaction coefficients for populations A and B.
        x (float, in [-1,1]):      Relative population of infected/parastitically-influenced individuals.
        Qxx (float):    Fluctuations of relative infection population. Only used in the case of applying the
                            external control of "x" as a quasi-dynamic quantity (dx_shift=True in Kz() below).
        solve (bool):   Boolean value to determine whether or not to solve upon initiation of the object.
        s_mode (none or str or int): None or "direct" or 1 --> First order solution method (i.e. dx=f*dt).
                                     "runge-kutta" or "rk" or 4 --> Runge-Kutta fourth order method.
    """
        
    def __init__(self, trng=(0,25), yzi=(0,0), nmax=100, mode="const", omg_frac=2/3, offset=0, phase=0, 
                     na=1, nb=1, dlt=0, gam=0, kapa=0.5, kapb=0.5, x=0, Qxx=0, solve=True, s_mode="rk"):
        mds = ["const"]
        mds += [ (par + "-" + func) for par in ("x","dlt","gam","kapa","kapb")
                 for func in ("osc","shift","step","smooth") ]
        self._MODES = tuple(mds)
        self._mode = None
        self._MODEPROPS = {"offset":offset, "phase":phase}
        self._pars = {"na":na, "nb":nb, "dlt":dlt, "gam":gam, "kapa":kapa, "kapb":kapb, "x":x, "Qxx":Qxx}
        self._real_tmin = 0
        self._efftrng = trng
        self._trng = (self._real_tmin, trng[1])
        self._nmax = nmax
        self._tstep = (trng[1] - self._real_tmin) / nmax
        self._tspace = np.linspace(self._real_tmin, trng[1], nmax+1)
##        self._tstep = (trng[1] - trng[0]) / nmax
##        self._tspace = np.linspace(trng[0], trng[1], nmax+1)
        self._yzi = yzi
        self._nmax = nmax
        self._na_const, self._nb_const = na, nb
        self._dlt_const, self._gam_const = dlt, gam
        self._kapa_const, self._kapb_const = kapa, kapb
        self._x_const, self._Qxx_const = x, Qxx
        self._omgfrac = omg_frac
        self.setMode(mode)
        if solve:
            self.Solve(s_mode=s_mode)
        else:
            self._sol, self._isfp = None, (None, None)
    #   Accessors (and a setter: just initial conditions & mode)
    def maxSteps(self):
        return self._nmax
    def trng(self, eff=False):
        if eff:
            return self._efftrng
        else:
            return self._trng
    def tstep(self):
        return self._tstep
    def tspace(self, tmin=None, tmax=None):
        if (tmin == None or tmin == 0) and (tmax == None or tmax == self.trng()[-1]):
            return self._tspace
        else:
            i = -1
            tsp = self._tspace
            if tmin > max(tsp) or tmin > tmax:
                print("ERROR: the entered minimum time %1.2f is greater than the maximum time %1.2f"
                      % (tmin, max(tsp)) )
                return ()
            else:
                imin = 0
                imax = -1
                for t in tsp:
                    i += 1
                    if tmin != None and t >= tmin and imin == 0:
                        imin = i
                    elif tmax != None and t > tmax and imax < imin:
                        imax = i
                return tsp[imin:imax]
                        
    def IC(self):
        return {"yi":self._yzi[0], "zi":self._yzi[1]}
    def setIC(self, yzi=(0,0)):
        self._yzi = yzi
        return
    def sol(self):
        return self._sol
    def ysol(self):
        return self.sol()[0]
    def zsol(self):
        return self.sol()[1]
    def isfp(self):
        return self._isfp
    def fp(self):
        if self.isfp()[0]:
            return (self.sol()[0][-1], self.sol()[1][-1])
    def allModes(self):
        return self._MODES
    def Mode(self):
        return self._mode
    def setMode(self, mode="const"):        #   NOTE: want generic function forms "-osc", etc. for each par
        md = mode.lower()
        if md not in self._MODES:
            print("\n\t" + "ERROR: given mode '%s' is not one of the supported modes." % mode)
            return
        self._mode = md
        tanh_mds = ("shift", "step", "smooth")
        if md.count("-") == 1:
            ind = md.index("-") + 1
        else:
            ind = 0
        mdname = md[ind:]
        if mdname == "osc":
            omg = self._na_const * (self._omgfrac)
            self._symfunc = self._MODEPROPS["offset"] + amp*sy.cos(omega*time + self._MODEPROPS["phase"])
            self._symfunc = self._symfunc.subs(omega, omg)
            self._numfunc = sy.lambdify((amp,time), self._symfunc, "numpy")
        elif mdname in tanh_mds:
            if mdname == "shift":
                k = 1/2
            elif mdname == "step":
                k = 12
            elif mdname == "smooth":
                k = 1/5
            midpt = (self._trng[1] - self._trng[0]) / 2
            self._symfunc = self._MODEPROPS["offset"] + self._MODEPROPS["phase"]*amp*sy.tanh(k*(time - midpt))
            self._numfunc = sy.lambdify((amp,time), self._symfunc, "numpy")
        if md[0] == "x":
            self._dx = self._symfunc.diff(time)
            self._dx = sy.lambdify((amp,time), self._dx, "numpy")
        else:
            self._dx = 0
        return
    def pars(self, par = None, wantPrint = False):
        if wantPrint:
            prt = lambda s,v: print("\n\t" + "%s = %1.3f" % (s, v) )
        else:
            prt = lambda s,v: 0
        if type(par) == str:
            if par in self._pars:
                prt(par, self._pars[par])
                return self._pars[par]
            else:
                print("\n\t" + "ERROR: no parameter found matching '%s'" % par)
                return
        else:
            if wantPrint:
                print("\n")
                for par in self._pars:
                    prt(par, self._pars[par])
                print("\n")
            return self._pars
    #   Parameter functions (for time dependence)
    def na(self, t):
        return self._na_const
    def nb(self, t):
        return self._nb_const
    def dlt(self, t):
        if self._mode[:3] == "dlt":
            return self._numfunc(self._dlt_const, t)
        else:
            return self._dlt_const
    def gam(self, t):
        if self._mode[:3] == "gam":
            return self._numfunc(self._gam_const, t)
        else:
            return self._gam_const
    def kapa(self, t):
        if self._mode[:4] == "kapa":
            return self._numfunc(self._kapa_const, t)
        else:
            return self._kapa_const
    def kapb(self, t):
        if self._mode[:4] == "kapb":
            return self._numfunc(self._kapb_const, t)
        else:
            return self._kapb_const
    def x(self, t):
        if self._mode[0] == "x":
            return self._numfunc(self._x_const, t)
        else:
            return self._x_const
    def dx(self, t):
        if self._mode[0] == "x":
            return self._dx(self._x_const, t)
        else:
            return self._dx
    def Qxx(self, t):
        return self._Qxx_const
    #   'Drift Coefficients' (functions)
    def Ky(self, y, z, t):
        res = self.na(t) * ( (1 + self.x(t)) * np.sinh(self.dlt(t) + self.kapa(t)*y)
                             - (y + z) * np.cosh(self.dlt(t) + self.kapa(t)*y) )
        res += self.nb(t) * ( ( 1 - self.x(t)) * np.sinh(self.gam(t) + self.kapb(t)*y)
                             - (y - z) * np.cosh(self.gam(t) + self.kapb(t)*y) )
        return res
    def Kz(self, y, z, t, dx_shift=False):
##        res = nc * (y - z) - ni * (y + z) * np.exp(-kapi*x)   #   full-blown: if x is dynamic as well
        res = self.na(t) * ( (1 + self.x(t)) * np.sinh(self.dlt(t) + self.kapa(t)*y)
                             - (y + z) * np.cosh(self.dlt(t) + self.kapa(t)*y) )
        res += - self.nb(t) * ( ( 1 - self.x(t)) * np.sinh(self.gam(t) + self.kapb(t)*y)
                                - (y - z) * np.cosh(self.gam(t) + self.kapb(t)*y) )
        if dx_shift and abs(self.x(t)) < 1:           #   reduced: if x is externally controlled/specified
            res += ( (self.Qxx(t) + self.dx(t)) * (y - z)/(1 - self.x(t))
                     - (self.Qxx(t) - self.dx(t))*(y + z)/(1 + self.x(t)) )
        return res
    #   Solving the (coupled, nonlinear) ODEs numerically
    def Solve(self, fp_thrs=1e-5, fp_pts=12, s_mode=None):
        [ysol, zsol] = [ [self.IC()[var]] for var in ("yi","zi")]      #   initiate solution lists (for y & z)
        diffs = []                  #   list of differences between consecutive phase points
        smcount = 0                 #   initiate count of consecutively small differences
        n = 0                       #   initiate number of time steps taken
        (t, tmax) = self.trng()
##        (t, tmax) = ( self._real_tmin, self.trng()[1] )
        nmax = self.maxSteps()
        tstep = self.tstep()
        dfunc0 = lambda y,z,t: (tstep * self.Ky(y,z,t) , tstep * self.Kz(y,z,t))
        if s_mode == None or s_mode == "direct" or s_mode == 1:
            dfunc = dfunc0
        elif s_mode == "runge-kutta" or s_mode == "rk" or s_mode == 4:
            def dfunc(y,z,t):
                (dy1, dz1) = dfunc0(y, z, t)
                (dy2, dz2) = dfunc0(y + dy1/2, z + dz1/2, t + tstep/2)
                (dy3, dz3) = dfunc0(y + dy2/2, z + dz2/2, t + tstep/2)
                (dy4, dz4) = dfunc0(y + dy3, z + dz3, t + tstep)
                dy = (dy1 + 2*dy2 + 2*dy3 + dy4)/6
                dz = (dz1 + 2*dz2 + 2*dz3 + dz4)/6
                return (dy, dz)
        else:
            print("\n" +
                  "ERROR: invalid solution mode '%s' given; choose either 'direct' or 'runge-kutta'." % s_mode
                  + "\n"*2)
        while t < tmax and n < nmax:
            ynow, znow = ysol[n], zsol[n]
##            dely = tstep * self.Ky(ynow,znow,t)
##            delz = tstep * self.Kz(ynow,znow,t)
            (dely, delz) = dfunc(ynow, znow, t)
            ynext , znext = ynow + dely , znow + delz
            ysol.append(ynext)
            zsol.append(znext)
            t += tstep
            n += 1
            # consecutive differences
            diffs.append( [dely , delz] )
            if max([abs(dely), abs(delz)]) < fp_thrs:
                if smcount == 0:
                    smi = n
                smcount += 1
            else:
                smcount = 0
        # threshold for convergence to fixed point
        n_thrs = nmax//fp_pts       #   minimum number of small consecutive diff.s to call "fixed"
        fixed_pt = False            #   initialize fixed point status
        nfixd = 0
        if smcount > n_thrs:
            nfixd = smi             #   minimum index (time-step) for the fixed point
            fixed_pt = True
        (self._sol, self._isfp) = ( (tuple(ysol), tuple(zsol)), (fixed_pt, nfixd) )
        return ( self._sol, self._isfp )
    #   Generic (sub) plot function
    def aPlot(self, fig, x, y, xrng, yrng, title = "", pos = 1, sub = "23", sty = "b-", bg = "w", color=None):
            if pos > int(sub[0])*int(sub[1]):
                print("\n\n\nWARNING: Plot '%s' exceeds subplot range.\n\n")
                return
            ax = fig.add_subplot(int(sub)*10 + pos, facecolor = bg, label = fig.number)
            ax.set_xlim( xrng[0], xrng[1] )
            ax.set_ylim( yrng[0], yrng[1] )
            ax.set_title(title)
            ax.grid()
            myTickParams(ax)
            if color != None:
                ax.plot(x, y, sty, color=color)
            else:
                ax.plot(x, y, sty)
            return ax
    #   Plot the solutions: y,z vs. t for direct solutions; z vs. y for 'phase diagram'; y vs. par(t)
    def solPlot(self, zoom={}, trng=None, show=True, save=None, title=False):
        """
        Plot the solutions y(t), z(t) directly, alongside phase plane (y vs. z).
        Parameter evolution f(t) is also shown (if applicable), with hysteresis plot (y vs. f).

        Attributes/args:
            zoom (dict): Each key correponds to horizontal/vertical axes, e.g. "ty" or "zy", etc.
                         Each value is the pair of ranges for those axes, as: ((xmin,xmax),(ymin,ymax)).
                         One or more subplots may be 'zoomed' in this way,
                             as long as the keys correspond to actual plots displayed for the given case.
            trng (2-tuple): Time-range to be displayed on the plots (a sub-interval of the trng in Soln()).
            show (boolean): Whether or not to show the plot (use False to show multiple plot figures).
            save (none or str): None -> Do not save the plot figure.
                                filename/path -> Save the figure, given a filename (in fixed subdirectory).
            title (boolean): Whether or not to show titles on individual plots (automatically generated).
        """
        if trng == None:
            trng = self.trng(eff=True)
        # styles
        STY1 = "b-"
        STY_pt = "kx"
        HEIGHT = 8
        WIDTH = 15
        BG = "w"
        # function & other shortcut for labeling
        def axlabel(ax, xlbl="", ylbl=""):
            ax.set_xlabel(xlbl)
            ax.set_ylabel(ylbl)
            return
        psyms = {"dlt":r"$\delta$", "gam":r"$\gamma$", "kapa":r"$\kappa_A$", "kapb":r"$\kappa_B$", "x":"x"}
        if title:
            partitle = ", ".join([ psyms[p] + "=%1.2f" % self._pars[p] for p in psyms ])
        else:
            partitle = ""
        # plotting (after defining spaces)
        sol = self.sol()
        tspace = self.tspace()
        tmin_i = 0
        tmax_i = 0
        efftrng = trng
        for i in range(len(tspace)):
            if tspace[i] >= efftrng[0]:
                tmin_i = i
                break
        for i in range(tmin_i,len(tspace)):
            if tspace[i] == efftrng[1]:
                tmax_i = i
                break
            elif tspace[i] > efftrng[1]:
                tmax_i = i - 1
                break
        if tmax_i == 0:
            tmax_i = len(tspace)-1
        tspace = tspace[tmin_i:tmax_i]
        tmin, tmax = min(tspace), max(tspace)
##        yarr, zarr = np.array(sol[0]), np.array(sol[1])
        yarr, zarr = np.array(sol[0][tmin_i:tmax_i]), np.array(sol[1][tmin_i:tmax_i])
        # plot ranges / zoom
        if len(zoom) > 0:
            if ( all([len(list(zoom.values())[i])==2 for i in range(len(zoom))])
                 and all([len(list(zoom.values())[i][j])==2 for i in range(len(zoom)) for j in [0,1]]) ):
                    zcheck = True
            else:
                print("\n\n" +
                      "ERROR: you must enter pairs of ranges ((xmin,xmax),(ymin,ymax))." + "\n"*3)
                zcheck = False
        else:
            zcheck = None
        rngs1246 = ((tmin, tmax), (-1, 1))
        rngs37 = ((-1,1),)*2
        rngs5 = lambda i: ( (1.5*self._pars[self._mode[:i]])*np.array([-1,1]), (-1,1) )
        rngs = {}
        PKeys = (("ty", "tz", "tf", "tx"), ("zy", "xy"), ("fy",))
        Drngs = (rngs1246, rngs37, rngs5)
        for ind in range(len(Drngs)):
            for pkey in PKeys[ind]:
                if zcheck and (pkey in zoom):
                    rngs.update({pkey:zoom[pkey]})
                else:
                    rngs.update({pkey:Drngs[ind]})
        # figure and main plots
        fig = plt.figure("Solutions of ODEs w/ phase portrait", (WIDTH, HEIGHT), facecolor = BG)
        ax1 = self.aPlot(fig, tspace, yarr, rngs["ty"][0], rngs["ty"][1], partitle, sty=STY1, bg=BG, pos=1)
        axlabel(ax1, "t (time)", "y solution (behavior)")
        ax2 = self.aPlot(fig, tspace, zarr, rngs["tz"][0], rngs["tz"][1], partitle, sty=STY1, bg=BG, pos=2)
        axlabel(ax2, "t (time)", "z solution (belonging)")
        ax3 = self.aPlot(fig, zarr, yarr, rngs["zy"][0], rngs["zy"][1], partitle, sty=STY1, bg=BG, pos=3)
        axlabel(ax3, "z solution (belonging)", "y solution (behavior)")
        # time-depedent parameter plots (with hysteresis)
        if self._mode[0] == "x":
            x_arr = self.x(tspace)
            ax6 = self.aPlot(fig, tspace, x_arr, rngs["tx"][0], rngs["tx"][1], partitle, sty=STY1, bg=BG, pos=5)
            axlabel(ax6, "t (time)", "given x fcn. (infection status)")
            ax7 = self.aPlot(fig, x_arr, yarr, rngs["xy"][0], rngs["xy"][1], partitle, sty=STY1, bg=BG, pos=6)
            axlabel(ax7, "given x (infection)", "y solution (behavior)")
        elif self._mode.count("-") == 1:
            i = self._mode.index("-")
            par = {"dlt":(self.dlt, r"$\delta$ (healthy opinion)"),
                     "gam":(self.gam, r"$\gamma$ (infected opinion)"),
                       "kapa":(self.kapa, r"$\kappa_A$ (healthy interaction)"),
                        "kapb":(self.kapb, r"$\kappa_B$ (infected interaction)")}
            par = par[self._mode[:i]]
            arr = par[0](tspace)
            ax4 = self.aPlot(fig, tspace, arr, rngs["tf"][0], rngs["tf"][1], partitle, sty=STY1, bg=BG, pos=4)
            axlabel(ax4, "t (time)", par[1])
            if "fy" in zoom:
                r5 = rngs["fy"]
            else:
                r5 = rngs["fy"](i)
            ax5 = self.aPlot(fig, arr, yarr, r5[0], r5[1],
                          partitle, sty=STY1, bg=BG, pos=5)
            axlabel(ax5, par[1], "y solution (behavior)")
        # fixed point markers
        (isfp, nfixd) = self.isfp()
        if isfp and (tmin_i <= nfixd <= tmax_i):
            ax1.plot([tspace[nfixd]], [yarr[nfixd]], STY_pt)
            ax2.plot([tspace[nfixd]], [zarr[nfixd]], STY_pt)
            ax3.plot([zarr[-1]], [yarr[-1]], STY_pt)
        # show / save
        fig.tight_layout()
        if type(save) == str:
            filepath = os.path.join(os.getcwd(), "socSolver_figs", "solns", save)
            os.makedirs(os.path.dirname(filepath), exist_ok=True)
            if os.path.isfile(filepath):
                print("ERROR: requested file name '%s' already exists is the directory '%s'."
                      % (save, os.path.dirname(filepath)) )
            else:
                plt.savefig(filepath)
        if show:
            plt.show()
            plt.close()
            return
        return fig
    #   Plot the Drift Coefficients (fcn.s of y & z), taking one variable at fixed point or initial point
    def driftPlot(self, y0=0, z0=0, t=0, show=True, save=None):
        """
        For plotting drift coefficients: Ky(y,z0) vs. y  &  Kz(y0,z) vs. z
        """
        # styles
        STY1 = "b-"
        STY_pt = "kx"
        SIZE = 7
        BG = "w"
        # arrange a 2D grid space (here: independent 1D axes)
        y_space = np.linspace(-1.,1.,50)
        z_space = y_space.copy()
        # now use the chosen point to plot drifts
        fig = plt.figure("Drift: behavior & 'belonging'", (SIZE, SIZE), facecolor = BG)
        ax1 = self.aPlot(fig, y_space, self.Ky(y_space,z0,t), (-1, 1), (-1.2, 1.2),
                   "y-drift Ky vs. y (x,z fixed)", sub="21", sty=STY1, bg=BG, pos=1)
        ax1.plot(y0, self.Ky(y0,z0,t), STY_pt)
        ax2 = self.aPlot(fig, z_space, self.Kz(y0,z_space,t), (-1, 1), (-1.5, 1.5),
                   "z-drift Kz vs. z (x,y fixed)", sub="21", sty=STY1, bg=BG, pos=2)
        ax2.plot(z0, self.Kz(y0,z0,t), STY_pt)
        # show / save
        fig.tight_layout()
        if type(save) == str:
            filepath = os.path.join(os.getcwd(), "socSolver_figs", "drifts", save)
            os.makedirs(os.path.dirname(filepath), exist_ok=True)
            if os.path.isfile(filepath):
                print("ERROR: requested file name '%s' already exists is the directory '%s'."
                      % (save, os.path.dirname(filepath)) )
            else:
                plt.savefig(filepath)
            plt.savefig(filepath)
        if show:
            plt.show(fig)
            plt.close()
            return
        return fig
    #   Plot the Drift Coefficients as surfaces in 3D (over both y & z)
    def drift3D(self, t=0, show=True, mode="surf", colorbar=-1, contour=False, save=None):
        """
        Plot drift coeffeicients as (surface) functions of two variables: Ky(y,z)  &  Kz(y,z).
        """
        # styles
        SIZE = 7
        BG = "w"
        CM = cm.viridis
        RStride, CStride = 5, 5
        # arrange a 2D grid space (axes)
        y_space = np.linspace(-1.,1.,50)
        z_space = y_space.copy()
        y,z = np.meshgrid(y_space,z_space)
        # make 3D plots
        fig = plt.figure("Drifts vs. (y,z)", (SIZE+2, SIZE), facecolor = BG)
        for i in (1,2):
            pos = int("12" + str(i))
            ax = fig.add_subplot(pos, projection='3d')
            drift = (self.Ky, self.Kz)[i-1]
            f = lambda var1,var2: drift(var1,var2,t)
            f = f(y,z)
            if mode == "surf":
                surf = ax.plot_surface(y, z, f, rstride=RStride, cstride=CStride, cmap=CM, alpha=0.7)
                if colorbar >= 0:
                    fig.colorbar(surf, shrink=colorbar)
            elif mode == "wire":
                wire = ax.plot_wireframe(y, z, f, rstride=RStride, cstride=CStride)
            if contour:
                cset = ax.contourf(y, z, f, zdir='z',cmap=cm.viridis,offset=-1)
            ax.set_xlabel("y (behavior)")
            ax.set_ylabel("z ('belonging')")
            ax.set_zlabel("drift %s" % ("Ky", "Kz")[i-1])
        # show
        fig.tight_layout()
        if type(save) == str:
            filepath = os.path.join(os.getcwd(), "socSolver_figs", "drift3D", save)
            os.makedirs(os.path.dirname(filepath), exist_ok=True)
            if os.path.isfile(filepath):
                print("ERROR: requested file name '%s' already exists is the directory '%s'."
                      % (save, os.path.dirname(filepath)) )
            else:
                plt.savefig(filepath)
        if show:
            plt.show(fig)
            plt.close()
            return
        return fig
        
#####   #####   #####   #####   #####   (end Soln def.)   #####   #####   #####   #####   #####  


##test = Soln(mode="const", x=0.4, dlt=0.2, gam=-0.2, kapa=1., kapb=1.)
##fig = test.drift3D(show=True, mode="surf",colorbar=0.6, contour=True)
##
##print(3*"\n" + "\t" + str(test.trng()) + "\n\t" + str(test.isfp()) #+ str(test.pars("dlt"))
##      + "\n\t" + str(test.Mode())
##      + "\n\t" + str(test.sol()[0][-5:]) + "\n\t" + str(test.sol()[1][-5:]) + "\n"*2)


#####   Simple Phase Evolution Plotting   #####
def phasearr(ics,kas,xs=0,ds=0,gs=0,kbfrac=1/2):
    res = []
    for ic in ics:
        for ka in kas:
            kb = kbfrac * ka
            a = Soln(trng=(0,40),yzi=ic,nmax=250,mode="const",dlt=ds,gam=gs,kapa=ka,kapb=kb,x=xs)
            res.append((ka, a.sol()[0][-1]))
    return res
def phaseplot(ics=((0.1,0.1),(-0.1,-0.1)), nk=150, kamax=2.5, kbfrac=1/2, d=0, g=0, x=0,
              sty1="-", sty2="-", color1="b", color2="b", show=True, save=None):
    kas = np.linspace(0,kamax,nk)
    parr1 = phasearr((ics[0],),kas,ds=d,gs=g,xs=x,kbfrac=kbfrac)
    parr2 = phasearr((ics[1],),kas,ds=d,gs=g,xs=x,kbfrac=kbfrac)
    y1, y2 = [p[1] for p in parr1], [p[1] for p in parr2]
    fig=plt.figure("Phase as function of "+r"$\kappa_A$",figsize=(5,5))
    ax=fig.add_subplot(111)
    ax.set_xlim(np.amin(kas),np.amax(kas))
    ax.set_ylim(-1,1)
    ax.grid()
    ax.plot(kas,y1,sty1, color=color1)
    line1 = mlines.Line2D([], [], linestyle=sty1, color=color1)
    ax.plot(kas,y2,sty2, color=color2)
    line2 = mlines.Line2D([], [], linestyle=sty2, color=color2)
    ax.set_xlabel(r"$\kappa_A$" + " (interaction)")
    ax.set_ylabel(r"$y^{\ast}$"+ " (behavior)")
    myTickParams(ax)
    if type(save) == str:
        filepath = os.path.join(os.getcwd(), "socSolver_figs", "1par-phases", save)
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        if os.path.isfile(filepath):
            print("ERROR: requested file name '%s' already exists is the directory '%s'."
                  % (save, os.path.dirname(filepath)) )
        else:
            plt.savefig(filepath)
    fig.tight_layout()
    if show:
        plt.show(fig)
        plt.close()
        return
    else:
        return (ax, line1, line2)
                
#####   #####   #####   #####   #####   #####


#####   Class Definition: Phase Diagram of behavior for a given pair of parameters   #####
class Phases:
    """
    Phase diagram for fixed-point solutions of 2-population, 2-opinion parasitic model.

    Attributes/args:
        yzi (2-tuple): Initial conditions, as (yi, zi).
        nmax (int): Maximum number of steps in each solution (from Soln).
        respts (int): Number of elements in each variable list; total_solutions=respts**2.
        makefp (boolean): Whether or not to make f.p. list upon initialization.
        pars (2-tuple of str): Chosen parameter variables (order is unimportant); e.g. ("x","dlt"), ("gam","dlt")
                --> possible pars: ["dlt", "gam", "kapa", "kapb", "x"]
    [the following parameter spec.s supply max. value, if that parameter is one selected for a plot variable]
        dlt, gam (float):          Internal opinions for populations A and B.
        kapa, kapb (float):        Indirect interaction coefficients for populations A and B.
        x (float, in [-1,1]):      Relative population of infected/parastitically-influenced individuals.
    [in-progress]
        parFrac (none or list): None -> Normal operation; each parameter is fixed (unless its chosen in plot).
                                [{"par":frac},"option"] -> The given parameter actually varies, with:
                                                        par=frac*chosen_par ("option" = "mult")
                                                        par=frac-chosen_par ("option" = "diff").
                                                e.g. [{"kapb":0.5},"mult"] -> kapb=0.5*kapa ;
                                                                            kapb changes with kapa.
    """
    def __init__(self, yzi=(0.1,-0.1), nmax=60, respts=10, makefp=True, pars=("x","dlt"),
                 dlt=0, gam=0, kapa=0.5, kapb=0.5, x=0, parFrac=None):
        ## parFrac=None  or  parFrac={"par":frac}  to indicate scaling of given parameter
        ##      e.g. to have kapb=0.5*kapa, set  parFrac={"kapb":0.5}
        self._yzi = yzi
        self._nmax = nmax
        self._respts = respts
        otherpars = ["dlt", "gam", "kapa", "kapb", "x"]
        ks, vs = otherpars, [dlt, gam, kapa, kapb, x]
        kv = [ [ [ks, vs][i][j] for i in range(2) ] for j in range(len(vs)) ]
        self._parvals = { k:v for (k,v) in kv }
        [ otherpars.remove(p) for p in pars ]
        self._fixed = tuple(otherpars)
        # shorthand function for creating solution objects, grab variables for x/y axes on colorplot
        (self._soln, self._pars) = self.chooseSoln(pars, self._parvals, yzi, nmax, parFrac)
        # define parameter spaces to be explored
        if self._pars[0][0] == "k":
            xmin = 0
        else:
            xmin = -1
        if self._pars[1][0] == "k":
            ymin = 0
        else:
            ymin = -1
        self._Xspace = np.linspace(xmin, 1, respts+1)
        self._Yspace = np.linspace(ymin, 1, respts+1)
        self._Xspace *= self._parvals[self._pars[0]]
        self._Yspace *= self._parvals[self._pars[1]]
        if makefp:
            self.makeFPlist()
        else:
            self._fplist = None
    #   create a solution function depending on pair of parameters chosen
    def chooseSoln(self, pars, parvals, yzi, nmax, parFrac=None):
        (Xpar, Ypar) = pars
        combs = list(combinations(["x", "kapb", "kapa", "gam", "dlt"], 2))
        if Xpar in combs[0] and Ypar in combs[0]:
            f = lambda x,kb: Soln(mode="const", yzi=yzi, nmax=nmax, solve=False, x=x, kapb=kb,
                                      kapa=parvals["kapa"], dlt=parvals["dlt"], gam=parvals["gam"])
            return (f, combs[0])
        elif Xpar in combs[1] and Ypar in combs[1]:
            if parFrac==None:
                f = lambda x,ka: Soln(mode="const", yzi=yzi, nmax=nmax, solve=False, x=x, kapa=ka,
                                      kapb=parvals["kapb"], dlt=parvals["dlt"], gam=parvals["gam"])
            elif "kapb" in parFrac.keys():
                f = lambda x,ka: Soln(mode="const", yzi=yzi, nmax=nmax, solve=False, x=x, kapa=ka,
                                      kapb=ka*parFrac["kapb"], dlt=parvals["dlt"], gam=parvals["gam"])
            return (f , combs[1])
        elif Xpar in combs[2] and Ypar in combs[2]:
            f = lambda x,g: Soln(mode="const", yzi=yzi, nmax=nmax, solve=False, x=x, gam=g,
                                      kapa=parvals["kapa"], kapb=parvals["kapb"], dlt=parvals["dlt"])
            return (f, combs[2])
        elif Xpar in combs[3] and Ypar in combs[3]:
            f = lambda x,d: Soln(mode="const", yzi=yzi, nmax=nmax, solve=False, x=x, dlt=d,
                                      kapa=parvals["kapa"], kapb=parvals["kapb"], gam=parvals["gam"])
            return (f, combs[3])
        elif Xpar in combs[4] and Ypar in combs[4]:
            f = lambda kb,ka: Soln(mode="const", yzi=yzi, nmax=nmax, solve=False, kapb=kb, kapa=ka,
                                      x=parvals["x"], dlt=parvals["dlt"], gam=parvals["gam"])
            return (f, combs[4])
        elif Xpar in combs[5] and Ypar in combs[5]:
            f = lambda kb,g: Soln(mode="const", yzi=yzi, nmax=nmax, solve=False, kapb=kb, gam=g,
                                      x=parvals["x"], kapa=parvals["kapa"], dlt=parvals["dlt"])
            return (f, combs[5])
        elif Xpar in combs[6] and Ypar in combs[6]:
            f = lambda kb,d: Soln(mode="const", yzi=yzi, nmax=nmax, solve=False, kapb=kb, dlt=d,
                                      x=parvals["x"], kapa=parvals["kapa"], gam=parvals["gam"])
            return (f, combs[6])
        elif Xpar in combs[7] and Ypar in combs[7]:
            f = lambda ka,g: Soln(mode="const", yzi=yzi, nmax=nmax, solve=False, kapa=ka, gam=g,
                                      x=parvals["x"], kapb=parvals["kapb"], dlt=parvals["dlt"])
            return (f, combs[7])
        elif Xpar in combs[8] and Ypar in combs[8]:
            if parFrac==None:
                f = lambda ka,d: Soln(mode="const", yzi=yzi, nmax=nmax, solve=False, kapa=ka, dlt=d,
                                      x=parvals["x"], kapb=parvals["kapb"], gam=parvals["gam"])
            elif "kapb" in parFrac[0].keys():
                # kapb=parFrac["kapb"]*ka
                if parFrac[1] == "mult":
                    f = lambda ka,d: Soln(mode="const", yzi=yzi, nmax=nmax, solve=False, kapa=ka, dlt=d,
                                      kapb=parFrac[0]["kapb"]*ka, x=parvals["x"], gam=parvals["gam"])
                # kapb=parFrac["kapb"]-ka
                elif parFrac[1] == "diff":
                    f = lambda ka,d: Soln(mode="const", yzi=yzi, nmax=nmax, solve=False, kapa=ka, dlt=d,
                                      kapb=parFrac[0]["kapb"]-ka, x=parvals["x"], gam=parvals["gam"])
                else:
                    print("\n\nError: improper form of 'parFrac'.\n\n")
                    f = None
            return (f, combs[8])
        elif Xpar in combs[9] and Ypar in combs[9]:
            f = lambda g,d: Soln(mode="const", yzi=yzi, nmax=nmax, solve=False, gam=g, dlt=d,
                                      kapa=parvals["kapa"], kapb=parvals["kapb"], x=parvals["x"])
            return (f, combs[9])
        else:
            print("ERROR: given pair ('%s', '%s') is not a valid combination." % self._pars)
            return
    #   generate solution & put f.p. in nested list -- solution object is deleted after storing the f.p.
    def makeFPlist(self, fp_thrs=1e-4, fp_pts=6, fp_warn=False):
        pars = self._pars
        fplist = []
        for ypt in self._Yspace:
            xlst = []
            for xpt in self._Xspace:
                s = self._soln(xpt, ypt)
                s.Solve(fp_thrs=fp_thrs, fp_pts=fp_pts)
                if s.isfp()[0]:
                    xlst.append( s.fp()[0] )
                else:
                    if fp_warn:
                        print(("WARNING: no fixed point found for pars:" +
                               ("(" + pars[0] + "=%1.2f, " + pars[1] + "=%1.2f).") % (xpt, ypt)))
                    xlst.append( s.sol()[0][-1] )
                del s
            fplist.append(xlst)
        self._fplist = tuple(fplist)
        return
    #   give the list of fixed points
    def getFPlist(self):
        return self._fplist
    #   make colormap plots of fixed points (just for y: behavior) for each pair of varying parameters
    def phasePlot(self, interp=None, cmap=cm.plasma, cmpadding=12, zoom=(),
                  show=True, save=None, title=False):
        """
        Plotting the 2D phase diagram of behavior (y) for the given parameters.
        Optional settings: ineterp(olation) type, colormap, padding to colormap, zoom, show, save, plot title
        
        some interpolation choices:
            "nearest"
            "bilinear"
            "gaussian"
            "bessel"
            
        some colormap choices:
            cm.viridis
            cm.plasma
            cm.jet
        """
        pars = self._pars
        # styles
        #STY_pt = "kx"
        HEIGHT = 5
        WIDTH = 5
        BG = "w"
        # labels etc.
        lbls = {"x":"x", "dlt":r"$\delta$", "gam":r"$\gamma$", "kapa":r"$\kappa_A$", "kapb":r"$\kappa_B$"}
        descs = {"x":" (healthy - infected)", "dlt":" (healthy opinion)", "gam":" (infected opinion)",
                 "kapa":" (healthy interaction)", "kapb":" (infected interaction)"}
        Xlbl, Xdesc = lbls[pars[0]], descs[pars[0]]
        Ylbl, Ydesc = lbls[pars[1]], descs[pars[1]]
        # initiate
        others_lbl = [ (s + " = %1.2f, " % self._parvals[s]) for s in self._fixed]
        others_lbl = "".join(others_lbl)[:-2]
        fig = plt.figure( ("Phase Diagram from ODE Solutions: " + others_lbl), (WIDTH, HEIGHT), facecolor=BG)
        # plot range / zoom
        if len(zoom) > 0:
            if ( len(zoom) == 2 ) and ( all([len(zoom[i])==2 for i in range(2)]) ):
                    zcheck = True
            else:
                print("\n\n" +
                      "ERROR: you must enter a pair of ranges ((xmin,xmax),(ymin,ymax))." + "\n"*3)
                zcheck = False
        else:
            zcheck = None
        Drng = ( np.amin(self._Xspace), np.amax(self._Xspace), np.amin(self._Yspace), np.amax(self._Yspace) )
        if zcheck:
            rng = zoom[0][:] + zoom[1][:]
        else:
            rng = Drng
        # plot the given parameter pairing
        y = self._fplist
        ax = fig.add_subplot(111, facecolor = BG, label = self.getFPlist()[-1])
        im = ax.imshow( y, cmap=cmap, origin='lower', interpolation=interp, vmin=-1., vmax=1., extent=rng )
##        divider = make_axes_locatable(ax)
##        cax = divider.append_axes("right", size="5%", pad=0.05)
##        cbar = fig.colorbar(im, cax=cax, format="%1.1f")
##        cax = fig.add_axes([0.8, 0.1, 0.05, 0.8])
##        cbar = fig.colorbar(im, ax=ax, format="%1.1f")
        cbar = plt.colorbar(im, format = "%1.1f", fraction=0.046, pad=0.04)   # magic colorbar numbers!
        cbar.ax.get_yaxis().labelpad = cmpadding
        cbar.ax.set_ylabel("y (behavior)",rotation=270)
        others_lbl = [ (lbls[s] + " = %1.2f, " % self._parvals[s]) for s in self._fixed]
        others_lbl = "".join(others_lbl)[:-2]
        if title:
            ax.set_title(others_lbl)
        ax.set_xlabel(Xlbl + Xdesc)
        ax.set_ylabel(Ylbl + Ydesc)
        myTickParams(ax, dr="inout")
        # show
        fig.tight_layout()
        if type(save) == str:
            filepath = os.path.join(os.getcwd(), "socSolver_figs", "phases", save)
            os.makedirs(os.path.dirname(filepath), exist_ok=True)
            if os.path.isfile(filepath):
                print("ERROR: requested file name '%s' already exists is the directory '%s'."
                      % (save, os.path.dirname(filepath)) )
            else:
                plt.savefig(filepath)
        if show:
            plt.show()
            plt.close()
            return
        return ax
#####   #####   #####   #####   #####   (end Phases def.)   #####   #####   #####   #####   #####


