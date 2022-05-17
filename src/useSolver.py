##  Using ODE solver objects to obtain phase diagrams, etc.
##  by Mike Phillips, 7/13/2018
##
##  With a useful ODE solver for our social dynamics equations in hand (see: socSolver.py),
##      we seek to apply it to construct phase diagrams, etc. -- using long-term solutions for various pars.

import socSolver as soc
#import itertools as it
from socSolver import np, sy, plt, cm

# shorthand function for creating solution objects -- here: only x, dlt, gam are allowed to vary
yzi = (0.1,-0.1)
kapa, kapb = 0.75, 0.10
def solnXDG(x=0, d=0.1, g=0.1):
    return soc.Soln(yzi=yzi, kapa=kapa, kapb=kapb, x=x, dlt=d, gam=g, mode="const", nmax=50, solve=False)

# define parameter spaces to be explored
respts = 10
xspace = np.linspace(-1, 1, respts+1)
dgmax = 0.8
dspace = np.linspace(-dgmax, dgmax, respts+1)
gspace = dspace.copy()

# generate solution & put fixed point in nested list -- solution object is deleted after storing the f.p.
fplist = []
for x in xspace:
    dlst = []
    for d in dspace:
        glst = []
        for g in gspace:
            pars = (x, d, g)
            s = solnXDG(x, d, g)
            s.Solve(fp_thrs=1e-3, fp_pts=6)
            if s.isfp()[0]:
                glst.append( s.fp()[0] )
            else:
                print("WARNING: no fixed point found for pars: (x=%1.2f, dlt=%1.2f, gam=%1.2f)." % pars)
                glst.append( s.sol()[0][-1] )
            if x != xspace[-1]:
                del s
        dlst.append(glst)
    print("x = %1.2f \t done" % x)
    fplist.append(dlst)

# make colormap plots of fixed points (just for y: behavior) for each pair of varying parameters
#STY1 = "b-"
#STY_pt = "kx"
HEIGHT = 8
WIDTH = 15
BG = "w"

x_lbl, x_desc = "x", ": infection"
d_lbl, d_desc = r"$\delta$", ": healthy opinion"
g_lbl, g_desc = r"$\gamma$", ": infected opinion"
x_ind, d_ind, g_ind = respts//4, respts//2, 3*respts//4

Xlbl, Xdesc = (x_lbl, x_lbl, g_lbl), (x_desc, x_desc, g_desc)
Ylbl, Ydesc = (d_lbl, g_lbl, d_lbl), (d_desc, g_desc, d_desc)
fixlbl, fixdesc, fixind = (g_lbl, d_lbl, x_lbl), (g_desc, d_desc, x_desc), (g_ind, d_ind, x_ind)
spaces = {x_lbl:xspace, d_lbl:dspace, g_lbl:gspace}
pos = 0

fig = plt.figure("Phase Diagram for Colletive Behavior from ODE Solutions", (WIDTH, HEIGHT), facecolor = BG)

for i_pair in range(len(Xlbl)):
    if i_pair == 0:
        y = [ [ fplist[i][j][g_ind] for i in range(respts+1) ] for j in range(respts+1) ]
    elif i_pair == 1:
        y = [ [ fplist[i][d_ind][j] for i in range(respts+1) ] for j in range(respts+1) ]
    elif i_pair == 2:
        y = [ [ fplist[x_ind][j][i] for i in range(respts+1) ] for j in range(respts+1) ]
    pos += 1
    ax = fig.add_subplot(220 + pos, facecolor = BG)
    ax = ax.imshow( y, cmap=cm.viridis, origin='lower', interpolation="bilinear",
                   extent=( np.amin(spaces[Xlbl[i_pair]]), np.amax(spaces[Xlbl[i_pair]]),
                           np.amin(spaces[Ylbl[i_pair]]), np.amax(spaces[Ylbl[i_pair]]) ) )
    cbar = fig.colorbar(ax)
    #cbar.ax.get_yaxis().labelpad=15
    cbar.ax.set_ylabel("y: Behavior",rotation=270)
    plt.title(g_lbl + (" = %1.2f " % spaces[fixlbl[i_pair]][fixind[i_pair]]) + fixdesc[i_pair])
    plt.xlabel(Xlbl[i_pair] + Xdesc[i_pair])
    plt.ylabel(Ylbl[i_pair] + Ydesc[i_pair])

fig.tight_layout()
plt.show(fig)




#######   Class Definition: Phase Diagrams of behavior for some parameters   #####
##class PhasesXDG:
##    def __init__(self, yzi=(0.1,-0.1), kapa=0.75, kapb=0.10, nmax=60, respts=10, dgmax=1.0, makefp=True):
##        self._yzi = yzi
##        self._kapa, self._kapb = kapa, kapb
##        self._nmax = nmax
##        self._respts = respts
##        # shorthand function for creating solution objects -- here: only x, dlt, gam are allowed to vary
##        self.solnXDG = lambda x,d,g: Soln(mode="const", yzi=yzi, kapa=kapa, kapb=kapb, nmax=nmax,
##                            solve=False, x=x, dlt=d, gam=g)
##        # define parameter spaces to be explored
##        self._xspace = np.linspace(-1, 1, respts+1)
##        self._dspace = np.linspace(-dgmax, dgmax, respts+1)
##        self._gspace = self._dspace.copy()
##        if makefp:
##            self.FPlist()
##        else:
##            self._fptupl = None
##    #   generate solution & put f.p. in nested list -- solution object is deleted after storing the f.p.
##    def FPlist(self, fp_thrs=1e-4, fp_pts=6, fp_warn=False):
##        fplist = []
##        for x in self._xspace:
##            dlst = []
##            for d in self._dspace:
##                glst = []
##                for g in self._gspace:
##                    pars = (x, d, g)
##                    s = self.solnXDG(x, d, g)
##                    s.Solve(fp_thrs=fp_thrs, fp_pts=fp_pts)
##                    if s.isfp()[0]:
##                        glst.append( s.fp()[0] )
##                    else:
##                        if fp_warn:
##                            print(("WARNING: no fixed point found for pars:" +
##                                   "(x=%1.2f, dlt=%1.2f, gam=%1.2f)." % pars))
##                        glst.append( s.sol()[0][-1] )
##                    del s
##                dlst.append(glst)
##            #print("x = %1.2f \t done" % x)
##            fplist.append(dlst)
##        self._fptupl = tuple(fplist)
##        return
##    #   make colormap plots of fixed points (just for y: behavior) for each pair of varying parameters
##    def phasePlot(self, interp="bilinear", wantfp=False, ind_fracs=(1/4, 1/2, 3/4), cmpad=10):
##        # styles
##        #STY_pt = "kx"
##        HEIGHT = 8
##        WIDTH = 15
##        BG = "w"
##        # labels etc.
##        x_lbl, x_desc = "x", ": infection"
##        d_lbl, d_desc = r"$\delta$", ": healthy opinion"
##        g_lbl, g_desc = r"$\gamma$", ": infected opinion"
##        (x_ind, d_ind, g_ind) = ( round(self._respts*frac) for frac in ind_fracs )
##        Xlbl, Xdesc = (x_lbl, x_lbl, g_lbl), (x_desc, x_desc, g_desc)
##        Ylbl, Ydesc = (d_lbl, g_lbl, d_lbl), (d_desc, g_desc, d_desc)
##        fixlbl, fixdesc, fixind = (g_lbl, d_lbl, x_lbl), (g_desc, d_desc, x_desc), (g_ind, d_ind, x_ind)
##        spaces = {x_lbl:self._xspace, d_lbl:self._dspace, g_lbl:self._gspace}
##        # initiate
##        pos = 0
##        fig = plt.figure( ("Phase Diagram for Colletive Behavior from ODE Solutions: "
##                           + ( (r"kapA = %1.2f, kapB = %1.2f") % (self._kapa, self._kapb) ) ),
##                          (WIDTH, HEIGHT), facecolor=BG)
##        # plot all parameter pairings
##        for i_pair in range(len(Xlbl)):
##            if i_pair == 0:
##                y = [ [ self._fptupl[i][j][g_ind] for i in range(len(spaces[x_lbl])) ]
##                      for j in range(len(spaces[d_lbl])) ]
##            elif i_pair == 1:
##                y = [ [ self._fptupl[i][d_ind][j] for i in range(len(spaces[x_lbl])) ]
##                      for j in range(len(spaces[g_lbl])) ]
##            elif i_pair == 2:
##                y = [ [ self._fptupl[x_ind][j][i] for i in range(len(spaces[g_lbl])) ]
##                      for j in range(len(spaces[d_lbl])) ]
##            pos += 1
##            ax = fig.add_subplot(220 + pos, facecolor = BG)
##            ax = ax.imshow( y, cmap=cm.viridis, origin='lower', interpolation=interp,
##                           extent=( np.amin(spaces[Xlbl[i_pair]]), np.amax(spaces[Xlbl[i_pair]]),
##                                   np.amin(spaces[Ylbl[i_pair]]), np.amax(spaces[Ylbl[i_pair]]) ) )
##            cbar = fig.colorbar(ax)
##            cbar.ax.get_yaxis().labelpad = cmpad
##            cbar.ax.set_ylabel("y: Behavior",rotation=270)
##            plt.title(fixlbl[i_pair] + (" = %1.2f " % spaces[fixlbl[i_pair]][fixind[i_pair]])
##                      + fixdesc[i_pair])
##            plt.xlabel(Xlbl[i_pair] + Xdesc[i_pair])
##            plt.ylabel(Ylbl[i_pair] + Ydesc[i_pair])
##        # show
##        fig.tight_layout()
##        plt.show(fig)
##        if wantfp:
##            return self._fptupl
##        else:
##            return
#######   #####   #####   #####   #####   (end PhasesXDG def.)   #####   #####   #####   #####   #####
