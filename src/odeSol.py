##  File for solving system of ODEs and plotting solutions (and phase portrait)

from meanValueDefs import (np, plt, Ky, Kz,
                           dlt_const, gam_const, kapa_const, kapb_const, x_const, dlt, gam, kapa, kapb, x, dx)

def odesol(yi,zi,tmin,tmax=25,nmax=120,thrs=10**-5):
    t = tmin        #   starting time
    ### Iterating to get apx. solutions to system of ODEs up to some max time
    #tmax = 30      #   max time for solutions
    #nmax = 120     #   max number of time-steps not to exceed (past t=0)
    delt = tmax/nmax    #   size of time-steps: (delta)t
    #delt = 0.5     #   alternatively: state the size of time-steps explicitly

    tspace = np.linspace(tmin,tmax,((tmax-tmin)//delt)+1)     #   initiated time space
    ysol = [yi]                 #   initiate solution list (for y)
    zsol = [zi]                 #   (for z)
    diffs = []                  #   list of differences between consecutive phase points
    smcount = 0                 #   initiate count of consecutively small differences
    n = 0                       #   initiate number of time steps taken
    while t < tmax and n < nmax:
        ynow, znow = ysol[n], zsol[n]
        dely = delt * Ky(ynow,znow,t)
        delz = delt * Kz(ynow,znow,t)
        ynext, znext = ynow + dely, znow + delz
        ysol.append(ynext)
        zsol.append(znext)
        t += delt
        n += 1
        # consecutive differences
        diffs.append( [dely , delz] )
        if max([abs(dely), abs(delz)]) < thrs:
            if smcount == 0:
                smi = n
            smcount += 1
        else:
            smcount = 0

    ### Threshold for convergence to fixed point
    n_thrs = nmax//12            #   minimum number of small consecutive diff.s to call "fixed"
    fixed_pt = False            #   initialize fixed point status
    if smcount > n_thrs:
        nfixd = smi
        fixed_pt = True

    ### Plot the solutions: y,z vs. t for direct solutions; z vs. y for 'phase diagram'
    # solution arrays
    yarr, zarr = np.array(ysol), np.array(zsol)
    # styles
    STY1 = "b-"
    STY_pt = "kx"
    HEIGHT = 8
    WIDTH = 15
    BG = "w"
    # plotting
    fig = plt.figure("Solutions of ODEs w/ phase diagram", (WIDTH, HEIGHT), facecolor = BG)
    def odePlot(x, y, xrng, yrng, title = "", pos = 1, sub = "23", sty = STY1):
        if pos > int(sub[0])*int(sub[1]):
            print("\n\n\nWARNING: Plot '%s' exceeds subplot range.\n\n")
            return
        ax = fig.add_subplot(int(sub)*10 + pos, facecolor = BG)
        ax.set_xlim( xrng[0], xrng[1] )
        ax.set_ylim( yrng[0], yrng[1] )
        ax.set_title(title)
        ax.grid()
        ax.plot(x, y, sty)
        return ax
    ax1 = odePlot(tspace, yarr, (tmin, tmax), (-1, 1), "y solution (behavior) vs. t", 1)
    ax2 = odePlot(tspace, zarr, (tmin, tmax), (-1, 1), "z solution (belonging) vs. t", 2)
    ax3 = odePlot(zarr, yarr, (-1, 1), (-1, 1), "y solution (behavior) vs. z solution (belonging)", 3)
    dlt_arr = dlt(tspace)
    x_arr = x(tspace)
    if type(dlt_arr) == type(tspace):   # and type(x_arr) != type(tspace):
        ax4 = odePlot(tspace, dlt_arr, (tmin, tmax), (-1, 1), r"$\delta$ (base opinion) vs. t", 4)
        ax5 = odePlot(dlt_arr, yarr, min(1.5*dlt_const,1)*np.array([-1,1]), (-1,1),
                      r"y solution (behavior) vs. $\delta$ (opinion)", 5)
    
    if type(x_arr) != type(tspace):
        x_arr = x_arr*np.ones(len(tspace))
    if type(dlt_arr) != type(tspace):
        ax6 = odePlot(tspace, x_arr, (tmin, tmax), (-1, 1), "given x fcn. (infection status) vs. t", 5)
        ax7 = odePlot(x_arr, yarr, (-1, 1), (-1, 1), "y solution (behavior) vs. given x (infection)", 6)
    
    if fixed_pt:
        ax1.plot([tspace[nfixd]], [ysol[nfixd]], STY_pt)
        ax2.plot([tspace[nfixd]], [zsol[nfixd]], STY_pt)
        ax3.plot([zsol[-1]], [ysol[-1]], STY_pt)
        res = (tspace, ysol, zsol, (ysol[-1],zsol[-1]))
    else:
        res = (tspace, yarr, zarr, ())
        

    fig.tight_layout()
    plt.show(fig)

    return res
