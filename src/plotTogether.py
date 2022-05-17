##  Using ODE solver objects to obtain combined plots
##  by Mike Phillips, 8/14/2018
##
##  With a useful ODE solver for our social dynamics equations in hand (see: socSolver.py),
##      we seek to apply it to construct combined plots with good colors, legend, etc.

import socSolver as soc
from socSolver import np, plt

select = "kablrg"

arrows = True

colors, stys = ["#0000FF", "#CC6600", "#336633"], ["-","--","-."]
##colors, stys = ["#CC6600", "#0000FF", "#336633"], ["--","-","-."]

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

if select == "trans":
    ##  Phase transition curve
    d = np.linspace(0.8,0,3)
    g = -d
    kamax = 2.5
    kbfrac = 1/2
    x = 0
    ics = ((0.5,0.5),(-0.5,-0.5))
    nk = 150
    axs = []
    for i in range(len(d)):
        axs.append( soc.phaseplot(ics=ics, nk=nk, kamax=kamax, kbfrac=kbfrac, d=d[i], g=g[i], x=x,
                                  sty1=stys[-1-i], sty2=stys[-1-i],
                                  color1=colors[-1-i], color2=colors[-1-i], show=False) )
    leglbl = [ (r"$\delta$" + "=%1.2f=-" + r"$\gamma$") % di for di in d]
    leglbl.reverse()
    lines = [ a[1] for a in axs ]
    lines.reverse()
    axs[-1][0].legend(handles=lines, labels=leglbl, loc="upper left", frameon=True, edgecolor="inherit")
    plt.tight_layout()
    plt.show()
    plt.close()
elif select == "dosc":
    ##  Oscillating opinion dlt
    d = 0.8
    g = -0.4
    ka, kb = 0.7, 0.7
    x = [0.4, 0,-0.4]
    sols = []
    for i in range(len(x)):
        sols.append( soc.Soln(trng=(0,50),yzi=(0,0),nmax=300,mode="dlt-osc",
                              dlt=d,gam=g,kapa=ka,kapb=kb,x=x[i],omg_frac=0.25,phase=-np.pi/2) ) 
    tsp = sols[0].tspace(6,31.5)
    imin,imax = list(sols[0].tspace()).index(tsp[0]), list(sols[0].tspace()).index(tsp[-1])+1
    fig = plt.figure("dlt hyst: cases 2-3", (6,6))
    axs = []
    for i in range(len(sols)):
        axs.append( sols[i].aPlot(fig, sols[i].dlt(tsp), sols[i].ysol()[imin:imax],
                                  (-0.85,0.85), (-1,1), sub="11", sty=stys[i], color=colors[i]) )
    axs[-1].set_xlabel(r"$\delta$ (healthy opinion)")
    axs[-1].set_ylabel("y solution (behavior)")
    xlbl = list(zip(["x="]*len(x),[str(round(k,2)) for k in x] ))
    xlbl = [ "".join(l) for l in xlbl ]
    axs[-1].legend(xlbl,loc="upper left",frameon=True,edgecolor="inherit")
##    if arrows:
        
    fig.tight_layout()
    plt.show()
    plt.close()
else:
    if select == "kabmed":
        ##  Medium ka+kb (oscillating x)
        d = 0
        g = -0.4
        x = 0.6
        ka = np.linspace(0.2,1.2,3)
        kb = 1.4-ka
        sols = []
        for i in range(len(ka)):
            sols.append( soc.Soln(trng=(0,50),yzi=(0,0),nmax=300,mode="x-osc",
                                  dlt=d,gam=g,kapa=ka[i],kapb=kb[i],x=x,omg_frac=0.25,phase=-np.pi/2) )
        tsp = sols[0].tspace(6,31.5)
        imin,imax = list(sols[0].tspace()).index(tsp[0]), list(sols[0].tspace()).index(tsp[-1])+1
        fig = plt.figure("x hyst: cases 4-6", (6,6))
        axs = []
        for i in range(len(sols)):
            axs.append( sols[i].aPlot(fig, sols[i].x(tsp), sols[i].ysol()[imin:imax],
                                      (-0.65, 0.65), (-1,0.2), sub="11", sty=stys[i], color=colors[i]) )
    elif select == "kablrg":
        ##  Large ka+kb (oscillating x)
        d = 0.4
        g = -0.4
        x = 0.6
        ka = np.linspace(1.1,2.0,3)
        kb = 2.2-ka
        sols = []
        for i in range(len(ka)):
            sols.append( soc.Soln(trng=(0,50),yzi=(0,0),nmax=300,mode="x-osc",
                                  dlt=d,gam=g,kapa=ka[i],kapb=kb[i],x=x,omg_frac=0.25,phase=-np.pi/2) )	
        tsp = sols[0].tspace(6,31.5)
        imin,imax = list(sols[0].tspace()).index(tsp[0]), list(sols[0].tspace()).index(tsp[-1])+1
        tsp = [tsp]*2 + [sols[-1].tspace(0,45)]
        imin = [imin]*2 + [list(sols[-1].tspace()).index(tsp[2][0])]
        imax = [imax]*2 + [list(sols[-1].tspace()).index(tsp[2][-1])+1]
        fig = plt.figure("x hyst: cases 7-9", (6,6))
        axs = []
        for i in range(len(sols)):
            axs.append( sols[i].aPlot(fig, sols[i].x(tsp[i]), sols[i].ysol()[imin[i]:imax[i]],
                                      (-0.65, 0.65), (-1,1), sub="11", sty=stys[i], color=colors[i]) )
    ##  Final Touches & Show
    kalbl = list(zip([r"$\kappa_A$="]*len(ka),[str(round(k,2)) for k in ka] ))
    kalbl = [ "".join(l) for l in kalbl ]
    kblbl = list(zip([r"$\kappa_B$="]*len(kb),[str(round(k,2)) for k in kb] ))
    kblbl = [ "".join(l) for l in kblbl ]
    leglbl = [ ", ".join([kalbl[i], kblbl[i]]) for i in range(len(ka)) ]
    axs[-1].set_xlabel("x (relative population)")
    axs[-1].set_ylabel("y solution (behavior)")
    axs[-1].legend(leglbl,loc="upper left",frameon=True,edgecolor="inherit")
    fig.tight_layout()
    plt.show()
    plt.close()

plt.close()
