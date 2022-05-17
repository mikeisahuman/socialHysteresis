##  File for plotting drift coefficients vs. its variable alone (the other at fixed pt., if possible)

from meanValueDefs import np, plt, Ky, Kz

def driftplots(y0,z0,t):
    #   Arrange a 3D grid space (axes)
    ##x = np.linspace(-1.,1.,50)
    y_space = np.linspace(-1.,1.,50)
    z_space = y_space.copy()

    ### plot x-drift (guides mean-value eq.)
    STY1 = "b-"
    STY_pt = "kx"
    SIZE = 7
    BG = "w"
    ##fig1 = plt.figure("Drift: isolated infection eq.", (SIZE, SIZE), facecolor = BG)
    ##ax1 = fig1.add_subplot(111, facecolor = BG, aspect = "equal")
    ##ax1.set_xlim( -1, 1 )
    ##ax1.set_ylim( -1, 1 )
    ##ax1.plot(x, Kx(x), STY1)
    ##ax1.plot(x0, Kx(x0), STY_pt)
    ##plt.show(fig1)
        
    # now use the x-equil. to plot the other drifts
    fig = plt.figure("Drift: attitude & 'belonging'", (SIZE, SIZE), facecolor = BG)
    ax21 = fig.add_subplot(211, facecolor = BG)
    ax21.set_xlim( -1, 1 )
    ax21.set_ylim( -1, 1 )
    ax21.set_title("y-drift Ky vs. y (x,z fixed)")  # at equil.)")
    ax21.plot(y_space, Ky(y_space,z0,t), STY1)
    ax21.plot(y0, Ky(y0,z0,t), STY_pt)
    ax22 = fig.add_subplot(212, facecolor = BG)
    ax22.set_xlim( -1, 1 )
    ax22.set_ylim( -1, 1 )
    ax22.set_title("z-drift Kz vs. z (x,y fixed)")  #  at equil.)")
    ax22.plot(z_space, Kz(y0,z_space,t), STY1)
    ax22.plot(z0, Ky(y0,z0,t), STY_pt)

    fig.tight_layout()
    plt.show(fig)

    return
