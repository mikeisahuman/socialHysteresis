##  File for "phase curves", i.e. fixed points vs. parameters

from meanValueDefs import np, optimize, plt, Ky, Kz, kapa

def phasecurves(t):
    #   The loop to make a phase diagram
    ##phase_list = []
    y_list = [] ;   z_list = []
    gam_list = np.linspace(-1.,1.,10)
    kapb_list = np.linspace(0.,kapa(t),10)
    for gam_const in gam_list:
        y_kap = [] ;    z_kap = []
        for kapb_const in kapb_list:
    ##        # find the equil. point (x eq. is self-contained)
    ##        x0 = optimize.brentq(Kx, -1., 1.)

            # use x-equil to optimize both equations for equilibrium y, z
            def Kyz(vecYZ):
                return [ Ky(vecYZ[0],vecYZ[1],t) , Kz(vecYZ[0],vecYZ[1],t) ]
            resYZ = optimize.root(Kyz, (0.,0.))
            [y0, z0] = resYZ.x
            if not resYZ.success:
                print("\n\tWARNING!\tThe root was not found successfully...\n\n")
                y_kap.append(0)
                z_kap.append(0)
                break
            else:
                y_kap.append(y0)
                z_kap.append(z0)
            ##        phase = (x0, y0, z0)
        if len(y_kap) > 0:
            y_list.append(y_kap)
            z_list.append(z_kap)
    ##        phase_list.append(phase)

    STY = "b-"
    STY_pt = "ko"
    SIZE = 8
    BG = "w"

    gam_ind = 3
    kap_ind = 2
    y_gam = [row[kap_ind] for row in y_list]
    z_gam = [row[kap_ind] for row in z_list]

    fig = plt.figure("Phase Evolution: fixed pts vary with pars", (SIZE, SIZE), facecolor = BG)
    ax1 = fig.add_subplot(221, facecolor = BG)
    ax1.set_xlim( min(kapb_list), max(kapb_list) )
    ax1.set_ylim( -1, 1 )
    ax1.set_title("y0 vs. kapb: gam = %3.3f" % gam_list[gam_ind])
    ax1.set_xlabel("kapb (interaction)")
    ax1.set_ylabel("y0 (attitude)")
    ax1.plot(kapb_list, y_list[gam_ind], STY)
    ax2 = fig.add_subplot(222, facecolor = BG)
    ax2.set_xlim( min(kapb_list), max(kapb_list) )
    ax2.set_ylim( -1, 1 )
    ax2.set_title("z0 vs. kapb: gam = %3.3f" % gam_list[gam_ind])
    ax2.set_xlabel("kapb (interaction)")
    ax2.set_ylabel("z0 (agreement)")
    ax2.plot(kapb_list, z_list[gam_ind], STY)
    ax3 = fig.add_subplot(223, facecolor = BG, aspect = "equal")
    ax3.set_xlim( min(gam_list), max(gam_list) )
    ax3.set_ylim( -1, 1 )
    ax3.set_title("y0 vs. gam: kapb = %3.3f" % kapb_list[kap_ind])
    ax3.set_xlabel("gamma (neg. shift)")
    ax3.set_ylabel("y0 (attitude)")
    ax3.plot(gam_list, y_gam, STY)
    ax4 = fig.add_subplot(224, facecolor = BG, aspect = "equal")
    ax4.set_xlim( min(gam_list), max(gam_list) )
    ax4.set_ylim( -1, 1 )
    ax4.set_title("z0 vs. gam: kapb = %3.3f" % kapb_list[kap_ind])
    ax4.set_xlabel("gamma (neg. shift)")
    ax4.set_ylabel("z0 (agreement)")
    ax4.plot(gam_list, z_gam, STY)

    fig.tight_layout()
    plt.show()

    return [y0,z0]
