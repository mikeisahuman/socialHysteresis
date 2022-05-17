##  File for plotting drift coefficients vs. both variables (3D surface/wire plots)

from meanValueDefs import np, plt, axes3d, cm, Ky, Kz

def drift3D(t):
    SIZE = 7
    BG = "w"
    fig = plt.figure("Drifts vs. (y,z)", (SIZE+2, SIZE), facecolor = BG)
    #   Arrange a 3D grid space (axes)
    ##x = np.linspace(-1.,1.,50)
    y_space = np.linspace(-1.,1.,50)
    z_space = y_space.copy()
    y,z = np.meshgrid(y_space,z_space)
    for i in (1,2):
        pos = "12" + str(i)
        pos = int(pos)
        ax = fig.add_subplot(pos, projection='3d')
        #ax = axes3d(fig)
        drift = (Ky, Kz)[i-1]
        f = lambda var1,var2: drift(var1,var2,t)
        f = f(y,z)
        #ax.plot_wireframe(y, z, f, rstride=5, cstride=5)
        surf = ax.plot_surface(y, z, f, rstride=4, cstride=4,cmap=cm.viridis,alpha=0.7)
        #fig.colorbar(surf,shrink=0.5)
        #cset = ax.contourf(y, z, f, zdir='z',cmap=cm.viridis,offset=-1)
        ax.set_xlabel('y (opinion)')
        ax.set_ylabel('z ("belonging")')
        ax.set_zlabel('drift %s' % ("Ky", "Kz")[i-1])
    fig.tight_layout()
    plt.show(fig)
