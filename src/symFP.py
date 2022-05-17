### (symbolic) / analytical computation for approximate fixed points

from meanValueDefs import (np, na, nb, dlt, gam, kapa, kapb, x, dx)

def SYMFP(t=0, wantprint=True, thrs=1e-10):
    nas = (na(t)*np.sinh(dlt(t)))
    nac = (na(t)*np.cosh(dlt(t)))
    nbs = (nb(t)*np.sinh(gam(t)))
    nbc = (nb(t)*np.cosh(gam(t)))
    ka , kb , X = kapa(t) , kapb(t) , x(t)
    if abs(1-X**2) > thrs:
        chi = 2*dx(t)/(1-(x(t)**2))
    else:
        chi = 0

    ##       ka:kapa(t), kb:kapb(t), X:x(t), chi:(2*dx(t)/(1-x(t)**2))

    ##import sympy as spy
    ##nas, nac, nbs, nbc, ka, kb = spy.symbols("nas nac nbs nbc ka kb", positive = True)
    ##X, chi, y, z = spy.symbols("X chi y z", real = True)

    #   y-pars
    As = nas - nbs
    Bs = nas + nbs
    Ec = nac*(ka*(1+X)-1) + nbc*(kb*(1-X)-1)
    Ac = nac - nbc
    Gs = nas*ka + nbs*kb
    Hs = nas*ka*(ka/2-1) - nbs*kb*(kb/2-1)

    #   z-pars (not already defined)
    Bc = nac + nbc
    Dc = nac*(ka*(1+X)-1) - nbc*(kb*(1-X)-1)
    Ks = nas*ka*(ka/2-1) + nbs*kb*(kb/2-1)
    Fs = nas*ka - nbs*kb

    #   quad. eq. : A y^2 + B y + C = 0
    A = Dc*Gs + Ac*Ks - Bc*Hs - Ec*Fs + (Gs - Hs)*chi
    B = Gs*Bs - Fs*As + Ac*Dc - Ec*Bc + (Gs*As - Fs*Bs)*X + (Ac - Ec)*chi
    C = Ac*Bs - As*(Bc + chi) + (Ac*As - Bs*(Bc + chi))*X

    #   fixed points
    if abs(A) > thrs:
        surd = B**2 - 4*A*C
        fpy = [ (-B + n * (surd**(0.5))) / (2*A) for n in (-1,1) ]
    elif abs(B) > thrs:
        fpy = [ -C / B ]
    else:
        return None
    fpz = [0]*len(fpy)
    ##fpy1 = (-B + spy.sqrt(surd))/(2*A)
    ##fpy2 = (-B - spy.sqrt(surd))/(2*A)
    ##fpy = [spy.simplify(fpy1), spy.simplify(fpy2)]
    ##fpz = [z,z]
    j = 0
    for y0 in fpy:
        if abs(Ac + Gs*y0) > thrs:      ##spy.simplify(Ac + Gs*y0) != 0:
            fpz[j] = (As + Bs*X + Ec*y0 + Hs*(y0**2)) / (Ac + Gs*y0)
            j += 1
            continue
        elif abs(Bc + chi + Fs*y0) > thrs:  ##spy.simplify(Bc + chi + Fs*y0) != 0:
            fpz[j] = (Bs + As*X + (Dc + chi)*y0 + Ks*(y0**2)) / (Bc + chi + Fs*y0)
            j += 1
            continue
        else:
            fpz[j] = 0
            j += 1
            continue
    ##fpz = [ spy.simplify(ex) for ex in fpz ]


    #   numbers with given parameters
    ##t = 0
    ##sub = {nas:(na(t)*spy.sinh(dlt(t))), nac:(na(t)*spy.cosh(dlt(t))),
    ##       nbs:(nb(t)*spy.sinh(gam(t))), nbc:(nb(t)*spy.cosh(gam(t))),
    ##       ka:kapa(t), kb:kapb(t), X:x(t), chi:(2*dx(t)/(1-x(t)**2))}
    ##ynum = [ ex.subs(sub) for ex in fpy ]
    ##znum = [ ex.subs(sub) for ex in fpz ]

    ynum, znum = [], []
    sol = []
    ##jac = []
    det = lambda mat: (mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0])
    tr = lambda mat: (mat[0][0] + mat[1][1])

    nullfunc = lambda x: ()
    if wantprint:
        pfunc = print
    else:
        pfunc = nullfunc
    pfunc("\n")
    for j in range(len(fpy)):
        jac = []
        if -1 <= fpy[j] <= 1 and -1 <= fpz[j] <= 1:
            ynum.append(fpy[j])
            znum.append(fpz[j])
            jac.append([Ec - Gs*znum[-1] + 2*Hs*ynum[-1], - Ac - Gs*ynum[-1]])
            jac.append([Dc + chi - Fs*znum[-1] + 2*Ks*ynum[-1], - Bc - chi - Fs*ynum[-1]])
            trjac = tr(jac)
            detjac = det(jac)
            if detjac >= 0:
                if abs(trjac) < thrs:
                    stabl = "center"
                elif trjac**2 < 4*detjac:
                    stabl = ["stable","unstable"][int(trjac>0)] + " spiral"
                elif trjac**2 >= 4*detjac:
                    stabl = ["stable","unstable"][int(trjac>0)] + " node"
            else:
                stabl = "unstable saddle"
            pfunc( ("\n\n\tysol%i : %1.3f" % (j+1, ynum[-1])) + ("\n\tzsol%i : %1.3f" % (j+1, znum[-1])) )
            pfunc("\tclassification : %s" % stabl)
            sol.append(((ynum[-1],znum[-1]), stabl))
        if j == 1 and len(ynum) == 0:
            pfunc("\n\t" + "No fixed points found.")
    pfunc("\n"*2)

    return tuple(sol)
