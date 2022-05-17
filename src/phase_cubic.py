##  Phase Diagram -- approx. from cubic f.p. eq.
##  by Mike Phillips, 9/24/2019
##
##  The transcendental fixed point equation is approximated (near y=0) as a cubic.
##  Phase transition curves are found by plotting the points where the cubic discriminant
##  is (approximately) zero; this is where a new root of the cubic equation emerges (or vanishes).

import numpy as np, matplotlib.pyplot as plt

# necessary setup pars
thrs = 1e-3     #   threshold for calling discriminant == 0
Npts = 250       #   number of points to use on each axis (total pts = Npts^2)
(xvar, yvar) = ("kb", "g")  #   chosen (x,y) variables for the phase diagram

minmax = {"x":(-1,1), "d":(-1,1), "g":(-1,1), "ka":(0,2.5), "kb":(0,2.5)}   #   min & max values dictionary
myvals = {"x":0, "d":0, "g":0, "ka":1, "kb":1}  #   values dictionary, when not chosen for the plot

# hyperbolic function shortcuts
tanh = np.tanh
sech = lambda x: 1/np.cosh(x)

# coefficients of cubic: A*x^3 + B*x^2 + C*x + D = 0
def A(x=0, d=0, g=0, ka=1, kb=1):
    res = (ka**3) * (2*(sech(d)**2)*(tanh(d))**2 + (sech(d))**4) * (1+x) / 3
    res += (kb**3) * (2*(sech(g)**2)*(tanh(g))**2 + (sech(g))**4) * (1-x) / 3
    return res

def B(x=0, d=0, g=0, ka=1, kb=1):
    res = - (ka**2) * (sech(d)**2) * tanh(d) * (1+x)
    res += - (kb**2) * (sech(g)**2) * tanh(g) * (1-x)
    return res

def C(x=0, d=0, g=0, ka=1, kb=1):
    res = ka * (sech(d)**2) * (1+x)
    res += kb * (sech(g)**2) * (1-x) - 2
    return res

def D(x=0, d=0, g=0, ka=1, kb=1):
    res = tanh(d) * (1+x)
    res += tanh(g) * (1-x)
    return res

# discriminant function
def disc(x=0, d=0, g=0, ka=1, kb=1):
    fA = A(x, d, g, ka, kb)
    fB = B(x, d, g, ka, kb)
    fC = C(x, d, g, ka, kb)
    fD = D(x, d, g, ka, kb)
    res = 18*fA*fB*fC*fD
    res += -4*(fB**3)*fD
    res += (fB**2)*(fC**2)
    res += -4*fA*(fC**3)
    res += -27*(fA**2)*(fD**2)
    return res

# lists to use (and function)
xminmax = minmax[xvar]
xlist = np.linspace(xminmax[0], xminmax[1], Npts)
yminmax = minmax[yvar]
ylist = np.linspace(yminmax[0], yminmax[1], Npts)
if (xvar, yvar) == ("ka","d"):
    func = lambda x,y: disc(x=myvals["x"], d=y, g=myvals["g"], ka=x, kb=myvals["kb"])
elif (xvar, yvar) == ("ka","kb"):
    func = lambda x,y: disc(x=myvals["x"], d=myvals["d"], g=myvals["g"], ka=x, kb=y)
elif (xvar, yvar) == ("ka","x"):
    func = lambda x,y: disc(x=y, d=myvals["d"], g=myvals["g"], ka=x, kb=myvals["kb"])
elif (xvar, yvar) == ("kb","g"):
    func = lambda x,y: disc(x=myvals["x"], d=myvals["d"], g=y, ka=myvals["ka"], kb=x)
elif (xvar, yvar) == ("kb","x"):
    func = lambda x,y: disc(x=y, d=myvals["d"], g=myvals["g"], ka=myvals["ka"], kb=x)
else:
    print("\n\nERROR: the given combination (x,y) = (%s, %s) is not supported.\n\n")
    1/0

# construct plot
STY = "k."
BG = "w"
WIDTH = 7
HEIGHT = 7
font = 14

fig = plt.figure("Phase Diagram: " + str(myvals), figsize=(WIDTH, HEIGHT), facecolor=BG)
ax = fig.add_subplot(111, facecolor=BG, xlim=xminmax, ylim=yminmax)
for x in xlist:
    for y in ylist:
        if abs(func(x,y)) < thrs:
##            print("\tPlotting a point...")
            ax.plot([x], [y], STY)
ax.set_xlabel(xvar, size=font)
ax.set_ylabel(yvar, size=font)
ax.grid(True)
ax.tick_params(size=font, labelsize=font)

# show
fig.tight_layout()
plt.show()
plt.close()

# done
del np, plt, A, B, C, D, disc
