##  Model Diagram
##  by Mike Phillips, 8/21/2018
##
##  A simple diagram is constructed to demonstrate populations and transitions.

import numpy as np
import matplotlib.pyplot as plt

originAxes = True
labelNeutral = False

colors = ["#339900", "#3366FF", "#FF6633", "#CC33CC"]

##  Plot font sizes
SMALL_SIZE = 10
MEDIUM_SIZE = 12
LARGE_SIZE = 14

plt.rc('font', size=LARGE_SIZE)        # controls default text sizes
plt.rc('axes', titlesize=LARGE_SIZE)    # fontsize of the axes title
plt.rc('axes', labelsize=LARGE_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)   # legend fontsize
plt.rc('figure', titlesize=LARGE_SIZE)  # fontsize of the figure title

##  Circle Function (quarter-circles centered on origin)
def makeCircle(quad=(1,1), r=1, npts=120):
    r = abs(r)
    if len(quad) == 2:
        if quad[0] == 0 and quad[1] == 0:
            print("\n" + "ERROR: must enter at least one nonzero value to select a quadrant." + "\n"*2)
    else:
        print("\n" + "ERROR: quadrant selection must be an iterable of length 2." + "\n"*2)
    xpm = np.sign(quad[0])
    ypm = np.sign(quad[1])
    x = np.linspace(0,xpm*r,npts)
    y = ypm*np.sqrt(r**2 - x**2)
    x, y = (list(x) + [0]), (list(y) + [0])
    return {"x":x, "y":y}

##  Annotation Function
def myAnnotate(ax, xycoord=(0,0), textcoord=(-1,-1), radius=0.5, alpha=0.5, size=20, fcolor="w", text=""):
    offsets = [0.04, 0.1]
    sx = np.sign(textcoord[0])
    sy = np.sign(textcoord[1])
    sgn = - sx*sy
    ax.annotate(text,
                xy=xycoord,
                xytext=textcoord,
                size=20,
                arrowprops=dict(arrowstyle="simple",
                                facecolor = fcolor, edgecolor="k", alpha=0.5,
                                connectionstyle="arc3, rad=" + str(sgn*radius) ) )
    if sx > 0:
        lx = "A"
        offx = offsets[0]
    else:
        lx = "B"
        offx = - offsets[0]
    if sy < 0:
        liy, lfy = "2", "1"
        offy = offsets[1]
    else:
        liy, lfy = "1", "2"
        offy = - offsets[1]
        offx *= 3
    linit, lfinal = (lx + liy), (lx + lfy)
    ax.annotate(r"$p_{" + lfinal + r" \leftarrow " + linit + r"}$",
                (xycoord[0] + offx, xycoord[1] + offy), ha="center", va="center", color=fcolor)

##  Figure and Plots
fig = plt.figure("Population Transition Diagram", figsize=(6,6))
ax = fig.add_subplot(111, xlim=(-1.4,1.4), ylim=(-1.2, 1.2))
ax.set_xlabel("infection status")
ax.set_ylabel("behavior")
ax.tick_params(axis="both", direction="in", bottom="on", top="on", left="on", right="on",
               labelbottom="on", labeltop="off", labelleft="on", labelright="off")

##box1 = {"x":[0,1,1,0], "y":[1,1,0,0]}
##ax.fill(box1["x"],box1["y"])

circles = []

quads = [ (i,j) for i in [1,-1] for j in [1,-1] ]
poplabels = [ ",".join([s,t]) for s in ["A","B"] for t in ["1","2"]]
poploc = 0.4243
for i in range(len(quads)):
    circles.append( makeCircle(quads[i]) )
    ax.fill(circles[i]["x"], circles[i]["y"], color=colors[i], alpha=0.4)
    ax.annotate(r"$n_{%s}$" % poplabels[i], (quads[i][0]*poploc,quads[i][1]*poploc), ha="center", va="center")

ann_x = [0.75]*2 + [0.95]*2 + [-0.75]*2 + [-0.95]*2
ann_y = [-0.75, 0.75]
ann_y += reversed(ann_y)
ann_y *= 2
coords = list(zip(ann_x,ann_y))
for i in range(len(coords)//2):
    myAnnotate(ax, xycoord=coords[2*i+1], textcoord=coords[2*i], fcolor=colors[i])

if originAxes:
    ax.plot([-1,1],[0,0],"k--")
    ax.plot([0,0],[-1,1],"k--")

tks = [-0.5, 0.5]
xlb = ["B" + "\n(infected)","A" + "\n(healthy)"]
ylb = ["2","1"]
if labelNeutral:
    tks += [0]
    xlb += [""]
    ylb += [""]
plt.xticks(tks, xlb)
plt.yticks(tks, ylb)
ax.set_aspect("equal")

plt.tight_layout()
plt.show()
plt.close()
