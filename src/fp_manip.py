## Check emergence of fixed points (graphical)
## by Mike Phillips, 9/21/2019
##
## The l.h.s. side of the fixed point equation is plotted alongside the r.h.s.
## Solutions are given by intersection points.

import numpy as np, matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, CheckButtons

# setup
tanh = np.tanh
(x, d, g, ka, kb) = (0, 0, 0, 1, 1)
var_step = None

def rhs(y, x, d, g, ka, kb):
    res = (1 + x) * tanh(d + ka*y)
    res += (1 - x) * tanh(g + kb*y)
    res *= 0.5
    return res

# plot
STY_l = "k-"
STY_r = "r-"
BG = "w"
WIDTH = 12
HEIGHT = 7.5
Npts = 150
font = 14

fig = plt.figure("Fixed Point Intersections", facecolor=BG, figsize=(WIDTH, HEIGHT))
ax = fig.add_subplot(111, facecolor=BG)
plt.subplots_adjust(right=0.7, top=0.95)
y = np.linspace(-1, 1, Npts)
l1 = ax.plot(y, y, STY_l)[0]
l2 = ax.plot(y, rhs(y, x, d, g, ka, kb), STY_r)[0]
ax.set_xlabel("y", size=font)
ax.set_ylabel("lhs & rhs", size=font)
ax.grid(True)
ax.tick_params(size=font, labelsize=font)
##ax.margins(x=0.5)

# manipulate --> sliders & buttons to control behavior
#   # axes: [left, bottom, width, height], facecolor
scolor = "silver"
w, h = 0.02, 0.65

ax_x = plt.axes([0.72, 0.3, w, h], facecolor=scolor)
ax_d = plt.axes([0.78, 0.3, w, h], facecolor=scolor)
ax_g = plt.axes([0.84, 0.3, w, h], facecolor=scolor)
ax_ka = plt.axes([0.90, 0.3, w, h], facecolor=scolor)
ax_kb = plt.axes([0.96, 0.3, w, h], facecolor=scolor)

#   # sliders: ax, label, valmin, valmax, valinit, valstep, orientation
sliders = {}
sliders.update({ "x" : Slider(ax=ax_x, label="x", valmin=-1, valmax=1, valinit=x,
                           valstep=var_step, orientation="vertical") })
sliders.update({ "d" : Slider(ax=ax_d, label="$\delta$", valmin=-1, valmax=1, valinit=d,
               valstep=var_step, orientation="vertical") })
sliders.update({ "g" : Slider(ax=ax_g, label="$\gamma$", valmin=-1, valmax=1, valinit=g,
               valstep=var_step, orientation="vertical") })
sliders.update({ "ka" : Slider(ax=ax_ka, label="$\kappa_A$", valmin=0, valmax=2.5, valinit=ka,
                valstep=var_step, orientation="vertical") })
sliders.update({ "kb" : Slider(ax=ax_kb, label="$\kappa_B$", valmin=0, valmax=2.5, valinit=kb,
                valstep=var_step, orientation="vertical") })

def update(val):
    (xp, dp, gp, kap, kbp) = ( s.val for s in sliders.values() )
    l2.set_ydata( rhs(y, xp, dp, gp, kap, kbp) )
    fig.canvas.draw_idle()

##for s in sliders.values():    # update all the same (no kb scaling)
##    s.on_changed(update)
for sname in ("x", "d", "g"):
    sliders[sname].on_changed(update)

justnow = False

def ka_update(val):
    global justnow
    is_scaled = check.get_status()[0]
    (xp, dp, gp, kap, kbp) = ( s.val for s in sliders.values() )
    if is_scaled and not justnow:
        justnow = True
        kbp = kap/2
        sliders["kb"].set_val(kbp)
        justnow = False
    l2.set_ydata( rhs(y, xp, dp, gp, kap, kbp) )
    fig.canvas.draw_idle()

sliders["ka"].on_changed(ka_update)

def kb_update(val):
    global justnow
    is_scaled = check.get_status()[0]
    (xp, dp, gp, kap, kbp) = ( s.val for s in sliders.values() )
    if is_scaled and not justnow:
        justnow = True
        kap = 2*kbp
        sliders["ka"].set_val(kap)
        justnow = False
    l2.set_ydata( rhs(y, xp, dp, gp, kap, kbp) )
    fig.canvas.draw_idle()

sliders["kb"].on_changed(kb_update)

#   # checkbutton: ax, labels, actives
ccolor = "lightgray"

ax_check = plt.axes([0.73, 0.1, 0.15, 0.07], facecolor=ccolor)

check = CheckButtons(ax=ax_check, labels=("autoscale $\kappa_B$",))

def scalekb(event):
    is_scaled = check.get_status()[0]
    if is_scaled:
        (xp, dp, gp, kap, kbp) = ( s.val for s in sliders.values() )
        kbp = kap/2
        sliders["kb"].set_val(kbp)
        l2.set_ydata( rhs(y, xp, dp, gp, kap, kbp) )
        fig.canvas.draw_idle()
        
check.on_clicked(scalekb)

# (custom checkbox properties)
ratio = 2* WIDTH / HEIGHT
rect = check.rectangles[0]
ht = rect.get_width() * ratio
rect.set_height(ht)
ry = rect.get_y() * 3/4
rect.set_y(ry)
for j in range(len(check.lines[0])):
    l = check.lines[0][j]
    ydat = l.get_ydata()
    xdat = l.get_xdata()
    if j==0:
        l.set_ydata((ry+ht, ry))
    else:
        l.set_ydata((ry, ry+ht))
        

#   # button: ax, label, image, color, hovercolor
bcolor = "palegreen"
hcolor = "lightblue"

ax_button = plt.axes([0.89, 0.1, 0.09, 0.05], facecolor=bcolor)

the_button = Button(ax=ax_button, label="Reset", color=bcolor, hovercolor=hcolor)

def reset(event):
    is_scaled = check.get_status()[0]
    if is_scaled:
        check.set_active(0)
    for s in sliders.values():
        s.reset()
        
the_button.on_clicked(reset)
    
# show
##fig.tight_layout()
plt.show()
plt.close()


del np, plt
