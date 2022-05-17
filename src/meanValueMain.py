## Social Dynamics - 2 populations (infected/non-infected) / 2 opinions (attitudes / behaviors)
## by Mike Phillips, 6/13/2018
##
## Exploring the drift coefficients and other properties appearing in the simple model.
## The form of transition rates / drift coefficients comes from the Master / Fokker-Planck equations
## Here we look at the case when the infected population is input as an explicit function of time.
##  The reason is that the infected fraction can remain nearly stationary when death processes are considered.
##  (If the base rate of infection is higher than that of curing, as expected since the rate of curing
##      can be nearly zero, then the stationary point for the infected fraction is close to unity, but
##      in reality new [uninfected] generations are born while older [infected] generations die quickly.)
## Taking this case, we look for solutions of mean value equations for the remaining two variables:
##  x -> infection (fixed/controlled), y -> opinion (variable of interest), z -> 'belonging'
##

##import numpy as np
##from scipy import optimize
##import matplotlib.pyplot as plt

from meanValueDefs import *

PHASES = False           #   if you want to see (simple) phase curves: y,z vs. parameters
DRIFTS = False           #   if you want to see plots of raw drifts: Ky, Kz vs. its own variable
DRIFT3D = False          #   if you want 3D (surface/wireframe) plots of drifts
FIXEDPTS = True         #   display fixed points / classification (from apx. nonlinear ODEs -> Jacobian)
ODESOL = True           #   if you want (apx.) solutions to system of ODEs


#   Initialize Time
tmin = 0
t = tmin
#   Initialize Fixed Point
[y0,z0] = [0,0]
#   Initialize Initial Values
yi, zi = 0.5, 0.5
#   Initialize Final / Threshold Values
tmax = 20
nmax = 120

if PHASES:
    import phaseCurves
    [y0,z0] = phaseCurves.phasecurves(t)

##if DRIFTS:
##    import driftPlots
##    driftPlots.driftplots(y0,z0,t)

if FIXEDPTS:
    import symFP
    fpts = symFP.SYMFP(t, True)     #   second argument is for printing values/classificiations

if DRIFT3D:
    import driftPlots3D
    driftPlots3D.drift3D(t)

if ODESOL:
    import odeSol
    (tspace, ysol, zsol, fp) = odeSol.odesol(yi,zi,tmin,tmax,nmax)

if DRIFTS and len(fp) > 0:
    import driftPlots
    driftPlots.driftplots(fp[0],fp[1],tmax)


    



