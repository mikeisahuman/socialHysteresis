### File for defining basic parameters and functions

import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d

#   List of parameters
# rates
##ni = 1 ; nc = ni/2      #   infection/cure base rates
na_const = 1. ; nb_const = na_const#/3     #   opinion base rates for each population
def na(t):
    return na_const
def nb(t):
    return nb_const
# shifts
dlt_const = 0. ; gam_const = 0.
def dlt(t):
    return dlt_const
    omg = na_const/2
    #return dlt_const*np.sin(omg*t)
def gam(t):
    return gam_const
# interactions
##kapi = 0.15             #   infection influence (feedback)
kapa_const = 0.5 ; kapb_const = kapa_const#/10
def kapa(t):
    return kapa_const
def kapb(t):
    return kapb_const

# externally-controlled healthy/infected fraction of the population, its derivative, and fluctuation term
x_const = 0.4
def x(t):
    #return x_const
    omg = na_const*(2/3)
    return x_const*np.cos(omg*t)
def dx(t):
    return 0
Qxx_const = 0
def Qxx(t):
    return Qxx_const

#   Define functional forms of drift coefficients
##def Kx(x):
##    return ( nc * (1 - x) - ni * (1 + x) * np.exp(-kapi*x) )

def Ky(y,z,t):
    res = na(t) * ( (1 + x(t)) * np.sinh(dlt(t) + kapa(t)*y) - (y + z) * np.cosh(dlt(t) + kapa(t)*y) )
    res += nb(t) * ( ( 1 - x(t)) * np.sinh(-gam(t) + kapb(t)*y) - (y - z) * np.cosh(-gam(t) + kapb(t)*y) )
    return res

def Kz(y,z,t):
##    res = nc * (y - z) - ni * (y + z) * np.exp(-kapi*x)
    res = na(t) * ( (1 + x(t)) * np.sinh(dlt(t) + kapa(t)*y) - (y + z) * np.cosh(dlt(t) + kapa(t)*y) )
    res += - nb(t) * ( ( 1 - x(t)) * np.sinh(-gam(t) + kapb(t)*y) - (y - z) * np.cosh(-gam(t) + kapb(t)*y) )
    if abs(x(t)) < 1:
        res += (Qxx(t) + dx(t)) * (y - z)/(1 - x(t)) - (Qxx(t) - dx(t))*(y + z)/(1 + x(t))
    return res

