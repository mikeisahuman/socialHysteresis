##  Implicit Phase Diagram
##  by Mike Phillips, 10/9/2019
##
##  Following from a second-order Taylor expansion of the transition equation, the zeros
##  of the resulting quadratic discriminant are plotted implicitly (kapa vs. kapb).
##  These zeros should give approximate phase boundaries.

from sympy import symbols, Eq, plot_implicit, tanh, sech

# plot variables
kb, ka = symbols("x y")

# parameters
(d, g, x) = (0.4, -0.4, -0.5)

# This form of discriminant is ONLY for x=0!    [really: disc = (discriminant)/4]
disc = Eq( d*(ka**4) + g*(kb**4) + 2*d*g*(ka**2)*(kb**2) - (d**2 - 1)*ka*(ka**3 + kb**3)
           - (g**2 - 1)*kb*(ka**3 + kb**3) - 2*(ka**3 + kb**3) , 0 )

# Better, more complete form                    [really: disc_comp = (complete discriminant)/4]
disc_comp = Eq( -(1+x)**2 * (ka**4) * (sech(d)**4) - (1-x)**2 * (kb**4) * (sech(g)**4)
                + 2*(1+x) * (ka**3) * (sech(d)**2) * (tanh(d)**2+1)
                + 2*(1-x) * (kb**3) * (sech(g)**2) * (tanh(g)**2+1)
                - (1-x**2)*ka*kb*(sech(d)**2)*(sech(g)**2) * ( (ka*tanh(d)-kb*tanh(g))**2 + ka**2 + kb**2), 0)

# sympy implicit plotting
p = plot_implicit(disc_comp, (kb, 0, 2), (ka, 0, 2))
