import numpy as np
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerLine2D
#from mpl_toolkits import mplot3d

import constants as c

class SchroedingersEq: 
    def __init__(self, mass, boundry, energy_level):
        self.m = mass
        self.a = boundry
        self.n = energy_level

    def solve(self, x):#, t):
        """
        i hbar d(psi(x, t))/dt = -hbar^2/2m d^2(psi(x, t))/dx^2 + V(x)psi(x, t)
        
        let psi(x, t) = u(x)T(t)
        => i hbar d(T(t))/dt = (-hbar^2/2m d^2(u(x))/dx^2 + V(x)u(x)) * T(t)/u(x)
        
        let E = (-hbar^2/2m d^2(u(x))/dx^2 + V(x)u(x))/u(x)
        => i hbar d(T(t))/dt = E T(t)
        => 1/T(t) d(T(t))/dt = -iE / hbar
        => ln|T(t)| = -iEt/hbar + C
        => T(t) = e^(-iEt/hbar) * e^C
        
        let A = e^C
        => T(t) = Ae^(-Et/hbar)

        E = (-hbar^2/2m d^2(u(x))/dx^2 + V(x)u(x))/u(x)
        => Eu(x) = -hbar^2/2m d^2(u(x))/dx^2 + V(x)u(x)
        => d^2(u(x))/dx^2 = 2m(E - V(x))u(x)/-hbar^2

        at V(x) = 0
        => d^2(u(x))/dx^2 = 2mEu(x)/-hbar^2
        
        let k  = sqrt(2mE)/hbar
        => d^2(u(x))/dx^2 = -k^2 u(x)
        => u(x) = sin(kx) or cos(kx)

        let u(-)(x) = sin(kx)
        let u(+)(x) = sin(kx)

        We will first look at u(-)(x)
        
        at x = 0, psi(-)(x, t) = 0 (i.e.: u(-)(x) = 0)
        => 0 = sin(ka) = sin(0) = 0
            => ka = 0
            => ka = n pi
        => u(-)(x) = sin(n pi x/a)
        To normalize:
        => u(-)(x) = sqrt(2/a) sin(n pi x/a)
        => E(-)(n) = n^2 pi^2 hbar^2/(2m a^2)

        psi(-)(x) = u(-)(x) * T(t)
        = sqrt(2/a) sin(n pi x/a) e^(-i n^2 pi^2 hbar t/(2ma^2))
        psi(+)(x) = u(+)(x) * T(t)
        = sqrt(2/a) sin((n - 0.5) pi x/a) e^(-i (n - 0.5)^2 pi^2 hbar t/(2ma^2))
        """

        self.E = [(self.n * c.pi * c.hbar) ** 2/(2 * self.m * self.a ** 2), ((self.n - 0.5) * c.pi * c.hbar) ** 2/(2 * self.m * self.a ** 2)]

        return [np.sqrt(2/self.a) * np.sin(self.n * c.pi * x/self.a), # * np.exp(-complex(0, 1) * (self.n * c.pi) ** 2 * c.hbar * t/(2*self.m*self.a**2)),
                np.sqrt(2/self.a) * np.cos((self.n - 0.5) * c.pi * x/self.a)] # * np.exp(-complex(0, 1) * ((self.n - 0.5) * c.pi) ** 2 * c.hbar * t/(2*self.m*self.a**2))]


if __name__ == "__main__":
    a = float(input("Boundry of potential barrier (a): "))
    eq = SchroedingersEq(float(input("Mass of particle: ")), a, float(input("Energy level of particle: ")))
    
    #_x = np.linspace(0, a)
    #_t = np.linspace(0, 10)
    #, t 
    #x = np.meshgrid(_x, _t)
    xa = np.arange(0, a)
    psi = eq.solve(xa)#, t)
    x = np.arange(a - 10, a + 10, 0.1)
    psi_full = eq.solve(x)
    
    print(f"Energy of particle (E(-)): {eq.E[0]}")
    print(f"Energy of particle (E(+)): {eq.E[1]}")
    
    #fig = plt.figure()
    #ax = plt.axes(projection='3d')

    #ax.plot_surface(x, t, psi[0], rstride=1, cstride=1, cmap='viridis', edgecolor='none')
    #ax.plot_surface(x, t, psi[1], rstride=1, cstride=1, cmap='viridis', edgecolor='none')

    #psi1, = plt.plot(xa, psi[0], label=u"\u03C8(-)(x)")
    #psi2, = plt.plot(xa, psi[1], label=u"\u03C8(+)(x)")
    psi_full1, = plt.plot(x, psi_full[0], label=u"Full \u03C8(-)(x)")
    psi_full2, = plt.plot(x, psi_full[1], label=u"Full \u03C8(+)(x)")
    psi_full_probability1, = plt.plot(x, psi_full[0] * np.conj(psi_full[0]), label=u"<\u03C8(-)(x)|\u03C8(-)(x)> = |\u03C8(-)(x)|^2")
    psi_full_probability2, = plt.plot(x, psi_full[1] * np.conj(psi_full[1]), label=u"<\u03C8(+)(x)|\u03C8(+)(x)> = |\u03C8(+)(x)|^2")
    potential1, = plt.plot([0, 0, 0], [-1, 0, 1], label="Potential barrier")
    potential2, = plt.plot([a, a, a], [-1, 0, 1], label="Potential barrier")
    #plt.legend(handler_map={psi1: HandlerLine2D(numpoints=4)})
    plt.legend(handler_map={psi_full1: HandlerLine2D(numpoints=4)})
    plt.show()
