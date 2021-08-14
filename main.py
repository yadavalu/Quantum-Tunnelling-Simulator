import numpy as np
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerLine2D
from tkinter import *

import constants as c


# Gaussian curve
"""
def V(x, sigma):
    return np.exp(-0.5 * (x/0.1) ** 2)/(sigma * np.sqrt(2 * np.pi))
"""

def V(x, location, barrier):
    return 0 if x < -location/2 or x > location/2 else barrier


class SchroedingersEq:
    """
    -hbar^2/2m d^2(psi)/dx^2 + V * psi = E * psi
    (H psi = E psi)
    """

    def __init__(self, barrier, length, mass):
        self.potential = barrier
        self.length = length
        self.mass = mass
        # TODO: E
        self.E = 1

    def solve(self):
        self.psi = []

        ## Area 1; V(x) = 0
        # -hbar/2pi d^2(psi(x))/dx^2 = Epsi(x)
        # d^2(psi(x))/dx^2 = 2mEpsi(x)/(-hbar^2)
        # k = sqrt(2mE)/hbar
        k = np.sqrt(2 * self.mass * self.E)/c.hbar
        # d^2(psi(x))/dx^2 = -k^2 psi(x)
        # psi(x) = c1 sin(kx) + c2 cos(kx)
        x = np.arange(-10, 10, 0.1)
        # TODO: c1, c2
        c1and31 = 1
        c1and32 = 1
        #self.psi1 = c1and31 * np.sin(k * x) + c1and32 * np.cos(k * x)
        #yield self.psi1

        ## TODO: 
        # ValueError: operands could not be broadcast together with shapes (0,) (200,)
        # Couldn't multiply k * x
        self.psi += c1and31 * np.sin(k * x) + c1and32 * np.cos(k * x)

        ## Area 2; V(x) = self.potential
        # -hbar/2pi d^2(psi(x))/dx^2 + V(x) = Epsi(x)
        # d^2(psi(x))/dx^2 = 2m(E - V(x))psi(x)/(-hbar^2)
        # k = sqrt(2m(E - V(x)))/hbar
        k = np.sqrt(2 * self.mass * (self.E - self.potential))/c.hbar
        # d^2(psi(x))/dx^2 = -k^2 psi(x)
        # psi(x) = c1 sin(kx) + c2 cos(kx)
        x = np.arrange(-10, 10, 0.1)
        # TODO: c1, c2
        c21 = 1
        c22 = 1
        #self.psi2 = c21 * np.sin(k * x) + c22 * np.cos(k * x)
        #yield self.psi2

        ## TODO: 
        # ValueError: operands could not be broadcast together with shapes (0,) (200,)
        # Couldn't multiply k * x
        self.psi += c21 * np.sin(k * x) + c22 * np.cos(k * x)

        ## Area 3; V(x) = 0
        #yield self.psi1
        self.psi += self.psi[0]

        return self.psi




class Potential(Frame):
    def __init__(self, parent):
        self.parent = parent
        super().__init__(self.parent)
        self.grid(row=0, column=0)

        self.sigma = 0

        #Label(self, text="Strength of potential barrier (standard deviation^-1 for Gaussian graph): ").grid(row=0, column=0)
        Label(self, text="Energy of potential barrier (in eV): ").grid(row=0, column=0)
        self.barrier_height = Entry(self)
        self.barrier_height.grid(row=0, column=1)
        
        Label(self, text="Length of potential barrier (in 10^-15 m): ").grid(row=1, column=0)
        self.barrier_length = Entry(self)
        self.barrier_length.grid(row=1, column=1)

        Label(self, text="Energy of particle (int eV): ").grid(row=2, column=0)
        self.energy = Entry(self)
        self.energy.grid(row=2, column=1)

        Label(self, text="Mass (in GeV/(c^2)): ").grid(row=3, column=0)
        self.mass = Entry(self)
        self.mass.bind("<Return>", self._plot_barrier)
        self.mass.grid(row=3, column=1)

        Button(self, text="OK", command=self._plot_barrier).grid(row=4, column=1)

    def _plot_barrier(self, *args):
        x = np.arange(-10, 10, 0.1)
        y = []
        self.length = float(self.barrier_length.get())
        self.height = float(self.barrier_height.get())

        for i in range(-100, 100, 1):
            y.append(V(i/10, self.length, self.height))

        potential, = plt.plot(x, y, label="Potential barrier")
        particle, = plt.plot(-10, float(self.energy.get()), "ro", label="Particle")
        plt.legend(handler_map={potential: HandlerLine2D(numpoints=4)})
        plt.show()

        self.psi = SchroedingersEq(y, self.length, float(self.mass.get())).solve()
        print(self.psi)
        plt.plot(x, self.psi[0])
        plt.plot(x, self.psi[1])
        plt.plot(x, self.psi[2])


if __name__ == "__main__":
    root = Tk()
    #root.geometry("480x200")
    root.title("Quantum tunnelling")
    Potential(root)
    root.mainloop()
