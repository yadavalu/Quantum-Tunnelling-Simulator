import numpy as np
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerLine2D
from tkinter import *


# Gaussian curve
"""
def V(x, sigma):
    return np.exp(-0.5 * (x/0.1) ** 2)/(sigma * np.sqrt(2 * np.pi))
"""

def V(x, location, barrier):
    return 0 if x < -location/2 or x > location/2 else barrier


class SchroedingersEq:

    h = 6.62607004 * 10 ** -34
    hbar = h/(2 * np.pi)

    """
    -hbar^2/2m d^2(psi)/dx^2 + V * psi = E * psi
    (H psi = E psi)
    """

    def __init__(self, barrier, length, mass):
        self.potential = barrier
        self.length = length
        self.mass = mass

    def solve(self):
        # psi(x) ln|psi(x)| - psi(x) = m/(hbar^2) * (2*l^2*V(x) - E*x^2)
        pass


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
        for i in range(-100, 100, 1):
            y.append(V(i/10, float(self.barrier_length.get()), float(self.barrier_height.get())))

        potential, = plt.plot(x, y, label="Potential barrier")
        particle, = plt.plot(-10, float(self.energy.get()), "ro", label="Particle")
        #plt.legend([potential, particle], ["Potential barrier", "Particle"])
        plt.legend(handler_map={potential: HandlerLine2D(numpoints=4)})
        plt.show()


if __name__ == "__main__":
    root = Tk()
    root.geometry("480x200")
    root.title("Quantum tunnelling")
    Potential(root)
    root.mainloop()
