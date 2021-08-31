import cplot
import numpy as np

def f(x):
    return np.exp(0 + 1j * x)


cplot.plot(f, (-20, 20), (-20, 20), 400).show()
