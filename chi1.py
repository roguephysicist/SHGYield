from numpy import linspace, loadtxt, column_stack
from numpy.random import randn
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt

FILE = "./sample_data/chi1.kk_xx_yy_zz_3107_25-nospin_scissor_0.5_Nc_8"
energy, real, imaginary = loadtxt(FILE, unpack=True, usecols=[0,1,2])
#print data
#data = column_stack(data)
#z = real + 1j * imaginary
energy_new = 2 * energy
s = InterpolatedUnivariateSpline(energy, real)
ys = s(energy_new)

plt.plot(energy,ys)
plt.show()
