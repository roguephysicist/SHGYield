from numpy import linspace, exp, loadtxt, savetxt, column_stack
from numpy.random import randn
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt

#x = linspace(0, 20, 2001)
#y = exp(-x**2) + randn(2001)/10
x, y = loadtxt('sample_data/chi1.kk_xx_yy_zz_3107_25-nospin_scissor_0.5_Nc_8', unpack=True, usecols=[0,1])
s = InterpolatedUnivariateSpline(x, y)

print type(s(x))
