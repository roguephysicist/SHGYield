"""
Quick epsilon generation for testing purposes.
"""
import numpy as np

CHI1B = "SiBulk-chi1-xx_yy_zz"
CHI1L = "SiH1x1-chi1-xx_yy_zz"
NORM = 1.2659296143

def epsload(in_file, norm):
    '''
    Reads calculated chi1 file, that is organized as
    Energy(1w) Re[chi_xx] Im[chi_xx] Re[chi_yy] Im[chi_yy] Re[chi_zz] Im[chi_zz].
    Converts to epsilon = 1 + 4*pi*chi1, and normalized chi1 according to
    PRB 92, 245308 (2015). Returns the averaged epsilon
    epsilon_{avg} = (epsilon^{xx} + epsilon^{yy} + epsilon^{zz})/3
    in a numpy array.
    '''
    global FREQ
    FREQ, rexx, imxx, reyy, imyy, rezz, imzz = np.loadtxt(in_file, unpack=True)
    real = (rexx + reyy + rezz)/3      # real average
    imag = (imxx + imyy + imzz)/3      # imag average
    coma = real + 1j * imag                     # complex average
    epsa = 1 + (4 * np.pi * coma * norm) # epsilon with normalization
    return epsa

EPSB = epsload(CHI1B, 1)
EPSL = epsload(CHI1L, NORM)

np.savetxt('SiBulk-epsilon', np.column_stack((FREQ, EPSB.real, EPSB.imag)),
           fmt=('%05.2f', '% .8e', ' % .8e'), delimiter='    ',
           header='hw(1w)  real                imag')
np.savetxt('SiH1x1-epsilon', np.column_stack((FREQ, EPSL.real, EPSL.imag)),
           fmt=('%05.2f', '% .8e', ' % .8e'), delimiter='    ',
           header='hw(1w)  real                imag')
