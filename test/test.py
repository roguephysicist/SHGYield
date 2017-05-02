#!/usr/bin/env python
import os
from subprocess import call
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

call(['python', '../shgyield.py', 'test.yml'])

print('Loading reference data.')
reference = np.loadtxt('data/reference.dat', unpack=True)
print('Calculating SHG yield for comparison.')
test = np.loadtxt('test-output.dat', unpack=True)
print(48 * '=')

color1 = '#268bd2'
color2 = '#dc322f'
plt.rcParams['font.family'] = 'serif'
plt.rcParams['mathtext.fontset'] = 'dejavuserif'

fig, axes = plt.subplots(2, 2, figsize=(10,6))

axes[0][0].text(0.05, 0.10, r'$\mathcal{R}_{pP}$', fontsize=16, transform = axes[0][0].transAxes)
axes[0][0].plot(2*reference[0], reference[1], label='Reference', color=color1, lw=1.5, ls='-')
axes[0][0].plot(2*test[0], test[1], label='Test', color=color2, lw=2.5, ls=':')
axes[0][0].set_xlim([2.5, 5])
axes[0][0].set_ylim([0, 3.0])
axes[0][0].grid(True)
axes[0][0].set_ylabel(r'$\mathcal{R}\, (10^{-20} \times \mathrm{cm}^{2}/\mathrm{W})$', fontsize=14)
axes[0][0].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
axes[0][0].legend(loc=2)

axes[0][1].text(0.05, 0.10, r'$\mathcal{R}_{pS}$', fontsize=16, transform = axes[0][1].transAxes)
axes[0][1].plot(2*reference[0], reference[2], label='Reference', color=color1, lw=1.5, ls='-')
axes[0][1].plot(2*test[0], test[2], label='Test', color=color2, lw=2.5, ls=':')
axes[0][1].set_xlim([2.5, 5])
axes[0][1].set_ylim([0, 0.4])
axes[0][1].grid(True)

axes[1][0].text(0.05, 0.10, r'$\mathcal{R}_{sP}$', fontsize=16, transform = axes[1][0].transAxes)
axes[1][0].plot(2*reference[0], reference[3], label='Reference', color=color1, lw=1.5, ls='-')
axes[1][0].plot(2*test[0], test[3], label='Test', color=color2, lw=2.5, ls=':')
axes[1][0].set_xlim([2.5, 5])
axes[1][0].set_ylim([0, 0.06])
axes[1][0].grid(True)
axes[1][0].set_xlabel(r'Two-photon energy (eV)', fontsize=12)
axes[1][0].set_ylabel(r'$\mathcal{R}\, (10^{-20} \times \mathrm{cm}^{2}/\mathrm{W})$', fontsize=14)

axes[1][1].text(0.05, 0.10, r'$\mathcal{R}_{sS}$', fontsize=16, transform = axes[1][1].transAxes)
axes[1][1].plot(2*reference[0], reference[4], label='Reference', color=color1, lw=1.5, ls='-')
axes[1][1].plot(2*test[0], test[4], label='Test', color=color2, lw=2.5, ls=':')
axes[1][1].set_xlim([2.5, 5])
axes[1][1].set_ylim([0, 0.1])
axes[1][1].grid(True)
axes[1][1].set_xlabel(r'Two-photon energy (eV)', fontsize=12)

fig.tight_layout()

print('Drawing plots for comparing results.')
print('"Reference" and "Test" curves should be identical.')
print('Close plot to finish test program.')
print(48 * '=')
plt.show();


try:
    os.remove('test-output.dat')
    print('Removing test data. Test complete!')
except OSError:
    pass
    
