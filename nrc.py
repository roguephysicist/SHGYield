import math
import random
from numpy import array, loadtxt, savetxt, genfromtxt, column_stack, logical_or, fabs, sqrt
from scipy import constants

########### user input ###########
#### XYZ file disorder
DISORDER_IN = "./sample_data/disorder/si_h_14.xyz" # full path to xyz file you want to disorder
BOND_LENGTH = 2.351 # bond length in angstrom for Si
ATOMS_DISORDERED = [2, 3, 4] # atoms you want to disorder
DISORDER_AMOUNT = [0.1, 0.1, 0.1] # amount you want to disorder each atom
REPEAT = 1 # number of times you want to repeat the disordering

#### Nonlinear reflection coefficient
NRC_PATH = "./sample_data/res/" # full path to 'res' folder for desired structure
K_POINTS = 1219 # kpoints of the case you want to use
ECUT = 20 # ecut of the case you want to use
# Polarization
ENTRY = "s" # 's' or 'p'
EXIT = "s" # 's' or 'p'
# Angles
THETA_DEG = 65 # angle of incidence in degrees
PHI_DEG = 30 # azimuthal angle
# Misc
ELEC_DENS = 50 # electronic density
########### end user input ###########

## Housekeeping
THETA_RAD = math.radians(THETA_DEG)
PHI_RAD = math.radians(PHI_DEG)

########### functions ###########
def disorder():
    count = 1
    while count <= REPEAT:
        SELECTED = [x - 1 for x in ATOMS_DISORDERED]
        BL_BOHR = BOND_LENGTH * constants.angstrom / constants.value("Bohr radius")
        XYZ = load_matrix("disorder")
        DATA = XYZ[SELECTED]
        NEW_XYZ = XYZ.copy()
        FINAL = zip(DATA, DISORDER_AMOUNT)
        for atom in range(0, len(FINAL)):
            polar = random.random() * constants.pi
            azimuthal = random.random() * 2 * constants.pi
            NEW_XYZ[SELECTED[atom], 0] = FINAL[atom][0][0] + (FINAL[atom][1] * (BL_BOHR / 2) * math.sin(polar) * math.cos(azimuthal))
            NEW_XYZ[SELECTED[atom], 1] = FINAL[atom][0][1] + (FINAL[atom][1] * (BL_BOHR / 2) * math.sin(polar) * math.sin(azimuthal))
            NEW_XYZ[SELECTED[atom], 2] = FINAL[atom][0][2] + (FINAL[atom][1] * (BL_BOHR / 2) * math.cos(polar))
        out = "./sample_data/disorder/si_h_14_mod_" + str(count).zfill(3) + ".xyz"
        save_matrix(out, NEW_XYZ)
        count += 1

def nonlinear_reflection_coefficient():
    e, chi1 = load_matrix("chi1")
    twoe, chi1twoe = load_matrix("twoe")
    R = rif_constants(e) * fabs(fresnel_vs(EXIT, twoe, chi1twoe) * fresnel_sb(EXIT, twoe, chi1twoe) * ((fresnel_vs(ENTRY, e, chi1) * fresnel_sb(ENTRY, e, chi1)) ** 2) * r_factors(ENTRY, EXIT, e, twoe, chi1, chi1twoe)) ** 2
    R = column_stack((e, R))
    out = NRC_PATH + "R" + ENTRY + EXIT + "_" + str(K_POINTS) + "_" + str(ECUT)
    save_matrix(out, R)

def rif_constants(energy):
    const = (32 * (constants.pi ** 3) * (energy ** 2)) / ((ELEC_DENS * constants.e ** 2) * (constants.c ** 3) * (math.cos(THETA_RAD) ** 2))
    return const

def fresnel_vs(polarization, energy, chi1):
    if polarization == "s":
        t = (2 * math.cos(THETA_RAD)) / (math.cos(THETA_RAD) + wave_vector(energy, chi1))
    elif polarization == "p":
        t = (2 * math.cos(THETA_RAD)) / (epsilon(chi1) * math.cos(THETA_RAD) + wave_vector(energy, chi1))
    return t

def fresnel_sb(polarization, energy, chi1):
    if polarization == "s":
        t = (2 * wave_vector(energy, chi1)) / (wave_vector(energy, chi1) + wave_vector(energy, chi1))
    elif polarization == "p":
        t = (2 * wave_vector(energy, chi1)) / (epsilon(chi1) * wave_vector(energy, chi1) + epsilon(chi1) * wave_vector(energy, chi1))
    return t

def r_factors(polar_in, polar_out, energy, twoe, chi1, chi1twoe):
    xxx = load_matrix("xxx")
    if polar_in == "p" and polar_out == "p":
        zzz = load_matrix("zzz")
        zxx = load_matrix("zxx")
        xxz = load_matrix("xxz")
        r = math.sin(THETA_RAD) * epsilon(chi1twoe) * (((math.sin(THETA_RAD) ** 2) * (epsilon(chi1) ** 2) * zzz) + (wave_vector(energy, chi1) ** 2) * (epsilon(chi1) ** 2) * zxx) + epsilon(chi1) * epsilon(chi1twoe) * wave_vector(energy, chi1) * wave_vector(twoe, chi1twoe) * (-2 * math.sin(THETA_RAD) * epsilon(chi1) * xxz + wave_vector(energy, chi1) * epsilon(chi1) * xxx * math.cos(3 * PHI_RAD))
    elif polar_in == "s" and polar_out == "p":
        zxx = load_matrix("zxx")
        r = math.sin(THETA_RAD) * epsilon(chi1twoe) * zxx - wave_vector(twoe, chi1twoe) * epsilon(chi1twoe) * xxx * math.cos(3 * PHI_RAD)
    elif polar_in == "p" and polar_out == "s":
        r = -(wave_vector(energy, chi1) ** 2) * (epsilon(chi1) ** 2) * xxx * math.sin(3 * PHI_RAD)
    elif polar_in == "s" and polar_out == "s":
        r = xxx * math.sin(3 * PHI_RAD)
    return r

def epsilon(chi1):
    epsilon = 1 + (4 * constants.pi * chi1)
    return epsilon

def wave_vector(energy, chi1):
    k = (energy / constants.c) * sqrt(epsilon(chi1) - (math.sin(THETA_RAD) ** 2))
    return k

def load_matrix(comp):
    if comp == "disorder":
        data = loadtxt(DISORDER_IN)
    elif comp == "chi1":
        f = NRC_PATH + "chi1.kk_xx_yy_zz_" + str(K_POINTS) + "_" + str(ECUT) + "-nospin_scissor_0_Nc_29"
        data = genfromtxt(f, skip_footer=1000, unpack=True, usecols=[0, 6])
    elif comp == "twoe":
        f = NRC_PATH + "chi1.kk_xx_yy_zz_" + str(K_POINTS) + "_" + str(ECUT) + "-nospin_scissor_0_Nc_29"
        new = column_stack(loadtxt(f, unpack=True, usecols=[0, 6]))
        wanted = array(2 * new[:, 0]).reshape(-1).tolist()
        data = column_stack(new[logical_or.reduce([new[:, 0] == x for x in wanted])])
    elif comp == "zzz" or comp == "zxx" or comp == "xxz" or comp == "xxx":
        f = NRC_PATH + "shgC.kk_" + comp + "_" + str(K_POINTS) + "_half-slab_" + str(ECUT) + "-nospin_scissor_0_Nc_29"
        data = genfromtxt(f, skip_footer=1000, unpack=True, usecols=[4])
    return data

def save_matrix(f, data):
    savetxt(f, data, fmt=('%5.14e'), delimiter='\t')

#disorder()
nonlinear_reflection_coefficient()