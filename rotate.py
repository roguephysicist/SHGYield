import sys
import yaml
from scipy import ndimage
import numpy as np

infile = sys.argv[1]
angle = sys.argv[2]
GAMMA = np.radians(float(angle))

def parse_input(infile):
    with open(infile) as data_file:
        params = yaml.load(data_file)
    return params

def shgload(infile):
    global freq
    freq, re1w, im1w, re2w, im2w = np.loadtxt(infile, unpack=True)
    real = broad(re1w + re2w, PARAM['chi2']['sigma'])
    imag = broad(im1w + im2w, PARAM['chi2']['sigma'])
    chi2 = real + 1j * imag
    return chi2

def broad(target, sigma):
    '''
    A function for applying Gaussian broadening on the final output data.
    '''
    data = ndimage.filters.gaussian_filter(target, sigma)
    return data


PARAM = parse_input(infile)    # Parses input file

CHI2 = {}
ALL_COMPONENTS = ['xxx', 'xyy', 'xzz', 'xyz', 'xxz', 'xxy',
                  'yxx', 'yyy', 'yzz', 'yyz', 'yxz', 'yxy',
                  'zxx', 'zyy', 'zzz', 'zyz', 'zxz', 'zxy']
for component in ALL_COMPONENTS:
    if component in PARAM['chi2']:
        value = PARAM['chi2'][component]
        shg = shgload(value)
        CHI2[component] = shg

CHI2R = {
    'XXX' : + np.sin(GAMMA)**3*CHI2['xxx'] + np.sin(GAMMA)*np.cos(GAMMA)**2*CHI2['xyy'] - 2*np.sin(GAMMA)**2*np.cos(GAMMA)*CHI2['xxy'] \
            - np.sin(GAMMA)**2*np.cos(GAMMA)*CHI2['yxx'] - np.cos(GAMMA)**3*CHI2['yyy'] + 2*np.sin(GAMMA)*np.cos(GAMMA)**2*CHI2['yxy'],
    'XYY' : + np.sin(GAMMA)*np.cos(GAMMA)**2*CHI2['xxx'] + np.sin(GAMMA)**3*CHI2['xyy'] + 2*np.sin(GAMMA)**2*np.cos(GAMMA)*CHI2['xxy'] \
            - np.cos(GAMMA)**3*CHI2['yxx'] - np.sin(GAMMA)**2*np.cos(GAMMA)*CHI2['yyy'] - 2*np.sin(GAMMA)*np.cos(GAMMA)**2*CHI2['yxy'],
    'XZZ' : + np.sin(GAMMA)*CHI2['xzz'] - np.cos(GAMMA)*CHI2['yzz'],
    'XYZ' : + np.sin(GAMMA)**2*CHI2['xyz'] + np.sin(GAMMA)*np.cos(GAMMA)*CHI2['xxz'] \
            - np.sin(GAMMA)*np.cos(GAMMA)*CHI2['yyz'] - np.cos(GAMMA)**2*CHI2['yxz'],
    'XXZ' : - np.sin(GAMMA)*np.cos(GAMMA)*CHI2['xyz'] + np.sin(GAMMA)**2*CHI2['xxz'] \
            + np.cos(GAMMA)**2*CHI2['yyz'] - np.sin(GAMMA)*np.cos(GAMMA)*CHI2['yxz'],
    'XXY' : + np.sin(GAMMA)**2*np.cos(GAMMA)*CHI2['xxx'] - np.sin(GAMMA)**2*np.cos(GAMMA)*CHI2['xyy'] + (np.sin(GAMMA)**3 - np.sin(GAMMA)*np.cos(GAMMA)**2)*CHI2['xxy'] \
            - np.sin(GAMMA)*np.cos(GAMMA)**2*CHI2['yxx'] - np.sin(GAMMA)*np.cos(GAMMA)**2*CHI2['yyy'] + (np.cos(GAMMA)**3 - np.sin(GAMMA)**2*np.cos(GAMMA))*CHI2['yxy'],
    'YXX' : + np.sin(GAMMA)**2*np.cos(GAMMA)*CHI2['xxx'] + np.cos(GAMMA)**3*CHI2['xyy'] - 2*np.sin(GAMMA)*np.cos(GAMMA)**2*CHI2['xxy'] \
            + np.sin(GAMMA)**3*CHI2['yxx'] + np.sin(GAMMA)*np.cos(GAMMA)**2*CHI2['yyy'] - 2*np.sin(GAMMA)**2*np.cos(GAMMA)*CHI2['yxy'],
    'YYY' : + np.cos(GAMMA)**3*CHI2['xxx'] + np.sin(GAMMA)**2*np.cos(GAMMA)*CHI2['xyy'] + 2*np.sin(GAMMA)*np.cos(GAMMA)**2*CHI2['xxy'] \
            + np.sin(GAMMA)*np.cos(GAMMA)**2*CHI2['yxx'] + np.sin(GAMMA)**3*CHI2['yyy'] + 2*np.sin(GAMMA)**2*np.cos(GAMMA)*CHI2['yxy'],
    'YZZ' : + np.cos(GAMMA)*CHI2['xzz'] + np.sin(GAMMA)*CHI2['yzz'],
    'YYZ' : + np.sin(GAMMA)*np.cos(GAMMA)*CHI2['xyz'] + np.cos(GAMMA)**2*CHI2['xxz'] \
            + np.sin(GAMMA)**2*CHI2['yyz'] + np.sin(GAMMA)*np.cos(GAMMA)*CHI2['yxz'],
    'YXZ' : - np.cos(GAMMA)**2*CHI2['xyz'] + np.sin(GAMMA)*np.cos(GAMMA)*CHI2['xxz'] \
            - np.sin(GAMMA)*np.cos(GAMMA)*CHI2['yyz'] + np.sin(GAMMA)**2*CHI2['yxz'],
    'YXY' : + np.sin(GAMMA)*np.cos(GAMMA)**2*CHI2['xxx'] - np.sin(GAMMA)*np.cos(GAMMA)**2*CHI2['xyy'] - (np.cos(GAMMA)**3 - np.sin(GAMMA)**2*np.cos(GAMMA))*CHI2['xxy'] \
            + np.sin(GAMMA)**2*np.cos(GAMMA)*CHI2['yxx'] - np.sin(GAMMA)**2*np.cos(GAMMA)*CHI2['yyy'] + (np.sin(GAMMA)**3 - np.sin(GAMMA)*np.cos(GAMMA)**2)*CHI2['yxy'],
    'ZXX' : + np.sin(GAMMA)**2*CHI2['zxx'] + np.cos(GAMMA)**2*CHI2['zyy'] - 2*np.sin(GAMMA)*np.cos(GAMMA)*CHI2['zxy'],
    'ZYY' : + np.cos(GAMMA)**2*CHI2['zxx'] + np.sin(GAMMA)**2*CHI2['zyy'] + 2*np.sin(GAMMA)*np.cos(GAMMA)*CHI2['zxy'],
    'ZZZ' : + CHI2['zzz'],
    'ZYZ' : + np.sin(GAMMA)*CHI2['zyz'] + np.cos(GAMMA)*CHI2['zxz'],
    'ZXZ' : - np.cos(GAMMA)*CHI2['zyz'] + np.sin(GAMMA)*CHI2['zxz'],
    'ZXY' : + np.sin(GAMMA)*np.cos(GAMMA)*CHI2['zxx'] - np.sin(GAMMA)*np.cos(GAMMA)*CHI2['zyy'] - np.cos(2*GAMMA)*CHI2['zxy']
}

np.savetxt('xxx_rot{:03d}.txt'.format(int(angle)), np.column_stack((freq, CHI2R['XXX'].real, CHI2R['XXX'].imag)), fmt=('%05.2f % .4e % .4e'))
np.savetxt('xyy_rot{:03d}.txt'.format(int(angle)), np.column_stack((freq, CHI2R['XYY'].real, CHI2R['XYY'].imag)), fmt=('%05.2f % .4e % .4e'))
np.savetxt('xzz_rot{:03d}.txt'.format(int(angle)), np.column_stack((freq, CHI2R['XZZ'].real, CHI2R['XZZ'].imag)), fmt=('%05.2f % .4e % .4e'))
np.savetxt('xyz_rot{:03d}.txt'.format(int(angle)), np.column_stack((freq, CHI2R['XYZ'].real, CHI2R['XYZ'].imag)), fmt=('%05.2f % .4e % .4e'))
np.savetxt('xxz_rot{:03d}.txt'.format(int(angle)), np.column_stack((freq, CHI2R['XXZ'].real, CHI2R['XXZ'].imag)), fmt=('%05.2f % .4e % .4e'))
np.savetxt('xxy_rot{:03d}.txt'.format(int(angle)), np.column_stack((freq, CHI2R['XXY'].real, CHI2R['XXY'].imag)), fmt=('%05.2f % .4e % .4e'))
np.savetxt('yxx_rot{:03d}.txt'.format(int(angle)), np.column_stack((freq, CHI2R['YXX'].real, CHI2R['YXX'].imag)), fmt=('%05.2f % .4e % .4e'))
np.savetxt('yyy_rot{:03d}.txt'.format(int(angle)), np.column_stack((freq, CHI2R['YYY'].real, CHI2R['YYY'].imag)), fmt=('%05.2f % .4e % .4e'))
np.savetxt('yzz_rot{:03d}.txt'.format(int(angle)), np.column_stack((freq, CHI2R['YZZ'].real, CHI2R['YZZ'].imag)), fmt=('%05.2f % .4e % .4e'))
np.savetxt('yyz_rot{:03d}.txt'.format(int(angle)), np.column_stack((freq, CHI2R['YYZ'].real, CHI2R['YYZ'].imag)), fmt=('%05.2f % .4e % .4e'))
np.savetxt('yxz_rot{:03d}.txt'.format(int(angle)), np.column_stack((freq, CHI2R['YXZ'].real, CHI2R['YXZ'].imag)), fmt=('%05.2f % .4e % .4e'))
np.savetxt('yxy_rot{:03d}.txt'.format(int(angle)), np.column_stack((freq, CHI2R['YXY'].real, CHI2R['YXY'].imag)), fmt=('%05.2f % .4e % .4e'))
np.savetxt('zxx_rot{:03d}.txt'.format(int(angle)), np.column_stack((freq, CHI2R['ZXX'].real, CHI2R['ZXX'].imag)), fmt=('%05.2f % .4e % .4e'))
np.savetxt('zyy_rot{:03d}.txt'.format(int(angle)), np.column_stack((freq, CHI2R['ZYY'].real, CHI2R['ZYY'].imag)), fmt=('%05.2f % .4e % .4e'))
np.savetxt('zzz_rot{:03d}.txt'.format(int(angle)), np.column_stack((freq, CHI2R['ZZZ'].real, CHI2R['ZZZ'].imag)), fmt=('%05.2f % .4e % .4e'))
np.savetxt('zyz_rot{:03d}.txt'.format(int(angle)), np.column_stack((freq, CHI2R['ZYZ'].real, CHI2R['ZYZ'].imag)), fmt=('%05.2f % .4e % .4e'))
np.savetxt('zxz_rot{:03d}.txt'.format(int(angle)), np.column_stack((freq, CHI2R['ZXZ'].real, CHI2R['ZXZ'].imag)), fmt=('%05.2f % .4e % .4e'))
np.savetxt('zxy_rot{:03d}.txt'.format(int(angle)), np.column_stack((freq, CHI2R['ZXY'].real, CHI2R['ZXY'].imag)), fmt=('%05.2f % .4e % .4e'))
