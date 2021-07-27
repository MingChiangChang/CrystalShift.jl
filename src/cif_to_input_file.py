# Module to turn cif file into input text file that
# Can be turned into CrystalPhase object with the constructor
# CrystalPhase(path_of_cif)
# File start with the crystal information 
# Crystal system, a=?, b=?, c=?, alpha=?, beta=?, gamma=?
# and follow by
# the information of the sticks as follows
# h k l intensity
# For each sticks

import glob
import os
from pathlib import Path

from pymatgen.io.cif import CifParser
import pymatgen.analysis.diffraction.xrd as pm_xrd
import numpy as np

def cif_to_input(cif_path, output_path, q_range, wvlen=0.15406):
    '''
    cif_to_input(cif_path, output_path)

    Main function that you should be interfacing with this modul 
    '''
    cif = CifParser(cif_path)
    cif_dict = cif.as_dict()
    with open(output_path / (os.path.basename(cif_path)[:-4]+'.csv'), 'w') as f:
        write_crystal_info(cif_dict, f)
        f.write('\n')
        return write_peaks_info(cif, f, q_range, wvlen)
        

def write_crystal_info(cif_dict, f):
    phase_name = _get_phase_name(cif_dict)
    info_dict = cif_dict[phase_name]
    f.write(phase_name)
    f.write(',')
    f.write(_get_crystal_system(info_dict))
    f.write(',')
    f.write(_get_lattice_parameters(info_dict))

def write_peaks_info(cif, f, q_range, wvlen):
    structure = cif.get_structures()[0]
    two_theta_range = q_to_two_theta(wvlen, *q_range)
    print(f"two theta range: {two_theta_range}")
    xrd_cal = pm_xrd.XRDCalculator() # Instantiate calculator class
    diff_pattern = xrd_cal.get_pattern(structure=structure,
                                       two_theta_range=two_theta_range)
    print(diff_pattern.hkls)
    return diff_pattern


def q_to_two_theta(wvlen, *args):
    two_thetas = []
    for q in args:
        print(q)
        print(wvlen*q/(4*np.pi))
        two_thetas.append(360*np.arcsin(wvlen*q/(4*np.pi))/np.pi)
    return tuple(two_thetas)

def _get_phase_name(cif_dict):
    return list(cif_dict.keys())[0]

def _get_lattice_parameters(info_dict):
    a, b, c = _get_cell_length(info_dict)
    alpha, beta, gamma = _get_cell_angle(info_dict)
    return f"{a},{b},{c},{alpha},{beta},{gamma}"

def _get_cell_length(info_dict):
    return (info_dict['_cell_length_a'],
            info_dict['_cell_length_b'],
            info_dict['_cell_length_c'])

def _get_cell_angle(info_dict):
    return (info_dict['_cell_angle_alpha'],
            info_dict['_cell_angle_beta'],
            info_dict['_cell_angle_gamma'])

def _get_crystal_system(info_dict):
    # TODO turn space group info into crystal system
    return info_dict['_symmetry_space_group_name_H-M']

#if __name__ == "__main__":

home = Path.home()
path = home / 'Desktop' / 'github' /\
        'Crystallography_based_shifting' / 'data' 
cif_path = list(path.glob('*computed.cif'))[0]
out_path = path 
print(cif_path)
diff = cif_to_input(cif_path, out_path, (0, 45))

