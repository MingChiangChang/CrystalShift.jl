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
import argparse

from pymatgen.io.cif import CifParser
from xrayutilities.materials.cif import CIFFile
from xrayutilities.materials.material import Crystal
from xrayutilities.simpack import PowderDiffraction
import numpy as np

def main():
    parser = get_parser()
    args = parser.parse_args()
    cifs = list(Path(args.cif).glob('*.cif'))
    if not args.outpath.endswith('.csv'):
        args.outpath += '.csv'
    cif_to_input(cifs, args.outpath, (float(args.qmin), float(args.qmax)), args.wvlen)


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--cif', required=True, help='A folder that contain all the CIF files')
    parser.add_argument('-w', '--wvlen', default=1.5406, help='Wavelength in angstrom of X-ray. Default to Cu Kα')
    parser.add_argument('-o', '--outpath', required=True, help='Output path for csv input file for CrystalShift')
    parser.add_argument('-qmin', '--qmin', default=10., help='Minimum Q value in nm-1')
    parser.add_argument('-qmax', '--qmax', default=80., help='Maximum Q value in nm-1')
    return parser



def cif_to_input(cif_paths, output_path, q_range, wvlen=1.5406, _type=float):
    '''
    cif_to_input(cif_path, output_path, q_range, output_name='sticks')

    Main function that you should be interfacing with this module
    '''
    with open(output_path, 'w') as f:
        for idx, cif_path in enumerate(cif_paths):
            print(cif_path)
            cif = CifParser(cif_path)
            lattice = CIFFile(cif_path).SGLattice()
            write_cif(f, idx, cif, lattice, _type, q_range, wvlen)

def write_cif(f, idx, cif, lattice, _type, q_range, wvlen):
    f.write(f'{idx},')
    write_crystal_info(f, cif, _type)
    write_peaks_info(f, lattice, q_range, wvlen)

def write_crystal_info(f, cif, _type):
    cif_dict = cif.as_dict()
    key = _get_key(cif_dict)
    info_dict = cif_dict[key]
    f.write(_get_phase_name(info_dict))
    f.write(',')
    f.write(_get_crystal_system(info_dict))
    f.write(',')
    f.write(_get_lattice_parameters(info_dict, _type))

def write_peaks_info(f, lattice, q_range, wvlen):
    crystal = Crystal('test', lattice)
    xrd = PowderDiffraction(crystal, wl=wvlen).data

    for i, peak in enumerate(xrd):
        q = xrd[peak]['qpos']*10
        I = xrd[peak]['r']
        print(i, q, I, flush=True)
        if q_range[0] < q < q_range[1] and I>0.0001:
            f.write(f'\n{peak[0]},{peak[1]},{peak[2]},{q},{I}')
    f.write('#\n')

def q_to_two_theta(wvlen, *args):
    two_thetas = []
    for q in args:
        two_thetas.append(360*np.arcsin(wvlen*q/(4*np.pi))/np.pi)
    return tuple(two_thetas)

def _get_key(cif_dict):
    return list(cif_dict.keys())[0]

def _get_phase_name(info_dict):
    if "_chemical_formula_structural" in info_dict.keys():
        phase_name = remove_blank(info_dict["_chemical_formula_structural"])
    elif "_chemical_formula_moiety" in info_dict.keys():
        phase_name = remove_blank(info_dict["_chemical_formula_moiety"])
    elif "_chemical_formula_sum" in info_dict.keys():
        phase_name = remove_blank(info_dict["_chemical_formula_sum"])
    else:
        print("cannot find phase name in cif")
    try:
        if "_symmetry_space_group_name_H-M" in info_dict:
            space_group = remove_blank(info_dict["_symmetry_space_group_name_H-M"])
        elif "_space_group_name_H-M_alt" in info_dict:
            space_group = remove_blank(info_dict["_space_group_name_H-M_alt"])
    except KeyError:
        print(f"No _space_group_name_H-M_alt for {phase_name}")
        space_group = ""
    return phase_name + "_" + space_group

def _get_lattice_parameters(info_dict, _type):
    a, b, c = _get_cell_length(info_dict)
    alpha, beta, gamma = _get_cell_angle(info_dict)
    print(a, b, c, alpha, beta, gamma)
    return (f"{_type(a)},{_type(b)},{_type(c)},"
            f"{_type(alpha)},{_type(beta)},{_type(gamma)}")

def _get_cell_length(info_dict):
    return (remove_parentheses(info_dict['_cell_length_a']),
            remove_parentheses(info_dict['_cell_length_b']),
            remove_parentheses(info_dict['_cell_length_c']))

def remove_parentheses(string):
    if '(' in string:
        return string[:string.index('(')]
    return string

def remove_blank(string):
    while ' ' in string:
        string = string.replace(' ', '')
    return string

def _get_cell_angle(info_dict):
    return (remove_parentheses(info_dict['_cell_angle_alpha']),
            remove_parentheses(info_dict['_cell_angle_beta']),
            remove_parentheses(info_dict['_cell_angle_gamma']))

def _get_crystal_system(info_dict):
    try:
        sg_num = int(info_dict['_space_group_IT_number'])
    except KeyError:
        print("No space group info in cif. Default to triclinic")
        return "triclinic"
    if sg_num in [1,2]:
        return "triclinic"
    elif 3 <= sg_num <= 15:
        return "monoclinic"
    elif 16 <= sg_num <= 74:
        return "orthohombic"
    elif 75 <= sg_num <= 142:
        return "tetragonal"
    elif 143 <= sg_num <= 167:
        return "trigonal"
    elif 168 <= sg_num <= 194:
        return "hexagonal"
    elif 195 <= sg_num <= 230:
        return "cubic"

if __name__ == "__main__":
    main()