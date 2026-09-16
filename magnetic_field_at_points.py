#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculate the magnetic field components at given (R, Z, phi) points.

Input parameters are read from input_params.dat (through input_params.py).
Point coordinates are read from a separate text file with three columns:

    R (m)   Z (m)   phi (rad)

Lines starting with '#' and empty lines are skipped.

Usage:
    python magnetic_field_at_points.py [coords_file] [output_file]

Defaults:
    coords_file = points_RZphi_example.dat
    output_file = output/magnetic_field_at_points.dat

The output file contains the input coordinates and the field components
in toroidal (Br, Bvartheta, Bvarphi) and cylindrical (BR, BZ, Bphi)
coordinates, together with the field magnitude B.
"""
import sys
import numpy as np

import input_params as ip
import magnetic_feild_calculation as mfc
from fun_time_independent_mh import R_Z_to_r_vartheta, field_perturbation_form

# Same choice of flags as in fun_time_independent_mh.py
flag_form = field_perturbation_form()
flag_toroildal = ip.flag_vartheta


def read_coordinates(file_path):
    """Read R, Z, phi columns from a text file."""
    data = np.loadtxt(file_path, comments='#', ndmin=2)
    if data.shape[1] != 3:
        raise ValueError(f"Expected 3 columns (R, Z, phi) in '{file_path}', got {data.shape[1]}.")
    return data[:, 0], data[:, 1], data[:, 2]


def magnetic_field_at_points(R_arr, Z_arr, phi_arr):
    """Calculate the magnetic field at the given (R, Z, phi) points."""
    n_points = len(R_arr)

    Br_arr = np.zeros(n_points)
    Bvartheta_arr = np.zeros(n_points)
    Bvarphi_arr = np.zeros(n_points)
    BR_arr = np.zeros(n_points)
    BZ_arr = np.zeros(n_points)
    B_arr = np.zeros(n_points)

    for i in range(n_points):
        r_temp, vartheta_temp = R_Z_to_r_vartheta(R_arr[i], Z_arr[i])

        Br_temp, Bvartheta_temp, Bvarphi_temp = \
            mfc.magnetic_field_components(r_temp, vartheta_temp, phi_arr[i], ip.B_0, ip.R_0,
                                          iota_a=ip.iota_a, iota_b=ip.iota_b,
                                          m=ip.m, m1=ip.m1, m2=ip.m2, m3=ip.m3,
                                          iota_res=ip.iota_res, iota_res1=ip.iota_res1,
                                          iota_res2=ip.iota_res2, iota_res3=ip.iota_res3,
                                          A=ip.A, A1=ip.A1, A2=ip.A2, A3=ip.A3, k=ip.k,
                                          A_arr=ip.A_arr[:ip.n_perturbation], m_arr=ip.m_arr[:ip.n_perturbation],
                                          iota_res_arr=ip.iota_res_arr[:ip.n_perturbation],
                                          flag_form=flag_form,
                                          flag_toroildal=flag_toroildal,
                                          flag_type_A_fun=ip.flag_type_A_fun)

        BR_temp = Br_temp * np.cos(vartheta_temp) - Bvartheta_temp * np.sin(vartheta_temp)
        BZ_temp = Br_temp * np.sin(vartheta_temp) + Bvartheta_temp * np.cos(vartheta_temp)

        Br_arr[i] = Br_temp
        Bvartheta_arr[i] = Bvartheta_temp
        Bvarphi_arr[i] = Bvarphi_temp
        BR_arr[i] = BR_temp
        BZ_arr[i] = BZ_temp
        B_arr[i] = np.sqrt(Br_temp**2 + Bvartheta_temp**2 + Bvarphi_temp**2)

    return Br_arr, Bvartheta_arr, Bvarphi_arr, BR_arr, BZ_arr, B_arr


def main():
    coords_file = sys.argv[1] if len(sys.argv) > 1 else 'points_RZphi_example.dat'
    output_file = sys.argv[2] if len(sys.argv) > 2 else 'output/magnetic_field_at_points.dat'

    R_arr, Z_arr, phi_arr = read_coordinates(coords_file)
    print(f"Read {len(R_arr)} points from '{coords_file}'")

    Br_arr, Bvartheta_arr, Bvarphi_arr, BR_arr, BZ_arr, B_arr = \
        magnetic_field_at_points(R_arr, Z_arr, phi_arr)

    header = (f"Magnetic field at points from '{coords_file}'\n"
              f"B_0 = {ip.B_0}, R_0 = {ip.R_0}, iota_a = {ip.iota_a}, iota_b = {ip.iota_b}, "
              f"flag_form = {flag_form}, flag_toroildal = {flag_toroildal}, "
              f"flag_type_A_fun = {ip.flag_type_A_fun}\n"
              f"{'R (m)':>13} {'Z (m)':>14} {'phi (rad)':>14} "
              f"{'Br (T)':>14} {'Bvartheta (T)':>14} {'Bvarphi (T)':>14} "
              f"{'BR (T)':>14} {'BZ (T)':>14} {'B (T)':>14}")

    data_out = np.column_stack((R_arr, Z_arr, phi_arr,
                                Br_arr, Bvartheta_arr, Bvarphi_arr,
                                BR_arr, BZ_arr, B_arr))
    np.savetxt(output_file, data_out, fmt='% .7e', header=header)
    print(f"Magnetic field written to '{output_file}'")

    for i in range(len(R_arr)):
        print(f"R = {R_arr[i]: .4f}, Z = {Z_arr[i]: .4f}, phi = {phi_arr[i]: .4f} -> "
              f"Br = {Br_arr[i]: .4e}, Bvartheta = {Bvartheta_arr[i]: .4e}, "
              f"Bvarphi = {Bvarphi_arr[i]: .4e}, |B| = {B_arr[i]: .4e}")


if __name__ == '__main__':
    main()
