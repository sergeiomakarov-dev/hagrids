#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create the Kisslinger surfaces of the simplified island divertor.

Three surfaces are created:

    IncidenceHamiltsym72    the target plate
    IncidenceHamiltsym72v2  the same plate, closed on the outside by a contour
    VesselHamilt_r70        the vessel

The target follows the island: over one field period its cross section is
rotated poloidally from -PHI_ROT to +PHI_ROT, and at the same time the plate
moves from RBOT at the two ends of the period to RTOP in the middle, where it
reaches deepest into the island.

All lengths are in cm and all angles in degrees, as the Kisslinger format
requires. A Kisslinger file has the structure

    # <name>
    <n_phi> <n_RZ> <nfold> <R_ref> <Z_ref>
    <phi>
    <R> <Z>                     n_RZ lines, repeated for every phi
    ...

Usage:
    python target_generation.py [output_dir]

Defaults:
    output_dir = output
"""
import os
import sys

import numpy as np

# Common
R0 = 590.0          # Major radius (cm)
NFOLD = 5           # Field periodicity
PHI_1 = 0.0         # First toroidal angle of the surfaces (deg)
PHI_2 = 72.0        # Last toroidal angle of the surfaces (deg)

# Target plate
RBOT = 65.25        # Plate radius at the two ends of the field period (cm)
ZBOT = 25.0         # Plate half height at the two ends of the field period (cm)
RTOP = 53.0         # Plate radius in the middle of the field period (cm)
ZTOP = 35.0         # Plate half height in the middle of the field period (cm)
RTAR_ADD = 8.0      # Radial length of the two side wings of the plate (cm)
RC_ADD = 4.0        # Radial offset of the outer closing contour (cm)
PHI_ROT = 36.0      # Poloidal rotation over half a field period (deg)
TR = 100            # Number of toroidal slices
N_CLOSURE = 50      # Number of points of the outer closing contour

# Vessel
R_VESSEL = 71.0     # Vessel minor radius (cm)
TR_VESSEL = 10      # Number of toroidal slices of the vessel
TP_VESSEL = 20      # Number of poloidal points of the vessel


def save_mesh(name, R, P, Z, output_dir, nfold=NFOLD, R_ref=0.0, Z_ref=0.0):
    """Write a mesh as a Kisslinger file."""
    file_path = os.path.join(output_dir, name)
    with open(file_path, 'w') as file:
        file.write(f"# {name}\n")
        file.write(f"{R.shape[0]:d} {R.shape[1]:d} {int(nfold):d} {R_ref:f} {Z_ref:f}\n")
        for j in range(R.shape[0]):
            file.write(f"{P[j, 0]:.5f}\n")
            for i in range(R.shape[1]):
                file.write(f"{R[j, i]:.5f} {Z[j, i]:.5f}\n")
    print(f"{name}: {R.shape[0]} toroidal x {R.shape[1]} (R, Z) points "
          f"written to '{file_path}'")


def plate_shape():
    """Radius and half height of the plate along the field period."""
    shape = np.cos(np.linspace(-np.pi / 2, np.pi / 2, TR))
    vecR = RBOT - (RBOT - RTOP) * shape
    vecZ = ZBOT - (ZBOT - ZTOP) * shape
    return vecR, vecZ


def rotate(vecRZ, phi_rot):
    """Rotate a cross section by phi_rot and shift it to the magnetic axis."""
    M = [[np.cos(phi_rot), -np.sin(phi_rot)],
         [np.sin(phi_rot), np.cos(phi_rot)]]
    vecRZ_rot = np.matmul(M, vecRZ)
    return vecRZ_rot[0] + R0, vecRZ_rot[1]


def target():
    """IncidenceHamiltsym72: side wing, plate, side wing."""
    vecR, vecZ = plate_shape()
    phi_rot = np.linspace(-PHI_ROT, PHI_ROT, TR) * np.pi / 180

    tmpR = np.empty((TR, 4))
    tmpZ = np.empty((TR, 4))
    for j in range(TR):
        vecRZ = [[vecR[j] + RTAR_ADD, vecR[j], vecR[j], vecR[j] + RTAR_ADD],
                 [-vecZ[j], -vecZ[j], vecZ[j], vecZ[j]]]
        tmpR[j], tmpZ[j] = rotate(vecRZ, phi_rot[j])

    tmpP = np.tile(np.linspace(PHI_1, PHI_2, TR), (4, 1)).T
    return tmpR, tmpP, tmpZ


def target_closed():
    """IncidenceHamiltsym72v2: the plate closed on the outside by a contour."""
    vecR, vecZ = plate_shape()
    phi_rot = np.linspace(-PHI_ROT, PHI_ROT, TR) * np.pi / 180

    num = 2 + N_CLOSURE
    tmpR = np.empty((TR, num))
    tmpZ = np.empty((TR, num))
    for j in range(TR):
        # 1e-3 keeps the corner of the plate away from arctan2(0, 0)
        phis = np.arctan2(vecZ[j] + 1e-3, vecR[j] + 1e-3)
        Rc = RBOT + RC_ADD
        add_R = Rc * np.cos(np.linspace(phis, 2 * np.pi - phis, N_CLOSURE))
        add_Z = Rc * np.sin(np.linspace(phis, 2 * np.pi - phis, N_CLOSURE))
        vecRZ = [np.append([vecR[j], vecR[j]], add_R),
                 np.append([-vecZ[j], vecZ[j]], add_Z)]
        tmpR[j], tmpZ[j] = rotate(vecRZ, phi_rot[j])

    tmpP = np.tile(np.linspace(PHI_1, PHI_2, TR), (num, 1)).T

    # Close the contour on itself
    tmpR = np.concatenate((tmpR, tmpR[:, :1]), axis=1)
    tmpZ = np.concatenate((tmpZ, tmpZ[:, :1]), axis=1)
    tmpP = np.concatenate((tmpP, tmpP[:, :1]), axis=1)
    return tmpR, tmpP, tmpZ


def vessel():
    """VesselHamilt_r70: a circular vessel of minor radius R_VESSEL."""
    vartheta = np.linspace(0, 2 * np.pi, TP_VESSEL)

    tmpR = np.tile(R0 + R_VESSEL * np.cos(vartheta), (TR_VESSEL, 1))
    tmpZ = np.tile(R_VESSEL * np.sin(vartheta), (TR_VESSEL, 1))
    tmpP = np.tile(np.linspace(PHI_1, PHI_2, TR_VESSEL), (TP_VESSEL, 1)).T
    return tmpR, tmpP, tmpZ


def main():
    output_dir = sys.argv[1] if len(sys.argv) > 1 else 'output'

    save_mesh('IncidenceHamiltsym72', *target(), output_dir)
    save_mesh('IncidenceHamiltsym72v2', *target_closed(), output_dir)
    save_mesh('VesselHamilt_r70', *vessel(), output_dir)


if __name__ == '__main__':
    main()
