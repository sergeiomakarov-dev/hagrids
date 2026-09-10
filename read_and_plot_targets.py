#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Read the Kisslinger surfaces written by target_generation.py and plot them.

Two figures are produced: a 3D view of the whole torus, where the field period
of every surface is repeated nfold times, and the poloidal cross sections of the
surfaces at a few toroidal angles.

All lengths are in cm and all angles in degrees, as in the Kisslinger files.

Usage:
    python read_and_plot_targets.py [output_dir]

Defaults:
    output_dir = output
"""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

SURFACES = ('IncidenceHamiltsym72', 'IncidenceHamiltsym72v2', 'VesselHamilt_r70')

COLORS = {'IncidenceHamiltsym72': 'tab:red',
          'IncidenceHamiltsym72v2': 'tab:blue',
          'VesselHamilt_r70': 'tab:gray'}

PHI_CUTS = (0.0, 18.0, 36.0, 54.0, 72.0)


def read_kisslinger(file_path):
    """Read a Kisslinger mesh. Returns R (cm), P (deg), Z (cm) and nfold."""
    with open(file_path) as file:
        lines = [line for line in file if line.strip()]

    header = lines[1].split()
    n_phi, n_rz, nfold = int(header[0]), int(header[1]), int(header[2])

    R = np.empty((n_phi, n_rz))
    Z = np.empty((n_phi, n_rz))
    P = np.empty((n_phi, n_rz))

    i_line = 2
    for j in range(n_phi):
        P[j, :] = float(lines[i_line])
        i_line += 1
        for i in range(n_rz):
            R[j, i], Z[j, i] = (float(value) for value in lines[i_line].split())
            i_line += 1

    return R, P, Z, nfold


def plot_torus(surfaces):
    """3D view of the surfaces, with every field period repeated nfold times.

    One panel per surface: drawn together, the closed surfaces hide the others.
    """
    figure = plt.figure(figsize=(5 * len(surfaces), 5.5))

    for i_surface, (name, (R, P, Z, nfold)) in enumerate(surfaces.items()):
        axes = figure.add_subplot(1, len(surfaces), i_surface + 1, projection='3d')

        for i_period in range(nfold):
            phi = np.deg2rad(P + i_period * 360.0 / nfold)
            X, Y = R * np.cos(phi), R * np.sin(phi)
            axes.plot_surface(X, Y, Z, color=COLORS[name], alpha=0.9,
                              linewidth=0, shade=True)

        axes.set_xlabel(r'X, cm', fontsize=11)
        axes.set_ylabel(r'Y, cm', fontsize=11)
        axes.set_zlabel(r'Z, cm', fontsize=11)
        axes.set_title(name, fontsize=13)
        axes.set_box_aspect((1, 1, 0.4))
        axes.view_init(elev=28, azim=-55)

    # tight_layout clips the z axis of a 3D panel against the next one
    figure.subplots_adjust(left=0.02, right=0.94, wspace=0.18)
    plt.show()


def plot_cross_sections(surfaces):
    """Poloidal cross sections of the surfaces at a few toroidal angles.

    The surfaces are drawn from the outermost inwards, so that the target stays
    visible where the closed surface runs along it.
    """
    figure, axes_list = plt.subplots(1, len(PHI_CUTS), figsize=(3.4 * len(PHI_CUTS), 4.6),
                                     sharex=True, sharey=True)

    widths = {'VesselHamilt_r70': 1.0, 'IncidenceHamiltsym72v2': 1.5,
              'IncidenceHamiltsym72': 3.0}

    for axes, phi_cut in zip(axes_list, PHI_CUTS):
        for name in ('VesselHamilt_r70', 'IncidenceHamiltsym72v2', 'IncidenceHamiltsym72'):
            R, P, Z, _ = surfaces[name]
            j = int(np.argmin(np.abs(P[:, 0] - phi_cut)))
            axes.plot(R[j], Z[j], color=COLORS[name], linewidth=widths[name],
                      label=name if axes is axes_list[0] else None)

        axes.set_aspect('equal')
        axes.set_title(r'$\varphi = {:.0f}^\circ$'.format(phi_cut), fontsize=13)
        axes.set_xlabel(r'R, cm', fontsize=12)

    axes_list[0].set_ylabel(r'Z, cm', fontsize=12)
    figure.legend(loc='lower center', ncol=len(SURFACES), fontsize=11)
    plt.tight_layout(rect=(0, 0.10, 1, 1))
    plt.show()


def main():
    output_dir = sys.argv[1] if len(sys.argv) > 1 else 'output'

    surfaces = {}
    for name in SURFACES:
        file_path = os.path.join(output_dir, name)
        surfaces[name] = read_kisslinger(file_path)
        R, P, _, nfold = surfaces[name]
        print(f"{name}: {R.shape[0]} toroidal x {R.shape[1]} (R, Z) points, "
              f"{nfold}-fold, phi from {P[0, 0]:.2f} to {P[-1, 0]:.2f} deg")

    plot_torus(surfaces)
    plot_cross_sections(surfaces)


if __name__ == '__main__':
    main()
