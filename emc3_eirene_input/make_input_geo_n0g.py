#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Write the input.geo and input.N0G files of EMC3-EIRENE for the grid written by
write_emc3_grid.py.

The grid resolution is read from the resolution file written by
write_emc3_grid.py: three numbers (radial, poloidal, toroidal) for the one zone
grid, or nine numbers (the same three for the core, for an island and for the PFR)
for the island aligned grid, which has 1 + m + 1 zones in the order core, islands,
PFR. The field periodicity and the number of islands m are read from
input_params.dat (through input_params.py).

input.geo holds the grid resolution, the non default and the non transparent
surfaces of every zone and the physical cell definition. input.N0G holds the non
transparent surfaces for the neutrals, the additional physical cells for the
neutrals and the neutral source, and points to the ADD_SF_N0 file with the
additional surfaces (see target_generation.py).

Usage, from the repository root where input_params.dat is:
    python emc3_eirene_input/make_input_geo_n0g.py [resolution_file] [output_dir]

Defaults:
    resolution_file = output/resolution7z.dat (flag_type_init = 3)
                      output/resolution1z.dat (flag_type_init = 4)
    output_dir      = output
"""
import os
import sys

# input_params lives in the repository root, one level up
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import input_params as ip

# Density (cm^-3), electron and ion temperature (eV) and Mach number of the
# additional physical cells for the neutrals
CORE_CELLS_PLASMA = ('3.5E13  1000.  1000.  0.00', '3.5E13  800.  800.  0.00')
VACUUM_CELLS_PLASMA = '1.0E+07  0.1  0.1   0.00'


def read_resolution(file_path):
    """Read the resolution file: a list of (n_rad, n_pol, n_tor) per zone type."""
    with open(file_path) as file:
        numbers = [int(value) for value in file.read().split()]
    if len(numbers) not in (3, 9):
        raise ValueError(f"Expected 3 or 9 numbers in '{file_path}', got {len(numbers)}.")
    return [tuple(numbers[i:i + 3]) for i in range(0, len(numbers), 3)]


def zone_resolutions(zone_types):
    """Resolution of every zone in the zone order core, islands, PFR."""
    if len(zone_types) == 1:
        return list(zone_types)
    core, island, pfr = zone_types
    return [core] + [island] * ip.m + [pfr]


def surface(index, zone, kind, range_1, range_2):
    """One surface definition: its index, zone and kind, and its extent."""
    return (f"{index:>12d}{zone:>12d}{kind:>12d}\n"
            f"{0:>12d}{range_1:>12d}{0:>12d}{range_2:>12d}\n")


def input_geo(zones):
    """The text of input.geo."""
    n_zones = len(zones)
    lines = [f"* geometry information for EMC3\n"
             f"* SYMMETRY {ip.nfp}\n"
             f"*--------------------------------\n"
             f"*** 1. grid resolution\n"
             f"*--------------------------------\n"
             f"* number of zones/blocks\n"
             f"{n_zones}\n"
             f"* number of radial, poloidal and toroidal grid points\n"]
    for n_rad, n_pol, n_tor in zones:
        lines.append(f"{n_rad:>11d}{n_pol:>11d}{n_tor:>11d}\n")

    lines.append("*--------------------------------\n"
                 "*** 2. surface definitions\n"
                 "*--------------------------------\n"
                 "*** 2.1 non default surface\n")

    # Radial: the outer boundary of the core and of the islands, the inner
    # boundary of the PFR. None for the one zone grid.
    lines.append("* radial\n")
    if n_zones == 1:
        lines.append(f"{0:>12d}\n")
    else:
        lines.append(f"{n_zones:>12d}\n")
        for i_zone, (n_rad, n_pol, n_tor) in enumerate(zones):
            index = 0 if i_zone == n_zones - 1 else n_rad - 1
            lines.append(surface(index, i_zone, 3, n_pol - 2, n_tor - 2))

    # Poloidal and toroidal: both ends of every zone
    lines.append("* poloidal\n")
    lines.append(f"{2 * n_zones:>12d}\n")
    for i_zone, (n_rad, n_pol, n_tor) in enumerate(zones):
        lines.append(surface(0, i_zone, 1, n_rad - 2, n_tor - 2))
        lines.append(surface(n_pol - 1, i_zone, 1, n_rad - 2, n_tor - 2))

    lines.append("* toroidal\n")
    lines.append(f"{2 * n_zones:>12d}\n")
    for i_zone, (n_rad, n_pol, n_tor) in enumerate(zones):
        lines.append(surface(0, i_zone, 3, n_rad - 2, n_pol - 2))
        lines.append(surface(n_tor - 1, i_zone, 3, n_rad - 2, n_pol - 2))

    # Non transparent radial surfaces: the core boundary a few cells inside the
    # first zone, the inner boundary of every island, and the vacuum boundary a
    # few cells inside the outer boundary of the last zone
    lines.append("*** 2.2 non transparent surface (Boundary condition must be defined)\n"
                 "* radial\n")
    if n_zones == 1:
        n_rad, n_pol, n_tor = zones[0]
        lines.append(f"{2:>12d}\n")
        lines.append("* 1: CORE BOUNDARY 0.50000E+00\n")
        lines.append(surface(3, 0, 1, n_pol - 2, n_tor - 2))
        lines.append("* 3: VACUUM BOUNDARY\n")
        lines.append(surface(n_rad - 4, 0, -1, n_pol - 2, n_tor - 2))
    else:
        lines.append(f"{n_zones:>12d}\n")
        for i_zone, (n_rad, n_pol, n_tor) in enumerate(zones):
            if i_zone == 0:
                lines.append("* 1: CORE BOUNDARY 0.50000E+00\n")
                lines.append(surface(3, i_zone, 1, n_pol - 2, n_tor - 2))
            elif i_zone < n_zones - 1:
                lines.append("* 2: CORE BOUNDARY\n")
                lines.append(surface(0, i_zone, 1, n_pol - 2, n_tor - 2))
            else:
                lines.append("* 3: VACUUM BOUNDARY\n")
                lines.append(surface(n_rad - 4, i_zone, -1, n_pol - 2, n_tor - 2))

    lines.append("* poloidal\n"
                 f"{0:>12d}\n"
                 "* toroidal\n"
                 f"{0:>12d}\n"
                 "*** 2.3 plate surface (Bohm Boundary condition)\n"
                 "* radial\n"
                 f"{-1:>12d}\n"
                 "* poloidal\n"
                 f"{-1:>12d}\n"
                 "* toroidal\n"
                 f"{-1:>12d}\n"
                 "*--------------------------------\n"
                 "*** 3. physical cell definition\n"
                 "*--------------------------------\n")
    if n_zones == 1:
        lines.append(" 1  6\n")
    else:
        lines.append("1\n" + ' '.join(['10'] * n_zones) + "\n")
    lines.append("* run cell check?\n"
                 "F\n")
    return ''.join(lines)


def neutral_cells(zone, rad_1, rad_2, pol_1, pol_2, pol_step, tor_2, plasma):
    """One block of additional physical cells for the neutrals."""
    header = f"{2:>11d}{1 if plasma in CORE_CELLS_PLASMA else 0:>6d}\n"
    cells = (f"{zone:>6d}{rad_1:>5d}{rad_2:>6d}{1:>6d}{pol_1:>6d}{pol_2:>6d}"
             f"{pol_step:>7d}{0:>5d}{tor_2:>6d}{1:>6d}\n")
    return header + cells + plasma + "\n"


def input_n0g(zones):
    """The text of input.N0G."""
    n_zones = len(zones)
    lines = ["******** additional geometry and parameters for EIRENE ****\n"
             "*--------------------------------\n"
             "*** 1. non-transparent surfaces for neutral particles\n"
             "*--------------------------------\n"
             "*  non-transparent surfaces with informations about\n"
             "*  this surface being defined in EIRENE. The surface\n"
             "*  number must be indicated here.\n"
             "* radial\n"]

    # The inner boundary of every zone but the last and the outer boundary of the
    # last; for the one zone grid both boundaries of the single zone
    if n_zones == 1:
        n_rad, n_pol, n_tor = zones[0]
        lines.append(f"{2:>13d}\n")
        lines.append(surface(0, 0, -3, n_pol - 2, n_tor - 2))
        lines.append(surface(n_rad - 1, 0, -2, n_pol - 2, n_tor - 2))
    else:
        lines.append(f"{n_zones:>13d}\n")
        for i_zone, (n_rad, n_pol, n_tor) in enumerate(zones):
            if i_zone < n_zones - 1:
                lines.append(surface(0, i_zone, -3, n_pol - 2, n_tor - 2))
            else:
                lines.append(surface(n_rad - 1, i_zone, -2, n_pol - 2, n_tor - 2))

    lines.append("* poloidal\n"
                 f"{0:>12d}\n"
                 "* toroidal\n"
                 f"{0:>12d}\n"
                 "*--------------------------------\n"
                 "*** 2. DEFINE ADDITIONAL PHYSICAL CELLS FOR NEUTRALS\n"
                 "*--------------------------------\n"
                 "*   ZONE  R1    R2    DR    P1    P2    DP    T1    T2    DT\n"
                 "* ne       Te      Ti        M\n")

    # Two plasma blocks in the first radial cells of the core, then one vacuum
    # block per zone (the whole zone for the one zone grid, every zone but the
    # core for the island aligned grid)
    n_rad, n_pol, n_tor = zones[0]
    vacuum_zones = range(n_zones) if n_zones == 1 else range(1, n_zones)
    lines.append(f"{2 + len(vacuum_zones):>11d}{70:>12d}\n")
    lines.append(neutral_cells(0, 0, 2, 0, n_pol - 1, n_pol - 1, n_tor - 1, CORE_CELLS_PLASMA[0]))
    lines.append(neutral_cells(0, 2, 3, 0, n_pol - 1, n_pol - 1, n_tor - 1, CORE_CELLS_PLASMA[1]))
    for i_zone in vacuum_zones:
        n_rad, n_pol, n_tor = zones[i_zone]
        lines.append(neutral_cells(i_zone, 0, n_rad - 1, 0, n_pol - 1, 1, n_tor - 1,
                                   VACUUM_CELLS_PLASMA))

    lines.append("*--------------------------------\n"
                 "*** 3. Neutral Source distribution\n"
                 "*--------------------------------\n"
                 "* N0S NS_PLACE  NSSIDE\n"
                 f"{0:>12d}{0:>12d}{1:>11d}\n"
                 "*--------------------------------\n"
                 "*** 4. Additional surfaces\n"
                 "*--------------------------------\n"
                 "./ADD_SF_N0\n"
                 "*--------------------------------\n"
                 "*** 5. Neutral gas diagnostics\n"
                 "*** 5.1. Particle and energy densities for atoms and molecules\n"
                 f"{0:>12d}{0:>12d}{0:>12d}{0:>12d}\n"
                 "*** 5.2. Flux and spectrum on a given surface\n"
                 f"{0:>12d}\n")
    return ''.join(lines)


def write_input_files(zones, output_dir='output'):
    """Write input.geo and input.N0G for the given zone resolutions. Returns their paths."""
    paths = []
    for name, text in (('input.geo', input_geo(zones)), ('input.N0G', input_n0g(zones))):
        paths.append(os.path.join(output_dir, name))
        with open(paths[-1], 'w') as file:
            file.write(text)
    return paths


def main():
    if len(sys.argv) > 1:
        resolution_file = sys.argv[1]
    elif ip.flag_type_init == 3:
        resolution_file = os.path.join('output', 'resolution7z.dat')
    else:
        resolution_file = os.path.join('output', 'resolution1z.dat')
    output_dir = sys.argv[2] if len(sys.argv) > 2 else 'output'

    zones = zone_resolutions(read_resolution(resolution_file))
    print(f"{len(zones)} zone(s) read from '{resolution_file}'")

    for path in write_input_files(zones, output_dir):
        print(f"written '{path}'")


if __name__ == '__main__':
    main()
