#!/usr/bin/env python3
"""Convert LAMMPS NEB data files (final_neb.*) into a single extended XYZ
trajectory readable by OVITO.

Usage:
    python neb_to_xyz.py [work_dir] [-o output.xyz] [-s species]

    work_dir   : directory containing final_neb.* files (default: current dir)
    -o         : output file name (default: neb_trajectory.xyz)
    -s         : species symbol for all atoms (default: W)
"""
import os, argparse, glob
import numpy as np


def parse_lammps_data(filepath):
    """Parse a LAMMPS data file (orthogonal or triclinic) and return cell + positions."""
    with open(filepath, 'r') as f:
        lines = f.readlines()

    nb_atoms = None
    xlo, xhi, ylo, yhi, zlo, zhi = 0, 0, 0, 0, 0, 0
    xy, xz, yz = 0.0, 0.0, 0.0
    atoms = []
    in_atoms = False

    for line in lines:
        tokens = line.split()
        if not tokens:
            if in_atoms and len(atoms) == nb_atoms:
                break
            continue

        if len(tokens) >= 2 and tokens[1] == 'atoms':
            nb_atoms = int(tokens[0])
        elif len(tokens) >= 4 and tokens[2] == 'xlo':
            xlo, xhi = float(tokens[0]), float(tokens[1])
        elif len(tokens) >= 4 and tokens[2] == 'ylo':
            ylo, yhi = float(tokens[0]), float(tokens[1])
        elif len(tokens) >= 4 and tokens[2] == 'zlo':
            zlo, zhi = float(tokens[0]), float(tokens[1])
        elif len(tokens) >= 3 and tokens[-3] == 'xy' and tokens[-2] == 'xz' and tokens[-1] == 'yz':
            xy, xz, yz = float(tokens[0]), float(tokens[1]), float(tokens[2])
        elif tokens[0] == 'Atoms':
            in_atoms = True
            continue
        elif in_atoms and len(tokens) >= 5:
            atom_id = int(tokens[0])
            x, y, z = float(tokens[2]), float(tokens[3]), float(tokens[4])
            atoms.append((atom_id, x, y, z))
            if len(atoms) == nb_atoms:
                break

    # Sort by atom ID for consistency
    atoms.sort(key=lambda a: a[0])

    # Build cell matrix (LAMMPS triclinic convention)
    # a = (xhi-xlo, 0, 0)
    # b = (xy, yhi-ylo, 0)
    # c = (xz, yz, zhi-zlo)
    cell = np.array([
        [xhi - xlo, 0.0, 0.0],
        [xy, yhi - ylo, 0.0],
        [xz, yz, zhi - zlo]
    ])

    positions = np.array([[a[1], a[2], a[3]] for a in atoms])
    return nb_atoms, cell, positions


def write_extended_xyz(outfile, frames, species):
    """Write all frames as extended XYZ (OVITO-compatible)."""
    with open(outfile, 'w') as f:
        for i, (nb_atoms, cell, positions) in enumerate(frames):
            f.write(f'{nb_atoms}\n')
            # Extended XYZ Lattice property (row-major, 9 values)
            lat = (f'Lattice="{cell[0,0]:.10f} {cell[0,1]:.10f} {cell[0,2]:.10f} '
                   f'{cell[1,0]:.10f} {cell[1,1]:.10f} {cell[1,2]:.10f} '
                   f'{cell[2,0]:.10f} {cell[2,1]:.10f} {cell[2,2]:.10f}" '
                   f'Properties=species:S:1:pos:R:3 '
                   f'pbc="T T T" '
                   f'Frame={i}\n')
            f.write(lat)
            for j in range(nb_atoms):
                f.write(f'{species} {positions[j,0]:.10f} {positions[j,1]:.10f} {positions[j,2]:.10f}\n')


def main():
    parser = argparse.ArgumentParser(description='Convert LAMMPS NEB data files to extended XYZ trajectory')
    parser.add_argument('work_dir', nargs='?', default='.', help='Directory with final_neb.* files')
    parser.add_argument('-o', '--output', default='neb_trajectory.xyz', help='Output XYZ file')
    parser.add_argument('-s', '--species', default='W', help='Atomic species symbol')
    args = parser.parse_args()

    pattern = os.path.join(args.work_dir, 'final_neb.*')
    files = glob.glob(pattern)
    if not files:
        print(f'No final_neb.* files found in {args.work_dir}')
        return

    # Sort numerically by replica index
    files.sort(key=lambda f: int(f.split('.')[-1]))

    frames = []
    for filepath in files:
        nb_atoms, cell, positions = parse_lammps_data(filepath)
        frames.append((nb_atoms, cell, positions))
        print(f'  Read {os.path.basename(filepath)}: {nb_atoms} atoms')

    outpath = os.path.join(args.work_dir, args.output)
    write_extended_xyz(outpath, frames, args.species)
    print(f'Wrote {len(frames)} frames to {outpath}')


if __name__ == '__main__':
    main()
