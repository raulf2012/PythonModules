"""
"""

#| - Import Modules
import os
import sys

import re
import numpy as np

from ase import Atoms
from ase.utils import basestring, reader, writer
from ase.io.utils import ImageIterator, ImageChunk

import ase

#  from ase.io.vasp import (
#      iread_vasp_out,
#      #  outcarchunks, _read_outcar_header, _read_outcar_frame,
#      #  _outcar_check_line, OUTCARChunk,
#      #  _OUTCAR_SCF_DELIM,
#      )

#  iread_vasp_out

#__|


#  Denotes end of Ionic step for OUTCAR reading
_OUTCAR_SCF_DELIM = 'FREE ENERGIE OF THE ION-ELECTRON SYSTEM'


@reader
def read_vasp_out(filename='OUTCAR', index=-1):
    """Import OUTCAR type file.

    Reads unitcell, atom positions, energies, and forces from the OUTCAR file
    and attempts to read constraints (if any) from CONTCAR/POSCAR, if present.
    """
    #| - read_vasp_out
    f = filename
    g = iread_vasp_out(f, index=index)
    # Code borrowed from formats.py:read
    if isinstance(index, (slice, basestring)):
        # Return list of atoms
        return list(g)
    else:
        # Return single atoms object
        return next(g)
    #__|

def _read_outcar_frame(lines, header_data):
    #| - _read_outcar_frame
    from ase.calculators.singlepoint import (SinglePointDFTCalculator,
                                             SinglePointKPoint)

    mag_x = None
    mag_y = None
    mag_z = None
    magmoms = None
    magmom = None
    stress = None
    efermi = None

    symbols = header_data['symbols']
    constraints = header_data['constraints']
    natoms = header_data['natoms']
    # nkpts = header_data['nkpts']
    nbands = header_data['nbands']
    kpt_weights = header_data['kpt_weights']
    ibzkpts = header_data['ibzkpts']

    atoms = Atoms(symbols=symbols, pbc=True, constraint=constraints)

    cl = _outcar_check_line     # Aliasing

    spinc = 0                   # Spin component
    kpts = []
    forces = np.zeros((natoms, 3))
    positions = np.zeros((natoms, 3))
    f_n = np.zeros(nbands)      # kpt occupations
    eps_n = np.zeros(nbands)    # kpt eigenvalues

    # Parse each atoms object
    for n, line in enumerate(lines):
        line = line.strip()
        if 'direct lattice vectors' in line:
            cell = []
            for i in range(3):
                parts = cl(lines[n + i + 1]).split()
                cell += [list(map(float, parts[0:3]))]
            atoms.set_cell(cell)
        elif 'magnetization (x)' in line:
            # Magnetization in both collinear and non-collinear
            nskip = 4           # Skip some lines
            mag_x = [float(cl(lines[n + i + nskip]).split()[-1])
                     for i in range(natoms)]

        # XXX: !!!Uncomment these lines when non-collinear spin is supported!!!
        # Remember to check that format fits!

        # elif 'magnetization (y)' in line:
        #     # Non-collinear spin
        #     nskip = 4           # Skip some lines
        #     mag_y = [float(cl(lines[n + i + nskip]).split()[-1])
        #              for i in range(natoms)]
        # elif 'magnetization (z)' in line:
        #     # Non-collinear spin
        #     nskip = 4           # Skip some lines
        #     mag_z = [float(cl(lines[n + i + nskip]).split()[-1])
        #              for i in range(natoms)]
        elif 'number of electron' in line:
            parts = cl(line).split()
            if len(parts) > 5 and parts[0].strip() != "NELECT":
                i = parts.index('magnetization') + 1
                magmom = parts[i:]
                if len(magmom) == 1:
                    # Collinear spin
                    magmom = float(magmom[0])
                # XXX: !!!Uncomment these lines when non-collinear spin is supported!!!
                # Remember to check that format fits!
                # else:
                #     # Non-collinear spin
                #     # Make a (3,) dim array
                #     magmom = np.array(list(map(float, magmom)))
        elif 'in kB ' in line:

            if "*" in line:
                stress = None
            else:
                stress = -np.asarray([float(a) for a in cl(line).split()[2:]])
                stress = stress[[0, 1, 2, 4, 5, 3]] * 1e-1 * ase.units.GPa

        elif 'POSITION          ' in line:
            nskip = 2
            for i in range(natoms):
                parts = list(map(float, cl(lines[n + i + nskip]).split()))
                positions[i] = parts[0:3]
                forces[i] = parts[3:6]
            atoms.set_positions(positions, apply_constraint=False)
        elif 'E-fermi :' in line:
            parts = line.split()

            if "*" in parts[2]:
                efermi = None
            else:
                efermi = float(parts[2])
        elif 'spin component' in line:
            # Update spin component for kpts
            # Make spin be in [0, 1], VASP writes 1 or 2
            tmp = int(line.split()[-1]) - 1
            if tmp < spinc:
                # if NWRITE=3, we write KPTS after every electronic step,
                # so we just reset it, since we went from spin=2 to spin=1
                # in the same ionic step.
                # XXX: Only read it at last electronic step
                kpts = []
            spinc = tmp
        elif 'k-point  ' in line:
            if 'plane waves' in line:
                # Can happen if we still have part of header
                continue
            # Parse all kpts and bands
            parts = line.split()
            ikpt = int(parts[1]) - 1  # Make kpt idx start from 0
            w = kpt_weights[ikpt]

            nskip = 2
            for i in range(nbands):
                parts = lines[n + i + nskip].split()
                eps_n[i] = float(parts[1])
                f_n[i] = float(parts[2])
            kpts.append(SinglePointKPoint(w, spinc, ikpt,
                                          eps_n=eps_n, f_n=f_n))
        elif _OUTCAR_SCF_DELIM in line:
            # Last section before next ionic step
            nskip = 2
            parts = cl(lines[n + nskip]).strip().split()
            energy_free = float(parts[4])  # Force consistent

            nskip = 4
            parts = cl(lines[n + nskip]).strip().split()
            energy_zero = float(parts[6])  # Extrapolated to 0 K

            # For debugging
            # assert len(kpts) == 0 or len(kpts) == (spinc + 1) * nkpts

            if mag_x is not None:
                if mag_y is not None:
                    # Non-collinear
                    assert len(mag_x) == len(mag_y) == len(mag_z)
                    magmoms = np.zeros((len(atoms), 3))
                    magmoms[:, 0] = mag_x
                    magmoms[:, 1] = mag_y
                    magmoms[:, 2] = mag_z
                else:
                    # Collinear
                    magmoms = np.array(mag_x)

            atoms.set_calculator(
                SinglePointDFTCalculator(atoms,
                                         energy=energy_zero,
                                         free_energy=energy_free,
                                         ibzkpts=ibzkpts,
                                         forces=forces,
                                         efermi=efermi,
                                         magmom=magmom,
                                         magmoms=magmoms,
                                         stress=stress))
            atoms.calc.name = 'vasp'
            atoms.calc.kpts = kpts
    return atoms
    #__|



#########################################################
def iread_vasp_out(filename, index=-1):
    """Import OUTCAR type file, as a generator."""
    #| - iread_vasp_out
    it = ImageIterator(outcarchunks)
    return it(filename, index=index)
    #__|

def outcarchunks(fd):
    #| -  outcarchunks
    # First we get header info
    header_data = _read_outcar_header(fd)

    while True:
        try:
            # Build chunk which contains 1 complete atoms object
            lines = []
            while True:
                line = next(fd)
                lines.append(line)
                if _OUTCAR_SCF_DELIM in line:
                    # Add 4 more lines to include energy
                    for _ in range(4):
                        lines.append(next(fd))
                    break
        except StopIteration:
            # End of file
            return
        yield OUTCARChunk(lines, header_data)
    #__|

def _read_outcar_header(fd):
    #| - _read_outcar_header

    # Get the directory of the OUTCAR we are reading
    wdir = os.path.dirname(fd.name)
    # Try and see if we can get constraints
    if os.path.isfile(os.path.join(wdir, 'CONTCAR')):
        constraints = read_vasp(os.path.join(wdir, 'CONTCAR')).constraints
    elif os.path.isfile(os.path.join(wdir, 'POSCAR')):
        constraints = read_vasp(os.path.join(wdir, 'POSCAR')).constraints
    else:
        constraints = None

    cl = _outcar_check_line     # Aliasing

    species = []
    natoms = 0
    species_num = []
    symbols = []
    nkpts = 0
    nbands = 0
    kpt_weights = []
    ibzkpts = []

    # Get atomic species
    for line in fd:
        line = line.strip()
        if 'POTCAR:' in line:
            temp = line.split()[2]
            for c in ['.', '_', '1']:
                if c in temp:
                    temp = temp[0:temp.find(c)]
            species += [temp]
        elif 'ions per type' in line:
            species = species[:len(species) // 2]
            parts = cl(line).split()
            ntypes = min(len(parts) - 4, len(species))
            for ispecies in range(ntypes):
                species_num += [int(parts[ispecies + 4])]
                natoms += species_num[-1]
                for iatom in range(species_num[-1]):
                    symbols += [species[ispecies]]
        elif 'NKPTS' in line:
            parts = cl(line).split()
            nkpts = int(parts[3])
            nbands = int(parts[-1])
        elif 'k-points in reciprocal lattice and weights' in line:
            # Get kpoint weights
            for _ in range(nkpts):
                parts = next(fd).strip().split()
                ibzkpts.append(list(map(float, parts[0:3])))
                kpt_weights.append(float(parts[-1]))

        elif 'Iteration' in line:
            # Start of SCF cycle
            header_data = dict(
                natoms=natoms,
                symbols=symbols,
                constraints=constraints,
                nkpts=nkpts,
                nbands=nbands,
                kpt_weights=np.array(kpt_weights),
                ibzkpts=np.array(ibzkpts))
            return header_data

    # Incomplete OUTCAR, we can't determine atoms
    raise IOError('Incomplete OUTCAR')
    #__|

def _outcar_check_line(line):
    """Auxiliary check line function for OUTCAR numeric formatting.
    See issue #179, https://gitlab.com/ase/ase/issues/179
    Only call in cases we need the numeric values
    """
    #| - _outcar_check_line
    if re.search('[0-9]-[0-9]', line):
        line = re.sub('([0-9])-([0-9])', r'\1 -\2', line)
    return line
    #__|

class OUTCARChunk(ImageChunk):
    """
    """
    #| - OUTCARChunk
    def __init__(self, lines, header_data):
        self.lines = lines
        self.header_data = header_data

    def build(self):
        return _read_outcar_frame(self.lines, self.header_data)
    #__|

