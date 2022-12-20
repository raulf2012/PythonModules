"""Methods to run VASP-Raman jobs."""

# | - IMPORT MODULES
import os
import sys
from math import sqrt

from ase import io
import numpy as np

from plotly.graph_objs import Scatter
import re
# __|

# | - Methods from vasp_raman script (Github)

def get_modes_from_OUTCAR(outcar_fh, nat=None, free_nat=None, path_i="."):
    # | - get_modes_from_OUTCAR
    if nat == None:

        name = "OUTCAR.phon"
        if os.path.isfile(path_i + "/" + name):
            file_name = name
        elif os.path.isfile(path_i + "/" + "OUTCAR"):
            file_name = "OUTCAR"
        else:
            file_name = None

        atoms = io.read(path_i + "/" + file_name)
        nat = atoms.get_number_of_atoms()

        # nat = vp.num_of_atoms_OUTCAR(outcar_fh)

    if free_nat == None:
        free_nat = nat

    eigvals = [ 0.0 for i in range(nat*3) ]
    eigvecs = [ 0.0 for i in range(nat*3) ]
    norms   = [ 0.0 for i in range(nat*3) ]

    outcar_fh.seek(0) # just in case
    while True:
        line = outcar_fh.readline()

        if not line:
            break

        if "Eigenvectors after division by SQRT(mass)" in line:
            outcar_fh.readline() # empty line
            outcar_fh.readline() # Eigenvectors and eigenvalues of the dynamical matrix
            outcar_fh.readline() # ----------------------------------------------------
            outcar_fh.readline() # empty line

            for i in range(free_nat*3): # all frequencies should be supplied, regardless of those requested to calculate

                outcar_fh.readline() # empty line
                line_i = outcar_fh.readline()
                p = re.search(r"^\s*(\d+).+?([\.\d]+) cm-1", line_i)

                eigvals[i] = float(p.group(2))

                outcar_fh.readline() # X         Y         Z           dx          dy          dz
                eigvec = []

                for j in range(nat):
                    tmp = outcar_fh.readline().split()

                    eigvec.append([ float(tmp[x]) for x in range(3,6) ])

                eigvecs[i] = eigvec
                norms[i] = sqrt( sum( [abs(x)**2 for sublist in eigvec for x in sublist] ) )

            return eigvals, eigvecs, norms

    print("[get_modes_from_OUTCAR]: ERROR Couldn't find 'Eigenvectors after",
        " division by SQRT(mass)'' in OUTCAR.",
        " Use 'NWRITE=3' in INCAR. Exiting...")
    sys.exit(1)
    # __|

def parse_env_params(params):
    # | - parse_env_params
    tmp = params.strip().split("_")
    if len(tmp) != 4:
        print("[parse_env_params]: ERROR there should be exactly four parameters")
        sys.exit(1)

    [first, last, nderiv, step_size] = [int(tmp[0]), int(tmp[1]), int(tmp[2]), float(tmp[3])]

    return first, last, nderiv, step_size
    # __|

def parse_poscar(poscar_fh, cons_atoms=0):
    # | - parse_poscar
    # modified subroutine from phonopy 1.8.3 (New BSD license)

    poscar_fh.seek(0) # just in case
    lines = poscar_fh.readlines()

    scale = float(lines[1])
    if scale < 0.0:
        print("[parse_poscar]: ERROR negative scale not implemented.")
        sys.exit(1)

    b = []
    for i in range(2, 5):
        b.append([float(x)*scale for x in lines[i].split()[:3]])

    vol = b[0][0]*b[1][1]*b[2][2] + b[1][0]*b[2][1]*b[0][2] + b[2][0]*b[0][1]*b[1][2] - \
          b[0][2]*b[1][1]*b[2][0] - b[2][1]*b[1][2]*b[0][0] - b[2][2]*b[0][1]*b[1][0]

    try:
        num_atoms = [int(x) for x in lines[5].split()]
        line_at = 6
    except ValueError:
        symbols = [x for x in lines[5].split()]
        num_atoms = [int(x) for x in lines[6].split()]
        line_at = 7

    nat = sum(num_atoms)  # TEMP
    free_nat = sum(num_atoms) - cons_atoms  # TEMP

    if lines[line_at][0].lower() == "s":
        line_at += 1

    if (lines[line_at][0].lower() == "c" or lines[line_at][0].lower() == "k"):
        is_scaled = False
    else:
        is_scaled = True

    line_at += 1

    positions = []
    for i in range(line_at, line_at + nat):
        pos = [float(x) for x in lines[i].split()[:3]]

        if is_scaled:
            pos = MAT_m_VEC(T(b), pos)

        positions.append(pos)

    poscar_header = "".join(lines[1:line_at-1]) # will add title and "Cartesian" later
    return nat, free_nat, vol, b, positions, poscar_header
    # __|

def MAT_m_VEC(m, v):
    # | - MAT_m_VEC
    p = [ 0.0 for i in range(len(v)) ]
    for i in range(len(m)):
        assert len(v) == len(m[i]), "Length of the matrix row is not equal to the length of the vector"
        p[i] = sum( [ m[i][j]*v[j] for j in range(len(v)) ] )
    return p
    # __|

def T(m):
    # | - T
    p = [[ m[i][j] for i in range(len( m[j] )) ] for j in range(len( m )) ]
    return p
    # __|

def get_epsilon_from_OUTCAR(outcar_fh):
    # | - get_epsilon_from_OUTCAR
    epsilon = []

    outcar_fh.seek(0) # just in case
    while True:
        line = outcar_fh.readline()
        if not line:
            break

        if "MACROSCOPIC STATIC DIELECTRIC TENSOR" in line:
            print(line)
            outcar_fh.readline()
            epsilon.append([float(x) for x in outcar_fh.readline().split()])
            epsilon.append([float(x) for x in outcar_fh.readline().split()])
            epsilon.append([float(x) for x in outcar_fh.readline().split()])
            return epsilon

    raise RuntimeError("[get_epsilon_from_OUTCAR]: ERROR Couldn't find dielectric tensor in OUTCAR")
    return 1
    # __|


# | - Substitute Functions for VTST

def parse_freqdat(freqdat_fh, nat):
    # | - parse_freqdat
    freqdat_fh.seek(0) # just in case
    #
    eigvals = [ 0.0 for i in range(nat*3) ]
    #
    for i in range(nat*3): # all frequencies should be supplied, regardless of requested to calculate
        tmp = freqdat_fh.readline().split()
        eigvals[i] = float(tmp[0])
    #
    return eigvals
    # __|

def parse_modesdat(modesdat_fh, nat):
    # | - parse_modesdat
    # from math import sqrt
    modesdat_fh.seek(0) # just in case
    #
    eigvecs = [ 0.0 for i in range(nat*3) ]
    norms =   [ 0.0 for i in range(nat*3) ]
    #
    for i in range(nat*3): # all frequencies should be supplied, regardless of requested to calculate
        eigvec = []
        for j in range(nat):
            tmp = modesdat_fh.readline().split()
            eigvec.append([ float(tmp[x]) for x in range(3) ])
        #
        modesdat_fh.readline().split() # empty line
        eigvecs[i] = eigvec
        norms[i] = sqrt( sum( [abs(x)**2 for sublist in eigvec for x in sublist] ) )
    #
    return eigvecs, norms
    # __|

# __|

# __|

def modes_list(num_modes, step, mode_0=1, modes_to_run="All"):
    """Returns list of strings representing number of modes to run per job

    Args:
        atoms:
        modes_to_run:
    """
    # | - modes_list


    # atoms = io.read("POSCAR.phon")
    # num_atoms = atoms.get_number_of_atoms()
    # max_modes = 3 * num_atoms


    if modes_to_run == "All":
        high_bound = mode_0 + num_modes - 1
    else:
        high_bound = mode_0 + modes_to_run - 1

    mode_range = [mode_0, high_bound]

    modes_list = []
    num_1 = mode_range[0]
    num_2 = 0

    while num_2 < mode_range[1] - step:
        num_2 = num_1 + step

        modes = str(num_1).zfill(2) + "-" + str(num_2).zfill(2)
        modes_list.append(modes)
        num_1 = num_2 + 1

    if num_2 != mode_range[1]:
        modes_rem = mode_range[1] - num_2
        last_modes = str(num_2 + 1).zfill(2) + "-" + str(num_2 + modes_rem)

        modes_list.append(last_modes)

    # for ind, modes in enumerate(modes_list):
    #     rev_entry = "Modes_" + modes + "_2_0.01"
    #     modes_list[ind] = rev_entry


    return(modes_list)
    # __|

def vasp_raman_input_file(modes, dir="."):
    """

    Args:
        modes: <type 'str'>
        dir: <type 'str'>
            Directory to write the vasp_raman input file to
    """
    # | - vasp_raman_input_file
    VASP_RAMAN_RUN = "python 3_vasp_raman_run_ibrion_n1.py"

    #| Changing "-" character in mode to "_"
    modes = modes.replace("-", "_")

    # Create VASP_RAMAN_RUN and VASP_RAMAN_PARAMS file
    VASP_RAMAN_PARAMS = modes + "_2_0.01"

    filename = "param_vasp_raman"

    with open(dir + "/" + filename, "w") as file:
        file.write(VASP_RAMAN_PARAMS + "\n" + VASP_RAMAN_RUN + "\n")
    # __|

def concatenate_mode_files():
    """
    """
    # | - concatenate_mode_filess
    tmp = 7

    # __|

def to_plot(hw, ab, gam=0.05, type="Gaussian", scaling=1.):
    """

    Args:
        hw:
        ab:
        gam:
        type:
        scaling:
    """
    # | - to_plot
    ab /= np.max(np.abs(ab), axis=0)

    ab = ab * scaling

    fmin = min(hw)
    fmax = max(hw)
    erange = np.arange(fmin - 40. * gam, fmax + 40. * gam, gam / 10.)  #np.arange(fmin-40*gam,fmax+40*gam,gam/10)
    spectrum = 0.0 * erange
    for i in range(len(hw)):
        if type=="Gaussian":

            term_1 = ab[i] * (2 * np.pi) ** (-0.5) / gam
            term_2 = np.exp(np.clip(-1.0 * (hw[i] - erange) ** 2/(2 * gam ** 2), -300, 300))
            spectrum += term_1 * term_2

            # spectrum += (2 * np.pi) ** (-0.5) / gam * np.exp(np.clip(-1.0 * (hw[i] - erange) ** 2/(2 * gam ** 2), -300, 300))

        elif type=="Lorentzian":
            spectrum += ab[i]*1/np.pi*gam/((hw[i]-erange)**2+gam**2)

    return erange, spectrum
    # __|

def proc_data(
    data,
    add_to_label="",
    scaling=1.,
    freq_shift=0.0,
    gauss_gamma=8.0,
    color_palette=None,
    ):
    """Process data.

    Args:
        data:
        add_to_label="":
        scaling=1.:
        freq_shift=0.0:
        gauss_gamma=8.0:
        color_palette=None:
    """
    # | - proc_data
    data_sets = data

    plot_data_list = []

    ind = 0
    for data_series in data_sets:
        ind += 1
        data_i_labels = add_to_label + data_series["label"]
        x_dat = data_series["data"]["frequency"].tolist()
        y_dat = data_series["data"]["activity"].tolist()

        x_dat = np.array(x_dat) + freq_shift

        # | - Creating Gaussian Curves From Frequencies and Activities
        x_dat, y_dat = to_plot(
            x_dat,
            y_dat,
            gam=gauss_gamma,
            scaling=scaling
            )
        # __|

        if color_palette is None:
            data_i = Scatter(
                x=x_dat,
                y=y_dat,
                fill="tozeroy",
                name=data_i_labels,
                mode="lines",
                )

        else:
            data_i = Scatter(
                x=x_dat,
                y=y_dat,
                fill="tozeroy",
                name=data_i_labels,
                mode="lines",
                marker={
                    "color": color_palette[ind - 1],
                    },
                )

        plot_data_list.append(data_i)

    return(plot_data_list)
    # __|
