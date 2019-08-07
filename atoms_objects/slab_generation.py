#!/usr/bin/env python

"""Cut slabs from bulk structures.

Author: Raul A. Flores
"""

#| - Import Modules
from ase.build import surface
# from ase import io

from catkit.gen.surface import SlabGenerator as SlabGenerator_catkit
# from ase.visualize import view

from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.core.surface import SlabGenerator as SlabGenerator_pymatgen
#__|

def cut_slab_ase(
    bulk,
    facet,
    layers=6,
    vacuum=8,
    ):
    """Cut slab from bulk using ASE.

    Args:
        bulk:
        facet:
    """
    #| - cut_slab_ase
    # facet = [1, 1, 0]
    # layers = 6
    # vacuum = 8

    out_file = "out_slab_ase.traj"

    # facet = [int(i) for i in facet]

    # atoms = io.read("qn.traj")  # Bulk traj file (no vacuum)
    # atoms = io.read("init.cif")  # Bulk traj file (no vacuum)

    slab = surface(bulk, facet, layers, vacuum)
    slab.set_pbc(True)
    # slab.write(out_file)

    return(slab)
    #__|

def cut_slab_pymatgen(
    bulk,
    facet,
    min_slab_size=18.,
    min_vacuum_size=10.,
    ):
    """Cut slab from bulk using pymatgen.

    Args:
        bulk:
        facet:
    """
    #| - cut_slab_pymatgen
    # atoms = io.read("init.traj")
    # atoms = io.read("init.cif")

    structure = AseAtomsAdaptor.get_structure(bulk)

    # ZnO=Structure.from_file("ZnO.cif",primitive=False)

    slab_pymatgen = SlabGenerator_pymatgen(
        structure,
        # [1, 1, 0],
        facet,
        min_slab_size,
        min_vacuum_size,
        lll_reduce=True,
        center_slab=True,
        max_normal_search=2,
        ).get_slab()

    # initial_structure, miller_index, min_slab_size,
    # min_vacuum_size, lll_reduce=False, center_slab=False,
    # in_unit_planes=False, primitive=True, max_normal_search=None,
    # reorient_lattice=True

    slab = AseAtomsAdaptor.get_atoms(slab_pymatgen)

    slab.center()

    # slab.write("out_slab_pymatgen.traj")

    # IStructure.to(Al55, "poscar", filename="POSCAR")

    return(slab)
    #__|

def cut_slab_catkit(
    bulk,
    facet,
    slab_thickness=6,
    vacuum=8.,
    ):
    """Cut slab from bulk using CatKit.

    Args:
        bulk:
        facet:
    """
    #| - cut_slab_catkit
    gen = SlabGenerator_catkit(
        bulk,
        miller_index=facet,
        layers=slab_thickness,
        vacuum=vacuum,
        fixed=2,
        layer_type='ang',
        # attach_graph=True,
        standardize_bulk=False,
        tol=1e-8
        )

    terminations = gen.get_unique_terminations()

    images = []
    for i, t in enumerate(terminations):
        images += [gen.get_slab(iterm=i, size=1)]

    # view(images)

    return(images)
    #__|
