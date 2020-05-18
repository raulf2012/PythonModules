#!/usr/bin/env python

"""Methods to add atoms/adsorbates to surfaces."""

# | - Import Modules
import numpy as np
import copy
import math

from operator import itemgetter
from ase.build import add_adsorbate
from misc_modules.numpy_methods import angle_between
# __|

def add_adsorbate_centered(active_element, slab, adsorbate, ads_height=2.5):
    """Add adsorbate to surface.

    Places an adsorbate above an active element/atom of a slab. If multiple
    active elements are present within a unit cell, the adsorbate will be placed
    above the atom closest to the center.

    Args:
        active_element: String
        slab: ASE atoms object
        adsorbate: ASE atoms object
        ads_height: Float
    """
    # | - add_adsorbate_centered
    center = (sum(slab.cell)) / 2
    act_metals = []
    for atom in slab:
        if atom.symbol == active_element:
            act_metals.append([
                atom.index,
                np.linalg.norm(atom.position - center),
                ])

    active_metal = slab[sorted(act_metals, key=itemgetter(1))[0][0]]
    ads_pos = (active_metal.position[0], active_metal.position[1])
    add_adsorbate(slab, adsorbate, ads_height, position=ads_pos)

    return slab
    # __|


def add_graphene_layer(
    slab,
    graphene_units=1,
    graph_surf_d=3.0,
    graph_bond_d_real=1.4237,
    ):
    """
    Add graphene layer above surface of hexagonal unit cell slab.

    Args:
        slab:
        graphene_units:
        graph_surf_d:
        graph_bond_d_real:
    """
    # | - add_graphene_layer
    slab = copy.deepcopy(slab)

    # num_graph_units = graphene_units + 1
    graphene_units = int(graphene_units)

    num_graph_units = int(graphene_units)
    ngu = num_graph_units

    num_graph_bond_lengths = (2. + 1.) * ngu

    v1 = slab.cell[0]
    v2 = slab.cell[1]

    angle = angle_between(v1, v2)
    angle = math.degrees(angle)
    # print("Angle between lattice vectors: " + str(angle))

    mag1 = np.linalg.norm(v1)
    mag2 = np.linalg.norm(v2)

    assert round(mag1) == round(mag2)

    graph_bond_d = mag1 / num_graph_bond_lengths

    xy_cell = slab.cell[0:2, 0:2]  # x and y components of unit cell
    x_unit_v = xy_cell[0]
    y_unit_v = xy_cell[1]

    x_unit_v = x_unit_v / np.linalg.norm(x_unit_v)
    y_unit_v = y_unit_v / np.linalg.norm(y_unit_v)


    # | - STRAIN
    tmp = mag1 / num_graph_bond_lengths
    strain = 100. * (graph_bond_d_real - tmp) / graph_bond_d_real
    print("Strain: " + str(strain))
    # __|

    patt_cnt_x = 0
    patt_cnt_y = 0
    C_pos_lst = []
    for y_ind in range(graphene_units * 3):
        patt_cnt_x = patt_cnt_y
        for x_ind in range(graphene_units * 3):

            if patt_cnt_x == 0 or patt_cnt_x == 1:

                pos_x = x_ind * graph_bond_d * x_unit_v
                pos_y = y_ind * graph_bond_d * y_unit_v

                pos_i = np.array(pos_x) + np.array(pos_y)
                pos_i = np.append(pos_i, 0.)

                C_pos_lst.append(pos_i)
                # atom_cent = (origin[0] + x_ind * graph_bond_d
                # - graph_bond_d * y_ind * y_unit_v[0], origin[1] + y_ind
                # * graph_bond_d)

                add_adsorbate(
                    slab,
                    "C",
                    graph_surf_d,
                    position=(pos_i[0], pos_i[1]),
                    )

                patt_cnt_x += 1

            elif patt_cnt_x == 2:
                patt_cnt_x = 0

        if patt_cnt_y == 0:
            patt_cnt_y = 2
            continue

        if patt_cnt_y == 2:
            patt_cnt_y = 1
            continue

        elif patt_cnt_y == 1:
            patt_cnt_y = 0
            continue

    C_pos_lst = np.array(C_pos_lst)

    return(slab)
    # __|
