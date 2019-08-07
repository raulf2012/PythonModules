#!/usr/bin/env python

"""Charge density class.

Author: Raul A. Flores
"""

#| - Import Modules
# import os
# import sys
import random
import pickle

import pandas as pd
import numpy as np
from ase.io.cube import read_cube
# , write_cube
# read_cube_data,

# import plotly as py
import plotly.graph_objs as go
#__|

class ChargeDensity(object):
    """docstring for ChargeDensity."""

    #| - ChargeDensity ********************************************************
    def __init__(self,
        cube_filename,
        master_data_filter=0.98,
        lower_bound_density_filter=0.025,

        wrap_atoms=True,
        working_dir=".",
        ):
        """Initialize ChargeDensity instance.

        Args:
            cube_filename:
            master_data_filter:
                Fraction of data to be removed randomly
            lower_bound_density_filter:
                Lower bound of normalized density value to be discarded
            working_dir:
        """
        #| - __init__

        #| - Define User Attributes
        self.cube_filename = cube_filename
        self.master_data_filter = master_data_filter
        self.lower_bound_density_filter = lower_bound_density_filter
        self.wrap_atoms = wrap_atoms
        self.working_dir = working_dir
        #__|

        (
            self.atoms,
            self.cd_data,
            self.origin,
            ) = self.__load_cube_file__()


        self.master_data_df = self.__process_data__()

        # self.__save_dataframe__()

        self.num_data_points = self.__number_of_data_points__()

        self.__filter_data__()
        self.__norm_electron_density__()
        self.__filter_low_density__()
        # self.__keep_only_edges__()
        #__|

    def __load_cube_file__(self):
        """Load charge density cube file."""
        #| - __load_cube_file__
        filename = self.cube_filename

        with open(filename, "r") as fle:
            data_master = read_cube(fle)

        atoms = data_master["atoms"]
        cd_data = data_master["data"]
        origin = data_master["origin"]

        if self.wrap_atoms:
            atoms.wrap(pbc=True)

        return(atoms, cd_data, origin)
        #__|

    def __process_data__(self):
        """Create dataframe from charge density data grid."""
        #| - __process_data__
        cd_data = self.cd_data
        atoms = self.atoms

        master_data_list = []
        for index, dens_i in np.ndenumerate(cd_data):
            data_i = {}
            data_i["x_ind"] = index[0]
            data_i["y_ind"] = index[1]
            data_i["z_ind"] = index[2]
            data_i["density"] = dens_i

            master_data_list.append(data_i)

        df = pd.DataFrame(master_data_list)

        df["x_coord_norm"] = df["x_ind"] / df["x_ind"].max()
        df["y_coord_norm"] = df["y_ind"] / df["y_ind"].max()
        df["z_coord_norm"] = df["z_ind"] / df["z_ind"].max()

        df["x_coord"] = df["x_coord_norm"] * atoms.cell[:, 0].sum()
        df["y_coord"] = df["y_coord_norm"] * atoms.cell[:, 1].sum()
        df["z_coord"] = df["z_coord_norm"] * atoms.cell[:, 2].sum()

        def multiply_by_unit_cell(row, atoms=None):
            coord = np.array([
                row["x_coord_norm"],
                row["y_coord_norm"],
                row["z_coord_norm"],
                ])

            scaled_cell_vectors = (atoms.cell.T * coord).T

            new_coord = np.array([
                scaled_cell_vectors[:, 0].sum(),
                scaled_cell_vectors[:, 1].sum(),
                scaled_cell_vectors[:, 2].sum(),
                ])

            # new_coord = atoms.cell.dot(coord)

            return(new_coord[0], new_coord[1], new_coord[2])

        (
            df["x_coord"],
            df["y_coord"],
            df["z_coord"],
            ) = zip(*df.apply(multiply_by_unit_cell, axis=1, atoms=atoms))

        return(df)
        #__|

    def __number_of_data_points__(self):
        """Return the number of individual density data points."""
        #| - __number_of_data_points__
        master_data_df = self.master_data_df

        num_data_points = len(master_data_df)

        return(num_data_points)
        #__|

    def __filter_data__(self):
        """Filter data randomly to decrease the data set size."""
        #| - __filter_data__
        master_data_df = self.master_data_df
        num_data_points = self.num_data_points
        master_data_filter = self.master_data_filter

        bool_list = []
        for data_i in range(num_data_points):
            rand_i = random.random()
            if rand_i > master_data_filter:
                bool_list.append(True)
            else:
                bool_list.append(False)

        df_filtered = master_data_df[bool_list]

        self.master_data_df = df_filtered
        #__|

    def __norm_electron_density__(self):
        """Normalize electron density from 0 to 1."""
        #| - __norm_electron_density__
        df = self.master_data_df

        max_density = df.density.max()
        df["norm_dens"] = df["density"] / max_density
        #__|

    def __filter_low_density__(self):
        """Filter low density entries from the data."""
        #| - __filter_low_density__
        df = self.master_data_df
        lower_bound_density_filter = self.lower_bound_density_filter

        df = df[df["norm_dens"] > lower_bound_density_filter]

        self.master_data_df = df
        #__|

    def __keep_only_edges__(self):
        """Only keep the outer surface points for clarity."""
        #| - __keep_only_edges__
        df = self.master_data_df

        df_a = df[df["x_ind"] == df.x_ind.max()]
        df_b = df[df["y_ind"] == df.y_ind.max()]
        df_c = df[df["z_ind"] == df.z_ind.max()]

        df_1 = df[df["x_ind"] == 0]
        df_2 = df[df["y_ind"] == 0]
        df_3 = df[df["z_ind"] == 0]

        df_surf = pd.concat(
            [
                df_a, df_b, df_c,
                df_1, df_2, df_3,
                ]
            )

        self.master_data_df = df_surf
        #__|

    def __save_dataframe__(self):
        """Save dataframe to pickle file.

        COMBAK impliment this to save time
        """
        #| - __save_dataframe__
        df = self.master_data_df
        working_dir = self.working_dir

        with open(working_dir + "/dataframe.pickle", "w") as fle:
            pickle.dump(df, fle)
        #__|

    def __load_dataframe__(self):
        """Load dataframe from pickle file.

        COMBAK and finish this
        """
        #| - __load_dataframe__
        tmp = 42

        #__|

    def create_charge_density_plotting_trace(self,
        opacity=0.4,
        size=4,
        ):
        """Create plotly trace from charge density distribution.

        Args:
            opacity:
            size:
                <float> or <str>
                Either a float to represent a constant size
                or
                "variable", to scale individual marker size with charge density

        """
        #| - create_charge_density_plotting_trace
        df = self.master_data_df

        if size == "variable":
            size = df["norm_dens"] * 13.

        trace1 = go.Scatter3d(

            x=df["x_coord"],
            y=df["y_coord"],
            z=df["z_coord"],

            # x=df["x_coord_norm"],
            # y=df["y_coord_norm"],
            # z=df["z_coord_norm"],

            # x=df["x_ind"],
            # y=df["y_ind"],
            # z=df["z_ind"],

            mode='markers',
            text=df["norm_dens"],
            opacity=opacity,
            marker=dict(
                size=size,
                color=df["norm_dens"],

                colorscale=[
                    [0., 'rgb(255, 200, 200, 0.1)'],
                    [1.0, 'rgb(255, 0, 0, 0.9)'],
                    ],
                )
            )

        return(trace1)
        #__|

    def create_unit_cell_plotting_trace(self,
        color="red",
        width=1.,
        ):
        """TEMP.

        Args:
            color
        """
        #| - create_unit_cell_plotting_trace
        atoms = self.atoms

        def unit_cell_leg_trace(
            point_0,
            point_1,
            color="red",
            width=1.,
            ):
            """Create trace for unit cell line.

            Args:
                point_0:
                point_1:
                color:
                width:
            """
            #| - unit_cell_leg_trace
            line_i = np.array([
                point_0,
                point_1,
                ])

            trace_i = go.Scatter3d(
                x=line_i[:, 0], y=line_i[:, 1], z=line_i[:, 2],
                mode="lines",
                line=dict(
                    color=color,
                    width=width,
                    )
                )

            return(trace_i)
            #__|

        line_trace_list = [

            #| - Origin to 3 adjacent corners
            unit_cell_leg_trace(
                np.array([0., 0., 0.]),
                atoms.cell[0],
                ),

            unit_cell_leg_trace(
                np.array([0., 0., 0.]),
                atoms.cell[1],
                ),

            unit_cell_leg_trace(
                np.array([0., 0., 0.]),
                atoms.cell[2],
                ),
            #__|

            #| - Farthest corner and 3 adjacent cornerse
            unit_cell_leg_trace(
                atoms.cell[0] + atoms.cell[1] + atoms.cell[2],
                atoms.cell[0] + atoms.cell[2]
                ),

            unit_cell_leg_trace(
                atoms.cell[0] + atoms.cell[1] + atoms.cell[2],
                atoms.cell[1] + atoms.cell[2]
                ),

            unit_cell_leg_trace(
                atoms.cell[0] + atoms.cell[1] + atoms.cell[2],
                atoms.cell[0] + atoms.cell[1]
                ),
            #__|

            #| - TEMP
            unit_cell_leg_trace(
                atoms.cell[1] + atoms.cell[2],
                atoms.cell[1],
                color=color,
                ),

            unit_cell_leg_trace(
                atoms.cell[0] + atoms.cell[2],
                atoms.cell[0],
                color=color,
                ),
            #__|

            #| - TEMP
            unit_cell_leg_trace(
                atoms.cell[2],
                atoms.cell[0] + atoms.cell[2],
                color=color,

                ),

            unit_cell_leg_trace(
                atoms.cell[2],
                atoms.cell[1] + atoms.cell[2],
                color=color,
                ),
            #__|

            #| - TEMP
            unit_cell_leg_trace(
                atoms.cell[0] + atoms.cell[1],
                atoms.cell[0],
                color=color,
                ),

            unit_cell_leg_trace(
                atoms.cell[0] + atoms.cell[1],
                atoms.cell[1],
                color=color,
                ),
            #__|

            ]

        return(line_trace_list)
        #__|

    def create_atoms_plotting_trace(self,
        size=12,
        ):
        """Create atom positions plotly trace.

        Args:
            size:
        """
        #| - create_atoms_plotting_trace
        atoms = self.atoms

        element_color_dict = {}
        element_color_dict["Fe"] = "orange"
        element_color_dict["C"] = "grey"
        element_color_dict["N"] = "blue"

        color_list = []
        for atom in atoms:
            elem = atom.symbol
            color_list.append(element_color_dict[elem])

        trace_i = go.Scatter3d(
            x=atoms.positions[:, 0],
            y=atoms.positions[:, 1],
            z=atoms.positions[:, 2],
            mode='markers',

            marker=dict(
                size=size,
                color=color_list,
                )
            )

        return(trace_i)
        #__|

    #__| **********************************************************************
