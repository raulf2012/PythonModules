#!/usr/bin/env python

"""Module to analyze and plot simulated XRD spectra from VESTA

# 31.0:       | (111) | Heighest peak
# 25.9:       | (110) | 2nd highest peak
# 35.9:       | (021) | 3rd highest peak
# 56.2 & 56.3 | (221) | 4th highest peak
# 33.75       | (002) | 5th highest peak
# 43.03       | (112) | 6th
# 54.7 & 54.8 | (202) | 7th
# 58.65       | (113) | 8th
# 41.7 & 41.8 | (200) | 9th

Author: Raul A. Flores
"""

# | - IMPORT MODULES
import os
import math

import pandas as pd
import numpy as np
import pickle

from sympy import Symbol
from sympy import *

from pathlib import Path

from tqdm import tqdm
import time


import chart_studio.plotly as py
import plotly.graph_objs as go
from plotly import io as pyio

from scipy.signal import find_peaks
# __|


class XRD_Spectra:
    """TEMP.
    """

    # | - XRD_Spectra ******************************************************
    _SAVE_FOLDER_NAME = ".xrd_save_dir"
    _SAVE_FILE_NAME = "xrd_save_state.pickle"
    _SAVE_ATTRIBUTES = ["df_peaks", "temp_0"]


    def __init__(self,
        reflections_table_path=None,
        theta_range=[0, 150],
        theta_spacing=0.1,
        load_from_saved_state=True,
        peak_broadening=0.04,
        ):
        """

        Args:
          reflections_table_path: str
            path + file_name of reflections table text file
          theta_range: list, len=2
            x-axis range of 2theta (should be 2theta not theta)
          theta_spacing: float
            Discritized spacing of x-axis (2theta)
          load_from_saved_state: Bool
            Whether to attempt to load attributes from file to save time
          peak_broadening: float
            Dictates the broadness of peaks (smaller is thinner)
        """
        # | - __init__
        self.temp_0 = "TEMP TEMP TEMP TEMP"

        # | - Setting Argument Instance Attributes
        self.reflections_table_path = reflections_table_path
        self.theta_range = theta_range
        self.theta_spacing = theta_spacing
        self.load_from_saved_state = load_from_saved_state
        self.peak_broadening = peak_broadening
        # __|

        # | - Initializing Internal Instance Attributes
        self.df_peaks = None
        # __|

        # Try to load previous saved state
        self.loaded_state_data = self.__load_state__()

        self.theta_array = self.__create_theta_array__()

        self.df = self.__read_reflections_table__()

        self.__process_reflections_table__()

        self.__create_lorentz_for_all_signals__()

        self.summed_lorentz = self.__create_summed_lorentz_function__()

        self.spectrum = self.__create_total_spectrum__()
        # __|


    def __read_reflections_table__(self):
        """
        """
        # | - __read_reflections_table__
        file_path_i = self.reflections_table_path
        # file_path_i = "/home/raulf2012/Dropbox/01_norskov/00_projects/columbite_iro2_oer/workflow/01_latt_const_opt/an_xrd_pattern/vesta_xrd_gen/optimized_bulk.txt"

        with open(file_path_i) as f:
            content = f.readlines()
        content = [x.strip() for x in content]

        list_of_rows = []
        for i_cnt, line_i in enumerate(content):
            if i_cnt == 0:
                header_row = line_i
                continue

            list_of_vals = [i for i in line_i.split(" ") if i is not ""]
            row_i = {
                "h": int(list_of_vals[0]),
                "k": int(list_of_vals[1]),
                "l": int(list_of_vals[2]),
                "d": float(list_of_vals[3]),
                "F_real": float(list_of_vals[4]),
                "F_imag": float(list_of_vals[5]),
                "F_abs": float(list_of_vals[6]),
                "2theta": float(list_of_vals[7]),
                "I": float(list_of_vals[8]),
                "M": float(list_of_vals[9]),
                "theta": float(list_of_vals[10]),
                "Phase": float(list_of_vals[11])}
            list_of_rows.append(row_i)

        df = pd.DataFrame(list_of_rows)

        return(df)

        # column_headers = [i for i in header_row.split(" ") if i is not ""]
        # df.columns = column_headers
        # df = pd.read_csv("../vesta_xrd_gen/reflections_table.csv")
        # df = df.sort_values("I", ascending=False)
        # __|


    def __process_reflections_table__(self):
        """
        """
        # | - __process_reflections_table__
        self.__create_simplified_facet_string__()

        # __|


    def __create_simplified_facet_string__(self):
        """
        """
        # | - __create_simplified_facet_string__
        df = self.df

        # Create string facet representation
        df["facet_mine"] = abs(df["h"]).astype("str") + abs(df["k"]).astype("str") + abs(df["l"]).astype("str")
        # __|


    def __create_lorentz_for_all_signals__(self):
        """Create lorentz function for all signals.

        Adds column to df
        """
        # | - __create_lorentz_for_all_signals__
        df = self.df
        peak_broadening = self.peak_broadening


        # Sympy addition of all functions
        x = Symbol("x")
        x0 = Symbol("x0")
        gamma = Symbol("gamma")
        peak_height = Symbol("peak_height")

        funct_list = []
        for i_ind, row_i in df.iterrows():
            gamma_i = peak_broadening
            intensity_i = row_i["I"]
            x0_i = row_i["2theta"]
            funct_i = Lorentz_i(x0, x, gamma, peak_height).subs({
                x0: x0_i,
                gamma: gamma_i,
                peak_height: intensity_i})
            funct_list.append(funct_i)

        df["function"] = funct_list
        # __|


    def __create_summed_lorentz_function__(self):
        """
        """
        # | - __create_summed_lorentz_function__
        df = self.df

        sum_funct = np.sum(df["function"].tolist())

        return(sum_funct)
        # __|


    def __create_theta_array__(self):
        """
        """
        # | - __create_theta_array__
        # Compute Signal
        theta_range = self.theta_range
        theta_spacing = self.theta_spacing

        min_theta = theta_range[0]
        max_theta = theta_range[1]
        # # x_range = np.arange(1, 150, 1.) #interval of points
        # min_theta = 1
        # max_theta = 80

        x_range = np.arange(min_theta, max_theta, theta_spacing)

        return(x_range)
        # __|


    def __create_total_spectrum__(self):
        """
        """
        # | - __create_total_spectrum__
        summed_lorentz = self.summed_lorentz
        # theta_range = self.theta_range
        theta_array = self.theta_array
        x = Symbol("x")

        a = lambdify((x, ), summed_lorentz)
        spectrum = a(np.array(theta_array))

        # Normalizing signal to 100
        spectrum = 100 * (spectrum / spectrum.max())

        return(spectrum)
        # __|


    def compute_peak_positions(self):
        """
        """
        # | - compute_peak_positions

        # | - Class Attributes
        spectrum = self.spectrum
        df = self.df
        theta_array = self.theta_array
        theta_array = self.theta_array
        theta_range = self.theta_range

        loaded_state_data = self.loaded_state_data
        # __|


        bool_0 = (loaded_state_data is not None)
        if bool_0:
            bool_1 = ("df_peaks" in list(loaded_state_data.keys()))
            if bool_1:
                print("Loaded 'df_peaks' from save state")
                df_peaks = loaded_state_data["df_peaks"]

        else:
            print("Computing 'df_peaks' from scratch")

            time.sleep(3)  # TEMP

            min_theta = theta_range[0]
            max_theta = theta_range[1]

            x = Symbol("x")

            # | -  Computing Peak Positions
            # Used for peak_finder to set width
            grid_to_theta = len(theta_array) / (max_theta - min_theta)


            peaks = find_peaks(
                spectrum,
                height=4,
            #     threshold=1,
                distance=grid_to_theta * 1.)

            peaks_x = [theta_array[i] for i in peaks[0]]
            # __|


            # | - Computing the Main Facets for the Prominant Peaks
            main_facets_list = []
            peaks_x_tqdm = tqdm(peaks_x)
            for peak_x_i in peaks_x_tqdm:
                peaks_x_tqdm.set_description("Processing %s" % peak_x_i)

                intensity_list = []
                for i_ind, row_i in df.iterrows():
                    funct_i = row_i["function"]

                    a = lambdify((x, ), funct_i)
                    signal_i = a(np.array([peak_x_i]))

                    if type(signal_i) is int:
                        intensity_list.append(float(signal_i))
                    elif type(signal_i).__module__ == np.__name__:
                        assert len(signal_i) == 1, "ISDFJIDSJfi"
                        intensity_list.append(signal_i[0])

                df[peak_x_i] = intensity_list

                highest_contributions = df.sort_values(
                    peak_x_i, ascending=False)[peak_x_i].unique()[0:1]

                major_facets = df[df[peak_x_i].isin(
                    highest_contributions)].sort_values(
                        peak_x_i, ascending=False)["facet_mine"].unique()

                main_facets_i = "_".join(major_facets)

                main_facets_list.append(main_facets_i)
            # __|

            df_peaks = pd.DataFrame()
            df_peaks["main_facets"] = main_facets_list
            df_peaks["2theta"] = peaks_x
            df_peaks["intensity"] = peaks[1]["peak_heights"]


        self.df_peaks = df_peaks


        # df.sort_values(peak_x_i, ascending=False)[peak_x_i].unique()[0:4]
        #
        # # df[df[peak_x_i].isin(highest_contributions)].sort_values(
        # #     peak_x_i, ascending=False)["facet_mine"]
        #
        # df[df[peak_x_i].isin(highest_contributions)].sort_values(
        #     peak_x_i, ascending=False)
    # __|


    def __create_save_dir__(self):
        """
        """
        # | - __create_save_dir__
        save_folder_name = self._SAVE_FOLDER_NAME
        if not os.path.exists(save_folder_name):
            os.makedirs(save_folder_name)
        # __|


    def __save_state__(self):
        """
        """
        # | - __save_states__
        save_folder_name = self._SAVE_FOLDER_NAME
        save_file_name = self._SAVE_FILE_NAME

        self.__create_save_dir__()

        save_attributes = dict()
        for _SAVE_ATTRIBUTE_i in self._SAVE_ATTRIBUTES:
            attr_i = getattr(self, _SAVE_ATTRIBUTE_i)

            save_attributes[_SAVE_ATTRIBUTE_i] = attr_i

        with open(os.path.join(save_folder_name, save_file_name), "wb") as fle:
            pickle.dump(save_attributes, fle)
        # __|


    def __load_state__(self):
        """
        """
        # | - __load_state__
        save_folder_name = self._SAVE_FOLDER_NAME
        save_file_name = self._SAVE_FILE_NAME
        load_from_saved_state = self.load_from_saved_state

        save_data = None
        if load_from_saved_state:
            path_i = os.path.join(save_folder_name, save_file_name)
            file_i = Path(path_i)
            if file_i.is_file():
                print("SUCCESS: Loaded save state data")
                with open(path_i, "rb") as fle:
                    save_data = pickle.load(fle)
        else:
            pass

        return(save_data)
        # __|

    # __| **********************************************************************


class XRD_Plot():
    """TEMP.
    """

    # | - XRD_Plot *************************************************************
    def __init__(self,
        XRD_Spectra,


        ):
        """
        """
        # | - __init__

        # | - Setting Instance Attributes
        self.XRD_Spectra = XRD_Spectra
        # __|

    def create_figure(self):
        """
        """
        # | - create_figure
        data = []

        trace_i = self.__trace_peak_labels__()
        data.append(trace_i)

        trace_i = self.__trace_spectra__()
        data.append(trace_i)

        return(data)
        # __|

    def __trace_peak_labels__(self):
        """
        """
        # | - __trace_peak_labels__
        df_peaks = self.XRD_Spectra.df_peaks

        trace_i = go.Scatter(
            x=df_peaks["2theta"],
            y=df_peaks["intensity"] + 5,
            text=df_peaks["main_facets"],
        #     mode="markers+text",
            mode="text",

        #     text=df_peaks["main_facets"],

            hoverinfo=None,
            hoverinfosrc=None,
            hoverlabel=None,
            hoveron=None,
            hovertemplate=None,
            hovertemplatesrc=None,
            hovertext=df_peaks["main_facets"],
            hovertextsrc=None,
            )

        return(trace_i)

        # data.append(trace_1)
        # __|


    def __trace_spectra__(self):
        """
        """
        # | - __trace_spectra__
        theta_array = self.XRD_Spectra.theta_array
        spectrum = self.XRD_Spectra.spectrum

        trace_i = go.Scatter(
            x=theta_array,
            y=spectrum,
            mode="lines",
            hoverinfo="skip",
            )

        return(trace_i)
        # __|

    def get_layout(self):
        """
        """
        # | - get_layout
        layout = go.Layout()

        layout.width = 18.4 * 37.795275591
        layout.height = 10 * 37.795275591

        # layout.width = 10 * 37.795275591
        # layout.height = 6 * 37.795275591

        layout.margin = go.layout.Margin(
            b=60,
            l=80,
            r=10,
            t=10,
            )

        layout.paper_bgcolor='rgba(240,240,240,0.8)'
        layout.plot_bgcolor='rgba(250,250,250,0.9)'
        # layout.plot_bgcolor='rgba(100,100,100,0.9)'

        layout.font = dict(
            family='Arial',
        #     size=18,
            color='black',
            )

        shared_dict = dict(
            title_font=dict(
                size=10.5 * (4/3),
        #         family='Courier',
        #         color='crimson',
                ),
            tickfont=dict(size=9 * (4/3)),
            )

        layout.xaxis = go.layout.XAxis(
            mirror=True,
            range=[20, 70],
        #     tickmode=,
            automargin=True,
            showgrid=False,


            ticks="outside",
            title="2θ",

            showticklabels=True,

            zeroline=False,
            zerolinecolor="black",

            showline=True,
            linecolor="black",

            **shared_dict,
            )

        layout.yaxis = go.layout.YAxis(
            automargin=True,
            mirror=True,
            range=[0, 115],
            showgrid=False,
            showline=True,
            showticklabels=True,
        #     tickmode=,
            linecolor="black",
            zerolinecolor="black",
            ticks="outside",
            title="Intensity (Arbitrary Units)",
            zeroline=False,

            **shared_dict,
            )

        return(layout)
        # __|

        # __|


    # __| **********************************************************************


# | - METHODS

def get_weighted_xrange(x_bounds, gamma, x0, min_step_size=1):
    """
    """
    # | - get_weighted_xrange
    x_i = x_bounds[0]
    x_range= [x_i]
    while x_i < x_bounds[-1]:
        pi = math.pi
        term_1 = (gamma ** 2) * ((x_i - x0) ** 2 + gamma ** 2) ** (-1)
        step_i = (0.5) * (term_1) ** -1

        if step_i > min_step_size:
            step_i = min_step_size
        x_new = x_i + step_i

        x_i = x_new
        x_range.append(x_i)

    return(x_range)
    # __|

def Lorentz_i(x0, x, gamma, peak_height):
    """
    """
    # | - Lorentz_i
    pi = math.pi
    term_0 = (1 / (pi * gamma))
    term_1 = (gamma ** 2) * ((x - x0) ** 2 + gamma ** 2) ** (-1)
#     y_i = term_0 * term_1
    y_i = peak_height * term_1

    return(y_i)
    # __|

def Lorentz_distr_i(x0, x_bounds, gamma, intensity):
    """
    """
    # | - Lorentz_distr_i
    x_range = get_weighted_xrange(x_bounds, gamma, x0, min_step_size=2)


    lorentz_distr = []
    for x_i in x_range:
        y_i = Lorentz_i(x0, x_i, gamma, intensity)
        lorentz_distr.append(y_i)

    return(x_range, lorentz_distr)
    # __|

# __|
