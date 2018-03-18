#| - IMPORT MODULES
import numpy as np
import matplotlib.patches as mpatches
import copy

import plotly
from plotly.graph_objs import Scatter, Layout, Figure
#__|

class ORR_Free_E_Plot:
    """
    Development Notes:
        1. Should we consider the case where the bulk energy is not 0, and we
        have to normalize all of the species energies by it?

        2. H2O2 methods <-------------------------------------------------------
    """
    #| - ORR_Free_E_Plot *******************************************************
    def __init__(self, free_energy_dict, system_properties=None):
        """
        """
        #| - __init__
        self.fe_dict = free_energy_dict
        self.sys_props = system_properties

        self.num_of_states = len(self.fe_dict) + 1  # bulk, OOH, O, OH, bulk

        self.rxn_mech_states = ["bulk", "ooh", "o", "oh", "bulk"]
        self.energy_lst = self.rxn_energy_lst()
        self.num_of_elec = range(self.num_of_states)[::-1]

        self.overpotential = self.calc_overpotential()[0]
        self.limiting_step = self.calc_overpotential()[1]

        self.ideal_energy = [4.92, 3.69, 2.46, 1.23, 0]

        self.energy_lst_h2o2 = self.rxn_energy_lst_h2o2()
        self.overpotential_h2o2 = self.calc_overpotential_h2o2()
        #__|

    def rxn_energy_lst_h2o2(self):
        """
        """
        #| - rxn_energy_lst_h2o2
        h2o2_e = 3.52

        free_energy_list = []
        for key in self.fe_dict:
            if key == "bulk" or key == key == "ooh":
                free_energy_list.append(self.fe_dict[key])
        # for state in self.rxn_mech_states:
            # free_energy_list.append(self.fe_dict[state])
        if len(free_energy_list) != 2:
            raise ValueError("Not the correct # of steps for H2O2")

        free_energy_list[0] += 4.92
        free_energy_list.append(3.52)
        # print(free_energy_list)
        return free_energy_list
        #__|

    def rxn_energy_lst(self):
        """
        Produces a list corresponding to the steps of ORR
        (1. O2, 2. *OOH, 3. *O, 4. *OH, 5. 2H2O)
        """
        #| - rxn_energy_lst
        free_energy_list = []
        for state in self.rxn_mech_states:
            free_energy_list.append(self.fe_dict[state])

        free_energy_list[0] += 4.92
        return free_energy_list
        #__|

    def apply_bias(self, bias, energy_list):
        """
        Applies a potential to every species in the 4 and 2-electron process
        and adjusts their free energies accordingly
        """
        #| - apply_bias
        mod_free_e_lst = []        # Free energy reaction path at applied bias
        for energy, elec in zip(energy_list, range(len(energy_list))[::-1]):
            mod_free_e_lst.append(energy - elec*bias)

        return(mod_free_e_lst)
        #__|

    def calc_overpotential(self):
        """
        # Returns the limiting overpotential for the given species and the limiting
        reaction step in the form of a list, species_A -> species_B is
        [species_A, species_B]
        """
        #| - calc_overpotential
        e_lst        = self.energy_lst
        rxn_spec    = self.rxn_mech_states

        overpotential_lst = []
        for energy_i in enumerate(self.energy_lst[:-1]):
            energy_i_plus1 = self.energy_lst[energy_i[0]+1]
            overpotential_i = 1.23 + energy_i_plus1 - energy_i[1]
            overpotential_lst.append(overpotential_i)

        overpotential = max(overpotential_lst)
        lim_step_index = overpotential_lst.index(overpotential)

        limiting_step = [rxn_spec[lim_step_index],rxn_spec[lim_step_index+1]]

        return [overpotential, limiting_step]
        #__|

    def calc_overpotential_h2o2(self):
        """
        """
        #| - calc_overpotential_h2o2
        return self.fe_dict["ooh"] - 4.22
        #__|

    def convert_to_plotting_list(self, energy_lst, spacing=0.5, step_size=1):
        """

        Args:
            energy_lst:
                # <type 'list'>
                # Modifies an energy list for plotting by repeating each entry
                # Ex. [4.92, 3.69, ... ] -> [4.92, 4.92, 3.69, 3.69, ... ]
        """
        #| - convert_to_plotting_list
        energy_dupl_lst = [energy_lst[i//2] for i in range(len(energy_lst)*2)]

        rxn_coord_steps = self.create_rxn_coord_array(len(energy_lst), spacing=spacing, step_size=step_size)
        out_list = [rxn_coord_steps, energy_dupl_lst]

        return(out_list)
        #__|

    def create_plotly_series(self,
        energy_lst,
        name="TEMP",
        group="group1",
        color="rgb(22, 96, 167)",
        plot_mode="all",
        ):
        """
        Args:
            energy_lst:
            name:
            group:
            color:
            plot_mode:
                "all"
                "states_only"
                "full_lines"
        """
        #| - create_plotly_series
        e_list = self.convert_to_plotting_list(energy_lst)

        # print(e_list)

        x_dat = e_list[0]
        y_dat = e_list[1]


        #| - Parameters
        if plot_mode == "all":
            show_leg_2 = False
        elif plot_mode == "states_only":
            show_leg_2 = False
        elif plot_mode == "full_lines":
            show_leg_2 = True
        #__|

        #| - Adding Breaks in Data
        new_x_dat = copy.copy(x_dat)
        new_y_dat = copy.copy(y_dat)

        cnt = 2
        for i_ind in range(len(x_dat) / 2 - 1):

            fill = new_x_dat[cnt - 1]

            new_x_dat.insert(cnt, fill)
            new_y_dat.insert(cnt, None)

            cnt += 3
        #__|

        #| - Plotly Scatter Plot
        data_1 = Scatter(
            x = new_x_dat,
            y = new_y_dat,
            legendgroup = group,
            showlegend=True,

            name = name,
            hoverinfo="none",  # TEMP - 180317
            connectgaps = False,
            line = dict(
                color = color,
                width = 6,
                ),
            mode = "lines",
            )

        data_2 = Scatter(
            x = new_x_dat,
            y = new_y_dat,
            legendgroup = group,
            name = name,
            connectgaps = True,
            showlegend=show_leg_2,
            hoverinfo="none",
            line = dict(
                color = color,
                width = 1,
                ),
            mode = "lines",
            )


        # print(new_x_dat)

        #| - Creating x-data in middle of states
        new_x_dat_tmp = copy.copy(new_x_dat)
        new_y_dat_tmp = copy.copy(new_y_dat)

        # print(new_y_dat_tmp)
        # print(y_dat)

        short_y = np.array(y_dat)[::2]
        # print(tmp)

        xdat = list(set(new_x_dat))
        xdat.sort()

        cnt = 0
        short_x = []
        for i_ind in range(len(xdat) / 2):
            short_x.append(xdat[cnt] + 0.5) # TEMP Replace 0.5 with variable
            cnt += 2

        #__|

        data_3 = Scatter(
            x = short_x,
            y = short_y,

            legendgroup = group,
            name = name,
            # connectgaps = False,
            showlegend=False,
            hoverinfo="y+name",

            marker= dict(
                size= 14,
                opacity= 0.,
                ),
            mode = "markers",
            )

        #__|

        if plot_mode == "all":
            data_lst = [data_1, data_2, data_3]
        elif plot_mode == "states_only":
            data_lst = [data_1, data_3]
        elif plot_mode == "full_lines":
            data_lst = [data_2, data_3]


        print(len(data_lst))
        return(data_lst)
        #__|

    def create_rxn_coord_array(self, rxn_steps, spacing=0, step_size=1):
        """
        Creates a reaction coordinate array ([0, 1, 1, 2, 2, 3]) for plotting

        Args:
            rxn_steps:
                # <type 'int'>
                # Number of steps in reaction coordinate including initial
                and final state.
                # Ex. A -> Intermediate -> C has 3 steps/states

            spacing:
                # <type 'float'>
                # Spacing inbetween the energy levels. The default of 0 creates a
                free energy diagram that looks like steps
        """
        #| - create_rxn_coord_array
        lst = []
        for i in range(1, rxn_steps):
            if i == 1:
                lst.append(step_size)
                # print(step_size)
                # print(spacing)
                lst.append(step_size+spacing)
            if i != 1:
                lst.append(lst[-1]+step_size)
                lst.append(lst[-2]+step_size+spacing)

        lst.insert(0,0)
        lst.append(lst[-1]+step_size)

        return(lst)
        #__|

    #__| **********************************************************************






















#| - __old__

# def process_orr_adsorbates_fed(free_e_dict):
#     """
#     Calculates the overpotential for the ORR reaction given the Free Energies
#     of the species. Also outputs the Free energy landscape at V=0 and V=V_eq, as
#     well as the x-coordinate data (just integer pairs defining the steps).
#
#     Output: {'overpotential': ###, 'ydata_nobias': [##,##,...], 'ydata_eq':
#     [##,##,...], 'xdata': [0,1,1,2,2,...]}
#
#     Args
#         free_e_lst: Dictionary containing the free energies of intermediates
#         format the ORR reaction. An example is shown below.
#         Ex.
#         {'Ni': {'bulk': 0, 'ooh': 4.87, 'o': 4.17, 'oh': 2.11}}
#     """
#     #| - process_orr_adsorbates_fed
#
#     #| - Appending Species Free Energies to Output Dict                     1111
#     metal = free_e_dict.keys()[0]
#     species_fe = free_e_dict[metal]
#     #__|                                                                    1111
#
#     fe_all = free_e_dict
#
#     steps = range(5)
#
#     #| - Constructing Reaction Coordinate - AUTOMATIC                       1111
#     i = 0
#     species_lst = ['bulk','ooh','o','oh','bulk']
#     free_e_lst = []
#     for species in species_lst:
#         free_e_lst.append(fe_all[metal][species])
#     # print(free_e_lst)
#     #__|                                                                    1111
#
#     #| - Normalized energies by bulk energy                                 1111
#     bulk_free_e =  free_e_lst[0]
#     for i in range(len(free_e_lst)):
#         free_e_lst[i] = free_e_lst[i] - bulk_free_e
#     #__|                                                                    1111
#
#     #| - Replaces String "M" with Metal Name                                1111
#     for j in range(len(species_lst)):
#         species_lst[j] = str(species_lst[j]).replace('M', metal)
#     #__|                                                                    1111
#
#     #| - Constructing Free Energy Diagram Data                              1111
#     xdata = []
#     ydata = []
#     for step in range(len(steps)):
#         xdata.append(steps[step])
#         xdata.append(steps[step]+1)
#
#         if steps[step] == 0:
#             ydata.append(free_e_lst[step] + 4.92)
#             ydata.append(free_e_lst[step] + 4.92)
#         else:
#             ydata.append(free_e_lst[step])
#             ydata.append(free_e_lst[step])
#     #__|                                                                    1111
#
#     #| - Finding Overpotential                                              1111
#     overpotential = 0
#     for i in range(len(ydata)-1):
#         if ydata[i+1] == ydata[i]:
#             None
#         else:
#             d = 1.23 + ydata[i+1] - ydata[i]
#             if d > overpotential:
#                 overpotential = d
#     #__|                                                                    1111
#
#     #| - Constructing FED at Equilibrium and Overpotential                  1111
#     ydatalim    = []
#     ydataeq        = []
#     for j in range(len(free_e_lst)):
#         if j == 0:
#
#             ydataeq.append(free_e_lst[j] - (4-j)*1.23 + 4*1.23)
#             ydataeq.append(free_e_lst[j] - (4-j)*1.23 + 4*1.23)
#
#             ydatalim.append(free_e_lst[j] - (4-j)*(1.23-overpotential) + 4*1.23)
#             ydatalim.append(free_e_lst[j] - (4-j)*(1.23-overpotential) + 4*1.23)
#
#         else:
#             ydataeq.append(free_e_lst[j] - (4-j)*1.23)
#             ydataeq.append(free_e_lst[j] - (4-j)*1.23)
#
#             ydatalim.append(free_e_lst[j] - (4-j)*(1.23-overpotential))
#             ydatalim.append(free_e_lst[j] - (4-j)*(1.23-overpotential))
#
#     #__|                                                                    1111
#
#     #| - Constructing FED for H2O2 (2-e process)                            1111
#     ydata_h2o2    = []
#     for j in range(3):
#         if j == 0:
#             ydata_h2o2.append(free_e_lst[j] - (4-j)*(0.7/2.) + 4*1.23)
#             ydata_h2o2.append(free_e_lst[j] - (4-j)*(0.7/2.) + 4*1.23)
#
#         elif j == 2:
#             ydata_h2o2.append(3.52)
#             ydata_h2o2.append(3.52)
#
#         else:
#             ydata_h2o2.append(free_e_lst[j] - (4-2.*j)*(0.7/2.))
#             ydata_h2o2.append(free_e_lst[j] - (4-2.*j)*(0.7/2.))
#
#     xdata_h2o2 = [ 0 , 1 , 1 , 2 , 2 , 3 ]
#     #__|                                                                    1111
#
#     return {"ydata_nobias":ydata, "ydata_eq":ydataeq, "ydata_h2o2":ydata_h2o2,
#      "xdata":xdata,    "overpotential":overpotential, "species_fe": species_fe,
#     "xdata_h2o2":xdata_h2o2}
#
#     #__|
#
#
# def create_rxn_coord_array(rxn_steps, spacing=0, step_size=1):
#     """
#     Creates a reaction coordinate array ([0, 1, 1, 2, 2, 3]) for plotting
#
#     Args:
#         rxn_steps:
#             # <type 'int'>
#             # Number of steps in reaction coordinate including initial
#             and final state.
#             # Ex. A -> Intermediate -> C has 3 steps/states
#
#         spacing:
#             # <type 'float'>
#             # Spacing inbetween the energy levels. The default of 0 creates a
#             free energy diagram that looks like steps
#     """
#     #| - create_rxn_coord_array
#     lst = []
#     for i in range(1, rxn_steps):
#         if i == 1:
#             lst.append(step_size)
#             lst.append(step_size+spacing)
#         if i != 1:
#             lst.append(lst[-1]+step_size)
#             lst.append(lst[-2]+step_size+spacing)
#
#     lst.insert(0,0)
#     lst.append(lst[-1]+step_size)
#     return lst
#     #__|
#
#
# def plot_fed_line(plot_obj, *args, **kwargs):
#     """
#
#     Args:
#         plot_obj:
#             # <>
#             #
#         xdata:
#             # <type 'list'>
#             #
#         ydata:
#             # <type 'list'>
#             #
#     """
#     #| - plot_fed_line
#     xdata = args[0]; ydata = args[1]
#
#     plot_obj.plot(xdata[2*0:2*0+2], ydata[2*0:2*0+2], linewidth=3,**kwargs)
#
#     if "label" in kwargs: del kwargs["label"]
#
#     kwargs_even        = copy.copy(kwargs)
#     kwargs_odd        = copy.copy(kwargs)
#
#     if "path_effects" in kwargs_odd: del kwargs_odd["path_effects"]
#
#     for i in range(len(xdata)/2):
#         tmp = plot_obj.plot(xdata[2*i:2*i+2], ydata[2*i:2*i+2], linewidth=3, solid_capstyle='round', **kwargs_even)
#         plot_obj.plot(xdata[2*i+1:2*i+3], ydata[2*i+1:2*i+3], linewidth=0.5,**kwargs_odd)
#         # print(dir(tmp))
#         # tmp.set_solid_capstyle('round')
#         # ax.margins(.2)
#     #__|

#__|
