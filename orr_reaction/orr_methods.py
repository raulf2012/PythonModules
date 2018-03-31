"""ORR energetics classes and methods."""

#| - IMPORT MODULES
import numpy as np
import copy

from plotly.graph_objs import Scatter
#__|

class ORR_Free_E_Plot:
    """ORR FED Class.

    Development Notes:
        1. Should we consider the case where the bulk energy is not 0, and we
        have to normalize all of the species energies by it?

        2. H2O2 methods <-------------------------------------------------------
    """

    #| - ORR_Free_E_Plot *******************************************************
    # FIXME | I'm deprecating free_energy_dict in favor of a more generalizable
    # free_energy_df
    #
    # ksjfkjsdofjisjlfkjsdkfjsijf
    #
    #
    #



    def __init__(self,
        free_energy_dict,
        free_energy_df=None,
        system_properties=None,
        ):
        """

        """
        #| - __init__
        print("#()* - 180330 - New branch")

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
        # h2o2_e = 3.52

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

        return(free_energy_list)
        #__|

    def rxn_energy_lst(self):
        """List corresponding to the steps of ORR.

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
        """Apply bias to free energies.

        Applies a potential to every species in the 4 and 2-electron process
        and adjusts their free energies accordingly
        """
        #| - apply_bias
        mod_free_e_lst = []        # Free energy reaction path at applied bias
        for energy, elec in zip(energy_list, range(len(energy_list))[::-1]):
            mod_free_e_lst.append(energy - elec * bias)

        return(mod_free_e_lst)
        #__|

    def calc_overpotential(self):
        """
        # Returns the limiting overpotential for the given species and the
        limiting reaction step in the form of a list, species_A -> species_B is
        [species_A, species_B]
        """
        #| - calc_overpotential
        # e_lst = self.energy_lst
        rxn_spec = self.rxn_mech_states

        overpotential_lst = []
        for energy_i in enumerate(self.energy_lst[:-1]):
            energy_i_plus1 = self.energy_lst[energy_i[0] + 1]
            overpotential_i = 1.23 + energy_i_plus1 - energy_i[1]
            overpotential_lst.append(overpotential_i)

        overpotential = max(overpotential_lst)
        lim_step_index = overpotential_lst.index(overpotential)

        limiting_step = [rxn_spec[lim_step_index], rxn_spec[lim_step_index + 1]]

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
        tmp_list = range(len(energy_lst) * 2)
        energy_dupl_lst = [energy_lst[i // 2] for i in tmp_list]

        rxn_coord_steps = self.create_rxn_coord_array(
            len(energy_lst),
            spacing=spacing,
            step_size=step_size,
            )
        out_list = [rxn_coord_steps, energy_dupl_lst]

        return(out_list)
        #__|

    def create_plotly_series(self,
        energy_lst,
        name="TEMP",
        group="group1",
        hover_text=None,
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
            x=new_x_dat,
            y=new_y_dat,
            legendgroup=group,
            showlegend=True,
            name=name,
            hoverinfo="none",  # TEMP - 180317
            # text=hover_text,

            connectgaps=False,
            line=dict(
                color=color,
                width=6,
                ),
            mode="lines",
            )

        data_2 = Scatter(
            x=new_x_dat,
            y=new_y_dat,
            legendgroup=group,
            name=name,
            connectgaps=True,
            showlegend=show_leg_2,
            hoverinfo="none",
            # text=hover_text,

            line=dict(
                color=color,
                width=1,
                ),
            mode="lines",
            )



        #| - Creating x-data in middle of states
        # new_x_dat_tmp = copy.copy(new_x_dat)
        # new_y_dat_tmp = copy.copy(new_y_dat)

        short_y = np.array(y_dat)[::2]

        xdat = list(set(new_x_dat))
        xdat.sort()

        cnt = 0
        short_x = []
        for i_ind in range(len(xdat) / 2):
            short_x.append(xdat[cnt] + 0.5)  # TEMP Replace 0.5 with variable
            cnt += 2
        #__|

        data_3 = Scatter(
            x=short_x,
            y=short_y,

            legendgroup=group,
            name=name,
            # connectgaps = False,
            showlegend=False,
            # hoverinfo="y+name",
            text=hover_text,
            marker=dict(
                size=14,
                opacity=0.,
                ),
            mode="markers",
            )
        #__|

        if plot_mode == "all":
            data_lst = [data_1, data_2, data_3]
        elif plot_mode == "states_only":
            data_lst = [data_1, data_3]
        elif plot_mode == "full_lines":
            data_lst = [data_2, data_3]

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
                # Spacing inbetween the energy levels. The default of 0 creates
                a free energy diagram that looks like steps
        """
        #| - create_rxn_coord_array
        lst = []
        for i in range(1, rxn_steps):
            if i == 1:
                lst.append(step_size)
                lst.append(step_size + spacing)
            if i != 1:
                lst.append(lst[-1] + step_size)
                lst.append(lst[-2] + step_size + spacing)

        lst.insert(0, 0)
        lst.append(lst[-1] + step_size)

        return(lst)
        #__|

    #__| **********************************************************************


def calc_ads_e(
    df_row,
    bare_raw_e,
    oxy_ref_e=-443.70964,
    hyd_ref_e=-16.46018,
    ):
    """Calculate adsorption energies from raw DFT energetics.

    TEMP
    """
    #| - calc_ads_e
    row = df_row
    bare_slab = bare_raw_e
    oxy_ref = oxy_ref_e
    hyd_ref = hyd_ref_e

    # ads_e_list = []
    # for index, row in df.iterrows():
    try:
        num_O = row["atom_type_num_dict"][0]["O"]
    except:
        num_O = 0

    try:
        num_H = row["atom_type_num_dict"][0]["H"]
    except:
        num_H = 0

    try:
        raw_e = row["elec_energy"]
        ads_e_i = raw_e - bare_slab - num_O * oxy_ref - num_H * hyd_ref
    except:
        ads_e_i = None

    return(ads_e_i)

    # ads_e_list.append(ads_e_i)
    # ads_e_list = np.array(ads_e_list)
    # df["ads_e"] = ads_e_list
    #__|

def lowest_e_path(tmp=42):
    """Find the lowest energy pathway FED.

    From a set of FE pathways corresponding to different sites, the lowest
    energy states will be selected to construct a new FED.
    """
    #| - lowest_e_path
    print(tmp)
    #__|
