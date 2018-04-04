"""ORR energetics classes and methods."""

#| - IMPORT MODULES
import copy
import numpy as np
import pandas as pd
from plotly.graph_objs import Scatter
#__|

class ORR_Free_E_Plot:
    """ORR FED Class.

    Development Notes:
        # TODO Should we consider the case where the bulk energy is not 0, and
        we have to normalize all of the species energies by it?
    """

    #| - ORR_Free_E_Plot *******************************************************

    def __init__(self,
        free_energy_df=None,
        system_properties=None,
        state_title="adsorbate",
        free_e_title="ads_e"
        ):
        """
        Input variables to class instance.

        Args:
            free_energy_df:
                Pandas dataframe containing the adsorbates as rows
                Required columns, adsorbate, free energy
            system_properties:
            state_title:
            free_e_title:
        """
        #| - __init__
        self.fe_df = free_energy_df
        self.sys_props = system_properties
        self.state_title = state_title
        self.fe_title = free_e_title

        self.rxn_mech_states = ["bulk", "ooh", "o", "oh", "bulk"]
        self.ideal_energy = [4.92, 3.69, 2.46, 1.23, 0]

        self.add_bulk_entry()

        self.fill_missing_data()

        self.num_of_states = len(self.fe_df) + 1  # bulk, OOH, O, OH, bulk
        self.energy_lst = self.rxn_energy_lst()
        self.num_of_elec = range(self.num_of_states)[::-1]
        self.overpotential = self.calc_overpotential()[0]
        self.limiting_step = self.calc_overpotential()[1]
        # self.ideal_energy = [4.92, 3.69, 2.46, 1.23, 0]
        self.energy_lst_h2o2 = self.rxn_energy_lst_h2o2()
        self.overpotential_h2o2 = self.calc_overpotential_h2o2()
        #__|

    def add_bulk_entry(self,
        bulk_e=0.0,
        ):
        """
        Append a row entry to data frame corresponding to bulk state.

        Args:
            bulk_e:
        """
        #| - add_bulk_entry
        df = self.fe_df
        bulk_df = pd.DataFrame([{
            "adsorbate": "bulk",
            "ads_e": bulk_e,
            }])

        df = df.append(bulk_df, ignore_index=True)

        self.fe_df = df
        #__|

    def rxn_energy_lst_h2o2(self):
        """Construct energy list of h2o2 FED."""
        #| - rxn_energy_lst_h2o2
        # h2o2_e = 3.52

        df = self.fe_df

        free_energy_list = []
        for index, row in df.iterrows():
            if row["adsorbate"] == "bulk" or row["adsorbate"] == "ooh":
                free_energy_list.append(row["ads_e"])

        # print(free_energy_list)
        # Checking length of energy list
        if len(free_energy_list) != 2:
            raise ValueError("Not the correct # of steps for H2O2")

        free_energy_list[0] += 4.92
        free_energy_list.append(3.52)

        return(free_energy_list)
        #__|

    def property_list(self, column_name):
        """
        General method to create a list from a column in the dataframe.

        The length of the list will correspond to the steps in the ORR
        mechanism.

        Args:
            column_name:
        """
        #| - property_list
        df = self.fe_df

        property_list = []
        for state in self.rxn_mech_states:
            tmp = df.loc[df[self.state_title] == state]
            tmp1 = tmp.iloc[0][column_name]
            property_list.append(tmp1)

        # free_energy_list[0] += 4.92

        return(property_list)
        #__|

    def fill_missing_data(self):
        """
        """
        #| - fill_missing_data
        df = self.fe_df
        # print(len(self.fe_df))
        df_missing_data = pd.DataFrame()
        for state in self.rxn_mech_states:
            df_state = df.loc[df[self.state_title] == state]

            #| - If df is missing state fill in row with NaN for energy
            if df_state.empty:
                df_state = pd.DataFrame([{
                    self.state_title: state,
                    self.fe_title: np.nan,
                    }])
                df_missing_data = df_missing_data.append(df_state)
            #__|

        # print(df_missing_data)
        self.fe_df = self.fe_df.append(df_missing_data)
        # print(len(self.fe_df))
        #__|

    def rxn_energy_lst(self):
        """List corresponding to the steps of ORR.

        (1. O2, 2. *OOH, 3. *O, 4. *OH, 5. 2H2O)
        """
        #| - rxn_energy_lst
        df = self.fe_df
        free_energy_list = []
        for state in self.rxn_mech_states:

            df_state = df.loc[df[self.state_title] == state]

            #| - If df is missing state fill in row with NaN for energy
            if df_state.empty:
                df_state = pd.DataFrame([{
                    self.state_title: state,
                    self.fe_title: np.nan,
                    }])
            #__|

            tmp1 = df_state.iloc[0][self.fe_title]
            free_energy_list.append(tmp1)

        free_energy_list[0] += 4.92

        return(free_energy_list)
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
        Calculate overpotential for 4e- process.

        Returns the limiting overpotential for the given species and the
        limiting reaction step in the form of a list, species_A -> species_B is
        [species_A, species_B]
        """
        #| - calc_overpotential
        rxn_spec = self.rxn_mech_states

        overpotential_lst = []
        for energy_i in enumerate(self.energy_lst[:-1]):
            energy_i_plus1 = self.energy_lst[energy_i[0] + 1]
            overpotential_i = 1.23 + energy_i_plus1 - energy_i[1]
            overpotential_lst.append(overpotential_i)

        overpotential = max(overpotential_lst)
        lim_step_index = overpotential_lst.index(overpotential)

        limiting_step = [rxn_spec[lim_step_index], rxn_spec[lim_step_index + 1]]
        out_list = [overpotential, limiting_step]

        return(out_list)
        #__|

    def calc_overpotential_h2o2(self):
        """
        Calculate overpotential for 2e- process.

        The overpotential for the 2e- process depends only on the energy of the
        *OOH intermediate
        """
        #| - calc_overpotential_h2o2
        df = self.fe_df
        ooh_row = df[df["adsorbate"] == "ooh"]
        ooh_ads_e = ooh_row.iloc[0]["ads_e"]

        op_4e = ooh_ads_e - 4.22

        return(op_4e)
        #__|

    def create_rxn_coord_array(self,
        rxn_steps,
        spacing=0,
        step_size=1,
        ):
        """
        Create a reaction coordinate array ([0, 1, 1, 2, 2, 3]) for plotting.

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

    #| - Plotting @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    def convert_to_plotting_list(self,
        energy_lst,
        spacing=0.5,
        step_size=1,
        ):
        """
        Repeat entries in energy list to conform to FED plot.

        Modifies an energy list for plotting by repeating each entry
        Ex. [4.92, 3.69, ... ] -> [4.92, 4.92, 3.69, 3.69, ... ]

        Args:
            energy_lst: <type 'list'>
            spacing:
            step_size:
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

    def plot_fed_series(
        self,
        bias=0.,
        properties=None,
        color_list=None,
        i_cnt=0,
        hover_text_col=None,
        plot_mode="all",
        ):
        """
        Process data for FED plot.

        Args:
            bias:
            properties:
            color_list:
            i_cnt:
            hover_text_col:

        #FIXME | This is  fairly rough as of right now
        """
        #| - plot_fed_series
        key = properties
        if type(key) == tuple:
            pass
        else:
            key = (key,)

        e_list = self.energy_lst
        e_list = self.apply_bias(bias, e_list)

        overpot_i = self.overpotential

        for n, i in enumerate(e_list):
            if np.isnan(i) is True:
                e_list[n] = None

        name_i = "_".join([str(i) for i in key]) + \
            " (OP: " + str(round(overpot_i, 2)) + ")"

        #| - Hover Text
        if hover_text_col is not None:
            hover_text_list = self.property_list(hover_text_col)
        else:
            hover_text_list = [np.nan for j_cnt in list(range(5))]
        #__|

        dat_lst = self.create_plotly_series(
            e_list,
            group=name_i,
            name=name_i,
            hover_text=hover_text_list,
            color=color_list[i_cnt - 1],
            plot_mode=plot_mode,
            )

        return(dat_lst)
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
        Create a plotly series for the current instance.

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

        if hover_text is None:
            hover_text = [np.nan for i_ind in range(5)]

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

        #| - Creating x-data in middle of states
        short_y = np.array(y_dat)[::2]

        xdat = list(set(new_x_dat))
        xdat.sort()

        cnt = 0
        short_x = []
        for i_ind in range(len(xdat) / 2):
            short_x.append(xdat[cnt] + 0.5)  # FIXME Replace 0.5 with variable
            cnt += 2
        #__|

        #| - Plotly Scatter Plot Instances

        #| - Thick horizontal state lines
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
        #__|

        #| - Full, thin line
        data_2 = Scatter(
            x=new_x_dat,
            y=new_y_dat,
            legendgroup=group,
            name=name,
            connectgaps=True,
            showlegend=show_leg_2,
            hoverinfo="none",
            text=hover_text,

            line=dict(
                color=color,
                width=1,
                ),
            mode="lines",
            )
        #__|

        #| - Points in middle of energy states
        data_3 = Scatter(
            x=short_x,
            y=short_y,
            legendgroup=group,
            name=name,
            showlegend=False,
            hoverinfo="y+text",
            text=hover_text,
            marker=dict(
                size=14,
                color=color,
                opacity=0.,
                ),
            mode="markers",
            )
        #__|

        #__|

        #| - Plot Mode (which data series to plot)
        if plot_mode == "all":
            data_lst = [data_1, data_2, data_3]
        elif plot_mode == "states_only":
            data_lst = [data_1, data_3]
        elif plot_mode == "full_lines":
            data_lst = [data_2, data_3]
        #__|

        return(data_lst)
        #__|


    #__| @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    #__| **********************************************************************


#| - MISC Methods

def calc_ads_e(
    df_row,
    bare_raw_e,
    oxy_ref_e=-443.70964,
    hyd_ref_e=-16.46018,
    ):
    """Calculate adsorption energies from raw DFT energetics.

    Default oxygen reference energy is based on water
    """
    #| - calc_ads_e
    row = df_row
    bare_slab = bare_raw_e
    oxy_ref = oxy_ref_e
    hyd_ref = hyd_ref_e

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
    #__|

def lowest_e_path(tmp=42):
    """Find the lowest energy pathway FED.

    From a set of FE pathways corresponding to different sites, the lowest
    energy states will be selected to construct a new FED.
    """
    #| - lowest_e_path
    print(tmp)
    #__|

#__|
