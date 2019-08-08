#!/usr/bin/env python

"""ORR/OER adsorbate energy scaling class and methods.

Author: Raul A. Flores
"""


#| - IMPORT MODULES
# import numpy as np
#
# import pandas as pd
# pd.options.mode.chained_assignment = None
#
# from sklearn.linear_model import LinearRegression
#
# # import plotly.plotly as py
# import plotly.graph_objs as go
#
#
# from orr_reaction.orr_series import ORR_Free_E_Series
#__|


class Adsorbate_Scaling:
    """ORR/OER adsorbates' dG/E of adsorption scaling.

    Development Notes:
    """

    #| - Adsorbate_Scaling ****************************************************

    def __init__(self,
        tmp=42
        ):
        """
        """
        #| - __init__
        self.tmp = tmp

        #__|

    def tmp_meth(self,
        ):
        """
        """
        #| - tmp
        tmp = 42
        return(tmp)
        #__|

    #__| **********************************************************************


def get_g_ooh(m_ooh, b_ooh, g_oh):
    """
    """
    #| - get_g_ooh
    g_ooh = m_ooh * g_oh + b_ooh
    return(g_ooh)
    #__|

def get_g_o(m_o, b_o, g_oh):
    """
    """
    #| - get_g_o
    g_o = m_o * g_oh + b_o
    return(g_o)
    #__|

def get_g_oh(m_oh, b_oh, g_oh):
    """
    """
    #| - get_g_oh
    g_oh = m_oh * g_oh + b_oh
    return(g_oh)
    #__|

def lim_U_i(
    g_oh=None,
    g_o_minus_g_oh=None,

    mech_step=None,  # 'o2_to_ooh', 'ooh_to_o', 'o_to_oh', 'oh_to_h2o'
    gas_molec_dict=None,
    scaling_dict=None,
    rxn_direction="forward",
    ):
    """
    Calculate the following mechanistic step of the OER/ORR:
      O2 + (H+ + e-) --> *OOH

    Args:
        g_oh:
        gas_molec_dict:
        scaling_dict:
        rxn_direction:
    """
    #| - lim_U_i

    #| - Checking Input Types
    if g_oh is None and g_o_minus_g_oh is None:
        raise ValueError("Need to provide either g_oh or g_o_minus_g_oh")

    if g_oh is not None and g_o_minus_g_oh is not None:
        raise ValueError("Please don't provide both g_oh and g_o_minus_g_oh")

    assert gas_molec_dict is not None, "Please provide gas_molec_dict"
    assert scaling_dict is not None, "Please provide the scaling_dict"
    assert mech_step is not None, "Please provide the step to calculate"
    #__|

    #| - linear fit and gas molecule data
    m_ooh = scaling_dict["ooh"]["m"]
    b_ooh = scaling_dict["ooh"]["b"]

    m_o = scaling_dict["o"]["m"]
    b_o = scaling_dict["o"]["b"]

    m_oh = scaling_dict["oh"]["m"]
    b_oh = scaling_dict["oh"]["b"]


    g_o2 = gas_molec_dict["o2"]
    g_h2 = gas_molec_dict["h2"]
    g_h2o = gas_molec_dict["h2o"]
    #__|

    if g_o_minus_g_oh is not None:
        """
        (G_O-G_OH) = m_o*G_OH + b_o - (G_OH)
        (G_O-G_OH) - b_o = G_OH*(m_o - 1)
        G_OH = [(G_O-G_OH) - b_o] / (m_o - 1)
        """
        g_oh = (g_o_minus_g_oh - b_o) / (m_o - 1)

    elif g_oh is not None:
        g_oh = g_oh

    #| - Calculating Limiting Potential for all legs
    if mech_step == "o2_to_ooh":
        lim_U_out = get_g_ooh(m_ooh, b_ooh, g_oh) + \
            - g_o2

    elif mech_step == "ooh_to_o":
        lim_U_out = get_g_o(m_o, b_o, g_oh) + \
            g_h2o + \
            - get_g_ooh(m_ooh, b_ooh, g_oh)

    elif mech_step == "o_to_oh":
        lim_U_out = get_g_oh(m_oh, b_oh, g_oh) + \
            - get_g_o(m_o, b_o, g_oh)

    elif mech_step == "oh_to_h2o":
        lim_U_out = g_h2o + \
            - get_g_oh(m_oh, b_oh, g_oh)
    else:
        raise ValueError("Woops, error here (9sdfijsd9)")
    #__|

    if rxn_direction == "forward":
        lim_U_out = - lim_U_out
    elif rxn_direction == "reverse":
        lim_U_out = + lim_U_out

    return(lim_U_out)
    #__|




#| - __old__

# def lim_U_o2_to_ooh(g_oh,
# gas_molec_dict, scaling_dict, rxn_direction="forward"):
#     """
#     Calculate the following mechanistic step of the OER/ORR:
#       O2 + (H+ + e-) --> *OOH
#
#     Args:
#         g_oh:
#         gas_molec_dict:
#         scaling_dict:
#         rxn_direction:
#     """
#     #| - lim_U_o2_to_ooh
#
#     #| - linear fit and gas molecule data
#     m_ooh = scaling_dict["ooh"]["m"]
#     b_ooh = scaling_dict["ooh"]["b"]
#
#     m_o = scaling_dict["o"]["m"]
#     b_o = scaling_dict["o"]["b"]
#
#     m_oh = scaling_dict["oh"]["m"]
#     b_oh = scaling_dict["oh"]["b"]
#
#
#     g_o2 = gas_molec_dict["o2"]
#     g_h2 = gas_molec_dict["h2"]
#     g_h2o = gas_molec_dict["h2o"]
#     #__|
#
#
#     if False:
#
#         g_oh = (TMP - b_o) / (m_o - 1)
#
#     lim_U_out = get_g_ooh(m_ooh, b_ooh, g_oh) - g_o2
#
#     if rxn_direction == "forward":
#         lim_U_out = - lim_U_out
#     elif rxn_direction == "reverse":
#         lim_U_out = + lim_U_out
#
#     return(lim_U_out)
#     #__|
#
# def lim_U_ooh_to_o(g_oh,
# gas_molec_dict, scaling_dict, rxn_direction="forward"):
#     """
#     Calculate the following mechanistic step of the OER/ORR:
#       *OOH + (H+ + e-) --> *O + H2O
#
#     Args:
#         g_oh:
#         gas_molec_dict:
#         scaling_dict:
#         rxn_direction:
#     """
#     #| - lim_U_ooh_to_o
#
#     #| - linear fit and gas molecule data
#     m_ooh = scaling_dict["ooh"]["m"]
#     b_ooh = scaling_dict["ooh"]["b"]
#
#     m_o = scaling_dict["o"]["m"]
#     b_o = scaling_dict["o"]["b"]
#
#     m_oh = scaling_dict["oh"]["m"]
#     b_oh = scaling_dict["oh"]["b"]
#
#
#     g_o2 = gas_molec_dict["o2"]
#     g_h2 = gas_molec_dict["h2"]
#     g_h2o = gas_molec_dict["h2o"]
#     #__|
#
#     lim_U_out = get_g_o(m_o,
# b_o, g_oh) + g_h2o - get_g_ooh(m_ooh, b_ooh, g_oh)
#
#     if rxn_direction == "forward":
#         lim_U_out = - lim_U_out
#     elif rxn_direction == "reverse":
#         lim_U_out = + lim_U_out
#
#     return(lim_U_out)
#     #__|
#
# def lim_U_o_to_oh(g_oh,
# gas_molec_dict, scaling_dict, rxn_direction="forward"):
#     """
#     Calculate the following mechanistic step of the OER/ORR:
#       O* + (H+ + e-) --> *OH
#
#     Args:
#         g_oh:
#         gas_molec_dict:
#         scaling_dict:
#         rxn_direction:
#     """
#     #| - lim_U_o_to_oh
#
#     #| - linear fit and gas molecule data
#     m_ooh = scaling_dict["ooh"]["m"]
#     b_ooh = scaling_dict["ooh"]["b"]
#
#     m_o = scaling_dict["o"]["m"]
#     b_o = scaling_dict["o"]["b"]
#
#     m_oh = scaling_dict["oh"]["m"]
#     b_oh = scaling_dict["oh"]["b"]
#
#
#     g_o2 = gas_molec_dict["o2"]
#     g_h2 = gas_molec_dict["h2"]
#     g_h2o = gas_molec_dict["h2o"]
#     #__|
#
#     lim_U_out = get_g_oh(m_oh, b_oh, g_oh) - get_g_o(m_o, b_o, g_oh)
#
#     if rxn_direction == "forward":
#         lim_U_out = - lim_U_out
#     elif rxn_direction == "reverse":
#         lim_U_out = + lim_U_out
#
#     return(lim_U_out)
#     #__|
#
# def lim_U_oh_to_h2o(g_oh,
# gas_molec_dict, scaling_dict, rxn_direction="forward"):
#     """
#     Calculate the following mechanistic step of the OER/ORR:
#       *OH + (H+ + e-) --> H2O
#
#     Args:
#         g_oh:
#         gas_molec_dict:
#         scaling_dict:
#         rxn_direction:
#     """
#     #| - lim_U_oh_to_h2o
#
#     #| - linear fit and gas molecule data
#     m_ooh = scaling_dict["ooh"]["m"]
#     b_ooh = scaling_dict["ooh"]["b"]
#
#     m_o = scaling_dict["o"]["m"]
#     b_o = scaling_dict["o"]["b"]
#
#     m_oh = scaling_dict["oh"]["m"]
#     b_oh = scaling_dict["oh"]["b"]
#
#
#     g_o2 = gas_molec_dict["o2"]
#     g_h2 = gas_molec_dict["h2"]
#     g_h2o = gas_molec_dict["h2o"]
#     #__|
#
#     lim_U_out = g_h2o - get_g_oh(m_oh, b_oh, g_oh)
#
#     if rxn_direction == "forward":
#         lim_U_out = - lim_U_out
#     elif rxn_direction == "reverse":
#         lim_U_out = + lim_U_out
#
#     return(lim_U_out)
#     #__|
#
#__|
