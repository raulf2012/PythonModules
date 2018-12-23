"""My custom pandas methods and protocols.


Author: Raul A. Flores
"""

__author__ = "Raul A. Flores"
__copyright__ = "NaN"
__version__ = "0.1"
__maintainer__ = "Raul A. Flores"
__email__ = "flores12@stanford.edu; raulf2012@gmail.com"
__date__ = "181115"

def reorder_df_columns(col_order_list, df):
    """
    """
    #| - show

    #| - __old__
    # col_order_list = [
    #     # Main system variables
    #     "bulk_system",
    #     "facet",
    #     "adsorbate",
    #     "coverage_type",
    #     "ooh_direction",
    #
    #     # Energetics
    #     "ads_e",
    #     "elec_energy",
    #
    #     # Magnetic Moments
    #     "magmoms",
    #     "total_magmom",
    #     "abs_magmom",
    #
    #     "path_short",
    #
    #     "name_i",
    #
    #     "max_force",
    #     "sum_force",
    #     "elem_num_dict",
    #
    #     "incar",
    #     "incar_parsed",
    #
    #     # Atom properties
    #     "init_atoms",
    #     "atoms_object",
    #     "N_atoms",
    #
    #     "dipole_correction",
    #     "u_correction",
    #
    #     # Low priority
    #     "job_type",
    #     "max_revision",
    #     "revision_number",
    #     "success",
    #     "coverage",
    #
    #
    #     "path",
    #     "name_i_2",
    #     "name_i_3",
    #
    #
    #     # Not needed columns
    #     "Job",
    #     "layers",
    #     ]
    #__|

    col_order_list.reverse()

    df_col_list = list(df)
    for col_i in col_order_list:

        if col_i in df_col_list:
            df_col_list.remove(col_i)
            df_col_list.insert(0, col_i)
        else:
            pass

    df = df[df_col_list]

    return(df)
    #__|
