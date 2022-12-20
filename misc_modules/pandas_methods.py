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
    returns(df)

    Returns:
        Pandas DataFrame: Returns the inputted dataframe with newly ordered colums
    """
    # | - show
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
    # __|

def drop_columns(df=None, columns=None, keep_or_drop="keep"):
    """
    """
    # | - drop_columns
    if keep_or_drop == "keep":
        cols_to_keep = columns

        columns_to_drop = []
        for col_i in df.columns:
            if col_i not in cols_to_keep:
                columns_to_drop.append(col_i)

        df_out = df.drop(columns=columns_to_drop)

    elif keep_or_drop == "drop":
        cols_to_drop = columns

        df_out = df.drop(columns=cols_to_drop)

    else:
        raise ValueError('BAD BAD BAD')

    return(df_out)
    # __|

def drop_nonunique_cols(df=None, ):
    """
    """
    #| - drop_nonunique_cols
    group = df

    nunique = group.apply(pd.Series.nunique)
    cols_to_drop = nunique[nunique == 1].index
    group = group.drop(cols_to_drop, axis=1)

    return(group)
    #__|
