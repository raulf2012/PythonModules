"""My custom pandas methods and protocols.


Author: Raul A. Flores
"""

__author__ = "Raul A. Flores"
__copyright__ = "NaN"
__version__ = "0.1"
__maintainer__ = "Raul A. Flores"
__email__ = "flores12@stanford.edu; raulf2012@gmail.com"
__date__ = "181115"

#| - Import modules
import os
import pickle
from pathlib import Path
#__|


class ProcessLargeDataFrame:
    """Convienence class to carry out the processing of a large dataframe in
    chunks.

    Will do a subset of the overall dataframe, save that chunk to file, and
    when you rerun the script it will read that data, and the new processed
    data will be added onto it.

    Example script:
        $git_repos/PythonModules/misc_modules/test/ProcessLargeDataFrame__example_script.py
    """

    #| - ProcessLargeDataFrame ************************************************
    _TEMP = "TEMP"


    def __init__(self,
        df=None,
        chunk_size=10,
        unique_keys=None,
        save_dir='./temp_save_dir',
        temp_data_keys=['index', 'data', ]
        ):
        """
        """
        #| - __init__

        #| - Setting Argument Instance Attributes
        self.df=df
        self.chunk_size=chunk_size
        self.unique_keys=unique_keys
        self.save_dir=save_dir
        self.temp_data_keys=temp_data_keys
        #__|

        #| - Initializing Internal Instance Attributes
        #__|

        if not os.path.exists(save_dir):
            os.makedirs(save_dir)

        temp_data = self._get_temp_data()
        self.temp_data = temp_data
        #__|


    def method_0(self):
        """
        """
        #| - method_0
        tmp = 42
        #__|

    def df_drop_rows(self):
        df = self.df
        chunk_size = self.chunk_size
        temp_data = self.temp_data

        # temp_data = PLDF._get_temp_data

        df = df.drop(temp_data['index'], inplace=False)

        if df.shape[0] < chunk_size:
            df = df
        else:
            df = df.sample(n=chunk_size)
        self.df = df
        # return(df)

    def save_temp_data(self):
        #| - save_temp_data
        temp_data = self.temp_data
        unique_keys = self.unique_keys
        save_dir = self.save_dir

        file_path = self._gen_filepath_temp_data(unique_keys, save_dir=save_dir)
        with open(file_path, 'wb') as fle:
            pickle.dump(temp_data, fle)
        #__|

    def _gen_filepath_temp_data(self, unique_keys, save_dir=None):
        """Generate the filepath to the temp datafile."""
        #| - _gen_filepath_temp_data
        file_name_prt = ''
        for key in unique_keys:
            file_name_prt += key[1] + '__'
        file_name_prt += 'temp_data.pickle'
        file_name = file_name_prt

        if save_dir is not None:
            file_path = os.path.join(save_dir, file_name)
            return(file_path)
        else:
            return(file_name)
        #__|

    def _get_temp_data(self):
        """Retrieve the index_list and feat_list objects saved to file.

        These objects contained the index label and feature data for
        the datasets and will be processed into the final feature dataframes.
        """
        #| - _get_temp_data
        unique_keys = self.unique_keys
        temp_data_keys = self.temp_data_keys
        save_dir = self.save_dir

        file_path = self._gen_filepath_temp_data(unique_keys, save_dir=save_dir)

        if Path(file_path).is_file():
            with open(file_path, 'rb') as fle:
                temp_data = pickle.load(fle)
        else:
            temp_data = dict()
            for key in temp_data_keys:
                temp_data[key] = []

        return(temp_data)
        #__|

    #__| **********************************************************************


# ---------------------------------------------------------


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
