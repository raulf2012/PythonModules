import pandas as pd
from tqdm import tqdm
#  from pandas_process_large_df import ProcessLargeDataFrame
from misc_modules.pandas_methods import ProcessLargeDataFrame

df = pd.read_csv('https://raw.githubusercontent.com/mwaskom/seaborn-data/master/iris.csv')

# The saved data will have to be keyed and unique
unique_keys = [
    ['key_0', 'TEMP0'],
    ['key_1', 'TEMP1'],
    ['key_2', 'TEMP2'],
    ]

PLDF = ProcessLargeDataFrame(
    df=df,
    chunk_size=10,
    unique_keys=unique_keys,
    save_dir='./temp_save_dir',
    temp_data_keys=['index', 'test_feature_0', 'test_feature_1', ]
    )

PLDF.df_drop_rows()

temp_data = PLDF.temp_data

df_ = PLDF.df

# ---------------------------------------------------------
iterator = tqdm(df_.iterrows(), total=df_.shape[0])
for i_cnt, (index, row) in enumerate(iterator):
    y_0 = (row.sepal_length) ** 2 + (row.petal_length) ** 2
    y_1 = (row.sepal_width) ** 2 + (row.petal_width) ** 2

    temp_data['index'].append(index)
    temp_data['test_feature_0'].append(y_0)
    temp_data['test_feature_1'].append(y_1)
# ---------------------------------------------------------


PLDF.save_temp_data()
