"""Methods to interface with mysql built on mysql-connector-python


I was working in the following dir when I started developing this:
    $dropbox/Career/Job Searches/2023_post_postdoc_jobsearch/Interviews/2024-04-10__Citrine_Informatics/Take Home Exam/code

"""


#| - Import modules
from configparser import ConfigParser

import pandas as pd

from mysql.connector import MySQLConnection, Error

#  from config import read_config
#  table_name = 'user'
#__|



def read_config(filename='app.ini', section='mysql'):
    """Read sql config file, obtained this online from somewhere."""
    #| - read_config
    # Create a ConfigParser object to handle INI file parsing
    config = ConfigParser()

    # Read the specified INI configuration file
    config.read(filename)

    # Initialize an empty dictionary to store configuration data
    data = {}

    # Check if the specified section exists in the INI file
    if config.has_section(section):
        # Retrieve all key-value pairs within the specified section
        items = config.items(section)

        # Populate the data dictionary with the key-value pairs
        for item in items:
            data[item[0]] = item[1]
    else:
        # Raise an exception if the specified section is not found
        raise Exception(f'{section} section not found in the {filename} file')

    # Return the populated data dictionary
    return data
    #__|


def return_table_as_df(config, table_name=None):
    """Return table from database as a pandas dataframe."""
    #| - return_table_as_df

    df = None
    conn = None
    cursor = None
    cursorB = None

    try:
        # Establish a connection to the MySQL database using the provided configuration
        conn = MySQLConnection(**config)


        # Get column name info
        cursorB = conn.cursor()

        cursorB.execute("SHOW COLUMNS FROM " + table_name)

        columns_ = cursorB.description

        columns = []
        for i in columns_:
            columns.append(i[0])

        df_cols = pd.DataFrame(cursorB.fetchall(), columns=columns)



        # Create a cursor to interact with the database
        cursor = conn.cursor()

        # Execute a SELECT query to retrieve all rows from the 'books' table
        cursor.execute("SELECT * FROM " + table_name)

        # Fetch all rows from the result set
        rows = cursor.fetchall()


        df = pd.DataFrame(rows, columns=df_cols.Field.tolist())


    except Error as e:
        # Print an error message if an error occurs during the execution of the query
        print(e)

    finally:
        # Close the cursor and connection in the 'finally' block to ensure it happens
        if conn is not None and conn.is_connected():
            #  cursor.close()
            conn.close()

            if cursor is not None:
                cursor.close()

    return(df)
    #__|




#| - __old__


# if __name__ == '__main__':
#     # Read the database configuration from the 'config' module
#     config = read_config()

#     # Call the function with the obtained configuration to execute the query
#     query_with_fetchall(config)




#  if __name__ == '__main__':
#      # Read the configuration from the default section ('mysql') in 'app.ini'
#      config = read_config()
#
#      # Display the obtained configuration
#      print(config)


#__|
