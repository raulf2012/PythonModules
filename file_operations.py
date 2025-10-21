"""Methods for carrying out file operations in python.

Author: Raul A. Flores
Date: 2022-11-07
"""

#| - Import modules
import os
import pandas as pd
from datetime import datetime
from shutil import copyfile

from pathlib import Path
#__|



def copy_files_to_dir(
    files_to_copy,
    folder_path,
    warning_instead_of_error=False,
    verbose=False,
    ):
    """Copy list of files to folder.

    Args:
      files_to_copy:
        list of files to transfer, can also include a list within this list, where the first entry is the source file location and the second entry is an alternative filename to give the file when it's moved to folder_path
      folder_path:
        Directory to move files into
    """
    #| - copy_files_to_dir
    for file in files_to_copy:

        if type(file) == list:
            filename_new = file[-1]
            filename = filename_new
            file = file[0]

        else:
            if '/' in file:
                filename = file.split('/')[-1]
            else:
                filename = file

        if verbose:
            print('')
            print(file)
            print(os.path.join(folder_path, filename))
        if not warning_instead_of_error:
            copyfile(file, os.path.join(folder_path, filename))
        else:
            my_file = Path(file)
            if my_file.is_file():
                copyfile(file, os.path.join(folder_path, filename))
            else:
                print("WARNING FILE WASN'T FOUND:", file)

    #__|

def count_file_lines_words_chars(file_path):
    """Counts the number of lines, words, and characters in a plain text file.

    https://stackoverflow.com/questions/41504428/find-the-number-of-characters-in-a-file-using-python
    """
    #| - count_file_lines_words_chars
    with open(file_path) as infile:
        words = 0
        characters = 0
        for lineno, line in enumerate(infile, 1):
            wordslist = line.split()
            words += len(wordslist)
            characters += sum(len(word) for word in wordslist)
    
    return lineno, words, characters
    #__|



def list_files_info(
    directory,
    recursive=False,
    ignore_dirs=None,
    write_to_file=True,
    ):
    """
    ignore_dirs: Ignores files within folders matching this

    Developed this on 2025-10-21 (RF)
    """
    #| - list_files_info


    def add_file_info(entry):
        """Expects a PosixPath object."""
        stats = entry.stat()
        file_data.append({
            'File Name': entry.name,
            'Full Path': entry.as_posix(),
            'Size (KB)': round(stats.st_size / 1024, 2),
            'Created': datetime.fromtimestamp(stats.st_ctime),
            'Modified': datetime.fromtimestamp(stats.st_mtime),
            'Accessed': datetime.fromtimestamp(stats.st_atime),
            'Extension': entry.suffix,
        })




    file_data = []

    if recursive:
        for root, _, files in os.walk(directory):
            for name in files:
                path = os.path.join(root, name)

                my_file = Path(path)

                in_ignore_dir = False
                for dir_i in my_file.parent.parts:
                    if dir_i in ignore_dirs:
                        in_ignore_dir = True

                if my_file.is_file() and not in_ignore_dir:
                    add_file_info(my_file)
    else:
        print('Deprecated functionality')
        # for entry in os.scandir(directory):
        #     if entry.is_file():
        #         add_file_info(entry)

    df = pd.DataFrame(file_data)

    if not df.empty:
        df = df.sort_values(by='Modified', ascending=False).reset_index(drop=True)

    if write_to_file:
        directory = './out_data'
        if not os.path.exists(directory):
            os.makedirs(directory)

        df.to_csv('./out_data/df_files_info.csv')

    return df

    #__|
