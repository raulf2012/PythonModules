"""Methods for carrying out file operations in python.

Author: Raul A. Flores
Date: 2022-11-07
"""

#| - Import modules
import os
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
