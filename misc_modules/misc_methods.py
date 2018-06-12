"""Miscellaneous and helpful methods for everday work."""

#| - IMPORT MODULES
import os
import sys
#__|


def even_spaced_range(start_finish, spacing):
    """
    Produces evenly spaced list between "start" and "finish" with the interval
    between each entry defined by "spacing".
    Includes start and finish at the beginning and end, respecively, regardless
    of whether finish is naturally reached.


    Args:
        start_finish: <type 'list'>
            List with start and end points of ENCUT_range
            ex. [1, 87]
        spacing:
    """
    #| - even_spaced_range
    out_list = []
    entry_i = start_finish[0]
    while entry_i < start_finish[1]:
        out_list.append(entry_i)
        entry_i += spacing
    out_list.append(start_finish[1])

    return(out_list)
    #__|


def merge_two_dicts(x, y):
    """
    """
    #| - merge_two_dicts
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None

    return z
    #__|



#| - File and Directory Management

def remove_file_from_all_folders(
    file_name_list,
    remove_files=False,
    root_dir=".",
    ):
    """Removes charge files (CHG and CHGCAR) from Raman DFT job folders.

    Args:
        file_name_list:
        remove_files: Whether to perform file deletion or not.
        root_dir: Starting directory, assumed to be current directory
    """
    #| - remove_file_from_all_folders
    for dirName, subdirList, fileList in os.walk(root_dir):
        for file_name in file_name_list:
            if file_name in fileList:
                print(dirName)
                print(subdirList)
                print(fileList)
                print(file_name)
                print("_________")

                if remove_files:
                    try:
                        os.remove(dirName + "/" + file_name)
                        print("File removed!!")
                    except:
                        pass
    #__|


#__|
