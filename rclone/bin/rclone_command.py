#!/usr/bin/env python

"""TEMP.

Author(s): Raul A. Flores
"""

#| - IMPORT MODULES
import os

import argparse
#__|

#| - Argument Parsing Setup
parser = argparse.ArgumentParser("rclone_command")
parser.add_argument("item", help="folder/file to sync", type=str)
parser.add_argument("--rclone_mode", help="Either sync or copy")
parser.add_argument(
    "--run_command",
    help="Whether to actually run the rclone command or not", 
    default=False,
    )
args = parser.parse_args()
#__|

#| - Input Variables
#  selected_obj = "dos_calc"
selected_obj = args.item


#  remote = "rclone_dropbox"
#  remote = "rclone_gdrive"
remote = "rclone_gdrive_stanford"

if args.rclone_mode:
    command_type = args.rclone_mode
else:
    command_type = "sync"
#__|

#| - Find computer system
compenv = os.environ["COMPENV"]
if compenv == "nersc":
    root_dir = "/global/cscratch1/sd/flores12"
elif compenv == "sher":
    tmp = 42
else:
    print("Couldn't determine system")
#__|

#| - Getting destintation path
root_dir_len = len(root_dir)

cwd = os.getcwd()
#  print("cwd:", cwd)

cwd_rel = cwd[root_dir_len + 1:]
#  print("cwd_rel:", cwd_rel)

selected_obj_rel_path = os.path.join(
    cwd_rel,
    selected_obj)

#  print("selected_obj_rel_path:", selected_obj_rel_path)
#__|


remote_name = os.environ[remote]

tmp = "" + \
    "rclone " + \
    command_type + \
    " " + \
    selected_obj + \
    " " + \
    remote_name + \
    ":" + \
    selected_obj_rel_path

print("")
print(tmp)

#  print(tmp.split(" ")
print("")
[print(i) for i in tmp.split(" ")]

if args.run_command:
    os.system(tmp)
else:
    print("COMMAND NOT RUN (set --run_command flag to True)")
