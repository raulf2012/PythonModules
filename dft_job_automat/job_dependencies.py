#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Methods and code to handle jobs that depend on one another."""

#| - Import Modules
import sys
import os
import shutil
import copy

from ase import io
import pandas as pd

from dft_job_automat.job_analysis import DFT_Jobs_Analysis
from dft_job_automat.job_manager import DFT_Jobs_Manager
#__|

#| - FUNCTIONS

def full_path_i(root_dir, step_dir_names, step, path):
    """Formats relative path to full path for a given calculation step
    """
    #| - full_path_i
    path_out = root_dir + "/" + step_dir_names[step - 1] + "/" + path

    return(path_out)
    #__|

def copyfiles_onestep_up(
    job_var_lst,
    step,
    JobsInstances_lst,
    files_lst=[],
    root_dir_files=None
    ):
    """
    """
    #| - copyfiles_onestep_up
    def copy_if_not_in_dest(source, dest_file):
        """
        """
        #| - copy_if_not_in_dest
        if not os.path.isfile(dest_file):

            shutil.copy(source, dest_file)
            source_fle = source.split("/")[-1]

            print("File " + str(source_fle) + " copied")
        else:
            pass

        #__|


    curr_step = JobsInstances_lst[step - 1]
    next_step = JobsInstances_lst[step]


    root_dir = curr_step.root_dir

    dir_curr = curr_step.var_lst_to_path(
        job_var_lst,
        job_rev="Auto",
        relative_path=False,
        )
    next_var_lst = copy.deepcopy(next_step.job_var_lst)

    #| - Finding jobs in step+1 whose properties match step's
    for curr_property in job_var_lst:
        prop_name = curr_property["property"]  # force-cutoff
        prop_value = curr_property["value"]    # 0.001

        for job in next_step.job_var_lst:
            path_i = next_step.var_lst_to_path(
                job,
                job_rev="Auto",
                relative_path=False,
                )

            # match_lst = []
            for next_property in job:
                next_prop_name = next_property["property"]
                next_prop_value = next_property["value"]

                if next_prop_name == prop_name:
                    if next_prop_value != prop_value:

                        if job in next_var_lst:
                            next_var_lst.remove(job)
                        else:
                            pass
    #__|

    for next_job in next_var_lst:
        path_i = next_step.var_lst_to_path(
            next_job,
            job_rev="Auto",
            relative_path=False,
            )

        for file_i in files_lst:

            #| - Copy Files to Directory in New Step
            if type(file_i) == str:
                curr_dir = dir_curr + \
                    curr_step.cluster.cluster.job_data_dir + "/" + file_i
                copy_if_not_in_dest(curr_dir, path_i + "/" + file_i)

            elif type(file_i) == list:
                curr_dir = dir_curr + \
                    curr_step.cluster.cluster.job_data_dir + "/" + file_i[0]
                copy_if_not_in_dest(curr_dir, path_i + "/" + file_i[1])
            #__|


        if root_dir_files is not None:
            root_dir = curr_step.root_dir
            source_dir = root_dir

            dest_dir = path_i

            for file_i in root_dir_files:
                #| - Copy Files from Root Dir to New Job Folder
                if type(file_i) == str:
                    copy_if_not_in_dest(
                        source_dir + "/" + file_i,
                        dest_dir + "/" + file_i,
                        )

                elif type(file_i) == list:
                    copy_if_not_in_dest(
                        source_dir + "/" + file_i[0],
                        dest_dir + "/" + file_i[1],
                        )
                #__|


    open(dir_curr + "/.FILES_COPIED", "w")
    # file = open(dir_curr + "/.FILES_COPIED", "w")

    #__|

def create_atoms_list(atoms_name, file_ext, root_dir):
    """
    """
    #| - create_atoms_list
    atoms_dict = {}
    for atom in atoms_name:
        atoms_i = io.read(root_dir + "/dir_atoms/" + atom + file_ext)
        atoms_dict[atom] = atoms_i

    return(atoms_dict)
    #__|

def create_level_entries_dict(tree_level_labels, tree_level_values):
    """
    """
    #| - create_level_entries_dict
    level_entries_dict = {}
    for index, variable in enumerate(tree_level_labels):
        level_entries_dict[variable] = tree_level_values[index]

    return(level_entries_dict)
    #__|

#| - __OLD__

def job_runnable(df, root_dir_beg, path_i):
    """
    """
    #| - job_runnable
    df["full_path"] = root_dir_beg + "/" + df["root_dir"] + "/" + df["path"] + \
        "_" + df["job_revision_number"].astype(str)

    index = df.index[df["full_path"] == path_i].tolist()
    df_i = df.iloc[index]
    max_rev = df_i["revision_number"].max()
    df_fin = df_i[df_i["revision_number"] == max_rev]
    assert len(df_fin) == 1

    if df_fin["job_state_2"].iloc[0] == "RUNNABLE":
        return(True)
    else:
        return(False)
    #__|

def job_failed(df, root_dir_beg, path_i):
    """
    """
    #| - job_failed
    df["full_path"] = root_dir_beg + "/" + df["root_dir"] + "/" + df["path"] + \
        "_" + df["job_revision_number"].astype(str)


    index = df.index[df["full_path"] == path_i].tolist()
    df_i = df.iloc[index]
    max_rev = df_i["revision_number"].max()
    df_fin = df_i[df_i["revision_number"] == max_rev]
    # assert len(index) == 1  # Checking that job path is unique
    assert len(df_fin) == 1


    # index = df.index[df["full_path"] == path_i].tolist()
    # assert len(index) == 1  # Checking that job path is unique

    # job_state = df.iloc[index[0]]["job_state_2"]

    # if job_state == "FAILED":
    # return(True)


    if df_fin["job_state_2"].iloc[0] == "FAILED":
        return(True)
    else:
        return(False)

    #__|

#__|

#__|

class DFT_Jobs_Workflow:
    """Summary line.

    Useful class to set up multiple DFT job in a directory structure.
    Must be initialized with tree_level and level_entries inputs
    """

    def __init__(self,
        atoms_prefix=".traj",
        atoms_list_names=["init"],
        model_names=None,

        tree_level_labels_list=None,
        tree_level_values_list=None,
        indiv_dir_lst_list=None,
        indiv_job_lst_list=None,

        setup_function=None,
        maint_function=None,
        number_of_steps=1,
        root_dir=".",
        run_jobs=False,
        ):
        """Initialize DFT_Jobs_Workflow instance.

        Args:
            atoms_prefix:
            atoms_list_names:
            model_names:

            tree_level_labels_list:
            tree_level_values_list:
            indiv_dir_lst_list:

            setup_function:
            maint_function:
            number_of_steps:
            root_dir:
            run_jobs:
        """
        #| - __init__
        self.mod_dir = "dir_models"
        self.atoms_dir = "dir_atoms"

        self.atoms_ext = atoms_prefix
        self.atoms_list_names = atoms_list_names

        self.num_steps = number_of_steps

        self.tree_level_labels_list = self.list_of_None_if_None(
            tree_level_labels_list,
            )

        self.tree_level_values_list = self.list_of_None_if_None(
            tree_level_values_list,
            )

        self.indiv_dir_lst_list = self.list_of_None_if_None(
            indiv_dir_lst_list,
            )

        self.indiv_job_lst_list = self.list_of_None_if_None(
            indiv_job_lst_list,
            )

        # self.model_names = model_names

        self.setup_function = setup_function
        self.maint_function = maint_function

        self.run_jobs = run_jobs

        self.root_dir = self.__set_cwd__(root_dir)
        # self.atoms_dict = create_atoms_list(
        #     atoms_list_names,
        #     atoms_prefix,
        #     self.root_dir,
        #     )

        self.step_dir_names = self.__set_step_dir_names__()
        self.model_names = self.__set_model_names__(model_names)

        self.jobs_an_list = self.__create_jobs_an__()

        self.__create_parent_dirs__()
        self.__prep_dir_sys__()

        self.jobs_man_list = self.__create_jobs_man__()
        self.__job_maint__()
        #__|

    def list_of_None_if_None(self, input):
        """Return list of 'None' of length == # of steps, if input is None
        """
        #| - list_of_None_if_None
        if input is None:
            none_list = [None for i in range(self.num_steps)]

            return(none_list)
        else:
            return(input)
        #__|

    def __set_cwd__(self, root_dir):
        """Set the working directory.

        Args:
            root_dir:
        """
        #| - __set_cwd__
        if root_dir == ".":
            root_dir_out = os.getcwd()
        else:
            root_dir_out = root_dir

        return(root_dir_out)
        #__|

    def __set_step_dir_names__(self):
        """
        """
        #| - __set_step_dir_names__
        number_of_steps = self.num_steps
        step_dir_names = []
        for step_i in range(number_of_steps):
            step_dir_name_i = str(step_i + 1) + "STEP"
            step_dir_names.append(step_dir_name_i)

        return(step_dir_names)
        #__|

    def __set_model_names__(self, model_names):
        """Return list of model script names.

        Args:
            model_names:
        """
        #| - __set_model_names__
        # model_names = self.model_names
        if model_names is None:
            model_names_list = []
            for step_i in range(self.num_steps):
                model_i = str(step_i) + "model.py"
                model_names_list.append(model_i)
        else:
            model_names_list = model_names

        return(model_names_list)
        #__|

    def __create_jobs_an__(self):
        """Create Jobs_Analysis instances for each step of workflow."""
        #| - __create_jobs_an__
        # print("PREPARING EXTENDED FOLDER SYSTEM")  #PERM_PRINT
        step_dir_names = self.step_dir_names
        master_root_dir = self.root_dir

        Jobs_Inst_list = []
        for step in range(len(step_dir_names)):
            step_num = step + 1

            print("Initializing Job Instance: " + str(step_num))  # PERM_PRINT

            # dir_struct_file = master_root_dir + "/" + step_dir_names[step] + \
            #     "/jobs_bin/dir_structure.json"

            level_labels_tmp = self.tree_level_labels_list[step]
            level_entries_tmp = self.tree_level_values_list[step]

            indiv_dir_lst_tmp = self.indiv_dir_lst_list[step]
            indiv_job_lst_tmp = self.indiv_job_lst_list[step]

            JobsAn = DFT_Jobs_Analysis(
                tree_level=level_labels_tmp,
                level_entries=level_entries_tmp,
                indiv_dir_lst=indiv_dir_lst_tmp,
                indiv_job_lst=indiv_job_lst_tmp,

                root_dir=master_root_dir,
                working_dir=step_dir_names[step],
                update_job_state=False,
                load_dataframe=False,
                )

            Jobs_Inst_list.append(JobsAn)

        return(Jobs_Inst_list)
        #__|

    def __create_jobs_man__(self):
        """Create Jobs_Manager instance(s)."""
        #| - __create_jobs_man__
        step_dir_names = self.step_dir_names
        master_root_dir = self.root_dir

        Jobs_Inst_list = []
        for step in range(len(step_dir_names)):
            # step_num = step + 1

            level_labels_tmp = self.tree_level_labels_list[step]
            level_entries_tmp = self.tree_level_values_list[step]

            indiv_dir_lst_tmp = self.indiv_dir_lst_list[step]
            indiv_job_lst_tmp = self.indiv_job_lst_list[step]

            Jobs = DFT_Jobs_Manager(


                tree_level=level_labels_tmp,
                level_entries=level_entries_tmp,

                skip_dirs_lst=None,
                indiv_dir_lst=indiv_dir_lst_tmp,  # <-----------------------------------------------
                indiv_job_lst=indiv_job_lst_tmp,

                root_dir=master_root_dir,
                working_dir=step_dir_names[step],

                update_job_state=False,
                load_dataframe=False,

                # tree_level=level_labels_tmp,
                # level_entries=level_entries_tmp,
                # working_dir=master_root_dir + "/" + step_dir_names[step],
                # load_dataframe=False,
                )

            Jobs_Inst_list.append(Jobs)

        return(Jobs_Inst_list)
        #__|

    def __create_parent_dirs__(self):
        """Create parent folders."""
        #| - __create_parent_dirs__
        step_dir_names = self.step_dir_names
        master_root_dir = self.root_dir

        for step in range(len(step_dir_names)):
            # step_num = step + 1

            step_folder = master_root_dir + "/" + step_dir_names[step]
            if not os.path.isdir(step_folder):
                os.makedirs(step_folder)
        #__|

    def __prep_dir_sys__(self):
        """
        """
        #| - __prep_dir_sys
        print("PREPARING EXTENDED FOLDER SYSTEM")  # PERM_PRINT
        step_dir_names = self.step_dir_names
        master_root_dir = self.root_dir

        self.__create_parent_dirs__()

        wf_vars = vars(self)
        for step in range(len(step_dir_names)):
            # step_num = step + 1
            JobsAn = self.jobs_an_list[step]

            print("Placing Initial Files in Folders | LOOP OVER JOBS")
            files_placed_file = master_root_dir + "/" + \
                step_dir_names[step] + "/.FILES_PLACED"

            # if not os.path.isfile(files_placed_file):
            if True:

                #| - Create Step Folder Structure
                JobsAn.create_dir_struct(create_first_rev_folder="True")
                #__|

                for Job_i in JobsAn.Job_list:
                    path_i = Job_i.full_path

                    job_i_params = Job_i.job_params

                    self.setup_function(step, path_i, job_i_params, wf_vars)

                # for job_i in JobsAn.job_var_lst:
                #     path_i = JobsAn.var_lst_to_path(
                #         job_i,
                #         job_rev="Auto",
                #         relative_path=False,
                #         )
                #
                #     #| - Job_i Parameters
                #     job_i_params = {}
                #     for variable in JobsAn.tree_level_labels:
                #         job_i_var_j = JobsAn.extract_prop_from_var_lst(
                #             job_i,
                #             variable,
                #             )
                #         job_i_params[variable] = job_i_var_j
                #     #__|
                #
                #     self.setup_function(step, path_i, job_i_params, wf_vars)

                file_name = master_root_dir + "/" + step_dir_names[step]
                file = open(file_name + "/.FILES_PLACED", "w")

            # file = open(master_root_dir + "/.FOLDERS_CREATED", "w")
            open(master_root_dir + "/.FOLDERS_CREATED", "w")
        #__|

    def __job_maint__(self):
        """Manage jobs after being submitted.

        Tries to figure out what state the job is in and acts according to
        user defined methods
        """
        #| - __job_maint__
        step_dir_names = self.step_dir_names
        # master_root_dir = self.root_dir

        for step in range(len(step_dir_names)):
            step_num = step + 1
            Jobs = self.jobs_man_list[step]

            # df = Jobs.data_frame
            tally = {"successes": 0, "failures": 0, "running": 0, "pending": 0}

            #| - PRINT
            print("")  # PERM_PRINT
            print("###########################################################")

            str_i = "########################"
            str_i += " STEP " + str(step_num) + " ###########################"
            # "######################## STEP " + str(step_num) +
            # " ###########################"
            print(str_i)
            print("###########################################################")
            print("Total Jobs: " + str(Jobs.num_jobs))  # PERM_PRINT
            #__|

            #| - LOOP OVER JOBS
            if self.run_jobs:

                tally = {
                    "successes": 0,
                    "failures": 0,
                    "running": 0,
                    "pending": 0,
                    }

                wf_vars = vars(self)

                for Job_i in Jobs.Job_list:
                    path_i = Job_i.full_path
                    job_i_params = Job_i.job_params

                    # Why is this being run again #COMBAK
                    self.setup_function(step, path_i, job_i_params, wf_vars)

                    tally = self.maint_function(
                        step,
                        # job_i,
                        path_i,
                        job_i_params,
                        wf_vars,
                        tally,
                        )

                #| - __old__
                # for job_i in Jobs.job_var_lst:
                #     path_i = Jobs.var_lst_to_path(
                #         job_i,
                #         job_rev="Auto",
                #         relative_path=False,
                #         )
                #
                #     #| - Job_i Parameters
                #     job_i_params = {}
                #     for variable in Jobs.tree_level_labels:
                #         job_i_param_j = Jobs.extract_prop_from_var_lst(
                #             job_i,
                #             variable,
                #             )
                #
                #         job_i_params[variable] = job_i_param_j
                #     #__|
                #
                #     tally = self.maint_function(
                #         step,
                #         job_i,
                #         job_i_params,
                #         wf_vars,
                #         tally,
                #         )
                #__|

                    # TODO Check that tally is being incremented by 1 only

                print(tally)

            print("")
            #__|

        #__|


#| - Reinitiating the Jobs Instances (Is this needed)
# print("")  #PERM_PRINT
# print("Reinitiating the Jobs Instances")  #PERM_PRINT
# Jobs_Inst_list = []
# for step in range(len(step_dir_names)):
#     step_num = step + 1
#
#     # print("Initializing Job Instance: " + str(step_num))
#     Jobs = DFT_Jobs_Analysis(
#         system="aws",
#         tree_level=tree_level_labels_list[step],
#         level_entries=level_entries_dict_list[step],
#         working_dir=master_root_dir + "/" + step_dir_names[step],
#         load_dataframe=False,
#         )
#
#     Jobs_Inst_list.append(Jobs)
#__|
