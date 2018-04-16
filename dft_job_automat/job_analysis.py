#!/usr/bin/env python

"""Class to analyse data using the DFT_Jobs_Setup class.


Development Notes:
    TODO Automaticall read README files contained within job folders
    TODO Delete the duplicate methods for job_status
"""

#| - Import Modules
import sys
import os
# import cPickle as pickle
import pickle
import copy
import glob

import pandas as pd
import numpy as np

from ase.visualize import view

# My Modules
from dft_job_automat.job_setup import DFT_Jobs_Setup
#__|

class DFT_Jobs_Analysis(DFT_Jobs_Setup):
    """Analysis methods for jobs in tree structure.

    # TODO path --> path_i
    # TODO Modify .FINISHED implementation

    Parent class to DFT_Jobs_Setup
    """

    #| - DFT_Jobs_Analysis ****************************************************

    #| - Class Variables
    finished_fle = ".FINISHED.new"
    #__|

    def __init__(self,
        system="sherlock",
        tree_level=None,
        level_entries=None,
        skip_dirs_lst=None,
        working_dir=".",
        update_job_state=True,
        load_dataframe=True,
        job_type_class=None,
        ):
        """Initialize DFT_Jobs_Analysis Instance.

        Args:
            system:
            tree_level:
            level_entries:
            skip_dirs_lst:
            working_dir:
            update_job_state:
                Updates job status for all jobs in ensemble
            load_dataframe:
            job_type_class:
        """
        #| - __init__
        DFT_Jobs_Setup.__init__(self,
            system=system,
            tree_level=tree_level,
            level_entries=level_entries,
            skip_dirs_lst=None,
            working_dir=working_dir,
            )

        #| - General Methods
        if update_job_state is True:
            print("update_job_state == True")
            self.add_data_column(
                self.job_state,
                column_name="job_state",
                )
            self.add_data_column(
                self.job_state_3,
                column_name="N/A",
                allow_failure=False,
                )
        #__|

        #| - Job Type Specific Methods
        # method = DFT_Methods().atom_type_num_dict
        # self.add_data_column(method, column_name="TEMP", allow_failure=False)

        if load_dataframe is True:
            self.__load_dataframe__()

        else:
            if job_type_class is not None:
                job_type_inst = job_type_class

                for method in job_type_inst.methods_to_run:
                    method_ref = getattr(job_type_inst, method)
                    self.method = method_ref
                    self.add_data_column(
                        method_ref,
                        column_name=method,
                        allow_failure=True,
                        )

                    # self.add_data_column(
                    #     method_ref,
                    #     column_name=method,
                    #     allow_failure=False,
                    #     )

                # TEMP
                self.__write_dataframe__()

        self.add_all_columns_from_file()
        #__|

        #__|

    def add_jobs_queue_data(self):
        """Add jobs queue data to jobs.csv file."""
        #| - add_jobs_queue_data
        # COMBAK Not finished
        #__|

    def job_queue_info(self, path_i):
        """
        Return row corresponding to job in path_i from the job queue dir.

        # FIXME Set up the jobs.csv system to work

        Args:
            path_i:
        """
        #| - job_queue_info
        jobs_file_path = self.job_queue_dir + "/jobs.csv"
        df = pd.read_csv(jobs_file_path)

        path_i = path_i[len(self.root_dir):]
        full_path = self.root_dir_short + path_i

        # Looking for job path in jobs.csv file
        index = df[df["job_path"] == full_path].index.tolist()[0]

        job_info = df.iloc[index].to_dict()

        return(job_info)
        #__|

    def __load_dataframe__(self):
        """Attempt to load dataframe."""
        #| - __load_dataframe__
        fle_name = self.root_dir + "/jobs_bin/job_dataframe.pickle"
        with open(fle_name, "rb") as fle:

            if sys.version_info.major > 2:
                df = pickle.load(fle, encoding="latin1")
                # NOTE Added encoding="latin1" for p36 support (180415 - RF)
            else:
                df = pickle.load(fle)

        self.data_frame = df
        #__|

    def __write_dataframe__(self):
        """
        Write dataframe to file in root_dir/jobs_bin folder.

        Writes dataframe in csv and pickle format. CSV is easily human readable
        """
        #| - __write_dataframe__
        df = self.data_frame
        df.to_csv(self.root_dir + "/jobs_bin/job_dataframe.csv", index=False)

        df_pickle_fle = self.root_dir + "/jobs_bin/job_dataframe.pickle"
        with open(df_pickle_fle, "wb") as fle:
            pickle.dump(df, fle)
        #__|

    def add_all_columns_from_file(self):
        """
        Add column to dataframe from files in data_columns folder.

        Files in data_columns folder must have the following format:
            3 columns, value | revision # | path
            Columns are separated by the "|" character
            File name must end with ".col" extension.
        """
        #| - add_all_columns_from_file
        dir_list = glob.glob(self.root_dir + "/jobs_bin/data_columns/*.col*")
        col_data_file_list = [dir.split("/")[-1] for dir in dir_list]

        for col_file in col_data_file_list:
            self.__add_data_column_from_file__(col_file)
        #__|

    def __add_data_column_from_file__(self,
        file_name,
        ):
        """Add data in "file_name" file to dataframe.

        Searches in jobs_bin/data_columns/
        Args:
            file_name:
        """
        #| - __add_data_column_from_file__

        #| - Extracting Column Name From File Name
        col_name = file_name.split(".")[0]
        #__|

        #| - Reading Column Data File - NEW
        column_file = self.root_dir + "/jobs_bin/data_columns/" + file_name
        with open(column_file, "r") as fle:
            content = fle.readlines()

        content = [x.strip().split("|") for x in content]

        content_new = []
        for line in content:
            line_new = {}
            line_new["value"] = line[0].strip()
            line_new["revision"] = line[1].strip()
            line_new["path"] = line[2].strip()

            content_new.append(line_new)
        #__|

        #| - Matching Dataframe with New Data Column
        df = self.data_frame

        df["full_path"] = df["path"].astype(str) + "_" + \
            df["revision_number"].astype(str)

        column_data_list = []
        for i_ind, (index, row) in enumerate(df.iterrows()):

            row_has_entry = False
            for entry in content_new:

                entry_fullpath = entry["path"] + "_" + entry["revision"]
                if entry_fullpath == row["full_path"]:
                    column_data_list.append(entry["value"])
                    row_has_entry = True
                    continue

            if not row_has_entry:
                column_data_list.append(np.nan)

        df.drop("full_path", axis=1)

        df[col_name] = column_data_list
        #__|

        #__|

    def add_data_column(self,
        function,
        column_name="new_column",
        revision="auto",
        allow_failure=True):
        """
        Add data column to data frame by iterating thourgh job folders.

        Args:
            function: <type 'function'>
                Set of operations that will be applied to indiviudal job
              folders. Must return a scalar quantity that will be added to
              data frame. The function should know "what to do" given only a
              path that correpsonds to a unique job folder (including revision).

                If the function returns a dict, then various columns will be
              added with the column names corresponding to the key name.

            column_name: <type 'str'>
                Name of new data column

            revision: <type 'str' or 'int'>
                The job revision from which data is scraped
                "auto" | Selects most recent revision folder
                "all"  | Adds all revisions to data frame (NOT WORKING)
                "previous | Second to last revision"
            allow_failure: <True or False>
                If True, a failed method call will result in NaN
        """
        #| - __add_data_coumn__
        new_data_col = []
        job_rev_lst = []

        print("Adding " + str(column_name))

        for entry in self.data_frame["variable_list"]:

            path = self.var_lst_to_path(
                entry,
                relative_path=False,
                job_rev="False",
                )

            #| - Picking Revision Number(s) To Query
            largest_rev = self.job_revision_number(entry)

            if revision == "auto":
                rev = [largest_rev]
            elif revision == "previous":
                if largest_rev == 1:
                    rev = [1]
                else:
                    rev = [largest_rev - 1]
            elif revision == "all":
                rev = range(self.job_revision_number(entry) + 1)
            else:
                rev = [revision]
            #__|

            for rev_num in rev:

                #| - Run Function
                path += "_" + str(rev_num)

                path = path + self.cluster.cluster.job_data_dir

                if allow_failure is True:
                    try:
                        out = function(path)
                    except:
                        out = np.nan
                else:
                    out = function(path)

                new_data_col.append(out)
                job_rev_lst.append(rev_num)
                #__|

        data_type_list = [type(x) for x in new_data_col]
        dict_in_list = any(item == dict for item in data_type_list)

        if dict_in_list:
            new_col = []
            for x in new_data_col:
                if pd.isnull(x) is True:
                    new_col.append({"NA": np.nan})
                else:
                    new_col.append(x)

            new_columns_df = pd.DataFrame(new_col)

            df1 = self.data_frame
            df2 = new_columns_df

            out_df = pd.concat([df1, df2], axis=1)
            self.data_frame = out_df

        else:
            self.data_frame[column_name] = new_data_col
        #__|

    def job_revisions(self, path):
        """Return number of revision folders for a given job.

        Args:
            path:
        """
        #| - job_revisions
        path = "/".join(path.split("/")[0:-1]) + "/"

        # Attempting to remove duplicate job folders (usually have spaces)
        dir_list = [x for x in os.walk(path).next()[1] if " " not in x]

        return(len(dir_list))
        #__|

    def job_state_file(self, path_i="."):
        """
        Return contents of '.QUEUESTATE' if present in the job directory.

        Args:
            path:
        """
        #| - job_state_file
        file_path = path_i + "/.QUEUESTATE"
        if os.path.isfile(file_path):
            with open(file_path, "r") as fle:
                job_state = fle.read().rstrip()
        else:
            job_state = None

        return(job_state)
        #__|


    #| - Data Frame Methods

    def create_data_sets(self, data_frame, free_variable):
        """
        Splinter data_frame into distinct data sets based on columns.

        Returns data sets from the data frame where the data is split by the
        job variables. One variable is excluded from this grouping and is
        treated as a "continous" variable to be plotted on the x-axis

        ex. Parameter sweep with varying k-points, pw-cutoff, and
        lattice-constant. The lattice-constant is designated as the free
        variable and so data sets are created for each combination of
        k-points and pw-cutoff where each set contains the full range of
        latt-const.
        """
        #| - create_data_sets
        # df = copy.deepcopy(self.data_frame)
        df = data_frame

        var_lst = copy.deepcopy(self.tree_level_labels)
        var_lst.remove(free_variable)

        df_unique_params = df[var_lst].drop_duplicates()
        indices = df_unique_params.index.values

        data_lst = []
        for index in indices:
            df_tmp = df_unique_params.ix[[index]]

            #| - Data Labels
            data_label_full = ""
            for column in df_tmp:
                col_i = df_tmp[column]

                col_name = col_i.name
                col_val = col_i.iloc[0]

                data_label_i = str(col_name) + ": " + str(col_val)

                data_label_full += data_label_i + " | "

            data_label_full = data_label_full[:-3]
            #__|

            i1 = df.set_index(var_lst).index
            i2 = df_tmp.set_index(var_lst).index

            df_i = df[i1.isin(i2)]

            data_lst.append({"label": data_label_full, "data": df_i})

        return(data_lst)
        #__|

    def filter_early_revisions(self, dataframe):
        """Remove all entries (rows) which aren't the highest revision number.

        Args:
            dataframe:
        """
        #| - filter_early_revisions
        max_rev = dataframe["revision_number"] == dataframe["max_revision"]
        data_series_maxrev = dataframe[max_rev]

        return(data_series_maxrev)
        #__|

    def view_atoms(self, ind):
        """
        View last image in atoms object in GUI.

        Args:
            ind:
                Index of dataframe corresponding to entry of interest.
        """
        #| - view_atoms
        df = self.data_frame

        path_i = df.iloc[ind]["path"]
        rev_num = df.iloc[ind]["revision_number"].astype(str)
        full_path = path_i + "_" + rev_num

        print(full_path)

        try:
            atoms = df.iloc[ind]["atoms_object"][-1]
            view(atoms)
        except:
            print("Couldn't read atoms object")
        #__|

    #__|

    #| - Query Job Status *****************************************************

    #| - __old__
    def job_state(self, path_i):
        """
        Return job state.

        Implentation is cluster dependent.

        Args:
            path_i
        """
        #| - job_state
        job_state = self.cluster.cluster.job_state(path_i=path_i)

        return(job_state)

        #| - OLD
        # if not os.path.isdir(path + "/simulation"):
        #     return("no_sim_folder")
        #
        # dir_cont = os.listdir(path + "/simulation")
        #
        # if "completed" in dir_cont:
        #     return("complete")
        # elif "out" not in dir_cont and "err" not in dir_cont:
        #     return("running")
        # else:
        #     return("error")
        #__|

        #__|

    def job_state_2(self, path):
        """
        Return job state.

        # COMBAK Deprecated *********************

        Implentation is cluster dependent.

        Args:
            path_i
        """
        #| - job_state_2
        #
        # #| - Formatting path Depending on Whether It is Full or Relative
        # if self.root_dir in path:
        #     ind = path.find(self.root_dir_short)
        #     # path = path[ind:]
        #     full_path = path[ind:]
        # else:
        #     full_path = self.root_dir_short + "/" + path
        # #__|
        #
        # #| - Finding job in jobs.csv file
        # jobs_file_path = self.job_queue_dir + "/jobs.csv"
        # df = pd.read_csv(jobs_file_path)
        #
        # job_in_csv = True
        # try:
        #     index = df[df["job_path"] == full_path].index.tolist()[0]
        # except:
        #     job_in_csv = False
        # #__|
        #
        # try:
        #
        #     #| - Attempting to read job_id from file
        #     with open(path + "/job_id") as fle:
        #         job_id = fle.read().rstrip()
        #     #__|
        #
        # except:
        #
        #     #| - Attempting to read job_id from jobs.csv by matching paths
        #     if job_in_csv:
        #         job_id = df.iloc[index]["job_id"]
        #     #__|
        #
        # job_queue_dict = AWS_Queues().job_info_batch(job_id)
        #
        # #| - Handling Case Where Job ID Is Not In Batch System
        # if job_queue_dict == "job not in batch system":
        #
        #     if job_in_csv:
        #         job_stat_from_file = df.iloc[index]["job_status"]
        #         job_queue_dict = {"job_status": job_stat_from_file}
        #
        #     elif os.path.isfile( path + "/.STATUS"):
        #         with open(stat_file, "w+") as fle:
        #             job_status = fle.readline().rstrip()
        #
        #             job_queue_dict = {"job_status": job_status}
        # #__|
        #
        # job_status = job_queue_dict["job_status"]
        #
        # #| - Writing Job Status to File
        # with open(path + "/.STATUS", "w+") as fle:
        #     fle.write(job_status + "\n")
        # #__|
        #
        # #| - Writing Job Status to Master Jobs Queue File
        # if job_in_csv:
        #     df.at[index, "job_status"] = job_status
        #     df.to_csv(jobs_file_path, index=False)
        # #__|
        #
        # return(job_status)
        #__|

    #__|

    def job_state_3(self, path_i):
        """
        Return job state of job_i.

        Args:
            path_i
        """
        #| - job_state_3
        out_dict = {
            "job_ready": self._job_ready(path_i),
            "job_pending": self._job_pending(path_i),
            "job_running": self._job_running(path_i),
            "job_succeeded": self._job_succeeded(path_i),
            "job_failed": self._job_failed(path_i),
            "job_submitted": self._job_submitted(path_i),
            }

        return(out_dict)
        #__|

    #| - OLD Methods That Use job_i

    def job_ready(self, job_i, require_READY_tag=True):
        """
        Return whether job_i is in READY state (Ready for submission).

        Args:
            job_i:
            require_READY_tag:
                Require a ".READY" file start job
        """
        #| - job_ready
        path_i = self.var_lst_to_path(
            job_i,
            job_rev="Auto",
            relative_path=False,
            )

        crit_0 = False
        if os.path.isfile(path_i + "/.READY"):
            crit_0 = True
        elif require_READY_tag is False:
            crit_0 = True

        crit_1 = False
        if not os.path.isfile(path_i + "/.SUBMITTED"):
            crit_1 = True

        #| - Having trouble with AWS .READY files not being copied over
        # if self.cluster.cluster_sys == "aws":
        #     crit_
        #
        #__|

        crit_list = [crit_0, crit_1]

        if all(crit is True for crit in crit_list):
            return(True)
        else:
            return(False)
        #__|

    def job_pending(self, job_i):
        """
        Check whether job_i is in PENDING state.

        Args:
            job_i:
        """
        #| - job_pending
        path_i = self.var_lst_to_path(job_i,
            job_rev="Auto",
            relative_path=False,
            )

        crit_0 = False
        job_state = self.cluster.cluster.job_state(path_i=path_i)
        if job_state == "PENDING":
            crit_0 = True

        crit_1 = False
        if self.job_state_file(path_i) == "PENDING":
            crit_1 = True

        crit_list = [crit_0, crit_1]

        if all(crit is True for crit in crit_list):
            return(True)
        else:
            return(False)


        # df["full_path"] = root_dir_beg + "/" + df["root_dir"] + "/" +
        # df["path"] + "_" + df["job_revision_number"].astype(str)
        #
        # index = df.index[df["full_path"] == path_i].tolist()
        # df_i = df.iloc[index]
        # max_rev = df_i["revision_number"].max()
        # df_fin = df_i[df_i["revision_number"] == max_rev]
        # assert len(df_fin) == 1
        #
        # if df_fin["job_state_2"].iloc[0] == "RUNNABLE":
        #     return(True)
        # else:
        #     return(False)
        #__|

    def job_running(self, job_i):
        """
        Check whether job_i is in RUNNING state.

        Args:
            job_i:
        """
        #| - job_running
        path_i = self.var_lst_to_path(job_i,
            job_rev="Auto",
            relative_path=False,
            )

        crit_0 = True
        if os.path.isfile(path_i + "/.READY"):
            crit_0 = True

        crit_1 = False
        if os.path.isfile(path_i + "/.SUBMITTED"):
            crit_1 = True

        crit_2 = False
        # if not os.path.isfile(path_i + "/.FINISHED"):
        if not os.path.isfile(path_i + "/" + DFT_Jobs_Analysis.finished_fle):
            crit_2 = True

        crit_3 = False
        job_state = self.cluster.cluster.job_state(path_i=path_i)
        if job_state == "RUNNING":
            crit_3 = True

        crit_list = [crit_0, crit_1, crit_2, crit_3]

        if all(crit is True for crit in crit_list):
            return(True)
        else:
            return(False)
        #__|

    def job_succeeded(self, job_i):
        """
        Check whether job_i is in SUCCEEDED state.

        Args:
            job_i:
        """
        #| - job_succeeded
        path_i = self.var_lst_to_path(job_i,
            job_rev="Auto",
            relative_path=False,
            )

        crit_0 = True
        if os.path.isfile(path_i + "/.READY"):
            crit_0 = True

        crit_1 = False
        if os.path.isfile(path_i + "/.SUBMITTED"):
            crit_1 = True

        # Checking for '.FINSISHED' file OR checking batch queue

        crit_2_1 = False
        fle_name = path_i + "/" + DFT_Jobs_Analysis.finished_fle
        if os.path.isfile(fle_name):
            with open(fle_name, "r") as fle:
                lines = [line.strip() for line in fle.readlines()]
                if "job_completed" in lines:
                    crit_2_1 = True

        #| - DELETE THIS
        # TEMP COMBAK FIXME Delete this after migration to new FINISHED file
        # format is done!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        fle_name = path_i + "/" + DFT_Jobs_Analysis.finished_fle
        if os.path.isfile(fle_name):
            crit_2_1 = True
        #__|


        crit_2_2 = False
        job_state = self.cluster.cluster.job_state(path_i=path_i)
        if job_state == "SUCCEEDED":
            crit_2_2 = True

        if crit_2_2 or crit_2_1:
            crit_2 = True
        else:
            crit_2 = False

        crit_list = [crit_0, crit_1, crit_2]
        if all(crit is True for crit in crit_list):
            return(True)
        else:
            return(False)
        #__|

    def job_failed(self, job_i):
        """
        Check whether job_i is in failed state.

        Args:
            job_i:
        """
        #| - job_failed
        path_i = self.var_lst_to_path(job_i,
            job_rev="Auto",
            relative_path=False,
            )

        crit_0 = False
        job_state = self.cluster.job_state(path_i=path_i)
        if job_state == "FAILED":
            crit_0 = True

        crit_1 = False
        if os.path.isfile(path_i + "/.SUBMITTED"):
            crit_1 = True

        crit_2 = False
        # if not os.path.isfile(path_i + "/.FINISHED"):
        if not os.path.isfile(path_i + "/." + DFT_Jobs_Analysis.finished_fle):
            crit_2 = True

        #| - Parsing Error File for "Error" (Sherlock only for now)
        if self.cluster.cluster_sys == "sherlock":
            error = self.parse_job_error_file(path_i)
            if error:
                crit_0 = True
        #__|

        crit_list = [crit_0, crit_1, crit_2]

        if all(crit is True for crit in crit_list):
            return(True)
        else:
            return(False)
        #__|

    def job_submitted(self, path_i):
        """
        Check whether job is submitted.

        Args:
            path_i
        """
        #| - job_submitted
        try:
            if os.path.isfile(path_i + "/.SUBMITTED"):
                return(True)
            else:
                return(False)
        except:
            return(False)
        #__|

    #__|

    #| - NEW Methods That Use path_i Instead of job_i (Create col in df!!)
    def _job_ready(self, path_i, require_READY_tag=True):
        """
        Return whether job_i is in READY state.

        Args:
            job_i:
            require_READY_tag:
                Require a ".READY" file start job
        """
        #| - job_ready
        # path_i = self.var_lst_to_path(job_i,
        #     ob_rev="Auto",
        #     relative_path=False,
        #     )

        crit_0 = False
        if os.path.isfile(path_i + "/.READY"):
            crit_0 = True
        elif require_READY_tag is False:
            crit_0 = True

        crit_1 = False
        if not os.path.isfile(path_i + "/.SUBMITTED"):
            crit_1 = True

        #| - Having trouble with AWS .READY files not being copied over
        # if self.cluster.cluster_sys == "aws":
        #     crit_
        #__|

        crit_list = [crit_0, crit_1]

        if all(crit is True for crit in crit_list):
            return(True)
        else:
            return(False)
        #__|

    def _job_pending(self, path_i):
        """
        Return whether job_i is in PENDING state.

        Args:
            path_i:
        """
        #| - job_pending
        # path_i = self.var_lst_to_path(job_i, job_rev="Auto",
        # relative_path=False)

        crit_0 = False
        job_state = self.cluster.cluster.job_state(path_i=path_i)
        if job_state == "PENDING":
            crit_0 = True

        crit_1 = False
        if self.job_state_file(path_i) == "PENDING":
            crit_1 = True

        crit_list = [crit_0, crit_1]

        if all(crit is True for crit in crit_list):
            return(True)
        else:
            return(False)
        #__|

    def _job_running(self, path_i):
        """
        Return whether job_i is in RUNNING state.

        Args:
            path_i:
        """
        #| - job_running
        crit_0 = True
        if os.path.isfile(path_i + "/.READY"):
            crit_0 = True

        crit_1 = False
        if os.path.isfile(path_i + "/.SUBMITTED"):
            crit_1 = True

        crit_2 = False
        if not os.path.isfile(path_i + "/.FINISHED"):
            crit_2 = True

        crit_3 = False
        job_state = self.cluster.cluster.job_state(path_i=path_i)
        if job_state == "RUNNING":
            crit_3 = True

        crit_list = [crit_0, crit_1, crit_2, crit_3]

        if all(crit is True for crit in crit_list):
            return(True)
        else:
            return(False)
        #__|

    def _job_succeeded(self, path_i):
        """
        Return whether job_i is in SUCCEEDED state.

        Args:
            path_i:
        """
        #| - job_succeeded
        crit_0 = True
        if os.path.isfile(path_i + "/.READY"):
            crit_0 = True

        crit_1 = False
        if os.path.isfile(path_i + "/.SUBMITTED"):
            crit_1 = True

        # Checking for '.FINSISHED' file OR checking batch queue

        crit_2_1 = False
        if os.path.isfile(path_i + "/.FINISHED"):
            crit_2_1 = True

        crit_2_2 = False
        job_state = self.cluster.cluster.job_state(path_i=path_i)
        if job_state == "SUCCEEDED":
            crit_2_2 = True

        if crit_2_2 or crit_2_1:
            crit_2 = True
        else:
            crit_2 = False

        crit_list = [crit_0, crit_1, crit_2]
        if all(crit is True for crit in crit_list):
            return(True)
        else:
            return(False)
        #__|

    def _job_failed(self, path_i):
        """
        Return whether job_i is in FAILED state.

        Args:
            path_i:
        """
        #| - job_failed
        crit_0 = False
        job_state = self.cluster.job_state(path_i=path_i)
        if job_state == "FAILED":
            crit_0 = True

        crit_1 = False
        if os.path.isfile(path_i + "/.SUBMITTED"):
            crit_1 = True

        crit_2 = False
        if not os.path.isfile(path_i + "/.FINISHED"):
            crit_2 = True


        #| - Parsing Error File for "Error" (Sherlock only for now)

        if self.cluster.cluster_sys == "sherlock":
            error = self.parse_job_error_file(path_i)
            if error:
                crit_0 = True

        # if error and self.cluster.cluster_sys == "sherlock":
        #     # print(error)
        #     crit_0 = True
        #__|

        crit_list = [crit_0, crit_1, crit_2]

        if all(crit is True for crit in crit_list):
            return(True)
        else:
            return(False)
        #__|

    def _job_submitted(self, path_i):
        """
        Return whether job_i is in SUBMITTED state.

        Args:
            path_i:
        """
        #| - job_submitted
        try:
            if os.path.isfile(path_i + "/.SUBMITTED"):
                return(True)
            else:
                return(False)
        except:
            return(False)
        #__|

    #__|

    def parse_job_error_file(self, path_i):
        """
        Read error file and searches for "Error" in last line.

        TODO This is very specific to QE jobs so should be moved.

        Args:
            path_i
        """
        #| - parse_job_error_file
        err_file = self.cluster.cluster.error_file
        err_file = os.path.join(path_i, err_file)

        error = False
        if os.path.isfile(err_file):
            with open(err_file) as fle:
                lines = [line.strip() for line in fle]

            lines = lines[-4:]

            for line in lines:

                if "KohnShamConvergenceError" in line:
                    print("KohnShamConvergenceError Occured!! - RF")
                    error = True
                    break
                elif "DUE TO TIME LIMIT" in line:
                    print("Job Reached Time Limit!! - RF")
                    error = True
                    break
                elif "RuntimeError: SCF Calculation failed":
                    print("RuntimeError: SCF Calculation failed - RF")
                    error = True
                    break
                else:
                    error = False
                    pass

        else:
            error = False

        return(error)
        #__|

    #__| **********************************************************************

#__| **************************************************************************
