"""Class to analyse data using the DFT_Jobs_Setup class."""

#| - Import Modules
import os
# import pickle
import cPickle as pickle
import subprocess
import copy

from ase import io
import pandas as pd
import numpy as np

# My Modules

# from aws.aws_class import AWS_Queues
from dft_job_automat.job_setup import DFT_Jobs_Setup
from dft_job_automat.job_types_classes.dft_methods import DFT_Methods

#__|

class DFT_Jobs_Analysis(DFT_Jobs_Setup):
    """Summary line.
    TEMP
    """

    def __init__(self,
        system="sherlock",
        tree_level=None,
        level_entries=None,
        skip_dirs_lst=None,
        working_dir=".",
        update_job_state=True,
        # update_job_state=False,  # TEMP
        load_dataframe=True,
        job_type_class=None,
        ):
        """TMP_docstring.
        TEMP TEMP

        Args:
            system:
            tree_level:
            level_entries:
            update_job_state:
                Updates job status for all jobs in ensemble
        """
        #| - __init__
        DFT_Jobs_Setup.__init__(self,
            system=system,
            tree_level=tree_level,
            level_entries=level_entries,
            skip_dirs_lst=None,
            working_dir=working_dir,
            )

        #| - TEMP
        # if load_dataframe == True:
        #     try:
        #         self.__load_dataframe__()
        #
        #         if update_job_state == True:
        #             self.add_data_column(self.job_state_2, column_name="job_state_2")
        #         else:
        #             pass
        #         print("")
        #
        #     except:
        #         self.add_data_column(self.job_revisions, column_name="job_revision_number")
        #         self.add_data_column(self.job_state, column_name="job_state")
        #         # self.add_data_column(self.elec_energy, column_name="elec_energy")
        #         self.add_data_column(self.max_force, column_name="max_force")
        #         self.add_data_column(self.job_submitted, column_name="job_submitted")
        #         self.add_data_column(self.job_queue_info, column_name="N/A")
        #
        #         self.add_data_column(self.job_state_2, column_name="job_state_2")
        #
        #         self.__write_dataframe__()
        #         print("")
        #
        # else:
        #     self.add_data_column(self.job_revisions, column_name="job_revision_number")
        #     self.add_data_column(self.job_state, column_name="job_state")
        #     # self.add_data_column(self.elec_energy, column_name="elec_energy")
        #     self.add_data_column(self.max_force, column_name="max_force")
        #     self.add_data_column(self.job_submitted, column_name="job_submitted")
        #     self.add_data_column(self.job_state_2, column_name="job_state_2")
        #     self.add_data_column(self.job_queue_info, column_name="N/A")
        #
        #     self.__write_dataframe__()
        #     print("")
        #__|


        #| - General Methods
        if update_job_state == True:
            print("update_job_state == True")
            self.add_data_column(self.job_state, column_name="job_state")
            self.add_data_column(self.job_state_3, column_name="N/A", allow_failure=False)

        else:
            pass
        #__|

        #| - Job Type Specific Methods

        # method = DFT_Methods().atom_type_num_dict
        # self.add_data_column(method, column_name="TEMP", allow_failure=False)


        if load_dataframe == True:
            self.__load_dataframe__()

        else:
            if job_type_class != None:
                job_type_inst = job_type_class

                for method in job_type_inst.methods_to_run:
                    method_ref = getattr(job_type_inst, method)
                    self.method = method_ref
                    self.add_data_column(method_ref, column_name=method, allow_failure=True)
                    # self.add_data_column(method_ref, column_name=method, allow_failure=False)

                # TEMP
                self.__write_dataframe__()
        #__|

        #__|

    def add_jobs_queue_data(self):
        """
        """
        #| - add_jobs_queue_data


        #__|

    def job_queue_info(self, path):
        """
        """
        #| - job_queue_info
        jobs_file_path = self.job_queue_dir + "/jobs.csv"
        df = pd.read_csv(jobs_file_path)

        path = path[len(self.root_dir):]
        full_path = self.root_dir_short + path

        # Looking for job path in jobs.csv file
        index = df[df["job_path"] == full_path].index.tolist()[0]

        job_info = df.iloc[index].to_dict()


        return(job_info)
        #__|

    def __load_dataframe__(self):
        """Attempts to loead
        """
        #| - __load_dataframe__
        # return(donkle)
        df = pd.read_csv(self.root_dir + "/jobs_bin/job_dataframe.csv")

        with open(self.root_dir + "/jobs_bin/job_dataframe.pickle", "rb") as  fle:
            df = pickle.load(fle)

        # df = pickle.load(open(self.root_dir + "/jobs_bin/job_dataframe.pickle", "rb"))
        # df = pickle.load(open(self.root_dir + "/jobs_bin/job_dataframe.pickle", "rb"))

        self.data_frame = df

        print("loaded dataframe")
        #__|

    def __write_dataframe__(self):
        """
        """
        #| - __write_dataframe__
        print("writing dataframe")
        df = self.data_frame
        df.to_csv(self.root_dir + "/jobs_bin/job_dataframe.csv", index=False)

        with open(self.root_dir + "/jobs_bin/job_dataframe.pickle", "wb") as fle:
            pickle.dump(df, fle)

        # pickle.dump(df, open(self.root_dir + "/jobs_bin/job_dataframe.pickle", "wb"))
        #__|

    def add_data_column(self,
        function,
        column_name="new_column",
        revision="auto",
        allow_failure=True):
        """

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

            path = self.var_lst_to_path(entry, relative_path=False, job_rev="False")

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

                if allow_failure == True:
                    try:
                        out = function(path)
                    except:
                        out = np.nan
                else:
                    out = function(path)

                # print("!@)#(!@*#(@!DFJ))")
                # print(out)
                new_data_col.append(out)
                job_rev_lst.append(rev_num)
                #__|

        data_type_list = [type(x) for x in new_data_col]
        dict_in_list = any(item == dict for item in data_type_list)

        if dict_in_list:
            new_col = []
            for x in new_data_col:
                if pd.isnull(x) == True:
                    new_col.append({"NA": np.nan})
                else:
                    new_col.append(x)

            new_columns_df = pd.DataFrame(new_col)
            # print(new_columns_df)

            df1 = self.data_frame
            # print(list(self.data_frame))
            df2 = new_columns_df

            out_df = pd.concat([df1, df2], axis=1)

            self.data_frame = out_df

            # print(list(self.data_frame))

            # out_df.to_csv("combined_df.csv")  # TEMP


        else:
            self.data_frame[column_name] = new_data_col


        #| - OLD
        # if type(new_data_col[0]) != dict:
        #     self.data_frame[column_name] = new_data_col
        #     # print("&&&&&&&&&&&&&&&")
        # else:
        #     # print("@@@@@@@@@@@@@@")
        #     # print("list of dicts")
        #     # print(type(new_data_col))
        #     print("")
        #     # print(new_data_col)
        #
        #     # print(new_data_col)
        #
        #     new_col = []
        #     for x in new_data_col:
        #         if pd.isnull(x) == True:
        #             new_col.append({"NA": np.nan})
        #             # new_col.append({"NA": "#############"})
        #         # if x == np.nan:
        #         else:
        #             new_col.append(x)
        #
        #     new_columns_df = pd.DataFrame(new_col)
        #     # new_data_col = [{"NA": np.nan} for x in new_data_col if x == np.nan]
        #
        #     # print(new_col)
        #     # return(new_col)  # TEMP
        #__|


        #__|

    def job_revisions(self, path):
        """
        """
        #| - job_revisions
        path = "/".join(path.split("/")[0:-1]) + "/"

        # Attempting to remove duplicate job folders (usually have spaces)
        dir_list = [x for x in os.walk(path).next()[1] if " " not in x]

        return(len(dir_list))
        #__|

    def job_state_file(self, path="."):
        """Returns contents of '.QUEUESTATE' if present in the job directory
        """
        #| - job_state_file
        file_path = path + "/.QUEUESTATE"
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
        Returns data sets from the data frame where the data is split by the
        job variables. One variable is excluded from this grouping and is
        treated as a "continous" variable to be plotted on the x-axis

        ex. Parameter sweep with varying k-points, pw-cutoff, and
        lattice-constant. The lattice-constant is designated as the free
        variable and so data sets are created for each combination of k-points
        and pw-cutoff where each set contains the full range of latt-const.
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

            data_lst.append({"label":data_label_full, "data":df_i})

        return(data_lst)
        #__|


    def filter_early_revisions(self, dataframe):
        """
        Removes all entries (rows) which aren't the highest revision number.
        """
        #| - filter_early_revisions
        max_rev = dataframe["revision_number"] == dataframe["max_revision"]
        data_series_maxrev = dataframe[max_rev]

        return(data_series_maxrev)
        #__|

#__|

#| - Query Job Status *********************************************************

    #| - __old__
    def job_state(self, path_i):
        """
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
        """
        #| - job_state_2

        #| - Formatting path Depending on Whether It is Full or Relative
        if self.root_dir in path:
            ind = path.find(self.root_dir_short)
            # path = path[ind:]
            full_path = path[ind:]
        else:
            full_path = self.root_dir_short + "/" + path
        #__|


        #| - Finding job in jobs.csv file
        jobs_file_path = self.job_queue_dir + "/jobs.csv"
        df = pd.read_csv(jobs_file_path)

        job_in_csv = True
        try:
            index = df[df["job_path"] == full_path].index.tolist()[0]
        except:
            job_in_csv = False
        #__|

        try:

            #| - Attempting to read job_id from file
            with open(path + "/job_id") as fle:
                job_id = fle.read().rstrip()
            #__|

        except:

            #| - Attempting to read job_id from jobs.csv by matching paths
            if job_in_csv:
                job_id = df.iloc[index]["job_id"]
            #__|

        job_queue_dict = AWS_Queues().job_info_batch(job_id)

        #| - Handling Case Where Job ID Is Not In Batch System
        if job_queue_dict == "job not in batch system":

            if job_in_csv:
                job_stat_from_file = df.iloc[index]["job_status"]
                job_queue_dict = {"job_status": job_stat_from_file}

            elif os.path.isfile( path + "/.STATUS"):
                with open(stat_file, "w+") as fle:
                    job_status = fle.readline().rstrip()

                    job_queue_dict = {"job_status": job_status}
        #__|

        job_status = job_queue_dict["job_status"]

        #| - Writing Job Status to File
        with open(path + "/.STATUS", "w+") as fle:
            fle.write(job_status + "\n")
        #__|

        #| - Writing Job Status to Master Jobs Queue File
        if job_in_csv:
            df.at[index, "job_status"] = job_status
            df.to_csv(jobs_file_path, index=False)
        #__|

        return(job_status)
        #__|

    #__|

    def job_state_3(self, path_i):
        """
        """
        #| - job_state_3
        # print("job_state_3")
        out_dict = {
            "job_ready": self._job_ready(path_i),
            "job_pending": self._job_pending(path_i),
            "job_running": self._job_running(path_i),
            "job_succeeded": self._job_succeeded(path_i),
            "job_failed": self._job_failed(path_i),
            "job_submitted": self._job_submitted(path_i),
            }

        # print(out_dict)
        return(out_dict)
        #__|

    #| - OLD Methods That Use job_i

    def job_ready(self, job_i, require_READY_tag=True):
        """

        Args:
            job_i:
            require_READY_tag:
                Require a ".READY" file start job
        """
        #| - job_ready
        path_i = self.var_lst_to_path(job_i, job_rev="Auto", relative_path=False)

        crit_0 = False
        if os.path.isfile(path_i + "/.READY"):
            crit_0 = True
        elif require_READY_tag == False:
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

        # print(crit_list)
        if all(crit == True for crit in crit_list):
            return(True)
        else:
            return(False)
        #__|

    def job_pending(self, job_i):
        """
        """
        #| - job_pending
        path_i = self.var_lst_to_path(job_i, job_rev="Auto", relative_path=False)

        crit_0 = False
        job_state = self.cluster.cluster.job_state(path_i=path_i)
        if job_state == "PENDING":
            crit_0 = True

        crit_1 = False
        if self.job_state_file(path_i) == "PENDING":
            crit_1 = True
            # print("Sucess!!! - TEMP")

        # print("JOB_PENDING")
        # print(crit_0)
        # print(crit_1)
        # print("____")

        crit_list = [crit_0, crit_1]

        if all(crit == True for crit in crit_list):
            return(True)
        else:
            return(False)


        # df["full_path"] = root_dir_beg + "/" + df["root_dir"] + "/" + df["path"] + "_" + df["job_revision_number"].astype(str)
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
        """
        #| - job_running
        path_i = self.var_lst_to_path(job_i, job_rev="Auto", relative_path=False)

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

        # print("$&^**")
        crit_list = [crit_0, crit_1, crit_2, crit_3]
        # print(crit_list)

        if all(crit == True for crit in crit_list):
            return(True)
        else:
            return(False)
        #__|

    def job_succeeded(self, job_i):
        """
        """
        #| - job_succeeded
        path_i = self.var_lst_to_path(job_i, job_rev="Auto", relative_path=False)

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
        if all(crit == True for crit in crit_list):
            return(True)
        else:
            return(False)
        #__|

    def job_failed(self, job_i):
        """
        """
        #| - job_failed
        path_i = self.var_lst_to_path(job_i, job_rev="Auto", relative_path=False)

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

        # print(crit_list)
        if all(crit == True for crit in crit_list):
            return(True)
        else:
            return(False)
        #__|

    def job_submitted(self, path):
        """
        """
        #| - job_submitted
        try:
            if os.path.isfile(path + "/.SUBMITTED"):
                return(True)
            else:
                return(False)
        except:
            return(False)
        #__|

    #__|

    #| - NEW Methods That Use path_i Instead of job_i (Can create col in df!!)
    def _job_ready(self, path_i, require_READY_tag=True):
        """

        Args:
            job_i:
            require_READY_tag:
                Require a ".READY" file start job
        """
        #| - job_ready
        # path_i = self.var_lst_to_path(job_i, job_rev="Auto", relative_path=False)

        crit_0 = False
        if os.path.isfile(path_i + "/.READY"):
            crit_0 = True
        elif require_READY_tag == False:
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

        # print(crit_list)
        if all(crit == True for crit in crit_list):
            return(True)
        else:
            return(False)
        #__|

    def _job_pending(self, path_i):
        """
        """
        #| - job_pending
        # path_i = self.var_lst_to_path(job_i, job_rev="Auto", relative_path=False)

        crit_0 = False
        job_state = self.cluster.cluster.job_state(path_i=path_i)
        if job_state == "PENDING":
            crit_0 = True

        crit_1 = False
        if self.job_state_file(path_i) == "PENDING":
            crit_1 = True
            # print("Sucess!!! - TEMP")

        # print("JOB_PENDING")
        # print(crit_0)
        # print(crit_1)
        # print("____")

        crit_list = [crit_0, crit_1]

        if all(crit == True for crit in crit_list):
            return(True)
        else:
            return(False)


        # df["full_path"] = root_dir_beg + "/" + df["root_dir"] + "/" + df["path"] + "_" + df["job_revision_number"].astype(str)
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

    def _job_running(self, path_i):
        """
        """
        #| - job_running
        # path_i = self.var_lst_to_path(job_i, job_rev="Auto", relative_path=False)

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

        # print("$&^**")
        crit_list = [crit_0, crit_1, crit_2, crit_3]
        # print(crit_list)

        if all(crit == True for crit in crit_list):
            return(True)
        else:
            return(False)
        #__|

    def _job_succeeded(self, path_i):
        """
        """
        #| - job_succeeded
        # path_i = self.var_lst_to_path(job_i, job_rev="Auto", relative_path=False)

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
        if all(crit == True for crit in crit_list):
            return(True)
        else:
            return(False)
        #__|

    def _job_failed(self, path_i):
        """
        """
        #| - job_failed
        # path_i = self.var_lst_to_path(job_i, job_rev="Auto", relative_path=False)

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

        # print(crit_list)
        if all(crit == True for crit in crit_list):
            return(True)
        else:
            return(False)
        #__|

    def _job_submitted(self, path_i):
        """
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

###############################################################################

    def parse_job_error_file(self, path_i):
        """
        Reads error file and searches for "Error" in last line
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

                else:
                    error = False
                    pass

        else:
            error = False

        # print(error)
        return(error)
        #__|

#__| ***************************************************************************
















#| - DEPRECATED METHODS
# def max_force(self, path):
#     """
#     """
#     #| - max_force
#     with open(path + "/simulation/qn.log", "r") as file:
#         max_f = file.readlines()[-1].split(" ")[-1]
#         return(float(max_f))
#     #__|

# def job_id_names(self):
#     """
#     """
#     #| - job_id_names
#     try:
#         f = open("job_id_names", "r")
#         job_id_data = pickle.load(f)
#         f.close()
#
#         id_list = job_id_data[0]
#         name_list = job_id_data[1]
#
#         self.data_frame["job_id"] = id_list
#         self.data_frame["job_name"] = name_list
#     except:
#         print("Couldn't parse job_id_names file")
#
#     #__|

#__|
