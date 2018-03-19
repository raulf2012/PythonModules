"""Module to carry out common file operations in jobs directory"""

#| - Import Modules
# from dft_job_automat.job_setup import DFT_Jobs_Setup
from dft_job_automat.job_analysis import DFT_Jobs_Analysis

from aws.aws_class import AWS_Queues

import os
# import sys
import shutil
# import subprocess
# import pickle
# import boto3
# import datetime
import pandas as pd
import filecmp
#__|


class DFT_Jobs_Manager(DFT_Jobs_Analysis):
    """Summary line.

    Manages job submission, resubmission, revision managmnet, job status, etc.
    """
    def __init__(self,
        system="sherlock",
        tree_level=None,
        level_entries=None,
        skip_dirs_lst=None,
        working_dir=".",
        update_job_state=True,
        load_dataframe=True,
        ):
        """TMP_docstring.

        TEMP TEMP
        """
        #| - __init__
        DFT_Jobs_Analysis.__init__(self,
            system=system,
            tree_level=tree_level,
            level_entries=level_entries,
            working_dir=working_dir,
            update_job_state=update_job_state,
            load_dataframe=load_dataframe,
            )
        #__|


    def restart_job(self,
        job_i,
        prev_rev_files_list,
        root_dir=".",
        file_list=None,
        sub_params=None,
        source_rev=None,
        from_simulation_folder=True,
        run_job=False,
        ):
        """
        # TODO | Make copy_if_not_in_dest a global function or something

        Args:
            prev_rev_files_list:
            file_list:
            job_i
        """
        #| - restart_job
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

        self.create_job_dir(job_i, revision="Auto")
        self.copy_files_from_last_revision(
            prev_rev_files_list,
            job_i,
            revisions="Auto",
            source_rev=source_rev,
            from_simulation_folder=from_simulation_folder,
            )

        if file_list is not None:
            source_dir = root_dir

            dest_dir = self.var_lst_to_path(
                job_i,
                job_rev="Auto",
                relative_path=False,
                )

            for file_i in file_list:
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

        path_i = self.var_lst_to_path(
            job_i,
            job_rev="Auto",
            relative_path=False,
            )

        if sub_params is None:
            params_dict = {
                "path_i": path_i
                }
        else:
            params_dict = sub_params
            params_dict["path_i"] = path_i

        if run_job:
            self.submit_job(**params_dict)
        else:
            pass
        #__|


    def copy_files_from_last_revision(self,
        files_list,
        job_i,
        revisions="Auto",
        source_rev=None,
        from_simulation_folder=True,
        ):
        """
        Args:
            revisions: <type 'str' or list>
                [source revision, destination revision]
            files_list:
            job_i:
            revisions:
            source_rev:
            from_simulation_folder:
                If system is AWS, decidedes whether the previous jobs' files
                are obtained from the /simulation folder or the origial job dir.
        """
        #| - copy_files_from_last_revision

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

        job_path = self.var_lst_to_path(
            job_i,
            job_rev=False,
            relative_path=False,
            )
        rev_num = self.job_revision_number(job_i)

        if revisions == "Auto":
            rev_dest = rev_num

            if source_rev is not None:
                rev_source = source_rev
            else:
                rev_source = rev_dest - 1

        else:
            rev_dest = revisions[1]
            rev_source = revisions[0]

        dest_dir = job_path + "_" + str(rev_dest)
        source_dir = job_path + "_" + str(rev_source)

        for file_i in files_list:

            #| - Copy Files to Directory in New Step

            if from_simulation_folder is True:
                data_d = self.cluster.cluster.job_data_dir
            else:
                data_d = ""

            if type(file_i) == str:
                copy_if_not_in_dest(
                    source_dir + data_d + "/" + file_i,
                    dest_dir + "/" + file_i,
                    )

            elif type(file_i) == list:
                copy_if_not_in_dest(
                    source_dir + data_d + "/" + file_i[0],
                    dest_dir + "/" + file_i[1],
                    )
            #__|

        #__|

    def submit_job(self, **kwargs):
        """
        """
        #| - submit_job
        # path_i = kwargs["path_i"]

        self.cluster.submit_job(**kwargs)
        #__|

    def remove_rev_folder(self, revision_number):
        """Removes revision job folder in all jobs directories
        TEMP
        """
        #| - remove_rev_folder
        print("Removing job revision folder " + str(revision_number))

        for job in self.job_var_lst:
            path = self.var_lst_to_path(job)

            shutil.rmtree(path + "_" + str(revision_number))

        #__|

    def restart_job_2(self, prev_rev_file_list=[], root_dir_file_list=[]):
        """
        """
        #| - restart_job_2
        # stop = False
        for index, row in self.data_frame.iterrows():

            # if row["job_state"] == "complete" or
            # row["job_state"] == "running":
            #     continue

            # if row["job_state"] == "error" or
            # row["job_state"] == "no_sim_folder":
            if row["job_state"] == "error":
                #| - TMP
                job = row["variable_list"]
                path = row["path"]
                rev_n = self.job_revision_number(job)

                first_path = path + "_1"
                prev_path = path + "_" + str(rev_n)
                new_path = path + "_" + str(rev_n + 1)

                # if stop == True:
                #     break
                # stop = True

                self.create_job_dir(job, revision="Auto")

                for file in prev_rev_file_list:

                    if "out.traj" in file:
                        dest_path = new_path + "/init.traj"
                    else:
                        dest_path = new_path + "/"

                    try:
                        shutil.copy(prev_path + "/" + file, dest_path)
                    except:
                        shutil.copy(first_path + "/" + file, dest_path)


                for file in root_dir_file_list:
                    file_path = self.root_dir + "/" + file
                    shutil.copy(file_path, new_path + "/")

                os.system("chmod 777 " + new_path + "/*")

                os.chdir(new_path)
                print("Submitting Folder: " + new_path)

                bash_command = "$HOME/matr.io/bin/trisub -q medium -c 4"
                os.system(bash_command)

                os.chdir(self.root_dir)
                #__|
            else:
                continue
        #__|

    #| - TEMP - Old Restart Job Method
    # def restart_job(self, prev_rev_file_list=[], root_dir_file_list=[],
    # revision="Auto"):
    #     """Restart job from previous run
    #
    #
    #         variable_lst: <type 'list'>
    #             Produced from DFT_Jobs_Setup.job_var_lst.
    #
    #     Args:
    #         model_file: <type 'str'>
    #             "Auto" |
    #         atoms_file: <type 'str'>
    #             Name of atoms object file from previous run
    #         revision: <type 'str' or 'int'>
    #             Job revision number from which the job will be restarted from.
    #             "Auto" | Restarts the job from the most recent job revision.
    #     """
    #     #| - restart_job
    #     for job in self.job_var_lst:
    #         path = self.var_lst_to_path(job)
    #
    #         rev_n = self.job_revision_number(job)
    #
    #         if revision == "Auto":
    #             prev_path = path + "_" + str(rev_n)
    #         else:
    #             prev_path = path + "_" + str(revision)
    #
    #         new_path = path + "_" + str(rev_n + 1)
    #
    #         self.create_job_dir(job, revision="Auto")
    #
    #         for file in prev_rev_file_list:
    #             # print(file)
    #             if "out.traj" in file:
    #                 dest_path = new_path + "/init.traj"
    #             else:
    #                 dest_path = new_path + "/"
    #
    #             shutil.copy(prev_path + "/" + file, dest_path)
    #
    #
    #         for file in root_dir_file_list:
    #             file_path = self.root_dir + "/" + file
    #             shutil.copy(file_path, new_path + "/")
    #
    #         os.system("chmod 777 " + new_path + "/*")
    #
    #     #__|
    #__|

    def copy_files_jd(self, file_list, variable_lst, revision="Auto"):
        """Copy files to job directory
        TEMP
        """
        #| - copy_files_jd
        path = self.var_lst_to_path(variable_lst)
        path += "_" + str(self.job_revision_number(variable_lst))

        for file in file_list:
            shutil.copyfile(self.root_dir + "/" + file, path + "/" + file)

        #__|

    def create_job_dir(self, variable_lst, revision="Auto"):
        """Creates indiviudal job folders as leaves within the dir tree

        TEMP
        """
        #| - create_job_dir
        path = self.var_lst_to_path(
            variable_lst,
            job_rev="False",
            relative_path=False,
            )

        # print(path)
        if revision == "Auto":
            rev = self.job_revision_number(variable_lst) + 1
            path += "_" + str(rev)

            if not os.path.exists(path):
                print("Creating revision folder " + str(rev))  # PRINT
                os.makedirs(path)

            else:
                path += "_" + str(revision)

                if not os.path.exists(path):
                    os.makedirs(path)

        #__|

    #                     path variable not needed/used DELETE #NOTE
    def submit_jobs(self, path=None, queue="medium", copy_PythonModules=True):
        """Submit all jobs within data folder

        Submits the most recent version number in each job folder

        1. Copies the current model.py into the "model_version" folder
        2. Goes through
        """
        #| - submit_jobs


        #| - OLD
        # def submit_folder_job(folder_dir, bash_command):
        #     """
        #     """
        #     #| - submit_folder_job
        #
        #
        #     #| - Submtting Job
        #     try:
        #         os.chdir(folder_dir)
        #
        #         if os.path.isfile(".SUBMITTED"):
        #             print("Directory already submitted, will be skipped")
        #             os.chdir(self.root_dir)
        #             return(None)
        #         else:
        #             print("tmp - submitting job")
        #             output = subprocess.check_output(bash_command, shell=True)
        #             # TEMP
        #             sub_time = datetime.datetime.now().isoformat()
        #
        #     except:
        #         os.chdir(self.root_dir)
        #         print("JOB SKIPPED: " + folder_dir)
        #         return(None)
        #     #__|
        #
        #
        #     #| - Parsing Submission for Job ID
        #     output = output.splitlines()
        #     for line in output:
        #         if "jobId" in line:
        #             lst = line.split('"')
        #             job_id_ind = (lst.index("jobId") + 2)
        #             jobId = lst[job_id_ind]
        #
        #         # This is not necessary anymore
        #         elif "jobName" in line:
        #             lst = line.split('"')
        #             job_name_ind = (lst.index("jobName") + 2)
        #             jobName = lst[job_name_ind]
        #     #__|
        #
        #
        #     file = open(".submitted", "w")
        #     file.close()
        #
        #     os.chdir(self.root_dir)
        #
        #     #| - Querying AWS For Job Info
        #     job_queue_dict = AWS_Queues(jobId).job_info_batch()
        #     job_queue_dict["submit_time"] = sub_time
        #
        #     jobs_file_path = self.job_queue_dir + "/jobs.csv"
        #
        #     df_new = pd.DataFrame([job_queue_dict])
        #     if os.path.isfile(jobs_file_path):
        #         df = pd.read_csv(jobs_file_path)
        #         df = df.append(df_new)
        #     else:
        #         df = df_new
        #
        #     df.to_csv(jobs_file_path, index=False)
        #     #__|
        #
        #     return job_queue_dict
        #
        #     #__|
        #__|


        #| - Create Models Directoy
        models_dir = "model_versions"
        if not os.path.exists(models_dir):
            os.makedirs(models_dir)

        def num_of_files(dir):
            num_files = len(os.listdir("model_versions"))
            return(num_files)

        lead_cnt = num_of_files("model_versions")
        rev_file_path_prev = "model_versions/" + \
            str((lead_cnt)).zfill(2) + "_model.py"
        rev_file_path = "model_versions/" + \
            str((lead_cnt + 1)).zfill(2) + "_model.py"

        if num_of_files("model_versions") == 0:
            shutil.copyfile("model.py", rev_file_path)

        elif not filecmp.cmp("model.py", rev_file_path_prev):
            shutil.copyfile("model.py", rev_file_path)
        #__|

        for job in self.job_var_lst:
            rev_num = self.job_revision_number(job)
            path_i = self.var_lst_to_path(job) + "_" + str(rev_num)
            # job_queue_dict = submit_folder_job(path_i, bash_command)
            # job_queue_dict = AWS_Queues().submit_job(path=path_i, queue=queue)
            AWS_Queues().submit_job(path=path_i, queue=queue)

        #__|

    def cancel_jobs(self, state="RUNNABLE"):
        """
        """
        #| - cancel_jobs


        # for job in self.job_var_lst:

        # full_path = self.root_dir_short + "/" + path

        jobs_file_path = self.job_queue_dir + "/jobs.csv"
        df = pd.read_csv(jobs_file_path)

        # df.loc[df["job_path"] == ]
        df_proj = df[df["job_path"].str.contains(self.root_dir_short)]

        return(df_proj)

        # index = df[df["job_path"] == full_path].index.tolist()[0]
        # queue = df.iloc[index]["job_queue"]
        # job_id = df.iloc[index]["job_id"]

        # job_queue_dict = AWS_Queues(job_id).job_info_batch()
        # job_status = job_queue_dict["job_status"]

        # df.at[index, "job_status"] = job_status
        # jobs_file_path = self.job_queue_dir + "/jobs.csv"
        # df.to_csv(jobs_file_path, index=False)

        # return(job_status)

        #__|

    def update_jobs_queue_file(self):
        """
        """
        #| - update_jobs_queue_file
        tmp = 42
        print(tmp)

        #__|



    #| - OUT OF VIEW | TEMP

    # def job_revision_number(self, variable_lst):
        """Returns the largest revision number for the job with the given
        variable list

        Args:
            variable_lst:
        """
        #| - job_revision_number
        # # path = "data/" + self.var_lst_to_path(variable_lst)
        # path = self.var_lst_to_path(variable_lst)
        # bash_comm = ""
        # os.chdir(path)
        #
        # dirs = filter(os.path.isdir, os.listdir(os.getcwd()))
        #
        # if self.system == "sherlock":
        #     num_jobs = len([dir for dir in dirs if "_jd" in dir])
        # elif self.system == "aws":
        #     num_jobs = len([dir for dir in dirs if dir[0] == "_" and
        #     dir[1].isdigit()])
        #
        # # print(job_dirs)
        # # os.system(bash_comm)
        # os.chdir(self.root_dir)
        #
        # return(num_jobs)
        #__|

        #__|
