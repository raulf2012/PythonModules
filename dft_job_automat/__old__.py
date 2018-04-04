#| - job_analysis

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

#__|

#| - job_setup

#| - __old__

#| - TEMP
# def __parse_cluster_type__(self):
#     """
#     """
#     #| - __parse_cluster_type__
#     clusters_dict = {
#         "aws": "AWSCluster",
#         "slac": "SLACCluster",
#         "sherlock": "SherlockCluster",
#         }
#
#     cluster_sys = os.environ.get("COMPENV")
#
#     if cluster_sys in clusters_dict:
#         package = "dft_job_automat.compute_env"
#         name = clusters_dict[cluster_sys]
#         cluster = getattr(__import__(package, fromlist=[name]), name)
#
#         self.cluster_sys = cluster_sys
#         self.cluster = cluster
#         # from dft_job_automat.compute_env import clusters_dict[cluster_sys]
#    return(cluster_sys)
#__|


#__|

#__|

#| - job_manager

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

#| - __old__

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

#__|
