#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Class for computer cluster operations, mainly batch.

Development Notes:
    TODO Modify .FINISHED implementation

    TODO Move the location of following command to base ComputerCluster class
        with open(".SUBMITTED", "w") as fle:
            fle.write("")

    TODO I don't think the child class attributes for error and out file are
    being respected by ComputerCluster

"""

#| - Import Modules
import os
import sys
import subprocess
import datetime
import re
import shutil
import copy
import time

import pandas as pd

import json

# My Modules
from misc_modules.misc_methods import merge_two_dicts
#__|

#| - Methods

def slurm_squeue_parse(
    job_id,
    path_i=None,
    queue_state_key="STAT",
    job_state_dict={
        "PD": "PENDING",
        "R": "RUNNING",
        "CF": "CONFIGURING",
        "SUCCEEDED": "SUCCEEDED",
        }

    ):
    """Parse slurm squeue command for job state.

    Args:
        job_id:
        path_i:
        queue_state_key:
    """
    #| - slurm_squeue_parse
    bash_comm = "squeue -j " + str(job_id)

    try:
        out = subprocess.check_output(
            bash_comm,
            shell=True,
            stderr=subprocess.STDOUT,
            )

        # 'JOBID PARTITION NAME USER ST TIME NODES NODELIST(REASON)'
        out = out.splitlines()
        out = out[1].split(" ")
        out = [i for i in out if i != '']

        data_dict = {
            "PARTITION": out[1],
            "STAT": out[4],
            # "CF":
            }

        if path_i is not None:
            with open(path_i + "/.QUEUESTATE", "w") as fle:
                fle.write(job_state_dict[data_dict[queue_state_key]])
                fle.write("\n")

    except subprocess.CalledProcessError:
        data_dict = None
        pass

    except:
        data_dict = None
        pass

    return(data_dict)
    #__|

#__|

################################################################################
class ComputerCluster():
    """Base class for interfacing with computing resources.

    Development Notes:
        TODO Create dummy cluster for WSL system
    """

    #| - ComputerCluster ******************************************************

    def __init__(self,
        ):
        """Initialize the base ComputerCluster class."""
        #| - __init__
        self.root_dir = os.getcwd()
        self.default_sub_params = self.default_submission_parameters()

        # self.username = self.__parse_username__()

        self.__parse_cluster_type__()
        #__|

    #| - __old__
    # def __parse_username__(self):
    #     """
    #     """
    #     #| - __parse_username__
    #     username = os.environ.get("USER")
    #
    #     cond_1 = False
    #     if type(username) == str:
    #         cond_1 = True
    #
    #     cond_2 = False
    #     if username is not None:
    #         cond_2 = True
    #
    #     print(cond_1)
    #     print(cond_2)
    #
    #     print(username)
    #     print("*******")
    #     return(username)
    #     #__|
    #__|

    def __parse_cluster_type__(self):
        """Parse for the current cluster system."""
        #| - __parse_cluster_type__
        clusters_dict = {
            "aws": "AWSCluster",
            "slac": "SLACCluster",
            "sherlock": "SherlockCluster",
            "nersc": "EdisonCluster",
            }

        cluster_sys = os.environ.get("COMPENV")

        if cluster_sys in clusters_dict:
            package = "dft_job_automat.compute_env"
            name = clusters_dict[cluster_sys]
            cluster = getattr(__import__(package, fromlist=[name]), name)

            self.cluster_sys = cluster_sys

            self.cluster = cluster(
                root_dir=self.root_dir,
                )
        else:
            pass
            self.cluster_sys = "dummy"
            # TODO Create instance to dummy cluster class

            self.cluster = DummyCluster(
                root_dir=self.root_dir,
                )

        home = os.environ.get("HOME")
        try:
            with open(home + "/.sys_files_rf/jobs_list_dir", "r") as fle:
                jobs_dir = fle.read().rstrip()
        except:
            jobs_dir = None

        self.jobs_list_dir = jobs_dir
        #__|

    def default_submission_parameters(self):
        """Global submission parameters for cluster jobs."""
        #| - default_submission_parameters
        def_params = {
            "path_i": ".",
            "job_name": "Default",
            "job_script": "model.py",
            "out_file": "job.out",
            "err_file": "job.err",
            }

        return(def_params)
        #__|

    def is_job_submitted(self, path_i="."):
        """Check if job has been submitted.

        Return TRUE if the job at 'path' has been submtted
        A job is considered submitted if there exists a '.submitted' file

        Args:
            path
        """
        #| - add_jobs_queue_data
        root_dir = os.getcwd()
        os.chdir(path_i)
        if os.path.isfile(".SUBMITTED"):
            print("Directory already submitted, will be skipped")
            os.chdir(root_dir)
            submitted = True
        else:
            os.chdir(root_dir)
            submitted = False

        return(submitted)
        #__|

    def job_state(self, path_i="."):
        """

        Args:
            path_i:
        """
        #| - job_state
        job_state = self.cluster.job_state(path_i=path_i)

        if job_state is None:
            try:
                with open(path_i + "/.QUEUESTATE", "r") as fle:
                    job_state = fle.read().rsplit()

                    len_js = len(job_state)
                    assert type(job_state) == list and len_js == 1, "error"
                    # if type(job_state) == list and len(job_state) == 1:
                    job_state = job_state[0]

            except:
                job_state = None
                pass

        return(job_state)
        #__|

    def job_info_batch(self, path_i="."):
        """

        Args:
            path_i:
        """
        #| - job_info_batch

        data_dict = self.cluster.job_info_batch(path_i=path_i)

        # jobs_file_path = self.jobs_list_dir


        # jobs_file_path = self.jobs_list_dir + "/jobs.csv"
        # df_new = pd.DataFrame([data_dict])
        # if os.path.isfile(jobs_file_path):
        #     df = pd.read_csv(jobs_file_path)
        #     df = df.append(df_new)
        # else:
        #     df = df_new
        # df.to_csv(jobs_file_path, index=False)


        return(data_dict)
        #__|

    def submit_job(self, **kwargs):
        """Call cluster specific job submission method.

        Notes:
            make sure to add 'path_i' key

        Args:
            **kwargs:
        """
        #| - submit_job
        kwargs = merge_two_dicts(self.default_sub_params, kwargs)

        #| - Checking if job has already been submitted
        if "path_i" in kwargs:
            path_i = kwargs["path_i"]
        else:
            path_i = "."

        if self.is_job_submitted(path_i=path_i):
            return(None)
        #__|

        #| - Writing Job Submission Parameters
        with open(".submission_params.json", "w") as fle:
            json.dump(kwargs, fle, indent=2, skipkeys=True)
        #__|

        self.cluster.submit_job_clust(**kwargs)

        #| - Writing Cluster System Info to File
        if "path_i" in kwargs:
            path_i = kwargs["path_i"]

            with open(".cluster_sys", "w") as fle:
                # fle.write(self.cluster_sys)
                fle.write(self.cluster_sys + "\n")
        #__|

        #__|

    #__| **********************************************************************

################################################################################


class EdisonCluster(ComputerCluster):
    """NERSC Edison computing cluster.


    I'll try to get this to work with Cori as well
    """

    #| - EdisonCluster ********************************************************
    def __init__(self, root_dir="."):
        """Initialize Sherlock cluster instance.

        Args:
            root_dir:
        """
        #| - __init__
        # 24 cores (edison) 32 cpus (cori)

        nersc_host = os.environ["NERSC_HOST"]
        if nersc_host == "cori":
            self.cores_per_node = 32  # <-- Cori Haswell
            # self.cores_per_node = 58  # <-- Cori KNL  #TODO Switch to these

        elif nersc_host == "edison":
            self.cores_per_node = 24


        # self.job_queue_dir = "/u/if/flores12/usr/bin"

        self.job_data_dir = ""
        self.root_dir = root_dir

        self.default_sub_params = self.default_submission_parameters()

        self.queues = self.__queue_types__()

        # self.job_queue_dir = "/u/if/flores12/usr/bin"

        self.job_state_keys = self.job_state_dict()
        self.job_queue_state_key = "STAT"  # COMBAK

        self.error_file = "job.err"
        self.out_file = "job.out"

        # self.aws_dir = os.environ["aws_sc"]
        # self.job_queue_dir = self.aws_dir + "/jobs_bin"
        # self.job_queue_state_key = "job_status"
        #__|

    def default_submission_parameters(self):
        """Defaul SLURM parameters for Sherlock cluster."""
        #| - default_submission_parameters

        def_params = {
            "queue": "regular",  # -p flag | regular, debug
            ##SBATCH -p regular
            #SBATCH -p debug
            "nodes": "10",
            #SBATCH -N 10 # 24 cpus per node on edison, 32 per node on cori haswell, 50? on knl

            "account": "m2997",  # -A flag
            #SBATCH -A  m2997

            "wall_time": "180",
            #SBATCH -t 00:30:00

            "priority": "scavenger",  # --qos -q flag

            ##SBATCH --qos=scavenger
            ##SBATCH --qos=premium
            ##SBATCH --qos=debug
            ##SBATCH --qos=regular

            "constraints": "haswell",
            ##SBATCH -C haswell #this is for cori haswell (old)
            ##SBATCH -C knl #new cori

            #SBATCH -e job.err#SBATCH -o job.out
            }

        return(def_params)
        #__|

    def submit_job_clust(self, **kwargs):
        """Submit job to sherlck.

        Args:
            **kwargs:
        """
        #| - submit_job
        time.sleep(1.5)

        #| - Merging Submission Parameters
        params = merge_two_dicts(self.default_sub_params, kwargs)

        path = params["path_i"]


        # Fixing debug flag specification
        if params["priority"] == "debug":
            params["queue"] = "debug"

        if params["queue"] == "debug":
            params["priority"] = "debug"
        #__|

        #| - Submit Job
        os.chdir(path)

        if params["job_name"] == "Default":
            params["job_name"] = os.getcwd()

        print("submitting job")
        os.system("chmod 777 *")
        # bash_command = "/u/if/flores12/bin/qv model.py"
        #__| **** TEMP

        #| - Create vasp_run script
        os.system("cd $SLURM_SUBMIT_DIR")
        os.system("export TMPDIR=$SLURM_SUBMIT_DIR")
        os.system("export VASP_SCRIPT=./run_vasp.py")
        os.system("echo import os > run_vasp.py")

        exitcode_line = "exitcode = os.system('srun -n " + \
            str(int(self.cores_per_node * int(params["nodes"]))) + \
            " /project/projectdirs/m2997/special_edison')"

        line_2 = 'echo ' + '"' + exitcode_line + '" >> run_vasp.py'
        os.system(line_2)  # on edison
        #__|

        #| - Bash Submisssion Command
        bash_command = "/usr/bin/sbatch "


        # The -q flag is being used in place of the -p flag
        # Only the -q needs to be defined
        bash_command += "-q " +                 str(params["queue"])       + " "
        # bash_command += "-p " +                 str(params["queue"])       + " "

        bash_command += "--nodes " +            str(params["nodes"])       + " "
        bash_command += "--time " +             str(params["wall_time"])   + " "

        # bash_command += "--qos " +              str(params["priority"])    + " "  # Didn't work

        bash_command += "--output " +           str(params["out_file"])    + " "
        bash_command += "--error " +            str(params["err_file"])    + " "
        bash_command += "-C haswell "

        bash_command += params["job_script"]

        print("Bash Submission Command:")
        print(bash_command)
        #__|

        try:
            output = subprocess.Popen(
                bash_command,
                stdout=subprocess.PIPE,
                shell=True,
                )
            sub_time = datetime.datetime.now().isoformat()
        # except subprocess.CalledProcessError, e:
        except subprocess.CalledProcessError as e:
            print("Ping stdout output:\n", e.output)

            os.chdir(self.root_dir)
            print("JOB SKIPPED: ")
            return(None)

        #| - Parsing Output
        # out, err = pickle.load(open("job_sub_output.pickle", "r"))

        try:
            # job_id = int(out_list[-1])
            out, err = output.communicate()
            out_copy = copy.deepcopy(out)
            out = out.strip()
            out_list = out.split(" ")
            job_id = int(out_list[-1])

        except:
            print("Couldn't parse for jobid")
            job_id = None

        # out = output.communicate()[0]
        # out_copy = copy.deepcopy(out)
        #
        # ind = out.find("job")
        # out = out[ind + 3:]
        #
        # jobid = re.sub("[^0-9]", "", out)
        #
        # try:
        #     jobid = int(jobid)
        #
        # except:
        #     print("Couldn't parse for jobid | !@!!")
        #     jobid = None
        #     pass
        #
        # if type(jobid) == int:
        #     jobid = jobid
        # else:
        #     jobid = None
        #__|

        #| - Writing Files
        with open(".SUBMITTED", "w") as fle:
            fle.write("\n")

        with open(".bash_comm", "w") as fle:
            fle.write(str(bash_command) + str("\n"))

        with open(".jobid", "w") as fle:
            fle.write(str(job_id) + str("\n"))

        if sys.version_info >= (3, 0):
            with open(".sub_out", "wb") as fle:
                fle.write(out_copy)

        else:
            with open(".sub_out", "w") as fle:
                fle.write(out_copy)
        #__|

        os.chdir(self.root_dir)

        #| - Save subprocess output for analysis
        # import pickle
        #
        # pickle.dump(
        #     output.communicate(),
        #     open("job_sub_output.pickle", "wb"),
        #     )
        # return(output)
        #__|

        # return(out, jobid)
        #__|

    def job_state_dict(self):
        """
        """
        #| - job_state_dict
        job_state_dict = {
            "PD": "PENDING",
            "R": "RUNNING",
            "CF": "CONFIGURING",
            "SUCCEEDED": "SUCCEEDED",

            # "FAILED": "FAILED",
            # "STARTING": "STARTING",
            # "RUNNABLE": "PENDING",
            # "SUBMITTED": "SUBMITTED",
            }

        return(job_state_dict)
        #__|

    def __queue_types__(self):
        """Queue types for Edison cluster
        """
        #| - __queue_types__
        queue_list = [
            "regular",
            "debug",
            "premium",
            ]

        return(queue_list)
        #__|

    def job_info_batch(self, job_id, path_i=None):
        """
        """
        #| - job_info_batch
        data_dict = slurm_squeue_parse(
            job_id,
            path_i=path_i,
            queue_state_key=self.job_queue_state_key,
            job_state_dict=self.job_state_keys,
            )

        return(data_dict)

        #| - __old__
        # bash_comm = "squeue -j " + str(job_id)
        #
        # try:
        #     out = subprocess.check_output(
        #         bash_comm,
        #         shell=True,
        #         stderr=subprocess.STDOUT,
        #         )
        #
        #     # 'JOBID PARTITION NAME USER ST TIME NODES NODELIST(REASON)'
        #     out = out.splitlines()
        #     out = out[1].split(" ")
        #     out = [i for i in out if i != '']
        #
        #     data_dict = {
        #         "PARTITION": out[1],
        #         "STAT": out[4],
        #         # "CF":
        #         }
        #
        #     if path_i is not None:
        #         key = self.job_queue_state_key
        #         with open(path_i + "/.QUEUESTATE", "w") as fle:
        #             fle.write(self.job_state_keys[data_dict[key]])
        #             fle.write("\n")
        #
        # except subprocess.CalledProcessError:
        #     data_dict = None
        #     pass
        #
        # except:
        #     data_dict = None
        #     pass
        #
        #
        # return(data_dict)
        #
        # bash_comm = "squeue -j " + str(job_id)
        #
        # try:
        #     out = subprocess.check_output(
        #         bash_comm,
        #         shell=True,
        #         stderr=subprocess.STDOUT,
        #         )
        #
        #     # 'JOBID PARTITION NAME USER ST TIME NODES NODELIST(REASON)'
        #     out = out.splitlines()
        #     out = out[1].split(" ")
        #     out = [i for i in out if i != '']
        #
        #     data_dict = {
        #         "PARTITION": out[1],
        #         "STAT": out[4],
        #         # "CF":
        #         }
        #
        #     if path_i is not None:
        #         key = self.job_queue_state_key
        #         with open(path_i + "/.QUEUESTATE", "w") as fle:
        #             fle.write(self.job_state_keys[data_dict[key]])
        #             fle.write("\n")
        #
        # except subprocess.CalledProcessError:
        #     data_dict = None
        #     pass
        #
        # except:
        #     print("tmp except final")
        #     data_dict = None
        #     pass
        #
        # # TEMP_PRINT
        # print(data_dict)
        #
        # return(data_dict)
        #__|

        #__|

    def completed_file(self, path_i="."):
        """Check whether ".FINISHED" file exists.

        Indicates that the job has gone to completion

        Args:
            path_i:
        """
        #| - completed_file
        completed_fle = False
        if os.path.exists(path_i + "/.FINISHED"):
            completed_fle = True

        return(completed_fle)
        #__|

    def job_state(self, path_i="."):
        """Return job state of path_i --> job_i.

        Args:
            path_i
        """
        #| - job_state
        job_id = self.get_jobid(path_i=path_i)

        # print("compute_env job_state job_id:")
        # print(job_id)


        job_state_out = None
        if job_id is not None:
            job_info = self.job_info_batch(job_id, path_i=path_i)

            # print("job_info **(***(*))")
            # print(job_info)

            if job_info is not None:
                key = self.job_queue_state_key
                if key in job_info:
                    job_state_out = job_info[key]
                    job_state_out = self.job_state_keys[job_state_out]

        #| - Checking for "completed" file indicating success
        completed_fle = self.completed_file(path_i=path_i)
        if completed_fle:
            job_state_out = self.job_state_keys["SUCCEEDED"]
        #__|

        return(job_state_out)
        #__|

    def get_jobid(self, path_i="."):
        """Return job ID of job_i.

        Args:
            path_i:
        """
        #| - get_jobid
        fileid_path = path_i + "/.jobid"
        if os.path.isfile(fileid_path):
            with open(path_i + "/.jobid") as fle:
                jobid = fle.read().strip()
        else:
            jobid = None

        return(jobid)
        #__|

    #__| **********************************************************************


class SLACCluster(ComputerCluster):
    """SLAC computing cluster."""

    #| - SLACCluster **********************************************************
    def __init__(self,
        root_dir=".",
        ):
        """
        """
        #| - __init__
        self.root_dir = root_dir
        self.default_sub_params = self.default_submission_parameters()
        self.job_data_dir = ""

        self.job_queue_dir = "/u/if/flores12/usr/bin"
        self.queues = self.__queue_types__()

        self.job_state_keys = self.job_state_dict()
        self.job_queue_state_key = "STAT"
        #__|

    def job_state_dict(self):
        """
        """
        #| - job_state_dict
        job_state_dict = {
            # "PENDING": "PEND",
            # "FINISHED": "DONE",
            # "RUNNING": "RUN",
            # "FAILED": "EXIT",
            "PEND": "PENDING",
            "DONE": "SUCCEEDED",
            "RUN": "RUNNING",
            "EXIT": "FAILED",
            "UNKWN": "UNKNOWN",
            }

        return(job_state_dict)
        #__|

    def __queue_types__(self):
        """Queue types for SLAC cluster
        """
        #| - __queue_types__
        queue_list = [
            "suncat-test",
            "suncat",
            "suncat2",
            "suncat2-xlong",
            "suncat3",
            "suncat3-xlong",
            "suncat-xlong",
            ]

        return(queue_list)
        #__|

    def default_submission_parameters(self):
        """
        """
        #| - default_submission_parameters
        def_params = {
            "queue": "suncat",
            "cpus": "8",
            "wall_time": "3000",
            # "memory": "6000",
            # "job_name":     "Default",
            "job_script": "model.py"
            }

        return(def_params)
        #__|

    def submit_job_clust(self, **kwargs):
        """FIX This should update data table

        ABSOLUTELY NEED TO DEFINE 'path_i' FOR THIS TO WORK!!!!!!!!!!!!!!!!!!

        Submits job to aws cluster. Copies PythonModules folder into
        job directory
        """
        #| - submit_job

        #| - Merging Submission Parameters
        params = merge_two_dicts(self.default_sub_params, kwargs)

        path = params["path_i"]
        #__|

        #| - Checking if job has already been submitted
        # if self.is_job_submitted():
        #     return(None)
        #__|

        #| - Submit Job *******************************************************

        os.chdir(path)

        if params["job_name"] == "Default":
            params["job_name"] = os.getcwd()

        print("submitting job")

        os.system("chmod 777 *")

        # bash_command = "/u/if/flores12/bin/qv model.py"

        bash_command = "/afs/slac/g/suncat/bin/dobsub "
        bash_command += "-q " + str(params["queue"]) + " "
        bash_command += "-n " + str(params["cpus"]) + " "
        bash_command += "-W " + str(params["wall_time"]) + " "
        bash_command += "-o " + str(params["out_file"]) + " "
        bash_command += "-e " + str(params["err_file"]) + " "

        # bash_command += "-o job.out "
        # bash_command += "-e job.err "

        # bash_command += "-M " + params["memory"] + " "
        bash_command += "-J " + params["job_name"] + " "
        bash_command += params["job_script"]

        #| - FIXME Python 2 --> 3
        print("Python version info")
        print(sys.version_info)
        if (sys.version_info > (3, 0)):
            try:
                output = subprocess.Popen(
                    bash_command,
                    stdout=subprocess.PIPE,
                    shell=True,
                    )
                sub_time = datetime.datetime.now().isoformat()


            except subprocess.CalledProcessError as e:
                print("Ping stdout output:\n", e.output)
                os.chdir(self.root_dir)
                print("JOB SKIPPED: ")
                return(None)

        else:
            try:
                output = subprocess.Popen(
                    bash_command,
                    stdout=subprocess.PIPE,
                    shell=True,
                    )
                sub_time = datetime.datetime.now().isoformat()

            # except subprocess.CalledProcessError, e:
            except subprocess.CalledProcessError as e:
                print("Ping stdout output:\n", e.output)
                # print "Ping stdout output:\n", e.output
                os.chdir(self.root_dir)
                print("JOB SKIPPED: ")
                return(None)

        #| - __old
        # try:
        #     output = subprocess.Popen(
        #         bash_command,
        #         stdout=subprocess.PIPE,
        #         shell=True,
        #         )
        #     sub_time = datetime.datetime.now().isoformat()
        #
        # # except subprocess.CalledProcessError, e:
        # except subprocess.CalledProcessError as e:
        #
        #     print "Ping stdout output:\n", e.output
        #
        #     os.chdir(self.root_dir)
        #     print("JOB SKIPPED: ")
        #     return(None)

        #__|

        #__|

        #__|

        #| - Parsing Output
        out = output.communicate()[0]
        ind = out.find("Job <")

        out1 = out[ind:]
        ind1 = out1.find("<")
        ind2 = out1.find(">")

        jobid = out1[ind1 + 1:ind2]

        if jobid.isdigit():
            jobid = int(jobid)
        else:
            jobid = None
        #__|

        #| - Writing Files
        with open(".SUBMITTED", "w") as fle:
            fle.write("\n")

        with open(".bash_comm", "w") as fle:
            fle.write(str(bash_command) + str("\n"))

        with open(".jobid", "w") as fle:
            fle.write(str(jobid) + str("\n"))

        with open(".sub_out", "w") as fle:
            fle.write(out)
        #__|

        os.chdir(self.root_dir)

        return(out, jobid)
        #__|

    def get_jobid(self, path_i="."):
        """
        """
        #| - get_jobid
        # path_i = "."

        if os.path.isfile(path_i + "/.jobid"):
            with open(path_i + "/.jobid") as fle:
                jobid = fle.read().strip()
        else:
            jobid = None

        return(jobid)
        #__|

    def job_info_batch(self, path_i="."):
        """
        """
        #| - job_info_batch
        job_id = self.get_jobid(path_i)

        #| - If No Job ID File Return None
        if job_id is None:
            return(None)
        #__|

        bash_comm = "/usr/local/bin/bjobs -w" + " " + job_id
        out = subprocess.check_output(
            bash_comm,
            shell=True,
            stderr=subprocess.STDOUT,
            )

        #| - Checking if Job Id Still in Batch System
        if "is not found" in out:
            print("Job ID no longer in batch system, or ID is wrong")
            return(None)
        #__|

        #| - Parsing bsub Output into Dict
        out = out.split()
        headers = out[0:8]
        data = out[8:]
        data_1 = data[0:7]

        data_dict = {
            "JOBID": data_1[0],
            "USER": data_1[1],
            "STAT": data_1[2],
            "QUEUE": data_1[3],
            "FROM_HOST": data_1[4],
            "EXEC_HOST": data_1[5],
            "JOB_NAME": data_1[6],
            }

        time = data[7:]
        time = "_".join(time)
        data_dict["SUBMIT_TIME"] = time
        #__|

        #| - bjob Command to Get Job Path From Job ID
        bash_comm_2 = "/usr/local/bin/bjobs -o" + " 'exec_cwd' " + job_id
        out2 = subprocess.check_output(bash_comm_2, shell=True)
        out2 = out2.split()
        data_dict["EXEC_CWD"] = out2[1]
        #__|

        # self.write_job_queue_state_file(path=path_i)

        #| - write_job_queue_state_file
        key = self.job_queue_state_key
        with open(path_i + "/.QUEUESTATE", "w") as fle:
            fle.write(self.job_state_keys[data_dict[key]])
            fle.write("\n")
        #__|

        return(data_dict)
        #__|

    def job_state(self, path_i="."):
        """Query job state.

        # FIXME This should update jobs table

        Args:
            path_i
        """
        #| - job_state
        job_info = self.job_info_batch(path_i=path_i)

        if job_info is not None:
            key = self.job_queue_state_key
            if key in job_info:
                job_state_out = job_info[key]
                job_state_out = self.job_state_keys[job_state_out]
        else:
            job_state_out = None

        return(job_state_out)
        #__|

    #__| **********************************************************************


class SherlockCluster(ComputerCluster):
    """Sherlock computing cluster."""

    #| - SherlockCluster ******************************************************
    def __init__(self, root_dir="."):
        """Initialize Sherlock cluster instance.

        Args:
            root_dir:
        """
        #| - __init__
        # self.job_queue_dir = "/u/if/flores12/usr/bin"

        self.job_data_dir = ""
        self.root_dir = root_dir

        self.username = self.__parse_username__()
        self.default_sub_params = self.default_submission_parameters()

        self.queues = self.__queue_types__()

        # self.job_queue_dir = "/u/if/flores12/usr/bin"

        self.job_state_keys = self.job_state_dict()
        self.job_queue_state_key = "STAT"


        self.error_file = "job.err"
        self.out_file = "job.out"

        # self.aws_dir = os.environ["aws_sc"]
        # self.job_queue_dir = self.aws_dir + "/jobs_bin"

        # self.job_queue_state_key = "job_status"
        #__|


    def __parse_username__(self):
        """
        """
        #| - __parse_username__
        username = os.environ.get("USER")

        cond_1 = False
        if type(username) == str:
            cond_1 = True

        cond_2 = False
        if username is not None:
            cond_2 = True

        return(username)
        #__|


    def default_submission_parameters(self):
        """Defaul SLURM parameters for Sherlock cluster."""
        #| - default_submission_parameters

        def_params = {
            "queue": "owners,iric,normal",  # -p flag
            "nodes": "1",  # --nodes
            "cpus": "16",  # --ntasks-per-node
            "memory": "4000",  # --mem-per-cpu
            "wall_time": "720",  # --time (720min -> 12hrs)
            "job_name": "Default",  # --job-name
            "priority": "normal",  # --qos
            # "email": "flores12@stanford.edu",  # --mail-user
            "email": self.username + "@stanford.edu",  # --mail-user
            "email_mess": "FAIL",  # --mail-type
            }

        return(def_params)
        #__|

    def submit_job_clust(self, **kwargs):
        """Submits job to sherlck.

        Args:
            **kwargs:
        """
        #| - submit_job
        time.sleep(1.5)

        #| - Merging Submission Parameters
        params = merge_two_dicts(self.default_sub_params, kwargs)

        path = params["path_i"]
        #__|

        #| - Submit Job
        os.chdir(path)

        if params["job_name"] == "Default":
            params["job_name"] = os.getcwd()

        print("submitting job")
        os.system("chmod 777 *")
        # bash_command = "/u/if/flores12/bin/qv model.py"
        #__| **** TEMP

        #| - Bash Submisssion Command
        bash_command = "/usr/bin/sbatch "

        bash_command += "-p " +                 str(params["queue"])       + " "
        bash_command += "--nodes " +            str(params["nodes"])       + " "
        bash_command += "--ntasks-per-node " +  str(params["cpus"])        + " "
        bash_command += "--mem-per-cpu " +      str(params["memory"])      + " "
        bash_command += "--time " +             str(params["wall_time"])   + " "
        bash_command += "--job-name " +         str(params["job_name"])    + " "
        bash_command += "--qos " +              str(params["priority"])    + " "
        bash_command += "--mail-user " +        str(params["email"])       + " "
        bash_command += "--mail-type " +        str(params["email_mess"])  + " "
        bash_command += "--output " +           str(params["out_file"])    + " "
        bash_command += "--error " +            str(params["err_file"])    + " "
        bash_command += "-C CPU_GEN:HSW "  # COMBAK Formalize this cpu architecture filter

        bash_command += params["job_script"]

        print("Bash Submission Command:")
        print(bash_command)
        #__|

        try:

            if sys.version_info[0] < 3:
                # raise Exception("Must be using Python 3")
                output = subprocess.Popen(
                    bash_command,
                    stdout=subprocess.PIPE,
                    shell=True,
                    # encoding="utf8",
                    )

            else:
                output = subprocess.Popen(
                    bash_command,
                    stdout=subprocess.PIPE,
                    shell=True,
                    encoding="utf8",
                    )

            sub_time = datetime.datetime.now().isoformat()
        # except subprocess.CalledProcessError, e:
        except subprocess.CalledProcessError as e:
            print("Ping stdout output:\n", e.output)

            os.chdir(self.root_dir)
            print("JOB SKIPPED: ")
            return(None)

        #| - Parsing Output
        out = output.communicate()[0]
        out_copy = copy.deepcopy(out)

        ind = out.find("job")

        out = out[ind + 3:]

        jobid = re.sub("[^0-9]", "", out)

        try:
            jobid = int(jobid)

        except:
            print("Couldn't parse for jobid | !@!!")
            jobid = None
            pass

        if type(jobid) == int:
            jobid = jobid
        else:
            jobid = None
        #__|

        #| - Writing Files
        with open(".SUBMITTED", "w") as fle:
            fle.write("\n")

        with open(".bash_comm", "w") as fle:
            fle.write(str(bash_command) + str("\n"))

        with open(".jobid", "w") as fle:
            fle.write(str(jobid) + str("\n"))

        with open(".sub_out", "w") as fle:
            fle.write(out_copy)

        #| - Writing Job Submission Parameters
        with open(".submission_params_2.json", "w") as fle:
            json.dump(params, fle, indent=2, skipkeys=True)
        #__|

        #__|

        os.chdir(self.root_dir)

        return(out, jobid)

        #__|

    def job_state_dict(self):
        """
        """
        #| - job_state_dict
        job_state_dict = {
            "PD": "PENDING",
            "R": "RUNNING",
            "CF": "CONFIGURING",
            "SUCCEEDED": "SUCCEEDED",

            # "FAILED": "FAILED",
            # "STARTING": "STARTING",
            # "RUNNABLE": "PENDING",
            # "SUBMITTED": "SUBMITTED",
            }

        return(job_state_dict)
        #__|

    def __queue_types__(self):
        """Queue types for SLAC cluster
        """
        #| - __queue_types__
        queue_list = [
            "owners",
            "iric",
            ]

        return(queue_list)
        #__|

    def get_jobid(self, path_i="."):
        """Return the job id.

        Args:
            path_i:
        """
        #| - get_jobid
        # # path_i = "."
        # fileid_path = path_i + "/.jobid"
        # # print(fileid_path)
        # if os.path.isfile(fileid_path):
        #     with open(path_i + "/.jobid") as fle:
        #         jobid = fle.read().strip()
        # else:
        #     jobid=None
        #
        # return(jobid)
        #__|

    def job_info_batch(self, job_id, path_i=None):
        """
        """
        #| - job_info_batch
        data_dict = slurm_squeue_parse(
            job_id,
            path_i=path_i,
            queue_state_key=self.job_queue_state_key,
            job_state_dict=self.job_state_keys,
            )

        return(data_dict)

        #| - __old__
        # bash_comm = "squeue -j " + str(job_id)
        #
        # try:
        #     out = subprocess.check_output(
        #         bash_comm,
        #         shell=True,
        #         stderr=subprocess.STDOUT,
        #         )
        #
        #     # 'JOBID PARTITION NAME USER ST TIME NODES NODELIST(REASON)'
        #     out = out.splitlines()
        #     out = out[1].split(" ")
        #     out = [i for i in out if i != '']
        #
        #     data_dict = {
        #         "PARTITION": out[1],
        #         "STAT": out[4],
        #         # "CF":
        #         }
        #
        #     if path_i is not None:
        #         key = self.job_queue_state_key
        #         with open(path_i + "/.QUEUESTATE", "w") as fle:
        #             fle.write(self.job_state_keys[data_dict[key]])
        #             fle.write("\n")
        #
        # except subprocess.CalledProcessError:
        #     data_dict = None
        #     pass
        #
        # except:
        #     data_dict = None
        #     pass
        #
        #
        # return(data_dict)
        #__|

        #__|

    def completed_file(self, path_i="."):
        """
        Check whether ".FINISHED" file exists.

        Indicates that the job has gone to completion

        Args:
            path_i:
        """
        #| - completed_file
        completed_fle = False
        if os.path.exists(path_i + "/.FINISHED"):
            completed_fle = True

        return(completed_fle)
        #__|

    def job_state(self, path_i="."):
        """
        Return job state of path_i --> job_i.

        Args:
            path_i
        """
        #| - job_state
        job_id = self.get_jobid(path_i=path_i)

        job_state_out = None
        if job_id is not None:
            job_info = self.job_info_batch(job_id, path_i=path_i)

            if job_info is not None:
                key = self.job_queue_state_key
                if key in job_info:
                    job_state_out = job_info[key]
                    job_state_out = self.job_state_keys[job_state_out]

        #| - Checking for "completed" file indicating success
        completed_fle = self.completed_file(path_i=path_i)
        if completed_fle:
            job_state_out = self.job_state_keys["SUCCEEDED"]
        #__|

        return(job_state_out)

        #__|

    def get_jobid(self, path_i="."):
        """
        Return job ID of job_i.

        Args:
            path_i:
        """
        #| - get_jobid
        fileid_path = path_i + "/.jobid"
        if os.path.isfile(fileid_path):
            with open(path_i + "/.jobid") as fle:
                jobid = fle.read().strip()
        else:
            jobid = None

        return(jobid)
        #__|

    #__| **********************************************************************


class AWSCluster(ComputerCluster):
    """AWS EC2 computing resource."""

    #| - AWSCluster ***********************************************************
    def __init__(self,
        root_dir=".",
        ):
        """
        """
        #| - __init__
        self.job_data_dir = "/simulation"
        self.root_dir = root_dir
        self.default_sub_params = self.default_submission_parameters()
        self.aws_dir = os.environ["aws_sc"]
        self.job_queue_dir = self.aws_dir + "/jobs_bin"
        self.job_state_keys = self.job_state_dict()
        self.queues = self.__queue_types__()
        self.job_queue_state_key = "job_status"
        self.error_file = "err"
        self.out_file = "out"
        #__|

    def default_submission_parameters(self):
        """
        """
        #| - default_submission_parameters
        def_params = {
            "queue": "medium",
            "cpus": "default",
            # "wall_time":    "2000",
            # "memory": "6000",
            # "job_name":     "Default",
            "job_script": "model.py",

            "copy_PythonModules": True,
            "copy_PythonPackages": True,
            }

        return(def_params)
        #__|

    def submit_job_clust(self, **kwargs):
        """Submit job to AWS cluster.

        Copies PythonModules and PythonPackages folder into job directory

        Args:
            kwargs
        """
        #| - submit_job_clust

        #| - Job Parameters
        params = merge_two_dicts(self.default_sub_params, kwargs)

        path = params["path_i"]
        copy_PythonModules = params["copy_PythonModules"]
        copy_PythonPackages = params["copy_PythonPackages"]
        cpus = params["cpus"]
        queue = params["queue"]
        #__|

        root_dir = os.getcwd()
        if path is None:
            path = root_dir

        #| - Checking if job has already been submitted
        os.chdir(path)
        if os.path.isfile(".SUBMITTED"):
            print("Directory already submitted, will be skipped")
            os.chdir(root_dir)
            return(None)  # <-------- SKIP JOB ---------------------------------
        else:
            os.chdir(root_dir)
        #__|

        self.__copy_pyth_mods_packs_to_job_dir__(
            path,
            copy_mods=copy_PythonModules,
            copy_packs=copy_PythonPackages,
            )

        #| - Submit Job
        # Args: path, root_dir, queue, cpus

        os.chdir(path)

        if os.path.isfile(".SUBMITTED"):
            print("Directory already submitted, will be skipped")
            os.chdir(root_dir)
            return(None)
        else:
            print("submitting job")
            aws_dir = os.environ["aws_sc"]

            if cpus == "default":
                bash_command = aws_dir + "/bin/trisub -q " + queue
                # bash_command = aws_dir + "/matr.io/bin/trisub -q " + queue
            else:
                bash_command = aws_dir + "/bin/trisub -c " + str(cpus) + \
                    " -q " + queue

            try:
                output = subprocess.check_output(bash_command, shell=True)
                sub_time = datetime.datetime.now().isoformat()
            # except subprocess.CalledProcessError, e:
            except subprocess.CalledProcessError as e:
                print("Ping stdout output:\n", e.output)

                os.chdir(root_dir)
                print("JOB SKIPPED: ")
                return(None)
        #__|

        os.system("chmod 777 " + path + "/*")
        os.system("chmod 777 " + path)

        #| - Parsing Submission for Job ID
        output = output.splitlines()
        for line in output:
            if "jobId" in line:
                lst = line.split('"')
                job_id_ind = (lst.index("jobId") + 2)
                jobId = lst[job_id_ind]

                file = open(".jobid", "w")
                file.write(jobId + "\n")
        #__|

        file = open(".SUBMITTED", "w")
        file.close()

        os.chdir(root_dir)

        #| - Querying AWS For Job Info
        job_queue_dict = self.job_info_batch(jobId)
        job_queue_dict["submit_time"] = sub_time

        jobs_file_path = self.job_queue_dir + "/jobs.csv"

        df_new = pd.DataFrame([job_queue_dict])
        if os.path.isfile(jobs_file_path):
            df = pd.read_csv(jobs_file_path)
            df = df.append(df_new)
        else:
            df = df_new

        df.to_csv(jobs_file_path, index=False)
        #__|

        return job_queue_dict
        #__|

    def job_state_dict(self):
        """
        """
        #| - job_state_dict
        job_state_dict = {
            "PENDING": "PENDING",
            "SUCCEEDED": "SUCCEEDED",
            "FAILED": "FAILED",
            "RUNNING": "RUNNING",
            "STARTING": "STARTING",

            # RUNNABLE in AWS really means pending
            # "RUNNABLE": "RUNNABLE",

            "RUNNABLE": "PENDING",
            "SUBMITTED": "SUBMITTED",
            }

        return(job_state_dict)
        #__|

    def __queue_types__(self):
        """Queue types for AWS cluster
        """
        #| - __queue_types__
        queue_list = [
            "test",
            "small",
            "medium",
            "large",
            ]

        return(queue_list)
        #__|

    def __copy_pyth_mods_packs_to_job_dir__(
        self,
        path_i,
        copy_mods=True,
        copy_packs=True,
        ):
        """
        """
        #| - __copy_pyth_mods_packs_to_job_dir__
        copy_PythonModules = copy_mods
        copy_PythonPackages = copy_packs

        #| - Copy PYTHONMODULES to Job Directory
        if copy_PythonModules:
            if os.path.isdir(path_i + "/PythonModules") is True:
                print("PythonModules already exists, erasing and recopying")
                shutil.rmtree(path_i + "/PythonModules")
                # py_mod = os.environ["python_modules"]  # Old pyth mods dir
                py_mod = os.environ["PYTHONMODULES"]
                shutil.copytree(py_mod, path_i + "/PythonModules")
            else:
                # py_mod = os.environ["python_modules"]
                py_mod = os.environ["PYTHONMODULES"]
                shutil.copytree(py_mod, path_i + "/PythonModules")
        #__|

        #| - Copy Python Packages to Job Directory
        if copy_PythonPackages:
            if os.path.isdir(path_i + "/PythonPackages") is True:
                print("PythonPackages already exists, erasing and recopying")
                shutil.rmtree(path_i + "/PythonPackages")
                py_pack = os.environ["PYTHONPACKAGES"]
                shutil.copytree(py_pack, path_i + "/PythonPackages")
            else:
                py_pack = os.environ["PYTHONPACKAGES"]
                shutil.copytree(py_pack, path_i + "/PythonPackages")
        #__|

        #__|

    def get_jobid(self, path_i="."):
        """
        """
        #| - get_jobid
        # path_i = "."
        fileid_path = path_i + "/.jobid"
        if os.path.isfile(fileid_path):
            with open(path_i + "/.jobid") as fle:
                jobid = fle.read().strip()
        else:
            jobid = None

        return(jobid)
        #__|

    def job_state(self, path_i="."):
        """
        """
        #| - job_state
        job_id = self.get_jobid(path_i=path_i)
        job_state_out = None
        if job_id is not None:
            job_info = self.job_info_batch(job_id)

            if job_info is not None:
                key = self.job_queue_state_key
                if key in job_info:
                    job_state_out = job_info[key]
                    job_state_out = self.job_state_keys[job_state_out]

        #| - Checking for Spot Termination
        spot_term = self.spot_terminated(path_i=path_i)
        if spot_term:
            job_state_out = self.job_state_keys["FAILED"]
        #__|

        #| - Checking for "completed" file indicating success
        completed_fle = self.completed_file(path_i=path_i)
        if completed_fle:
            job_state_out = self.job_state_keys["SUCCEEDED"]
        #__|

        if os.path.isfile(path_i + "/.QUEUESTATE"):
            bash_comm = "chmod 777 " + path_i + "/.QUEUESTATE"
            os.system(bash_comm)

        with open(path_i + "/.QUEUESTATE", "w") as fle:
            fle.write(str(job_state_out))
            fle.write("\n")

        return(job_state_out)
        #__|

    def spot_terminated(self, path_i="."):
        """
        """
        #| - spot_terminated
        symlink_dir = path_i + "/simulation"

        spot_terminated = False
        if os.path.isdir(symlink_dir):

            if os.path.exists(symlink_dir + "/spotTerminated"):
                spot_terminated = True
                print("Spot Termination")

        return(spot_terminated)
        #__|

    def completed_file(self, path_i="."):
        """
        """
        #| - completed_file
        symlink_dir = path_i + "/simulation"

        completed_fle = False
        if os.path.isdir(symlink_dir):
            if os.path.exists(symlink_dir + "/completed"):
                completed_fle = True

        return(completed_fle)
        #__|

    def job_info_batch(self, job_id):
        """
        """
        #| - job_info_batch
        import boto3
        batch = boto3.client("batch")
        job_descriptions = batch.describe_jobs(jobs=[job_id])

        #| - Checking if Job is in AWS Batch
        if len(job_descriptions["jobs"]) == 0:
            return("job not in batch system")
        else:
            job_info = job_descriptions["jobs"][0]
        #__|

        job_status = job_info["status"]
        job_path = job_info["parameters"]["model"]
        job_ram = job_info["parameters"]["ram"]
        job_cpus = job_info["parameters"]["cpus"]
        job_name = job_info["jobName"]

        job_queue_str = job_info["jobQueue"]
        if "small" in job_queue_str:
            job_queue = "small"
        elif "medium" in job_queue_str:
            job_queue = "medium"
        elif "large" in job_queue_str:
            job_queue = "large"
        elif "test" in job_queue_str:
            job_queue = "test"

        job_queue_dict = {"job_status": job_status, "job_path": job_path,
                          "job_id": job_id, "job_ram": job_ram,
                          "job_queue": job_queue, "job_cpus": job_cpus,
                          "job_name": job_name}

        return(job_queue_dict)
        #__|

    #__| **********************************************************************


class DummyCluster(ComputerCluster):
    """Placeholder class for when current cluster isn't supported."""
    #| - DummyCluster
    def __init__(self, root_dir="."):
        """
        """
        #| - __init__
        # self.job_queue_dir = "/u/if/flores12/usr/bin"

        self.job_data_dir = ""
        self.root_dir = root_dir

        # self.default_sub_params = self.default_submission_parameters()
        # self.queues = self.__queue_types__()


        # self.job_queue_dir = "/u/if/flores12/usr/bin"
        # self.job_state_keys = self.job_state_dict()
        # self.job_queue_state_key = "STAT"

        self.error_file = "job.err"
        self.out_file = "job.out"

        # self.aws_dir = os.environ["aws_sc"]
        # self.job_queue_dir = self.aws_dir + "/jobs_bin"

        # self.job_queue_state_key = "job_status"
        #__|


    def submit_job_clust(self, **kwargs):
        """just a placeholder function.

        Args:
            **kwargs:
        """
        #| - submit_job
        print("submit_job_clust | DummyCluster")
        print("Nothing happens!!")
        #__|


    #| - Not Needed

    # def default_submission_parameters(self):
    #     """Defaul SLURM parameters for
    #     """
    #     #| - default_submission_parameters
    #     def_params = {
    #         "queue": "owners,iric,normal",  # -p flag
    #         "nodes": "1",  # --nodes
    #         "cpus": "16",  # --ntasks-per-node
    #         "memory": "4000",  # --mem-per-cpu
    #         "wall_time": "720",  # --time (720min -> 12hrs)
    #         "job_name": "Default",  # --job-name
    #         "priority": "normal",  # --qos
    #         "email": "flores12@stanford.edu",  # --mail-user
    #         "email_mess": "FAIL",  # --mail-type
    #         }
    #
    #     return(def_params)
    #     #__|
    #
    # def submit_job_clust(self, **kwargs):
    #     """
    #     Submits job to sherlck.
    #     """
    #     #| - submit_job
    #     time.sleep(1.5)
    #
    #     #| - Merging Submission Parameters
    #     params = merge_two_dicts(self.default_sub_params, kwargs)
    #
    #     path = params["path_i"]
    #     #__|
    #
    #     #| - Submit Job
    #     os.chdir(path)
    #
    #     if params["job_name"] == "Default":
    #         params["job_name"] = os.getcwd()
    #
    #     print("submitting job")
    #     os.system("chmod 777 *")
    #     # bash_command = "/u/if/flores12/bin/qv model.py"
    #     #__| **** TEMP
    #
    #     #| - Bash Submisssion Command
    #     bash_command = "/usr/bin/sbatch "
    #
    #     bash_command += "-p " +                 str(params["queue"])       + " "
    #     bash_command += "--nodes " +            str(params["nodes"])       + " "
    #     bash_command += "--ntasks-per-node " +  str(params["cpus"])        + " "
    #     bash_command += "--mem-per-cpu " +      str(params["memory"])      + " "
    #     bash_command += "--time " +             str(params["wall_time"])   + " "
    #     bash_command += "--job-name " +         str(params["job_name"])    + " "
    #     bash_command += "--qos " +              str(params["priority"])    + " "
    #     bash_command += "--mail-user " +        str(params["email"])       + " "
    #     bash_command += "--mail-type " +        str(params["email_mess"])  + " "
    #     bash_command += "--output " +           str(params["out_file"])    + " "
    #     bash_command += "--error " +            str(params["err_file"])    + " "
    #     bash_command += "-C CPU_GEN:HSW "  # COMBAK Formalize this cpu architecture filter
    #
    #     bash_command += params["job_script"]
    #
    #     print("Bash Submission Command:")
    #     print(bash_command)
    #     #__|
    #
    #     try:
    #         output = subprocess.Popen(
    #             bash_command,
    #             stdout=subprocess.PIPE,
    #             shell=True,
    #             )
    #         sub_time = datetime.datetime.now().isoformat()
    #     # except subprocess.CalledProcessError, e:
    #     except subprocess.CalledProcessError as e:
    #         print("Ping stdout output:\n", e.output)
    #
    #         os.chdir(self.root_dir)
    #         print("JOB SKIPPED: ")
    #         return(None)
    #
    #     #| - Parsing Output
    #     out = output.communicate()[0]
    #     out_copy = copy.deepcopy(out)
    #
    #     ind = out.find("job")
    #     out = out[ind + 3:]
    #
    #     jobid = re.sub("[^0-9]", "", out)
    #
    #     try:
    #         jobid = int(jobid)
    #
    #     except:
    #         print("Couldn't parse for jobid | !@!!")
    #         jobid = None
    #         pass
    #
    #     if type(jobid) == int:
    #         jobid = jobid
    #     else:
    #         jobid = None
    #     #__|
    #
    #     #| - Writing Files
    #     with open(".SUBMITTED", "w") as fle:
    #         fle.write("\n")
    #
    #     with open(".bash_comm", "w") as fle:
    #         fle.write(str(bash_command) + str("\n"))
    #
    #     with open(".jobid", "w") as fle:
    #         fle.write(str(jobid) + str("\n"))
    #
    #     with open(".sub_out", "w") as fle:
    #         fle.write(out_copy)
    #     #__|
    #
    #     os.chdir(self.root_dir)
    #
    #     return(out, jobid)
    #
    #     #__|
    #
    # def job_state_dict(self):
    #     """
    #     """
    #     #| - job_state_dict
    #     job_state_dict = {
    #         "PD": "PENDING",
    #         "R": "RUNNING",
    #         "CF": "CONFIGURING",
    #         "SUCCEEDED": "SUCCEEDED",
    #
    #         # "FAILED": "FAILED",
    #         # "STARTING": "STARTING",
    #         # "RUNNABLE": "PENDING",
    #         # "SUBMITTED": "SUBMITTED",
    #         }
    #
    #     return(job_state_dict)
    #     #__|
    #
    # def __queue_types__(self):
    #     """Queue types for SLAC cluster
    #     """
    #     #| - __queue_types__
    #     queue_list = [
    #         "owners",
    #         "iric",
    #         ]
    #
    #     return(queue_list)
    #     #__|
    #
    # def get_jobid(self, path_i="."):
    #     """
    #     """
    #     #| - get_jobid
    #     # # path_i = "."
    #     # fileid_path = path_i + "/.jobid"
    #     # # print(fileid_path)
    #     # if os.path.isfile(fileid_path):
    #     #     with open(path_i + "/.jobid") as fle:
    #     #         jobid = fle.read().strip()
    #     # else:
    #     #     jobid=None
    #     #
    #     # return(jobid)
    #     #__|
    #
    # def job_info_batch(self, job_id, path_i=None):
    #     """
    #     """
    #     #| - job_info_batch
    #     bash_comm = "squeue -j " + str(job_id)
    #
    #     try:
    #         out = subprocess.check_output(
    #             bash_comm,
    #             shell=True,
    #             stderr=subprocess.STDOUT,
    #             )
    #
    #         # 'JOBID PARTITION NAME USER ST TIME NODES NODELIST(REASON)'
    #         out = out.splitlines()
    #         out = out[1].split(" ")
    #         out = [i for i in out if i != '']
    #
    #         data_dict = {
    #             "PARTITION": out[1],
    #             "STAT": out[4],
    #             # "CF":
    #             }
    #
    #         if path_i is not None:
    #             key = self.job_queue_state_key
    #             with open(path_i + "/.QUEUESTATE", "w") as fle:
    #                 fle.write(self.job_state_keys[data_dict[key]])
    #                 fle.write("\n")
    #
    #     except subprocess.CalledProcessError:
    #         data_dict = None
    #         pass
    #
    #     except:
    #         data_dict = None
    #         pass
    #
    #
    #     return(data_dict)
    #
    #     #__|
    #
    # def completed_file(self, path_i="."):
    #     """
    #     Check whether ".FINISHED" file exists.
    #
    #     Indicates that the job has gone to completion
    #
    #     Args:
    #         path_i:
    #     """
    #     #| - completed_file
    #     completed_fle = False
    #     if os.path.exists(path_i + "/.FINISHED"):
    #         completed_fle = True
    #
    #     return(completed_fle)
    #     #__|
    #
    # def job_state(self, path_i="."):
    #     """
    #     Return job state of path_i --> job_i.
    #
    #     Args:
    #         path_i
    #     """
    #     #| - job_state
    #     job_id = self.get_jobid(path_i=path_i)
    #
    #     job_state_out = None
    #     if job_id is not None:
    #         job_info = self.job_info_batch(job_id, path_i=path_i)
    #
    #         if job_info is not None:
    #             key = self.job_queue_state_key
    #             if key in job_info:
    #                 job_state_out = job_info[key]
    #                 job_state_out = self.job_state_keys[job_state_out]
    #
    #     #| - Checking for "completed" file indicating success
    #     completed_fle = self.completed_file(path_i=path_i)
    #     if completed_fle:
    #         job_state_out = self.job_state_keys["SUCCEEDED"]
    #     #__|
    #
    #     return(job_state_out)
    #
    #     #__|
    #
    # def get_jobid(self, path_i="."):
    #     """
    #     Return job ID of job_i.
    #
    #     Args:
    #         path_i:
    #     """
    #     #| - get_jobid
    #     fileid_path = path_i + "/.jobid"
    #     if os.path.isfile(fileid_path):
    #         with open(path_i + "/.jobid") as fle:
    #             jobid = fle.read().strip()
    #     else:
    #         jobid = None
    #
    #     return(jobid)
    #     #__|

    #__|

    #__|
