#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Class for computer cluster operations, mainly batch.

"""

# | - Import Modules
import os
import sys
import subprocess
import datetime
import copy
import time

from shutil import copyfile

# My Modules
from misc_modules.misc_methods import merge_two_dicts

from dft_job_automat.compute_env import ComputerCluster
from dft_job_automat.compute_env import slurm_squeue_parse
# __|



class NERSC_Cluster(ComputerCluster):
    """NERSC computing cluster.

    Two main architectures, knl and haswell
    """

    # | - NERSC_Cluster ********************************************************
    def __init__(self,
        root_dir=".",
        ):
        """Initialize Sherlock cluster instance.

        Args:
            root_dir:
        """
        # | - __init__
        # print(80 * "#")
        # print(4 * "jfsd78yf7823w78yhs81")
        # print(80 * "#")

        # nersc_host = os.environ["NERSC_HOST"]

        self.cores_per_node = 32  # <-- Cori Haswell

        self.default_sub_params = self.default_submission_parameters()


        self.job_data_dir = ""
        self.root_dir = root_dir
        self.queues = self.__queue_types__()
        self.job_state_keys = self.job_state_dict()
        self.job_queue_state_key = "STAT"  # COMBAK
        self.error_file = "job.err"
        self.out_file = "job.out"
        # __|

    def submit_job_clust(self, **kwargs):
        """Submit job to sherlck.

        Args:
            **kwargs:
        """
        # | - submit_job
        time.sleep(1.5)

        # | - Merging Submission Parameters
        params = merge_two_dicts(self.default_sub_params, kwargs)

        path = params["path_i"]

        # Fixing debug flag specification
        if params["priority"] == "debug":
            params["queue"] = "debug"

        if params["queue"] == "debug":
            params["priority"] = "debug"

        if params["priority"] == "scavenger":
            params["queue"] = "regular"
            params["qos"] = "scavenger"
            # params["queue"] = "regular"


            if params["wall_time"] < 180:

                params["min_time"] = params["wall_time"] - 30
            else:
                params["min_time"] = 180
            # --time-min=01:30:00

        # __|

        self.__import_arch_class__(params)

        self.arch_inst.__make_run_vasp_script__(params)

        params = merge_two_dicts(params, self.arch_inst.sbatch_params)

        # | - Write submission script

        knl_sub_script = """#!/bin/bash
# KNL submission bash script

module load vasp-tpc/5.4.4-knl

# cd $SLURM_SUBMIT_DIR
# export TMPDIR=$SLURM_SUBMIT_DIR

export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

export VASP_SCRIPT=./run_vasp.py
export VASP_PP_PATH=/project/projectdirs/m2997/vasp-psp/pseudo52

python ./model.py"""




        cori_sub_script = """#!/bin/bash
# Cori submission bash script

module load vasp-tpc/5.4.4-hsw

# cd $SLURM_SUBMIT_DIR
# export TMPDIR=$SLURM_SUBMIT_DIR

export VASP_SCRIPT=./run_vasp.py
export VASP_PP_PATH=/project/projectdirs/m2997/vasp-psp/pseudo52

python ./model.py"""

        if params["architecture"] == "knl":
            sub_script = knl_sub_script
        elif params["architecture"] == "haswell":
            sub_script = cori_sub_script
        else:
            # TEMP
            print("IJDIFJIJDSIFJs8df8sdd")
            sub_script = None

        with open(os.path.join(path, 'submit_job.pbs'), 'w') as the_file:
            the_file.write(sub_script)
        # __|


        # | - Submit Job
        os.chdir(path)
        if params["job_name"] == "Default":
            params["job_name"] = os.getcwd()

        print("submitting job")
        os.system("chmod 777 *")
        # bash_command = "/u/if/flores12/bin/qv model.py"
        # __| **** TEMP

        # | - Bash Submisssion Command
        # The -q flag is being used in place of the -p flag
        # Only the -q needs to be defined

        sbatch_to_param_key_dict = {
            "constraints": "-C",
            "out_file": "--output",
            "err_file": "--error",
            "wall_time": "--time",
            "nodes": "--nodes",
            "queue": "-q",
            "tasks-per-node": "--tasks-per-node",
            "account": "-A",
            "qos": "--qos",
            "min_time": "--time-min",
            }


        bash_command = "/usr/bin/sbatch "
        for key, value in params.items():
            if key in sbatch_to_param_key_dict:
                sbatch_key = sbatch_to_param_key_dict[key]
                bash_command += sbatch_key + " "
                bash_command += str(value) + " "
            else:
                tmp = 42

                # print("kjfs8usd9fusd09")
                # print(key, value)

        # bash_command += params["job_script"]
        bash_command += "submit_job.pbs"

        print("Bash Submission Command:")
        print(bash_command)
        # __|

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

        # | - Parsing Output
        out, err = output.communicate()
        out_copy = copy.deepcopy(out)

        try:
            out = out.strip()
            out_list = out.decode().split(" ")
            job_id = int(out_list[-1])

        except:
            print("Couldn't parse for jobid")
            job_id = None
        # __|

        # | - Writing Files
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
        # __|

        os.chdir(self.root_dir)
        # __|



    def default_submission_parameters(self):
        """Defaul SLURM parameters for Sherlock cluster."""
        # | - default_submission_parameters

        def_params = {
            "queue": "regular",  # -p flag | regular, debug
            ##SBATCH -p regular
            #SBATCH -p debug
            "nodes": "10",
            # 24 cpus per node on edison, 32 per node on cori haswell
            #SBATCH -N 10

            "account": "m2997",  # -A flag
            #SBATCH -A  m2997

            "wall_time": "180",
            #SBATCH -t 00:30:00


            "priority": "regular",  # --qos -q flag
            # "priority": "scavenger",  # --qos -q flag

            ##SBATCH --qos=scavenger
            ##SBATCH --qos=premium
            ##SBATCH --qos=debug
            ##SBATCH --qos=regular

            "constraints": "haswell",

            ##SBATCH -C haswell #this is for cori haswell (old)
            ##SBATCH -C knl #new cori
            #SBATCH -e job.err#SBATCH -o job.out

            "architecture": "knl",
            }

        return(def_params)
        # __|


    def __import_arch_class__(self, params):
        """
        """
        # | - __import_arch_class__
        arch = params["architecture"]

        if arch == "knl":
            arch_inst = KNL_Arch()
        elif arch == "haswell":
            arch_inst = Haswell_Arch()


        self.arch_inst = arch_inst
        # __|


    # | - out of sight

    def job_state_dict(self):
        """
        """
        # | - job_state_dict
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
        # __|


    def __queue_types__(self):
        """Queue types for Edison cluster
        """
        # | - __queue_types__
        queue_list = [
            "regular",
            "debug",
            "premium",
            ]

        return(queue_list)
        # __|

    def job_info_batch(self, job_id, path_i=None):
        """
        """
        # | - job_info_batch
        data_dict = slurm_squeue_parse(
            job_id,
            path_i=path_i,
            queue_state_key=self.job_queue_state_key,
            job_state_dict=self.job_state_keys,
            )

        return(data_dict)

        # | - __old__
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
        # __|

        # __|

    def completed_file(self, path_i="."):
        """Check whether ".FINISHED" file exists.

        Indicates that the job has gone to completion

        Args:
            path_i:
        """
        # | - completed_file
        completed_fle = False
        if os.path.exists(path_i + "/.FINISHED"):
            completed_fle = True

        return(completed_fle)
        # __|

    def job_state(self, path_i="."):
        """Return job state of path_i --> job_i.

        Args:
            path_i
        """
        # | - job_state
        job_id = self.get_jobid(path_i=path_i)

        job_state_out = None
        if job_id is not None:
            job_info = self.job_info_batch(job_id, path_i=path_i)

            if job_info is not None:
                key = self.job_queue_state_key
                if key in job_info:
                    job_state_out = job_info[key]
                    job_state_out = self.job_state_keys[job_state_out]

        # | - Checking for "completed" file indicating success
        completed_fle = self.completed_file(path_i=path_i)
        if completed_fle:
            job_state_out = self.job_state_keys["SUCCEEDED"]
        # __|

        return(job_state_out)
        # __|

    def get_jobid(self, path_i="."):
        """Return job ID of job_i.

        Args:
            path_i:
        """
        # | - get_jobid
        fileid_path = path_i + "/.jobid"
        if os.path.isfile(fileid_path):
            with open(path_i + "/.jobid") as fle:
                jobid = fle.read().strip()
        else:
            jobid = None

        return(jobid)
        # __|

    # __|

    # __| **********************************************************************





class KNL_Arch():
    """
    """

    # | - KNL_Arch *************************************************************
    def __init__(self):
        """Initialize Sherlock cluster instance.

        Args:
            root_dir:
        """
        # | - __init__
        self.vasp_module = "vasp-tpc/5.4.4-knl"

        # Actually 68
        self.cores_per_node = 66

        self.__module_load_vasp__()

        self.sbatch_params = self.__self_sbatch_settings__()

        os.putenv("OMP_NUM_THREADS", "1")
        os.putenv("OMP_PLACES", "threads")
        os.putenv("OMP_PROC_BIND", "spread")

        os.putenv("TMPDIR", "$SLURM_SUBMIT_DIR")
        os.putenv("VASP_SCRIPT", "./run_vasp.py")
        # __|

    def __make_run_vasp_script__(self, params):
        """
        """
        # | - __make_run_vasp_script__
        os.system("echo import os > run_vasp.py")
        exitcode_line = "exitcode = os.system('srun -n" + " " + \
            str(int(self.cores_per_node * int(params["nodes"]))) + " " + \
            "-c 4 --cpu_bind=cores" + " " + \
            "vasp_std')"
        line_2 = 'echo ' + '"' + exitcode_line + '" >> run_vasp.py'
        os.system(line_2)  # on edison


        if params["path_i"] == ".":
            pass

        else:
            copyfile(
                "run_vasp.py",
                os.path.join(params["path_i"], "run_vasp.py"),
                )

        # __|

    def __self_sbatch_settings__(self):
        """
        """
        # | - __self_sbatch_settings__
        def_params = {
            "constraints": "knl",
            "tasks-per-node": 66,
            }

        return(def_params)
        # __|

    def __module_load_vasp__(self):
        """
        THIS DOESN'T WORK
        """
        # | - __module_load_vasp__
        bash_comm = "module load " + self.vasp_module

        # TEMP
        print("fjjisdjfjisdfu8h3w8whgiafd6gwgundssoijgg")
        print(bash_comm)
        os.system(bash_comm)
        # __|

    # __| **********************************************************************




class Haswell_Arch():
    """
    """

    # | - Haswell_Arch *********************************************************
    def __init__(self):
        """Initialize Sherlock cluster instance.

        Args:
            root_dir:
        """
        # | - __init__
        self.vasp_module = "vasp-tpc/5.4.4-hsw"

        # Actually 68
        self.cores_per_node = 32

        self.sbatch_params = self.__self_sbatch_settings__()
        # __|

    def __self_sbatch_settings__(self):
        """
        """
        # | - __self_sbatch_settings__
        def_params = {
            "constraints": "haswell",
            # "tasks-per-node": 66,
            }

        return(def_params)
        # __|

    def __make_run_vasp_script__(self, params):
        """
        """
        # | - __make_run_vasp_script__
        os.system("echo import os > run_vasp.py")
        exitcode_line = "exitcode = os.system('srun -n" + " " + \
            str(int(self.cores_per_node * int(params["nodes"]))) + " " + \
            "vasp_std')"
        line_2 = 'echo ' + '"' + exitcode_line + '" >> run_vasp.py'
        os.system(line_2)  # on edison


        # __|

    # __| **********************************************************************
