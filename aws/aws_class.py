"""Class for performing AWS operations."""

#| - Import Modules
from dft_job_automat.job_setup import DFT_Jobs_Setup
import os
import errno
import subprocess
from copy import deepcopy
import json
import datetime
import shutil
import pandas as pd

import boto3
#__|


def force_symlink(file1, file2):
    #| - force_symlink
    try:
        os.symlink(file1, file2)
    except OSError, e:
        if e.errno == errno.EEXIST:
            os.unlink(file2)
            os.symlink(file1, file2)
    #__|


class AWS_Queues():
    """
    """
    #| - AWS_Class *****************************************************************
    def __init__(self):
        """
        """
        #| - __init__
        try:
            self.aws_dir = os.environ["TRI_PATH"]
            self.job_queue_dir = self.aws_dir + "/bin/job_queues"
        except:
            pass

        #__|


    def job_info_batch(self, job_id):
        """
        """
        #| - job_info_batch
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

    def cancel_job(self, job_id, reason="N/A"):
        """
        """
        #| - cancel_job
        bash_command = "aws batch cancel-job --job-id "
        bash_command += job_id + " --reason " + reason

        print("Cancelling job | " + job_id)
        run_bash = subprocess.check_output(bash_command, shell=True)

        #__|

    def cancel_jobs(self, job_id_list):
        """Cancels all jobs in a given job id list
        """
        #| - cancel_jobs
        for job in job_id_list:
            self.cancel_job(job)
        #__|

    def list_jobs(self, queue="small"):
        """
        Returns all jobs for a particular queue.
        Included jobs in SUBMITTED, PENDING, RUNNABLE, STARTING, RUNNING, and
        FAILED states.
        """
        #| - list_jobs
        batch = boto3.client("batch")
        job_status_opt = ["SUBMITTED", "PENDING", "RUNNABLE",
                          "STARTING", "RUNNING", "SUCCEEDED", "FAILED"]

        job_id_list = []
        for status in job_status_opt:
            all_jobs = batch.list_jobs(jobQueue=queue, jobStatus=status)  # TEMP
            job_ids = [i["jobId"] for i in all_jobs["jobSummaryList"]]

            #| - Checking that queue is being used
            print(str(len(job_ids)) + " jobs in the " + status + " queue")
            #
            # if len(job_ids) < 1:
            #     print("No jobs in the " + status + " queue")
            #__|


            #| - Retreiving Job ID's From AWS
            job_ids_split = [job_ids[x:x+100] for x in xrange(0, len(job_ids), 100)]
            for j in range(len(job_ids_split)):
                job_descriptions = batch.describe_jobs(jobs=job_ids_split[j])
                for i in job_descriptions["jobs"]:
                    job_id_list.append({
                    "job_id": i["jobId"],
                    "job_name": i["jobName"],
                    "model": i["parameters"]["model"],
                    "job_state": status,
                    })
            #__|

        return(job_id_list)

        #__|

    def submit_job(self, path=None, queue="test", cpus="default", copy_PythonModules=True):
        """
        Submits job to aws cluster. Copies PythonModules folder into
        job directory
        """
        #| - submit_job
        root_dir = os.getcwd()
        if path == None:
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


        #| - Copy PYTHONMODULES to Job Directory
        if os.path.isdir(path + "/PythonModules") == True:
            print("PythonModules already exists, erasing and recopying")
            shutil.rmtree(path + "/PythonModules")
            py_mod = os.environ["python_modules"]
            shutil.copytree(py_mod, path + "/PythonModules")
        else:
            py_mod = os.environ["python_modules"]
            shutil.copytree(py_mod, path + "/PythonModules")
        #__|


        #| - Submit Job
        os.chdir(path)

        if os.path.isfile(".sSUBMITTED"):
            print("Directory already submitted, will be skipped")
            os.chdir(root_dir)
            return(None)
        else:
            print("submitting job")
            aws_dir = os.environ["aws_dir"]


            if cpus == "default":
                bash_command = aws_dir + "/matr.io/bin/trisub -q " + queue
            else:
                bash_command = aws_dir + "/matr.io/bin/trisub -c " + str(cpus) + " -q " + queue


            # bash_command = aws_dir + "/matr.io/bin/trisub " + "-c" +  + "-q " + queue


            try:
                output = subprocess.check_output(bash_command, shell=True)
                sub_time = datetime.datetime.now().isoformat()
            except subprocess.CalledProcessError, e:
                print "Ping stdout output:\n", e.output

                os.chdir(root_dir)
                print("JOB SKIPPED: ")
                return(None)

        #__|


        #| - Parsing Submission for Job ID
        output = output.splitlines()
        for line in output:
            if "jobId" in line:
                lst = line.split('"')
                job_id_ind = (lst.index("jobId") + 2)
                jobId = lst[job_id_ind]
        #__|

        file = open(".submitted", "w")
        file.close()

        file = open("job_id", "w")
        file.write(jobId + "\n")

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

    #__| ***************************************************************************


class AWS_Class_tmp():
    """Summary line.

    TEMP
    """
    #| - AWS_Class_tmp *************************************************************
    def __init__(self):
        """TMP_docstring.

        TEMP TEMP
        """
        #| - __init__
        self.tmp = "TMP AWS Class Atribute!!!!!"

        try:
            self.TRI_PATH = os.environ["TRI_PATH"]
        except:
            pass
        #__|


    def create_symlinks(self):
        """
        Attempts to create symlinks in all subfolders whose names starts with _<num>
        """
        #| - create_symlinks
        # root_dir = '.'
        root_dir = os.getcwd()
        for dir_name, subdirList, fileList in os.walk(root_dir):
            dir_name_i = dir_name.split("/")[-1]

            if dir_name_i[0] == "_" and dir_name_i[1].isdigit() == True:
                # print(dir_name)
                dest = dir_name + "/simulation"

                source = deepcopy(dir_name)
                source = source.replace("/model/", "/simulation/", 1)

                print(source)
                print(dest)
                print("____")
                force_symlink(source, dest)

            # print("Found directory: %s" % dirName)
            # for fname in fileList:
                # print('\t%s' % fname)
        #__|


    def create_symlink(self, path=None):
        """
        """
        #| - create_symlink

        # def force_symlink(file1, file2):
        #     #| - force_symlink
        #     try:
        #         os.symlink(file1, file2)
        #     except OSError, e:
        #         if e.errno == errno.EEXIST:
        #             os.unlink(file2)
        #             os.symlink(file1, file2)
        #     #__|

        if path == None:
            path = os.getcwd()
        else:
            pass

        dest = path + "/simulation"
        print(dest)

        source = deepcopy(path)
        source = source.replace("/model/", "/simulation/", 1)
        force_symlink(source, dest)

        #__|

    #__| ***************************************************************************

# NOT USEFUL FOR THIS TO BE CHILD CLASS -- FIX!!!!!!!!!!!!
class AWS_Class(DFT_Jobs_Setup):
    """Summary line.

    TEMP
    """
    #| - AWS_Class *****************************************************************
    def __init__(self, system="sherlock"):
        """TMP_docstring.

        TEMP TEMP
        """
        #| - __init__
        DFT_Jobs_Setup.__init__(self, system=system)
        self.tmp = "TMP AWS Class Atribute!!!!!"

        try:
            self.TRI_PATH = os.environ["TRI_PATH"]
        except:
            pass
        #__|


    def create_symlinks(self):
        """Creates symlinks for every job revision in job directory.

        TEMP TEMP
        """
        #| - create_symlinks

        def force_symlink(file1, file2):
            #| - force_symlink
            try:
                os.symlink(file1, file2)
            except OSError, e:
                if e.errno == errno.EEXIST:
                    os.unlink(file2)
                    os.symlink(file1, file2)
            #__|

        for job in self.job_var_lst:
            path = self.var_lst_to_path(job)
            path = self.root_dir + "/" + path

            for folder in os.walk(path).next()[1]:
                path_i = path + folder
                dest = path_i + "/simulation"

                source = deepcopy(path_i)
                source = source.replace("/model/", "/simulation/", 1)

                force_symlink(source, dest)
        #__|

    #__| ***************************************************************************





#| - DEPRECATED METHODS

# def list_jobs(self, status="RUNNING", queue="small"):
#     """
#     """
#     #| - list_jobs
#     bash_command = "aws batch list-jobs --job-status {} --job-queue {} > aws_list_jobs.json".format(status, queue)
#     run_bash = subprocess.check_output(bash_command, shell=True)  #TEMP
#     data = json.load(open("aws_list_jobs.json"))
#
#     job_id_list = []
#     for job in data["jobSummaryList"]:
#         job_id_list.append(job["jobId"])
#
#     return(job_id_list)
#
#     #__|


#__|
