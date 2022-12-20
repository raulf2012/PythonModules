# #########################################################
# SLURM python methods.
# #########################################################


#| - Import modules
import os
import numpy as np
#__|



def get_sbatch_command(
    walltime=None,
    nodes=None,
    job_path='.',
    save_sbatch_to_file=False,
    sbatch_filename='sbatch.sh',
    submit_script_filename='submit.sh',
    ):
    """
    Args:
      walltime
        - Walltime in minutes
    """
    #| - get_sbatch_command
    hours = int(np.floor(walltime / 60))
    minutes = int(np.floor(walltime - 60 * hours))

    slurm_key_vals = {
        'walltime': walltime,
        'nodes': nodes,
        }

    slurm_keys = {
        'walltime': 'time',
        'nodes': 'nodes',
        }



    # #####################################################
    # Constructing sbatch command

    #  sbatch_str = 'sbatch ' + \
    #      '--time=' + str(hours).zfill(2) + ':' + str(minutes).zfill(2) + ':00' + ' ' + \
    #      '--nodes=' + str(nodes) + ' ' + \
    #      submit_script_filename + \
    #      '\n'

    sbatch_str = 'sbatch '
    for key, val in slurm_key_vals.items():

        if key == 'walltime':
            walltime = val
            hours = int(np.floor(walltime / 60))
            minutes = int(np.floor(walltime - 60 * hours))

            sbatch_seg = '' + \
                '--' + slurm_keys[key] + \
                '=' + str(hours).zfill(2) + \
                ':' + str(minutes).zfill(2) + ':00' + ' '

        elif key == 'TEMP':
            tmp = 42
        else:
            sbatch_seg = '' + \
                '--' + slurm_keys[key] + \
                '=' + str(val) + ' '

        sbatch_str += sbatch_seg

    sbatch_str += submit_script_filename
    # #####################################################


    if save_sbatch_to_file:
        with open(os.path.join(job_path, sbatch_filename), 'w') as file:
            file.write(sbatch_str)

    return(sbatch_str)
    #__|

