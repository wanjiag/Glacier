#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Run fMRIPrep."""

## Info----------------------------------------------------------------------------
##
## Author: Zhifang Ye
##
## Email: zhifang.ye.fghm@gmail.com
##
## Date Created: 2019-08-26
##
## Last Updated: 2019-10-12
##
## Notes:
## Edited by Wanjia Guo 02/25/2022

from pathlib import Path
import argparse
import pandas as pd
from time import time

root_path = '/home/wanjiag/projects/GLACIER/'

# Get Slurm parameters
parser = argparse.ArgumentParser(description='Slurm parameters.')
parser.add_argument(
    '--nodes', '-N', action='store', default=1, help='number of nodes on which to run'
)
parser.add_argument(
    '--cpus-per-task',
    '-c',
    action='store',
    default=16,
    type=int,
    help='number of cpus required per task'
)
parser.add_argument(
    '--mem',
    action='store',
    default='156G',
    #type=int,
    help='minimum amount of real memory'
)
parser.add_argument(
    '--time',
    '-t',
    action='store',
    default='0-23:59:59',
    help='time limit (dd-hh:mm:ss)'
)
parser.add_argument('--log-dir', action='store', default=None, help='log dir')
parser.add_argument(
    '--account',
    action='store',
    default='kuhl_lab',
    help='charge job to specified account'
)
parser.add_argument(
    '--partition', action='store', default='kuhl', help='partition requested'
)
parser.add_argument(
    '--mail-user',
    action='store',
    default='wanjiag@uoregon.edu',
    help='who to send email notification for job state changes'
)
parser.add_argument(
    '--mail-type',
    action='store',
    default='END',
    help='notify on state change: BEGIN, END, FAIL or ALL'
)
args = parser.parse_args()

# Directories
base_dir = Path(root_path)
deriv_dir = base_dir.joinpath('derivatives')
fmriprep_dir = deriv_dir.joinpath('fmriprep')
sbatch_dir = deriv_dir.joinpath('scripts/fmriprep', 'sbatch')
sbatch_dir.mkdir(exist_ok=True, parents=True)

# Get subjects list
sub_info = pd.read_csv(deriv_dir.joinpath('scripts','fmriprep','participants.tsv'), delimiter='\t')
sub_list = [i.replace('sub-', '') for i in sub_info['participant_id'].tolist()]

# Job information
job_name = 'fMRIPrep'
if args.log_dir is None:
    log_dir = deriv_dir.joinpath('scripts', 'fmriprep','log')
    log_dir.mkdir(exist_ok=True, parents=True)

# Create batch file for submission
cmd_file = sbatch_dir.joinpath('run_fMRIPrep_{}.sbatch'.format(time()))
cmd_str = ''
for sub_id in sub_list:
    cmd_str += (
        'sbatch '
        f'--nodes={args.nodes} '
        f'--cpus-per-task={args.cpus_per_task} '
        f'--mem={args.mem} '
        f'--time={args.time} '
        f'--account={args.account} '
        f'--partition={args.partition} '
        f'--mail-user={args.mail_user} '
        f'--mail-type={args.mail_type} '
        f'--job-name={job_name} '
        f'--output={log_dir}/%x_sub-{sub_id}_%j.log '
        '--wrap='
        '\"module load singularity && '
        'singularity run -e '
	'--bind /projects/kuhl_lab/wanjiag:/projects/kuhl_lab/wanjiag '
        '-B /gpfs/projects/kuhl_lab/wanjiag/GLACIER:/data '
        '-B /gpfs/projects/kuhl_lab/wanjiag/GLACIER/derivatives/fmriprep:/output '
        '-B /gpfs/projects/kuhl_lab/wanjiag/GLACIER/derivatives/fmriprep_tmp:/work '
        '/projects/kuhl_lab/shared/fmriprep-v23.2.0.simg '
        '/data /output participant '
        f'--participant_label {sub_id} '
        f'--nthreads {args.cpus_per_task} --omp-nthreads 16 --mem_mb {args.mem} '
        '--output-spaces MNI152NLin2009cAsym T1w fsaverage6 '
        '--fs-license-file /projects/kuhl_lab/wanjiag/NEUDIF/scripts/license.txt '
        f'-w /work --stop-on-first-crash --notrack\"\n'
    )
cmd_file.write_text(cmd_str)
