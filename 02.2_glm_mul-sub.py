#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pathlib import Path
import argparse
import pandas as pd
from time import time

root_path = '/home/wanjiag/projects/GLACIER/'

# Get Slurm parameters
parser = argparse.ArgumentParser(description='Slurm parameters.')
parser.add_argument(
    '--partition', action='store', default='kuhl', help='partition requested'
)
parser.add_argument(
    '--time',
    '-t',
    action='store',
    default='0-23:59:59',
    help='time limit (dd-hh:mm:ss)'
)
parser.add_argument(
    '--nodes', '-N', action='store', default=1, help='number of nodes on which to run'
)
parser.add_argument(
    '--ntasks', action='store', default=1, help='Number of tasks to be launched per Node'
)
parser.add_argument(
    '--account',
    action='store',
    default='kuhl_lab',
    help='charge job to specified account'
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
parser.add_argument('--log-dir', action='store', default=None, help='log dir')

args = parser.parse_args()

# Directories
base_dir = Path(root_path)
deriv_dir = base_dir.joinpath('derivatives')
sbatch_dir = deriv_dir.joinpath('scripts', 'glm_sbatch')
sbatch_dir.mkdir(exist_ok=True, parents=True)

# Get subjects list
sub_info = pd.read_csv(deriv_dir.joinpath('scripts','glm_participants.tsv'), delimiter='\t')
sub_list = sub_info['participant_id'].tolist()

# Job information
job_name = 'glm'
if args.log_dir is None:
    log_dir = deriv_dir.joinpath('scripts', 'logs')
    log_dir.mkdir(exist_ok=True, parents=True)

# Create batch file for submission
cmd_file = sbatch_dir.joinpath('run_GLM_{}.sbatch'.format(time()))
cmd_str = ''
for sub_id in sub_list:
    cmd_str += (
        'sbatch '
        f'--partition={args.partition} '
        f'--job-name={job_name} '
        f'--output={log_dir}/%x_sub-{sub_id}_%j.log '
        f'--error={log_dir}/%x_sub-{sub_id}_%j.err '
        f'--time={args.time} '
        f'--nodes={args.nodes} '
        f'--ntasks-per-node={args.ntasks} '
        f'--account={args.account} '
        f'--mail-user={args.mail_user} '
        f'--mail-type={args.mail_type} '      
        
        '--wrap='
        '\"module load miniconda && '
        'conda activate /home/wanjiag/projects/environments/glacier_env && '
        f'python3.10 /gpfs/projects/kuhl_lab/wanjiag/GLACIER/derivatives/scripts/python_code/glm.py --subs {sub_id}\"\n'
    )
cmd_file.write_text(cmd_str)
