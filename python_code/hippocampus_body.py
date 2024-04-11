LEFT = {'sub-GLACIER01':[33,8],
        'sub-GLACIER02':[31,10],
        'sub-GLACIER03':[32,10],
        'sub-GLACIER04':[33,8],
        'sub-GLACIER06':[31,8],
        'sub-GLACIER07':[31,10],
        'sub-GLACIER08':[34,10],
        'sub-GLACIER10':[31,10],
        'sub-GLACIER11':[32,9],
        'sub-GLACIER12':[32,11],
        'sub-GLACIER13':[25,10],
        'sub-GLACIER14':[32,10],
        'sub-GLACIER16':[33,9],
        'sub-GLACIER18':[34,10]}

RIGHT = {'sub-GLACIER01':[33,8],
         'sub-GLACIER02':[31,10],
         'sub-GLACIER03':[32,10],
         'sub-GLACIER04':[34,7],
         'sub-GLACIER06':[31,8],
         'sub-GLACIER07':[33,8],
         'sub-GLACIER08':[35,9],
         'sub-GLACIER10':[31,10],
         'sub-GLACIER11':[31,10],
         'sub-GLACIER12':[34,9],
         'sub-GLACIER13':[26,9],
         'sub-GLACIER14':[33,9],
         'sub-GLACIER16':[32,10],
         'sub-GLACIER18':[33,11]}

from glob import glob
import os 
from os.path import join as opj
import subprocess

from typing import Union, List
from shlex import split
import sys

def run_cmd(
    cmd: Union[str, list[str]],
    print_output: bool = True,
    shell: bool = False,
    check: bool = True,
    **kwargs,
) -> subprocess.CompletedProcess:
    """Executes command in Shell.

    Args:
        cmd: Command to be executed in external shell. It could be a
            string or a list of command parts (see subprocess function
            'run' for details).
        print_output: If true, print out the shell outputs.
        shell: If true, the command will be executed through the shell
            (see subprocess doc for details).
        check: If check is true, and the process exits with a non-zero
            exit code, a CalledProcessError exception will be raised.
            Attributes of that exception hold the arguments, the exit
            code, and stdout and stderr if they were captured.
        **kwargs: Additional keyword arguments pass to function 'run'.

    Returns:
        A subprocess.CompletedProcess object.
    """
    print(cmd, flush=True)
    try:
        if shell:
            if isinstance(cmd, list):
                cmd = " ".join(cmd)
            res = subprocess.run(
                cmd, shell=True, capture_output=True, check=check, encoding="utf-8", **kwargs
            )
        else:
            if isinstance(cmd, str):
                cmd = split(cmd)
            res = subprocess.run(cmd, capture_output=True, check=check, encoding="utf-8", **kwargs)
        if print_output:
            if res.stdout != "":
                print(res.stdout.rstrip("\n"), flush=True)
            if res.stderr != "":
                print(res.stderr, flush=True)
    except subprocess.CalledProcessError as e:
        print(e.stdout)
        print(e.stderr)
        sys.exit(1)
    return res

def sh(script):
    os.system("bash -c '%s'" % script)
    
def hippo_body(left, right, t2_file, subfield, subnum, output_dir, temp_dir):
    
    empty = opj(temp_dir, 'empty_t2.nii.gz')
    cmd = f'fslmaths {t2_file} -sub {t2_file} {empty}'
    run_cmd(cmd)
    
    # left
    
    left_body = opj(temp_dir, f'{subfield}_left-body.nii.gz')
    cmd = f'fslroi {left} {left_body} 0 -1 0 -1 {LEFT[subnum][0]} {LEFT[subnum][1]}'
    run_cmd(cmd)
    
    front = LEFT[subnum][0]-1
    end = LEFT[subnum][0]+LEFT[subnum][1]
    end_length = 65 - end + 1
    
    front_empty = opj(temp_dir, f'{subfield}_left-front-empty.nii.gz')
    cmd = f'fslroi {empty} {front_empty} 0 -1 0 -1 0 {front}'
    run_cmd(cmd)    
    
    end_empty = opj(temp_dir, f'{subfield}_left-end-empty.nii.gz')
    cmd = f'fslroi {empty} {end_empty} 0 -1 0 -1 {end} {end_length}'
    run_cmd(cmd)    
    
    left_t2 = opj(temp_dir, f'{subfield}_left-body_t2.nii.gz')
    cmd = f'fslmerge -z {left_t2} {front_empty} {left_body} {end_empty}'
    run_cmd(cmd)        
    
    # right
    
    right_body = opj(temp_dir, f'{subfield}_right-body.nii.gz')
    cmd = f'fslroi {right} {right_body} 0 -1 0 -1 {RIGHT[subnum][0]} {RIGHT[subnum][1]}'
    run_cmd(cmd)
    
    front = RIGHT[subnum][0]-1
    end = RIGHT[subnum][0]+RIGHT[subnum][1]   
    
    front_empty = opj(temp_dir, f'{subfield}_right-front-empty.nii.gz')
    cmd = f'fslroi {empty} {front_empty} 0 -1 0 -1 0 {front}'
    run_cmd(cmd)  
    
    end_empty = opj(temp_dir, f'{subfield}_right-end-empty.nii.gz')
    cmd = f'fslroi {empty} {end_empty} 0 -1 0 -1 {end} {end_length}'
    run_cmd(cmd)    
    
    right_t2 = opj(temp_dir, f'{subfield}_right-body_t2.nii.gz')
    cmd = f'fslmerge -z {right_t2} {front_empty} {right_body} {end_empty}'
    run_cmd(cmd)
    
    # combined roi
    combined_file = opj(temp_dir, f'{subfield}-body_T2w.nii.gz')
    cmd = f'fslmaths {left_t2} -add {right_t2} -bin {combined_file}'
    run_cmd(cmd)
    
    return combined_file
