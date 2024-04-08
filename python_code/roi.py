from glob import glob
import os 
from os.path import join as opj
import subprocess

from typing import Union, List
from shlex import split
import sys

print(sys.version)

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
    

epi_mask_threshold = 0.5
singularity_prefix = 'singularity exec --bind /projects/kuhl_lab/wanjiag/GLACIER/derivatives:/projects/kuhl_lab/wanjiag/GLACIER/derivatives /gpfs/projects/kuhl_lab/shared/fmriprep-v23.2.0.simg'

######## Running ########

derivative_dir = '/projects/kuhl_lab/wanjiag/GLACIER/derivatives/'
glm_base_dir = opj(derivative_dir, 'glm/')
fmriprep_base_dir = opj(derivative_dir, 'fmriprep/')
roi_base_dir = opj(derivative_dir, 'roi/')
automatic_detecting_subjects = True

if automatic_detecting_subjects:
    f_list = glob(os.path.join(glm_base_dir, '*sub-GLACIER*/'))
    subs = list(map(lambda f: f[len(os.path.commonpath(f_list))+1:-1], f_list))
    subs.sort()
    
    processed_list = glob(os.path.join(roi_base_dir, '*sub-GLACIER*/'))
    if len(processed_list) == 0:
        todo_subs = subs
    else:
        processed_subs = [x.split('/')[-2] for x in processed_list]
        todo_subs = [x for x in subs if x not in processed_subs]
    
    print(todo_subs)
    
for sub in todo_subs:
    output_dir = opj(roi_base_dir, sub)

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    fmriprep_dir = opj(fmriprep_base_dir, sub)
    glm_dir = opj(glm_base_dir, sub)
    
    print(f'--------------------{sub}-------------------')
    print(output_dir)
    
    # Getting brain mask & transformation matrix & fmriprep functional file
    brain_mask = opj(glm_dir, 'temporary_files', f'{sub}_space-T1w_desc-brain_intersect_mask.nii.gz')
    h5 = opj(fmriprep_dir, f'anat/{sub}_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5')
    func_file_list = [x for x in glob(opj(fmriprep_dir, 'func', f'{sub}_task-glacier_run-*_space-T1w_desc-preproc_bold.nii.gz'))] 
    func_file_list.sort()
    if len(func_file_list) != 8:
        print(f'--------------------Some thing is wrong with {sub} functional runs-------------------')
        break

    '''
    # Finding and coverting aparc into nifti file.
    mgz_file = opj(fmriprep_base_dir, 'sourcedata', 'freesurfer', sub, 'mri', 'aparc.a2009s+aseg.mgz')
    aparc_2009 = opj(output_dir, f'{sub}_aparc.a2009s+aseg.nii.gz')
    print(mgz_file)
    print(aparc_2009)
    cmd = f'mri_convert {mgz_file} {aparc_2009}'
    run_cmd(cmd)
    
    # Finding and copying aparcaseg file into ROIs folder
    aparc_file = os.path.join(fmriprep_base_dir, 'sourcedata', 'freesurfer', sub, 'mri','aparc+aseg.mgz')
    aparc = os.path.join(output_dir,  f'{sub}_aparc+aseg.nii.gz')
    print(aparc_file)
    print(aparc)
    cmd = f'mri_convert {aparc_file} {aparc}'
    run_cmd(cmd)
    '''
    
    # Calculate All runs average as a reference image
    mean_out_file = opj(output_dir, f'{sub}_space-T1w_desc-preproc_bold_mean_all.nii.gz')
    if not os.path.exists(mean_out_file): 
        add_string = ''
        mean_func_list = []
        for func_file in func_file_list:
            # calculate Tmean for each functional file
            out_file = opj(output_dir, f'{os.path.basename(func_file).split(".")[0]}_mean.nii.gz')
            cmd = f'fslmaths {func_file} -Tmean {out_file}'
            run_cmd(cmd)
            mean_func_list.append(out_file)

        # calculate mean across all functional runs
        cmd = f'fslmaths {mean_func_list[0]} '
        add_string = ''
        for mean_file in mean_func_list[1:]:
                add_string += f'-add {mean_file} '
        mask_string = f'-mas {brain_mask}'
        div_string = f'-div {len(mean_func_list)}'
        run_cmd(f'{cmd} {add_string} {div_string} {mask_string} {mean_out_file}')

        # Remove temporary files
        for mean_func in mean_func_list:
            sh(f'rm {mean_func}')
            
    # PPA from MNI space
    print('--------------------PPA-------------------')
    ppa_mni = '/home/wanjiag/projects/GLACIER/derivatives/roi/mni/ppa/ppa.nii.gz'
    ppa_out = opj(output_dir, 'ppa_mni-2-epi.nii.gz')
    cmd = f'{singularity_prefix} antsApplyTransforms -d 3 -i {ppa_mni} -r {mean_out_file} -t {h5} -f 0 -o {ppa_out}'
    run_cmd(cmd)
    # threshold & bin
    ppa_final_out = opj(output_dir, 'ppa_mni-2-epi_thr-0.5_masked_bin.nii.gz')
    cmd = f'fslmaths {ppa_out} -thr {epi_mask_threshold} -mas {brain_mask} -bin {ppa_final_out}'
    run_cmd(cmd)
    # remove temporary file
    sh(f'rm {ppa_out}')

    # EVC from MNI space
    print('--------------------EVC-------------------')
    ev_mni_path = '/home/wanjiag/projects/GLACIER/derivatives/roi/mni/visual_cortex/subj_vol_all'
    ev_files = ['perc_VTPM_vol_roi1_lh.nii.gz',
                'perc_VTPM_vol_roi1_rh.nii.gz',
                'perc_VTPM_vol_roi2_lh.nii.gz',
                'perc_VTPM_vol_roi2_rh.nii.gz']

    ev_file_threshold = 50
    ev_tmp = []
    ev_output = []

    for ev_file in ev_files:
        at_out_file = opj(output_dir, '{}-2-epi.nii.gz'.format(ev_file.split('.')[0]))
        cmd = f'{singularity_prefix} antsApplyTransforms -d 3 -i {opj(ev_mni_path,ev_file)} -r {mean_out_file} -t {h5} -f 0 -o {at_out_file}'
        run_cmd(cmd)
        ev_tmp.append(at_out_file)

        trh_out_file = opj(output_dir, '{}-2-epi_thr-{}.nii.gz'.format(ev_file.split('.')[0], ev_file_threshold))

        cmd = f'fslmaths {at_out_file} -thr {ev_file_threshold} -bin {trh_out_file}'
        run_cmd(cmd)

        ev_output.append(trh_out_file)


    # calculate mean across all functional runs
    cmd = f'fslmaths {ev_output[0]} '
    add_string = ''
    ev_out = opj(output_dir, f'evc-2-epi_thr-{ev_file_threshold}_masked_bin.nii.gz')
    for ev_file in ev_output[1:]:
            add_string += f'-add {ev_file} '
    mask_string = f'-mas {brain_mask}'
    run_cmd(f'{cmd} {add_string} {mask_string} -bin {ev_out}')

    # Remove temporary files
    for ev_file in ev_output:
        sh(f'rm {ev_file}')

    for ev_file in ev_tmp:
        sh(f'rm {ev_file}')