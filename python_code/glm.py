from os.path import join as opj
import glob
import re
from glm_utils import *
import argparse

parser = argparse.ArgumentParser(description="Parameters.")
parser.add_argument(
    "-s",
    "--subs",
    action="store",
    nargs="*",
    help=(
        "One or more subject identifiers (e.g., sub-GLACIER01)."
        "If this is omitted, using a pre-defined subject_list."
    ),
)
args = parser.parse_args()
# yapf: enable

confounds_list = [
    "trans_x", "trans_y", "trans_z",
    "rot_x", "rot_y", "rot_z",
    "framewise_displacement",
    'a_comp_cor_00', 'a_comp_cor_01','a_comp_cor_02','a_comp_cor_03','a_comp_cor_04','a_comp_cor_05', 
    'csf']
repetition_time = 2.
fd_thresh = 0.5
std_dvars_thresh = 1.5
num_TRs = 168
# Onset time offset (adjust for slice timing correction)
# see https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dTshift.html
# and fMRIPrep implementation https://github.com/nipreps/fmriprep/blob/7f68f590d3abfb7800067a88d74c4f18ac512c29/fmriprep/workflows/bold/stc.py#L100
# tzero = np.round(first + frac * (last - first), 3)
# first and last timing could be found in SliceTiming field in raw data json file
stim_times_subtract = np.round(0 + 0.5 * (1.9 - 0), 3)

derivative_dir = '/projects/kuhl_lab/wanjiag/GLACIER/derivatives/'
fmriprep_dir = opj(derivative_dir, 'fmriprep/')
behav_dir = opj(derivative_dir, 'behavior/')

for sub in args.subs:
    sub_id = re.findall(r'\d', sub)
    sub_id = "".join(list(map(str, sub_id)))

    behav_file_list =  [x for x in glob.glob(opj(behav_dir, f'sub{sub_id}',f'sub{sub_id}_exposure*_behav_*.csv'))]
    behav_file_list.sort()
    timing_file_list = [x for x in glob.glob(opj(behav_dir, f'sub{sub_id}', f'sub{sub_id}_exposure*_timing_*.csv'))]
    timing_file_list.sort()

    bold_dir=opj(fmriprep_dir, f'{sub}/func/')
    out_dir = opj(derivative_dir, f'glm/{sub}')

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    tmp_dir = opj(out_dir, 'temporary_files')
    if not os.path.isdir(tmp_dir):
        os.makedirs(tmp_dir)

    confound_dir = opj(out_dir, 'confound_files')
    if not os.path.isdir(confound_dir):
        os.makedirs(confound_dir)

    xmat_dir = opj(out_dir, 'xmat_files')
    if not os.path.isdir(xmat_dir):
        os.makedirs(xmat_dir)

    bucket_dir = opj(out_dir, 'bucket_files')
    if not os.path.isdir(bucket_dir):
        os.makedirs(bucket_dir)
        
    glm_dir = opj(out_dir, 'glm_files')
    if not os.path.isdir(glm_dir):
        os.makedirs(glm_dir)

    #trl_mod_dir = Path(tmp_dir).joinpath("bucket_files")
    #trl_mod_dir.mkdir(exist_ok=True, parents=True)

    func_file_list = [x for x in glob.glob(opj(bold_dir, f'{sub}_task-glacier_run-*_space-T1w_desc-preproc_bold.nii.gz'))] 
    func_file_list.sort()

    # making general mask
    mask_file_list = [x for x in glob.glob(opj(bold_dir, f'{sub}_task-glacier_run-*_space-T1w_desc-brain_mask.nii.gz'))]
    mask_file_list.sort()

    mask_intersect = opj(tmp_dir, f'{sub}_space-T1w_desc-brain_intersect_mask.nii.gz')
    if not os.path.exists(mask_intersect): 
        intersect_masks(mask_file_list, mask_intersect)
        
    # fmriprep confound files
    confound_file_list = [x for x in glob.glob(opj(bold_dir, f'{sub}_task-glacier_run-*_desc-confounds_timeseries.tsv'))]
    confound_file_list.sort()

    confounds = []
    for f in confound_file_list:
        confounds_run = pd.read_csv(f, sep="\t")

        # add useful information
        confounds_run["sub"] = sub
        confounds_run["run_id"] = int(f.split('/')[-1].split('_')[2].split('-')[-1])
        col_list = confounds_run.columns.tolist()
        confounds_run = confounds_run.loc[:, col_list]
        confounds.append(confounds_run)

    try:
        confounds = pd.concat(confounds).reset_index(drop=True)
    except InvalidIndexError as e:
        print(
            "This error is likely due to the different number of *_comp_cor columns in "
            "different runs.\nIf that's the case, read each run separately or set "
            "'exclude_compcor' to True."
        )
        raise e

    confounds_cur_list = [
                "trans_x",
                "trans_y",
                "trans_z",
                "rot_x",
                "rot_y",
                "rot_z",
                "framewise_displacement",
                "std_dvars",
                "dvars",
                "rmsd",
                "global_signal",
                "csf",
                "white_matter",
                "csf_wm",
                'a_comp_cor_00', 'a_comp_cor_01','a_comp_cor_02','a_comp_cor_03','a_comp_cor_04','a_comp_cor_05',
            ]
    if confounds_cur_list:
        confounds = confounds.loc[:, ["sub", "run_id"] + confounds_cur_list]

    # Make motion confounds regressor file
    confounds_file = make_confounds_regressor(df = confounds,
        out_dir = Path(confound_dir),
        demean = True,
        split_into_runs = True,
        confounds_list = confounds_list)

    for f in confounds_file:
        # check and remove allzero confounds column
        _, bad_col_idx = remove_allzero_column(f)
        if len(bad_col_idx) != 0:
            print(f"All zero column: {' '.join(list(np.array(f)[bad_col_idx]))}")

    # Make good TR file
    goodtr_file, _ = make_good_tr_regressor(
        confounds,
        Path(confound_dir),
        fd_thresh=fd_thresh,
        std_dvars_thresh=std_dvars_thresh,
        censor_prev_tr=False,
    )
    
    for func_id, func_file in enumerate(func_file_list):

        # Calculate run length
        run_length = calc_run_length(Path(func_file), cifti=True)[0]
        if run_length != num_TRs:
            break
            print(f'==============={func_file} has incorrect number of TRs, please check raw files===============')

        # Scaling functional file
        fname_prefix = '_'.join(func_file.split('/')[-1].split('_')[:4])
        scaled_file = scale_func_image(
                func_file, Path(tmp_dir, f"{fname_prefix}_scaled.nii.gz"), mask_file = mask_intersect, cifti=False
            )
        
        timing_file = timing_file_list[func_id]
        timing = pd.read_csv(timing_file)

        behav_file = behav_file_list[func_id]
        behav = pd.read_csv(behav_file)

        curr_confounds_file = confounds_file[func_id]
        curr_goodtr_file = goodtr_file[func_id]

        print(func_file)
        print(timing_file)
        print(behav_file)
        print(curr_confounds_file)
        print(curr_goodtr_file)
        
        run_id = int(re.findall(r'\d', os.path.basename(scaled_file).split('_')[2])[0])
        behav_id = int(re.findall(r'\d', os.path.basename(behav_file).split('_')[1])[0]) + 1 #behav start from 0
        timing_id = int(re.findall(r'\d', os.path.basename(timing_file).split('_')[1])[0]) + 1 #behav start from 0
        confound_id = int(re.findall(r'\d', os.path.basename(curr_confounds_file).split('_')[0])[0])
        goodtr_id = int(re.findall(r'\d', os.path.basename(curr_goodtr_file).split('_')[0])[0])

        if run_id != behav_id or run_id != timing_id or run_id != confound_id or run_id != goodtr_id:
            print(f'==============={func_file} has incorrect files: index inconsistent===============')
            print(f'run_id:{run_id}, behav_id:{behav_id}, timing_id:{timing_id}, confound_id:{confound_id}, goodtr_id:{goodtr_id}')
            break
            
        # Making design matrix
        design_matrix = timing.loc[(timing['stimli'] != 'fixation') & (timing['stimli'] != 'resp-6')].reset_index(drop=True)
        design_matrix = design_matrix.drop_duplicates(subset=['stimli']).loc[:,['stimli', 'design_onset']]
        trl_list = design_matrix["stimli"].values
        time_list = design_matrix['design_onset'].values
        
        for trl_id in trl_list:
            if 'lure' in trl_id:
                continue
                print(f"\n########## Skipped fitting model for {run_id} trial: {trl_id} ##########",
                    flush=True,)
            print(f"\n########## Fitting model for {run_id} trial: {trl_id} ##########",
                    flush=True,)
            # Make main regressors onset time file
            onset_file, _ = make_singletrial_onset_time(
                design_matrix, trl_id, Path(xmat_dir), prefix=fname_prefix
            )
            # Make model design matrix
            xmat_file = make_singletrial_design_matrix(
                onset_file,
                trl_id,
                Path(xmat_dir),
                run_length,
                repetition_time,
                confounds_file=curr_confounds_file,
                goodtr_file=curr_goodtr_file,
                stim_times_subtract=stim_times_subtract,
                prefix=fname_prefix,
            )
            
            trial_prefix = f"{fname_prefix}_trial-{trl_id}"
            # Fit model use AFNI's 3dREMLfit
            bucket_file = fit_3dREMLfit_cifti_separate(
                xmat_file,
                bucket_dir,
                volume_file=scaled_file,
                volume_mask_file=mask_intersect,
                prefix=trial_prefix,
            )
            
            output_files = extract_sub_bucket_3dTcat(
                input_bucket_file = bucket_file,
                out_dir = Path(glm_dir),
                prefix = trial_prefix,
            )
            
            print(output_files)