from hippocampus_body import *
import argparse

parser = argparse.ArgumentParser(description="Parameters.")

parser.add_argument(
    "-s",
    "--subs",
    action="store",
    nargs="*",
    default=None,
    help=(
        "One or more subject identifiers (e.g., sub-GLACIER01)."
        "If subject number is entered, will overwrite existing files."
        "If this is omitted, using a pre-defined automatic detection based on folder names."
    ),
)

parser.add_argument('--visual', '-v', action='store_true',
                   help = ("extracting visual attention region rois including EVC and PPA (converting from mni space)."))

parser.add_argument('--ashs', '-a', action='store_true',
                   help = ("extracting hippocampus subfields from ASHS"))

args = parser.parse_args()

print(args, flush = True)

if not args.visual and not args.ashs:
    print('No ROI was generated. Please add -v or -a argument to specify which ROIs should be generated')
    sys.exit()

epi_mask_threshold = 0.5
singularity_prefix = 'singularity exec --bind /projects/kuhl_lab/wanjiag/GLACIER/derivatives:/projects/kuhl_lab/wanjiag/GLACIER/derivatives /gpfs/projects/kuhl_lab/shared/fmriprep-v23.2.0.simg'

subfields = {'ca1':1,
             'ca23dg':[2,4] #choosing 2 to 4, including ca2, dg, and ca3
            }

######## Running ########

derivative_dir = '/projects/kuhl_lab/wanjiag/GLACIER/derivatives/'
glm_base_dir = opj(derivative_dir, 'glm/')
fmriprep_base_dir = opj(derivative_dir, 'fmriprep/')

roi_base_dir = opj(derivative_dir, 'roi/')
ashs_base_dir = opj(roi_base_dir, 'ASHS/')

if args.subs:
    todo_subs = args.subs
else:
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

print(todo_subs, flush = True)


for sub in todo_subs:
    output_dir = opj(roi_base_dir, sub)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
        
    fmriprep_dir = opj(fmriprep_base_dir, sub)
    glm_dir = opj(glm_base_dir, sub)
    
    print(f'--------------------{sub}-------------------', flush=True)
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
            
    if args.visual: 
        # PPA from MNI space
        print('--------------------PPA-------------------', flush=True)
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
        print('--------------------EVC-------------------', flush=True)
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

        # adding all evc files together
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

    if args.ashs: 
        # ASHS subfields
        ashs_output_dir = f'/home/wanjiag/projects/GLACIER/derivatives/roi/ASHS/{sub}'

        if os.path.exists(ashs_output_dir):
            print('--------------------ASHS-------------------', flush=True)
            left = opj(ashs_output_dir, f'final/{sub}_left_lfseg_corr_nogray.nii.gz')
            right = opj(ashs_output_dir, f'final/{sub}_right_lfseg_corr_nogray.nii.gz')

            temp_dir = opj(output_dir, 'temp')
            if not os.path.isdir(temp_dir):
                os.makedirs(temp_dir)

            t1_file = opj(ashs_output_dir, 'mprage.nii.gz')
            t2_file = opj(ashs_output_dir, 'tse.nii.gz')

            output_mat = opj(temp_dir, 't2-to-t1.mat')
            if not os.path.exists(output_mat): 
                cmd = f'flirt -in {t2_file} -ref {t1_file} -dof 6 -cost mutualinfo -omat {output_mat}'
                run_cmd(cmd)

            #ref_mat = opj(ashs_output_dir, 'flirt_t2_to_t1', 'flirt_t2_to_t1.mat') <-- didn't work with flirt, had to re-register with flirt first

            for i in subfields:
                if i == 'ca1':
                    op_string = f'-thr {subfields[i]} -uthr {subfields[i]} -bin'
                if i == 'ca23dg':
                    op_string = f'-thr {subfields[i][0]} -uthr {subfields[i][1]} -bin'

                # left roi
                left_out_file = opj(temp_dir, f'{i}-left_T2w.nii.gz')
                cmd = f'fslmaths {left} {op_string} {left_out_file}'
                run_cmd(cmd)

                # right roi
                right_out_file = opj(temp_dir, f'{i}-right_T2w.nii.gz')
                cmd = f'fslmaths {right} {op_string} {right_out_file}'
                run_cmd(cmd)

                # combined roi
                combined_file = opj(temp_dir, f'{i}_T2w.nii.gz')
                cmd = f'fslmaths {left_out_file} -add {right_out_file} -bin {combined_file}'
                run_cmd(cmd)

                # t2 to t1 space
                combined_file_t1 = opj(temp_dir, f'{i}_T2w-2-t1.nii.gz')
                cmd = f'flirt -ref {t1_file} -in {combined_file} -applyxfm -init {output_mat} -out {combined_file_t1}'
                run_cmd(cmd)

                # t1 to epi space
                combined_file_epi = opj(temp_dir, f'{i}_t1-2-epi.nii.gz')
                cmd = f'flirt -in {combined_file_t1} -ref {mean_out_file} -applyxfm -usesqform -o {combined_file_epi}'
                run_cmd(cmd)

                # epi space mask + threshold
                combined_file_epi_thr = opj(output_dir, f'{i}_t1-2-epi_masked_thr-{epi_mask_threshold}.nii.gz')
                cmd = f'fslmaths {combined_file_epi} -thr {epi_mask_threshold} -mas {brain_mask} -bin {combined_file_epi_thr}'
                run_cmd(cmd)

                print('--------------------ASHS BODY-------------------', flush=True)
                # Hippocampus body
                body_file = hippo_body(left_out_file, right_out_file, t2_file, i, sub, output_dir, temp_dir)

                # t2 to t1 space
                combined_file_t1 = opj(temp_dir, f'{i}-body_T2w-2-t1.nii.gz')
                cmd = f'flirt -ref {t1_file} -in {body_file} -applyxfm -init {output_mat} -out {combined_file_t1}'
                run_cmd(cmd)

                # t1 to epi space
                combined_file_epi = opj(temp_dir, f'{i}-body_t1-2-epi.nii.gz')
                cmd = f'flirt -in {combined_file_t1} -ref {mean_out_file} -applyxfm -usesqform -o {combined_file_epi}'
                run_cmd(cmd)

                # epi space mask + threshold
                combined_file_epi_thr = opj(output_dir, f'{i}-body_t1-2-epi_masked_thr-{epi_mask_threshold}.nii.gz')
                cmd = f'fslmaths {combined_file_epi} -thr {epi_mask_threshold} -mas {brain_mask} -bin {combined_file_epi_thr}'
                run_cmd(cmd)