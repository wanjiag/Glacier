from nilearn.input_data import NiftiMasker
import os
from os.path import join as opj
from glob import glob
import pandas as pd
import numpy as np

derivative_dir = '/projects/kuhl_lab/wanjiag/GLACIER/derivatives/'
roi_base_dir = opj(derivative_dir, 'roi')
csv_base_dir = opj(derivative_dir, 'csv')
glm_base_dir = opj(derivative_dir, 'glm')

rois = ['evc-2-epi_thr-50_masked_bin', 'ppa_mni-2-epi_thr-0.5_masked_bin',
        'ca1_t1-2-epi_masked_thr-0.5', 'ca1-body_t1-2-epi_masked_thr-0.5',
        'ca23dg_t1-2-epi_masked_thr-0.5', 'ca23dg-body_t1-2-epi_masked_thr-0.5']
#hippo_subfields = ['ashs/body/ca1-body_thre_0.5_masked', 'ashs/body/ca23dg-body_thre_0.5_masked',
#                   'ashs/whole/ca1_thre_0.5_masked', 'ashs/whole/ca23dg_thre_0.5_masked']

automatic_detecting_subjects = True
if automatic_detecting_subjects:
    f_list = glob(os.path.join(roi_base_dir, '*sub-GLACIER*/'))
    subs = list(map(lambda f: f[len(os.path.commonpath(f_list))+1:-1], f_list))
    subs.sort()
    
    processed_list = glob(os.path.join(csv_base_dir, '*sub-GLACIER*/'))
    if len(processed_list) == 0:
        todo_subs = subs
    else:
        processed_subs = [x.split('/')[-2] for x in processed_list]
        todo_subs = [x for x in subs if x not in processed_subs]
    
    print(todo_subs)
    
for sub in todo_subs:
    
    output_dir = opj(csv_base_dir, sub)
    print(f'---------------------------{sub}---------------------------')
    
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    
        glm_folder = opj(glm_base_dir, sub, 'glm_files')
        func_files = glob(f'{glm_folder}/*tstat*')
        
    for roi in rois:
        print(f'---------------------------{roi}---------------------------')
        
        region_mask = opj(roi_base_dir, f'{sub}/{roi}.nii.gz')
        print(roi, region_mask)
        masker = NiftiMasker(region_mask)
        
        output_csv = pd.DataFrame()
        
        for file_name in func_files:
            print(f'---------------------------{file_name}---------------------------')
            
            if 'lure' in file_name:
                continue
            
            run_id = os.path.basename(file_name).split('_')[2]
            trial = os.path.basename(file_name).split('_')[4].split('-')[-1]
                
            region_data = masker.fit_transform(file_name)
            
            region_data = pd.DataFrame(region_data)
            region_data["sub"] = sub
            region_data["roi"] = roi
            region_data["run"] = run_id
            region_data["trial"] = trial
            
            output_csv = pd.concat([output_csv, region_data], ignore_index=True)
            
        output_file = opj(output_dir, f'{sub}_{roi}.csv')
        output_csv.to_csv(output_file, index=True)