import pydicom
import sys
from os.path import join as opj
from os import rename
from os import path
from pathlib import Path
from re import sub
from glob import glob
import simplejson as json

import mrpyconvert

import pathlib

repo = pathlib.Path('/projects/lcni/dcm/kuhl_lab/Kuhl/')
glacier_path = repo / 'glacier'
all_bids_data = glob(opj(glacier_path, 'GLACIER*'))

output_path = '/projects/kuhl_lab/wanjiag/GLACIER'
processed_files = glob(opj(output_path, 'sub-GLACIER*'))
processed_ids = [x.split('/')[5].split('-')[1] for x in processed_files]

to_process_sub_path = []
to_process_sub_id = []

for d in all_bids_data:
    subid = d.split('/')[7].split('_')[0]
    # Special case:
    if subid not in processed_ids:
        if subid == 'GLACIER00':
            continue
        
        to_process_sub_id.append(subid)

        if subid == 'GLACIER01':
            print('GLACIER01: /home/wanjiag/projects/GLACIER/derivatives/special_dcm/GLACIER01_20240130_094638')
            to_process_sub_path.append('/home/wanjiag/projects/GLACIER/derivatives/special_dcm/GLACIER01_20240130_094638')
        elif subid == 'GLACIER03':
            print('GLACIER03: /home/wanjiag/projects/GLACIER/derivatives/special_dcm/GLACIER03_20240131_135222')
            to_process_sub_path.append('/home/wanjiag/projects/GLACIER/derivatives/special_dcm/GLACIER03_20240131_135222')
        elif subid == 'GLACIER08':
            print('GLACIER08: /home/wanjiag/projects/GLACIER/derivatives/special_dcm/GLACIER08_20240223_154711')
            to_process_sub_path.append('/home/wanjiag/projects/GLACIER/derivatives/special_dcm/GLACIER08_20240223_154711')
        elif subid == 'GLACIER10':
            print('GLACIER10: /home/wanjiag/projects/GLACIER/derivatives/special_dcm/GLACIER10_20240215_145541')
            to_process_sub_path.append('/home/wanjiag/projects/GLACIER/derivatives/special_dcm/GLACIER10_20240215_145541')
        else:
            to_process_sub_path.append(d)

# Edit fmap metadata
def write_metadata(json_file, intended_list):
    """Write IntendedFor field to json metadata.

    Parameters
    ----------
    json_file : os.PathLike
        Metadata json file.
    intended_list : list[str]
        Intended file list.

    """
    # Add field
    json_file.chmod(0o644)
    with json_file.open("r") as f:
        data = json.load(f)
    with json_file.open("w") as f:
        data["IntendedFor"] = intended_list
        json.dump(data, f, indent=2)
    json_file.chmod(0o444)
    
for sub_p in to_process_sub_path:
    print(f"Processing sub-{sub_p}...")
    
    converter = mrpyconvert.Converter()
    converter.set_bids_path(output_path)
    
    converter.inspect(sub_p)
    converter.add_dicoms(sub_p)
    
    converter.add_entry('mprage', search='mprage_p2', datatype='anat', suffix='T1w')
    converter.add_entry('t2', search='t2_tse', datatype='anat', suffix='T2w')
    converter.add_entry('se_epi_ap', search='se_epi_ap', datatype='fmap', suffix='epi', chain = {'dir':'AP'},json_fields={'B0FieldIdentifier': 'fieldmap'})
    converter.add_entry('se_epi_pa', search='se_epi_pa', datatype='fmap', suffix='epi', chain = {'dir':'PA'},json_fields={'B0FieldIdentifier': 'fieldmap'})
    
    
    converter.add_entry('task1', search='EPI_1', datatype='func', suffix='bold', chain={'task':'glacier', 'run':1},json_fields={'TaskName': 'glacier', 'B0FieldSource': 'fieldmap'})
    converter.add_entry('task2', search='EPI_2', datatype='func', suffix='bold', chain={'task':'glacier', 'run':2},json_fields={'TaskName': 'glacier', 'B0FieldSource': 'fieldmap'})
    converter.add_entry('task3', search='EPI_3', datatype='func', suffix='bold', chain={'task':'glacier', 'run':3},json_fields={'TaskName': 'glacier', 'B0FieldSource': 'fieldmap'})
    converter.add_entry('task4', search='EPI_4', datatype='func', suffix='bold', chain={'task':'glacier', 'run':4},json_fields={'TaskName': 'glacier', 'B0FieldSource': 'fieldmap'})
    converter.add_entry('task5', search='EPI_5', datatype='func', suffix='bold', chain={'task':'glacier', 'run':5},json_fields={'TaskName': 'glacier', 'B0FieldSource': 'fieldmap'})
    converter.add_entry('task6', search='EPI_6', datatype='func', suffix='bold', chain={'task':'glacier', 'run':6},json_fields={'TaskName': 'glacier', 'B0FieldSource': 'fieldmap'})
    converter.add_entry('task7', search='EPI_7', datatype='func', suffix='bold', chain={'task':'glacier', 'run':7},json_fields={'TaskName': 'glacier', 'B0FieldSource': 'fieldmap'})
    converter.add_entry('task8', search='EPI_8', datatype='func', suffix='bold', chain={'task':'glacier', 'run':8},json_fields={'TaskName': 'glacier', 'B0FieldSource': 'fieldmap'})
    
    converter.add_entry('sbref', search='EPI_.*_SBRef', datatype='func', suffix='sbref', autorun=True, chain={'task':'glacier'},json_fields={'TaskName': 'glacier', 'B0FieldSource': 'fieldmap'})

    converter.convert()


### adding intendedfor parameter for field map file
for sub_id in to_process_sub_id:
    print('Adding fmap intendedto parameter')
    print(sub_id)

    f_list = list(
        sorted(glob(opj(output_path,f"sub-{sub_id}",'func', '*_bold.nii.gz')))
    )
    file_names = [f.split('/')[-1] for f in f_list]
    intended_list = [f"func/{f}" for f in file_names]
    
    if path.exists(opj(output_path, f"sub-{sub_id}", "fmap")):
        json_file_list = list(
            sorted(glob(opj(output_path,f"sub-{sub_id}",'fmap', '*.json')))
        )
        print(json_file_list)
        for json_file in json_file_list:
            json_file = Path(json_file)
            write_metadata(json_file, intended_list)
        print(f"Completed sub-{sub_id} ...")
