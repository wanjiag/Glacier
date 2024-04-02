import pandas as pd
import numpy as np
import os
from os import PathLike
from pathlib import Path
from typing import Optional, Union
import subprocess
import tempfile
from shlex import split
import sys

afni_prefix = '/home/wanjiag/abin/'

wb_command_prefix = '/gpfs/projects/kuhl_lab/wanjiag/packages/workbench/bin_linux64/'

def make_confounds_regressor(
    df: pd.DataFrame,
    out_dir: PathLike,
    demean: bool = True,
    split_into_runs: bool = True,
    split_into_pad_runs: bool = False,
    confounds_list: list[str] = [
    "trans_x", "trans_y", "trans_z",
    "rot_x", "rot_y", "rot_z",
    "framewise_displacement",],
) -> list[Path]:
    """Makes confounds regressors file.

    This function makes 1D txt file that could be read by AFNI's
    programs. It also has the ability to make run-specific padded files.
    This is very useful when building a design matrix contains all runs.

    Args:
        df: Confounds dataframe.
        out_dir: Directory to store output file.
        demean: If true, remove mean value from each column.
        split_into_pad_runs: If true, make run-specific confounds files
            with same length as input df. Values in rows doesn't belong
            to the current run are filled with 0.
        confounds_list: Confounds names include in the output file.
            Every specified confound should present in the df.
        prefix: Filename prefix of the output file. If it's None, the
            default filename is confounds.1D (or {run_id}_confounds.1D).

    Returns:
        A confounds regressor file which could be used in AFNI's
        3dDeconvolve program.
        If 'split_into_pad_runs' is true, returning a list of filenames
        corresponds to each run in the df.

    Raises:
        ValueError: Less than 2 runs in df if 'split_into_pad_runs' is
            true.
    """

    print(f"Confounds regressor: {', '.join(confounds_list)}.")

    # Get run list if split_into_pad_runs
    if split_into_runs or split_into_pad_runs:
        run_list = df["run_id"].unique().tolist()
        if len(run_list) < 2:
            raise ValueError("There should be at least 2 runs if 'split_into_pad_runs' is true.")
    # Mean-center confounds regressors
    if demean:
        if split_into_runs or split_into_pad_runs:
            confounds = (
                df.loc[:, ["run_id"] + confounds_list]
                .groupby(by=["run_id"], sort=False)
                .transform(lambda x: x - x.mean())
            )
            confounds = confounds.fillna(0)
            print("Mean center all regressors within each run.")
        else:
            confounds = (df.loc[:, confounds_list] - df.loc[:, confounds_list].mean()).fillna(0)
            print("Mean center all regressors.")
    # Or not
    else:
        confounds = df.loc[:, confounds_list].fillna(0)
    # Convert confounds regressors for per run regression
    if split_into_runs:
        confounds_split = dict()
        for run_id in run_list:
            confounds_split[run_id] = confounds.loc[
                df.run_id == run_id, :
            ].to_numpy()
    if split_into_pad_runs:     
        confounds_split = dict()
        for run_id in run_list:
            confounds_split[run_id] = np.zeros((df.shape[0], len(confounds_list)))
            confounds_split[run_id][df.run_id == run_id, :] = confounds.loc[
                df.run_id == run_id, :
            ].to_numpy()
    # Write confounds regressors to file
    fname = out_dir.joinpath("confounds.1D")
    confounds_file = [fname]
    np.savetxt(fname, confounds, fmt="%.6f")
    if split_into_runs or split_into_pad_runs:
        confounds_file = []
        for run_id in run_list:
            fname = out_dir.joinpath(f"{run_id}_confounds.1D")
            confounds_file.append(fname)
            np.savetxt(fname, confounds_split[run_id], fmt="%.6f")

    return confounds_file


def remove_allzero_column(confounds_file: PathLike) -> tuple[Path, list]:
    """Removes all zero column in confounds regressor file.

    This functions overwrites the input confounds regressor file. It is
    useful when the head motions are very small in some direction. In
    such case, the regressors will contain only zeros under a give float
    precision, which could be problematic for GLM programs.

    Args:
        confounds_file: Confounds regressor file.

    Returns:
        A tuple (ConfoundsFile, Index), where ConfoundsFile is the
        filename of the input confounds regressor file, and the Index is
        the index of columns only have zeros.
    """

    confounds = np.loadtxt(confounds_file)
    ncol = confounds.shape[1]
    # Check each column
    sel = [True] * ncol
    for i in range(confounds.shape[1]):
        if np.allclose(confounds[:, i], 0):
            sel[i] = False
    # Remove column in place
    if np.sum(sel) != ncol:
        confounds = confounds[:, sel]
        print(f"WARNING: Removing {ncol-np.sum(sel)} all zero column!")
        np.savetxt(confounds_file, confounds, fmt="%.6f")
        # get bad column index
        allzero_column_index = list(np.arange(ncol)[np.invert(sel)])
    else:
        allzero_column_index = []
    return confounds_file, allzero_column_index


def make_good_tr_regressor(
    df: pd.DataFrame,
    out_dir: os.PathLike,
    fd_thresh: Optional[float] = 0.5,
    std_dvars_thresh: Optional[float] = 1.5,
    split_into_runs: bool = True,
    enorm_thresh: Optional[float] = None,
    extra_censor: Optional[pd.DataFrame] = None,
    censor_prev_tr: bool = False,
    dry_run: bool = False,
) -> tuple[Path, Path]:
    """Makes good TR regressors file based on motion parameters.

    Args:
        df: Confounds dataframe.
        out_dir: Directory to store output file.
        fd_thresh: Framewise displacement threshold. TRs exceed this are
            marked as bad.
        std_dvars_thresh: Standardized DVARS threshold. TRs exceed this
            are marked as bad.
        enorm_thresh: Eucilidean norm threshold. TRs exceed this are
            marked as bad.
        extra_censor: Extra censor information. It should be a numpy
            array contrain only 1s and 0s. The 1s mark the outliers and
            the 0s mark the valid TRs. The length of the array should
            match the number of rows of df. Note, the argument
            'censor_prev_tr' does not apply to this extra censor input.
        censor_prev_tr: If true, also mark the the time point before a
            bad TR as bad.
        dry_run: If true, only print out censored TR information,
            instead of writing output files.

    Returns:
        A tuple (good_tr_file, outlier_info_file), where good_tr_file is
        the file marks good TR with 1s, and outlier_info_file is file
        contains detailed information used in outlier selection process.
    """

    # Get run length list if ignore_first_volume_per_run
    if "run_id" not in df.columns:
        raise KeyError("Column 'run_id' not found in the input dataframe.")
    run_list = df["run_id"].unique().tolist()
    run_lengths = []
    for run_id in run_list:
        run_lengths.append(df.loc[df.run_id == run_id, :].shape[0])
    # Create good tr file for timepoint censor
    good_tr = np.ones(df.shape[0], dtype=np.int16)
    outlier = df[[]].copy()
    # Censor TR based on Framewise displacement (L1 norm)
    if fd_thresh:
        print(f"Framewise Displacement threshold: {fd_thresh}")
        fd = df["framewise_displacement"].to_numpy()
        fd = np.nan_to_num(fd, nan=0)
        # Set first volume of each run to 0
        for i in range(1, len(run_lengths)):
            fd[np.sum(run_lengths[:i])] = 0
        outlier["fd"] = fd
        outlier["fd_outlier"] = np.where(outlier["fd"] > fd_thresh, 1, 0)
        good_tr = good_tr * np.where(outlier["fd_outlier"] == 1, 0, 1)
    # Censor TR based on standardized DVARS
    if std_dvars_thresh:
        print(f"Standards DVARS threshold: {std_dvars_thresh}")
        std_dvars = df["std_dvars"].to_numpy()
        std_dvars = np.nan_to_num(std_dvars, nan=0)
        # Set first volume of each run to 0
        for i in range(1, len(run_lengths)):
            std_dvars[np.sum(run_lengths[:i])] = 0
        outlier["std_dvars"] = std_dvars
        outlier["std_dvars_outlier"] = np.where(outlier["std_dvars"] > std_dvars_thresh, 1, 0)
        good_tr = good_tr * np.where(outlier["std_dvars_outlier"] == 1, 0, 1)
    '''
    # Censor TR based on Euclidean Norm (L2 norm)
    if enorm_thresh:
        print(f"Euclidean Norm threshold: {enorm_thresh}")
        enorm = calc_motion_enorm(df)
        # Set first volume of each run to 0
        for i in range(1, len(run_lengths)):
            enorm[np.sum(run_lengths[:i])] = 0
        outlier["enorm"] = enorm
        outlier["enorm_outlier"] = np.where(outlier["enorm"] > enorm_thresh, 1, 0)
        good_tr = good_tr * np.where(outlier["enorm_outlier"] == 1, 0, 1)
    # Also censor previous TR when a TR is marked as bad
    if censor_prev_tr:
        good_tr[:-1] = good_tr[:-1] * good_tr[1:]
    # Extra censor from external source
    if extra_censor is not None:
        if (
            (not isinstance(extra_censor, np.ndarray))
            or (extra_censor.ndim != 1)
            or (extra_censor.shape[0] != df.shape[0])
        ):
            raise ValueError(
                "Argument 'extra_censor' should be a 1 dimension ndarray with only 1s and 0s."
            )
        print("Calculate bad TR based on extra censor data ...")
        good_tr *= np.where(extra_censor == 0, 1, 0)
    '''
    
    # Write good tr and motion censor info to file
    good_tr_file = out_dir.joinpath("goodtr.1D")
    outlier_file = out_dir.joinpath("censor_info.csv")
    if not dry_run:
        good_tr = good_tr.T.astype(np.int16)
        np.savetxt(good_tr_file, good_tr, fmt="%i")
        outlier.to_csv(outlier_file, index=False)
    # Print out useful information
    n_censor = np.sum(good_tr == 0)
    pct_censor = np.mean(good_tr == 0) * 100
    print(f"Total censored TR number: {n_censor}({pct_censor:.2f}%)", flush=True)

    #return good_tr_file, outlier_file

    # Convert confounds regressors for per run regression
    if split_into_runs:
        good_tr_split = dict()
        good_tr_df = pd.DataFrame(good_tr)
        for run_id in run_list:
            good_tr_split[run_id] = good_tr_df.loc[
                df.run_id == run_id, :
            ].to_numpy()
    # Write confounds regressors to file
    if split_into_runs:
        goodtr_files = []
        for run_id in run_list:
            fname = out_dir.joinpath(f"{run_id}_goodtr.1D")
            goodtr_files.append(fname)
            np.savetxt(fname, good_tr_split[run_id], fmt="%.6f")
    
    return goodtr_files, outlier_file

####################### 

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

def intersect_masks(in_file: list[PathLike], out_file:PathLike):
    
    in_file = "run-*".join(in_file[0].split('run-1'))
    run_cmd(
            f"{afni_prefix}3dmask_tool -input {in_file} -inter "
            f"-prefix {out_file}"
        )
    

def scale_func_image(
    in_file: PathLike,
    out_file: PathLike,
    mask_file: Optional[PathLike] = None,
    cifti: bool = False,
) -> Path:
    """Scales a image in range of 0 to 200.

    By scaling BOLD data in range [0, 200] with mean equals to 100, the
    beta derived from GLM roughly represents the percentage of signal
    change (all regressors are scaled to unit size). See AFNI documents
    for details.

    Args:
        in_file: A functional image file.
        out_file: Output functional file.
        mask_file: A mask file applies to input image during scaling. It
            only enables when input file is a NIFTI or GIFTI file.
        cifti: If true, treat input file as a CIFTI file and scale it
            with wb_command.

    Returns:
        A functional image file.
    """

    with tempfile.TemporaryDirectory() as tmp_dir:

        if cifti:
            # Calculate temporal mean image
            mean_file = Path(tmp_dir, "mean.dscalar.nii")
            run_cmd(f"wb_command -disable-provenance -cifti-reduce {in_file} MEAN {mean_file}")
            # Scale func image in to 0-200 range
            scale_expr = f"min(200,(a/b)*100)*(a>0)*(b>0)"
            run_cmd(
                f"wb_command -disable-provenance -cifti-math '{scale_expr}' "
                f"{out_file} -var a {in_file} -var b {mean_file} -select 1 1 -repeat",
                print_output=False,
            )
        else:
            # Calculate temporal mean image
            if Path(in_file).suffix.endswith("gii"):
                mean_file = Path(tmp_dir, "mean.shape.gii")
            else:
                mean_file = Path(tmp_dir, "mean.nii.gz")
            run_cmd(f"{afni_prefix}3dTstat -mean -prefix {mean_file} {in_file}")
            # Scale func image in to 0-200 range
            if mask_file:
                scale_expr = f"c*min(200,(a/b)*100)*step(a)*step(b)"
                run_cmd(
                    f"{afni_prefix}3dcalc -a {in_file} -b {mean_file} -c {mask_file} "
                    f"-expr '{scale_expr}' -prefix {out_file}"
                )
            else:
                scale_expr = f"min(200,(a/b)*100)*step(a)*step(b)"
                run_cmd(
                    f"{afni_prefix}3dcalc -a {in_file} -b {mean_file} "
                    f"-expr '{scale_expr}' -prefix {out_file}"
                )

    return Path(out_file)


def calc_run_length(in_file: Union[PathLike, list[PathLike]], cifti: bool = False) -> list[int]:
    """
    Calculates functional image length (number of time points).

    Args:
        in_file: A single or a list of functional image files.
        cifti: If true, treat input file as a CIFTI file.

    Returns:
        A list of lengths (number of timepoints) of the input files.
    """

    if not isinstance(in_file, list):
        in_file = [in_file]
    run_length = []
    for f in in_file:
        # CIFTI or GIFTI file
        if cifti or (Path(f).suffix.endswith("gii")):
            cmd = [f"{wb_command_prefix}wb_command", "-file-information", f, "-only-number-of-maps"]
        # NIFTI file
        else:
            cmd = ["fslnvols", f]
        ntp = int(run_cmd(cmd, print_output=False).stdout)
        run_length.append(ntp)
    return run_length


# single trial helper functions

def make_singletrial_onset_time(df, trl_id, out_dir, prefix=""):
    """Makes singletrial model event onset time file."""

    prefix = prefix if prefix.endswith("_") else f"{prefix}_"
    # Basic info
    regressor_list = ["Target", "Other"]
    regressor = dict()
    for r_id in regressor_list:
        regressor[r_id] = ""
    # Regressor: Stim, target
    sel = df.loc[df.stimli == trl_id, "design_onset"]
    regressor["Target"] += " ".join([f"{i:.3f}" for i in sel.to_list()])
    # Regressor: Other
    sel = df.loc[df.stimli != trl_id, "design_onset"]
    regressor["Other"] += " ".join([f"{i:.3f}" for i in sel.to_list()])

    # Write regressors to file
    out_list = []
    for regressor_name, onset_time in regressor.items():
        fname = out_dir.joinpath(f"{prefix}trl{trl_id}_ev-{regressor_name}.1D")
        fname.write_text(onset_time)
        out_list.append(fname)

    return out_list, regressor


def make_singletrial_design_matrix(
    onset_file,
    trl_id,
    out_dir,
    n_timepoint,
    repetition_time,
    confounds_file=None,
    goodtr_file=None,
    stim_times_subtract=None,
    prefix: str = "",
):
    """Makes singletrial model design matrix."""

    prefix = prefix if prefix.endswith("_") else f"{prefix}_"
    # Xmat file
    xmat_file = out_dir.joinpath(f"{prefix}trl{trl_id}_xmat.1D")
    xmat_nocensor_file = out_dir.joinpath(f"{prefix}trl{trl_id}_nocensor_xmat.1D")
    xmat_stim_file = out_dir.joinpath(f"{prefix}trl{trl_id}_stim_xmat.1D")
    xmat_plot_file = out_dir.joinpath(f"{prefix}trl{trl_id}_xmat.png")
    # Make design matrix using AFNI's 3dDeconvolve
    cmd = (
        f"{afni_prefix}3dDeconvolve -DAFNI_USE_ERROR_FILE=NO -nodata {n_timepoint} {repetition_time} "
        f"-local_times -polort A "
        f"-num_stimts 2 "
        f"-stim_times 1 {onset_file[0]} 'SPMG1(6)' -stim_label 1 Target "
        f"-stim_times 2 {onset_file[1]} 'SPMG1(6)' -stim_label 2 Other "
        f"-x1D {xmat_file} -x1D_uncensored {xmat_nocensor_file} -x1D_stop "
    )
    if confounds_file:
        cmd += f"-ortvec {confounds_file} confounds "
    if goodtr_file:
        cmd += f"-censor {goodtr_file} "
    if stim_times_subtract:
        cmd += f"-stim_times_subtract {stim_times_subtract} "
    print(cmd)
    _ = run_cmd(cmd)
    # Extract stim only matrix
    cmd = f"{afni_prefix}1d_tool.py -infile {xmat_nocensor_file} -write_xstim {xmat_stim_file} -overwrite"
    _ = run_cmd(cmd)

    return xmat_file


def fit_3dREMLfit_cifti_separate(
    design_matrix_file: PathLike,
    out_dir: PathLike,
    left_surf_file: Optional[PathLike] = None,
    right_surf_file: Optional[PathLike] = None,
    volume_file: Optional[PathLike] = None,
    volume_mask_file: Optional[PathLike] = None,
    left_roi_file: Optional[PathLike] = None,
    right_roi_file: Optional[PathLike] = None,
    volume_label_file: Optional[PathLike] = None,
    prefix: Optional[str] = None,
    extra_3dremlfit_args: str = "-tout -rout -noFDR -nobout",
    debug: bool = False,
) -> Path:
    """Fits GLM with AFNI's 3dREMLfit on CIFTI file.

    This function fits models on each part of the CIFTI file, which are
    specified separately. It is useful when fitting multiple models on
    the same data, for example, singletrial responses estimation.

    Args:
        design_matrix_file: Design matrix file in AFNI format.
        out_dir: Directory to store output files.
        left_surf_file: Left surface GIFTI file. Optional.
        right_surf_file: Right surface GIFTI file. Optional.
        volume_file: Volume NIFTI file. Optional.
        volume_mask_file: Volume mask file of volume_file. Optional.
        left_roi_file: Left surface mask file. Optional.
        right_roi_file: Right surface mask file. Optional.
        volume_label_file: Volume structure label file. This file is
            required if the input CIFTI file has volume part.
        prefix: The output filename prefix (before .dscalar.nii).
            If None, use default names.
        extra_3dremlfit_args: Extra arguments pass to AFNI's 3dREMLFit
            program.
        debug: If true, save intermediate files to fitted_bucket folder
            inside the out_dir.

    Returns:
        A CIFTI file dscalar contains all outputs from 3dREMLFit.

    Raises:
        ValueError: None of the input file is specified.
        ValueError: Input file's format is incorrect.
        ValueError: The volume_label_file is None when volume_file is
            specified.
    """
    
    # Parse prefix
    prefix = "fitted" if prefix is None else f"{prefix}_fitted"

    # Fit GLM using AFNI's 3dREMLfit
    mask = f"-mask {volume_mask_file}" if volume_mask_file is not None else ""
    out_file = Path(out_dir, f"{prefix}_volume_bucket.nii.gz")
    cmd = (
        f"{afni_prefix}3dREMLfit -input {volume_file} {mask} -matrix {design_matrix_file} "
        f"-Rbuck {out_file} {extra_3dremlfit_args} "
    )
    print(cmd)
    run_cmd(cmd, cwd=out_dir)

    return out_file
