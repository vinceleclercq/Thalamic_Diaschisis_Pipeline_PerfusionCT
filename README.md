# Thalamic diaschisis CT perfusion pipeline

## Overview

This repository contains the MATLAB code used to:

1. convert anatomical CT and CT-perfusion Tmax DICOM series to NIfTI;
2. coregister Tmax maps to anatomical CT;
3. normalize the maps to CTseg template space;
4. create subject-level validity masks and group maps;
5. resample the AAL atlas onto the normalized Tmax grid;
6. extract bilateral thalamic Tmax measurements.

The pipeline was developed for the study:

> **Thalamic Tmax Asymmetry, Acute EEG Abnormalities, and Seizure Outcomes After Middle Cerebral Artery Ischemic Stroke**

No patient-level clinical or imaging data are included.

## Repository structure

```text
thalamic-diaschisis-pipeline/
├── run_all.m
├── PerfCT_Reg_MNI_001.m
├── Extract_ROI_Thalamus_001.m
├── config/
│   └── example_config.m
├── src/
│   ├── load_pipeline_config.m
│   ├── setup_environment.m
│   ├── convert_dicom_to_nifti.m
│   ├── register_tmax_to_template.m
│   ├── create_group_maps.m
│   ├── extract_thalamic_tmax.m
│   └── validate_pipeline_outputs.m
├── docs/
│   ├── METHODS_TEXT.md
│   ├── OUTPUTS.md
│   └── ZENODO_RELEASE_CHECKLIST.md
├── CITATION.cff
├── .zenodo.json
├── .gitignore
├── CHANGELOG.md
└── LICENSE
```

## Requirements

- MATLAB R2021a or later recommended
- SPM12
- CTseg, including:
  - `spm_CTseg.m`
  - `mu_CTseg.nii`
- FieldTrip
- AAL atlas distributed with FieldTrip:
  - `template/atlas/aal/ROI_MNI_V4.nii`

## Expected input structure

```text
CT_perf/
├── Sub001/
│   ├── ANAT/
│   │   └── [anatomical CT DICOM files]
│   └── Perf_T/
│       └── [Tmax DICOM files]
├── Sub002/
│   ├── ANAT/
│   └── Perf_T/
└── ...
```

Subject directories must match the pattern configured in `cfg.subject_pattern`, which defaults to `Sub*`.

## Installation

1. Download or clone this repository.
2. Copy:

```text
config/example_config.m
```

to:

```text
config/local_config.m
```

3. Edit `config/local_config.m` and set the local paths to:
   - the patient-data root directory;
   - SPM12;
   - CTseg;
   - FieldTrip.

`local_config.m` is ignored by Git and must not be uploaded.

## Running the complete pipeline

In MATLAB:

```matlab
cd('PATH_TO_REPOSITORY')
run_all
```

The three processing stages are executed sequentially:

```matlab
convert_dicom_to_nifti(cfg)
register_tmax_to_template(cfg)
extract_thalamic_tmax(cfg)
```

A final validation report is then generated.

## Running the historical entry points

The two original script names are retained as wrappers:

```matlab
PerfCT_Reg_MNI_001
Extract_ROI_Thalamus_001
```

`PerfCT_Reg_MNI_001` performs DICOM conversion, registration, normalization, and group-map creation.

`Extract_ROI_Thalamus_001` performs bilateral thalamic ROI extraction.

## Main outputs

By default, outputs are written to:

```text
CT_perf/Group_results/
```

The principal outputs are:

```text
processing_summary.csv
dicom_conversion_summary.csv
Sub001_wfTmax.nii
Sub001_mask.nii
group_sum_wfTmax.nii
group_count_valid.nii
group_mean_wfTmax.nii
group_coverage_fraction.nii
group_subjects.txt
thalamus_CT_values.csv
thalamus_CT_results.mat
pipeline_validation_report.txt
```

See `docs/OUTPUTS.md` for details.

## Quality control

Automated checks do not replace visual review. Before analysis, inspect:

- anatomical CT–Tmax coregistration;
- left–right orientation;
- normalization to CTseg template space;
- coverage of both thalami;
- subject masks;
- the processing summary for errors;
- agreement of the final thalamic values with the original analysis dataset.

The pipeline should not be released as version 1.0.0 until it reproduces the study values for all included subjects.

## Important reproducibility choices

- Tmax is treated as a contrast-arrival-delay measure, not as absolute cerebral blood flow.
- The native validity mask is defined using a configurable lower threshold, set by default to `Tmax > -100`, matching the original processing script.
- ROI extraction uses both the propagated validity mask and finite values.
- To reproduce the original study, non-positive Tmax values are excluded by default during ROI extraction.
- AAL labels are resampled using nearest-neighbour interpolation.
- Continuous Tmax data are resampled using trilinear interpolation.

## Intended use

Research and reproducibility only. The code has not been validated for clinical decision-making.
