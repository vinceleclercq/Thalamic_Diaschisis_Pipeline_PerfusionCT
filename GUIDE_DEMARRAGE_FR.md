# Suggested Methods text

## CT perfusion processing and thalamic measurements

Admission non-contrast CT images and vendor-generated CT perfusion Tmax maps were processed in MATLAB [release; MathWorks, Natick, MA, USA] using two custom processing modules (`PerfCT_Reg_MNI_001`, version 1.0.0, and `Extract_ROI_Thalamus_001`, version 1.0.0). The archived code is available at [GitHub URL] and Zenodo [DOI]. Image processing relied on SPM12 [revision], CTseg [version/commit], FieldTrip [release], and the Automated Anatomical Labeling atlas (AAL, ROI_MNI_V4).

DICOM series were converted to NIfTI format. The anatomical CT was processed with CTseg to obtain a deformation field between native CT and template space. Tmax maps were converted to floating-point format, and a native validity mask was generated to distinguish acquired perfusion voxels from background or invalid values. When required, the initial Tmax affine matrix was reinitialized using the orientation and centre of the anatomical CT. Tmax maps were then rigidly coregistered to the anatomical CT using normalized mutual information, resliced to the anatomical CT grid with trilinear interpolation, and transformed to CTseg template space. Validity masks were propagated using nearest-neighbour interpolation. Processing outputs underwent visual and numerical quality control.

The AAL atlas was resampled onto the normalized Tmax grid using nearest-neighbour interpolation to preserve discrete anatomical labels. Left and right thalamic regions of interest were identified from the corresponding AAL labels. Mean thalamic Tmax was calculated from finite positive voxels within the propagated validity mask, reproducing the prespecified study pipeline. Thalamic values were relabelled as ipsilateral or contralateral according to stroke laterality. The signed thalamic asymmetry index was calculated as:

\[
AI_{signed} =
\frac{Tmax_{ipsilateral}-Tmax_{contralateral}}
{(Tmax_{ipsilateral}+Tmax_{contralateral})/2}
\]

Positive values indicated greater ipsilateral contrast-arrival delay.

## Code availability statement

The MATLAB scripts used for DICOM conversion, CT–Tmax coregistration, spatial normalization, validity-mask propagation, and bilateral thalamic Tmax extraction are publicly available at [GitHub URL] and archived on Zenodo (version 1.0.0; DOI: [TO BE ADDED]). No patient-level imaging or clinical data are included in the repository. Access to the underlying data is restricted by ethics and privacy requirements.
