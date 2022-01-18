# MIRDCalculator v1.2.0
Python application to calculate internal dosimetry following MIRD formulation

This python-based application calculates doses in Nuclear Medicine applications.
MIRDCalculator takes a SPECT or PET image and computes the dose values based on the voxel-based S-values published by [Lanconelli et al.] [1].
Supported radionuclides in this database are 89Sr, 90Y, 131I, 153Sm, 177Lu, 186Re and 188Re.
Dose is calculated using the voxel sizes in the activity map and then it is interpolated to the CT grid.
The interpolation method uses the inverse of the distance between each node in the CT grid and the 8 closest nodes in the SPECT/PET grid to weight each contribution.

## Use (distributed version):
Download the MIRDCalculator-1.2.0.tar.gz 
Install using pip via:
  
  `pip install MIRDCalculator-1.2.0.tar.gz`

You can import two modules: MIRD and DICOM_RT (for dealing with DICOM features)
To calculate MIRD dose and store it in DICOM RTDOSE format you need:
  
  `from MIRD import MIRDCalculator`
  
  `MIRDCalculator.GetMIRDDoseInDICOM(basepath, nameDicom, radionuclide, tissue, norm, unit, accum, countThreshold, ct_path, nm_path)`

The parameters of this method are:
* (string) basepath: path to the directory in which CT and NM studies are located. It is assumed the structure for both is basepath/CT/ and basepath/NM/, respectively.
* (string) nameDicom: name of the RTDOSE Dicom file with the calculated dose. It is stored in the folder basepath
* (string) radionuclide: name of the radionuclide considered ('89Sr', '90Y', '131I', '153Sm', '177Lu', '186Re' ot '188Re')
* [OPTIONAL - default is 'Soft] (string) tissue: specifies what tissue is used for calculation (either 'Soft' or 'Bone' in [1])
* [OPTIONAL - default is True] (boolean) norm: if selected, the activity map counts are normalized to 1, so that 1 MBq is assumed as the overall uptake
* [OPTIONAL - default is 'Gy/mCi'] (string) unit: can be 'mGy/MBq', 'Gy/MBq', 'mGy/MBq' or 'Gy/mCi'. If accumulated activity is not selected, then the units are mGy/(MBq s) and so on.
* [OPTIONAL - default is True] (boolean) accum: if selected, accumulated activity is calculated assuming complete decay at the image positions
* [OPTIONAL - default is 0] (int) countThreshold: threshold of counts in the activity map to be considered in the calculation (for speed reasons)
* [OPTIONAL - in case you don't use basepath] (string) ctpath: path to the directory in which CT study is located
* [OPTIONAL - in case you don't use basepath] (string) nmpath: path to the directory in which the activity map image (in DICOM format) is located

### References:
[1]: Lanconelli N, Pacilio M, Meo S Lo, Botta F, Dia A Di, Aroche L A T, Pérez M A C and Cremonesi M 2012 A free database of radionuclide voxel S values for the dosimetry of nonuniform activity distributions Phys. Med. Biol. 57 517–33

## For developers:
Compile and build changes using the .sh script. Make changes according to your directories
