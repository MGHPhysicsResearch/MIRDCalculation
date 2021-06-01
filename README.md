# MIRDCalculator v1.0
Python application to calculate internal dosimetry following MIRD formulation

This python-based application calculates doses in Nuclear Medicine applications.
MIRDCalculator takes a SPECT or PET image and computes the dose values based on the voxel-based S-values published by [Lanconelli et al.][1].
Supported radionuclides in this database are 89Sr, 90Y, 131I, 153Sm, 177Lu, 186Re and 188Re.
Dose is calculated using the voxel sizes in the activity map and then it is interpolated to the CT grid.
The interpolation method uses the inverse of the distance between each node in the CT grid and the 8 closest nodes in the SPECT/PET grid to weight each contribution.

*Requirement*: Python3 installed with pydicom (https://github.com/pydicom/pydicom)

##Use:
Download the complete package; install Python3 and pydicom. 
In the 'main.py' file the following parameters are specified:
* (string) ctpath: path to the directory in which CT study is located
* (string) nmpath: path to the directory in which the activity map image (in DICOM format) is located
* (string) tissue: choose between 'Soft' and 'Bone' to use the corresponding S-values in [1]
* (boolean) norm: if selected, the activity map counts are normalized to 1, so that 1 MBq is assumed as the overall uptake
* (float) unit: can be 1 (results in mGy/MBq), Gy (results in Gy/MBq), 1/mCi (results in (mGy/MBq) or Gy/mCi (results in Gy/mCi)
* (boolean) accum: if selected, accumulated activity is calculated assuming complete decay at the image positions
* (int) countThreshold: threshold of counts in the activity map to be considered in the calculation (for speed reasons)
* (string) nameDcm: name of the RTDOSE Dicom file with the calculated dose

###References:
[1]: Lanconelli N, Pacilio M, Meo S Lo, Botta F, Dia A Di, Aroche L A T, Pérez M A C and Cremonesi M 2012 A free database of radionuclide voxel S values for the dosimetry of nonuniform activity distributions Phys. Med. Biol. 57 517–33
