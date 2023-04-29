from DICOM_RT import DicomPatient as dcmpat

nm_folder = '/Users/ai925/workspace/LiverPatient11/NM/'
ct_folder = '/Users/ai925/workspace/LiverPatient11/CT/'

new_nm_folder = '/Users/ai925/workspace/LiverPatient11/NM_masked/'

# 1. Load patient and calculate body mask
ctPatient = dcmpat.PatientCT(ct_folder)
ctPatient.CreateBodyMask(write=True)

# 2. Load PET and apply mask
nmPatient = dcmpat.Patient3DActivity(nm_folder)
nmPatient.ApplyValueOutsideMask(ctPatient.bodyMask, 0, ctPatient)

# 3. Save new dicom files
nmPatient.WriteWithNewImg3D(new_nm_folder)