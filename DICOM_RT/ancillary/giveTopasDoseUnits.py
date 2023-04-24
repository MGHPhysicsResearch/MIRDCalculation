from DICOM_RT import DicomPatient as dcmpat

### USER PARAMETERS ###
# 1. Path to DICOM files
basePath = '/Users/ai925/Dropbox (Partners HealthCare)/RPT Project/workspace/Maitz-Data/TOPAS files/Y90_dogs/Bailey/Images/DICOM_01/0000/'
doseFile = 'DoseOnCTGrid.dcm'

# 2. Simulation parameters
nHistories = 1e4
desiredUnit = 'Gy/GBq'

# 3. Load patient and dose
ctPatient = dcmpat.PatientCT(basePath + "/CT_corrected")
ctPatient.LoadRTDose(basePath + doseFile, 'Dose', None, desiredUnit, nHistories)

# 4. Save new dicom file
ctPatient.WriteRTDose(basePath + doseFile)