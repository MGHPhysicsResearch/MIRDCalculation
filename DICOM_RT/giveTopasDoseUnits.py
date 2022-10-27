from DICOM_RT import DicomPatient as dcmpat

### USER PARAMETERS ###
# 1. Path to DICOM files
basePath = '/Users/ai925/Dropbox (Partners HealthCare)/RPT Project/BronchialSIR/Patient7/'
doseFile = 'DoseOnCTGrid.dcm'

# 2. Simulation parameters
nHistories = 1e7
desiredUnit = 'Gy/GBq'

# 3. Load patient and dose
ctPatient = dcmpat.PatientCT(basePath + "/CT")
ctPatient.LoadRTDose(basePath + doseFile, 'Dose', None, desiredUnit, nHistories)

# 4. Save new dicom file
ctPatient.WriteRTDose(basePath + doseFile)