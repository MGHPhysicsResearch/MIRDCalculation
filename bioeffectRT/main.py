#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 1 13:49:21 2022

@authors: mjlindsey, alejandrobertolet
"""
from EUBEDCalculator import *
import numpy as np

### USER PARAMETERS ###
# 1. Path to DICOM files
basePath = '/Users/ai925/Dropbox (Partners HealthCare)/RPT Project/workspace/BronchialSIR/Patient3/'

# 2. RTDOSE filename
doseFile = 'DoseOnCTGrid.dcm'

# 3. Radionuclide used and (optional) number of histories if RTDose file was obtained from MC/TOPAS,  and desired unit for output
radionuclide = 'Y90'
nHistories = 1e7
unit = 'Gy/GBq'

# 4. Tumor location to indicate parameters to use
site = "Lung"

# 5. Select EQDXs to calculate (X=0 -> BED). Calculate also EQD_2Gy by default
X = [2]

# 6. Metrics to consider
metrics = ['EUEQDX', 'MeanEQDX', 'EUEQDX', 'MeanEQDX','EUEQDX', 'MeanEQDX']
structures = ['LEFT TUMOR', 'LEFT TUMOR', 'RIGHT LUNG', 'RIGHT LUNG', 'LEFT LUNG', 'LEFT LUNG']

# Main script
calc = EUBEDCalculator(basePath, doseFile, radionuclide, unit, nHistories, site, rtstructpath='/RTSTRUCT_Corrected/')
# Option to get new structure such as lung - tumor (requires knowledge of the structure names for each patient)
#calc.ctPatient.addNewBooleanStructure('subtraction', 'Lung', ['Left-sided tumor', 'Right-sided tumor 1', 'Right-sided tumor 2', 'Right-sided tumor 3'])
calc.GetPredictiveActivityCurves(metrics, structures, X)
calc.PlotPredictiveActivityCurves(basePath)

# Evaluate results and calculate EUBED
calc.ShowDVHs()
