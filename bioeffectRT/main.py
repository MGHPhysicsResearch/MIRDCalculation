#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 1 13:49:21 2022

@authors: mjlindsey, alejandrobertolet
"""
from EUBEDCalculator import *

### USER PARAMETERS ###
# 1. Path to DICOM files
basePath = '/Users/ai925/Dropbox (Partners HealthCare)/RPT Project/BronchialSIR/Patient4/'

# 2. RTDOSE filename
doseFile = 'DoseOnCTGrid.dcm'

# 3. Radionuclide used and (optional) number of histories if RTDose file was obtained from MC/TOPAS,  and desired unit for output
radionuclide = 'Y90'
nHistories = 1e7
unit = 'Gy/GBq'

# 4. Tumor location to indicate parameters to use
site = "Lung"

# Main script
calc = EUBEDCalculator(basePath, doseFile, radionuclide, unit, nHistories, site)
calc.CalculateBED()
calc.WriteDICOMRTBED()
#calc.EUBED(RoiList, createFile)
#calc.EUD(RoiList, createFile)