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

# 3. Units for BED RTDose file
unit = "Gy/10MBq"

# 4. Tumor location to indicate parameters to use
site = "Lung"

# 5. ROIs for EUBED calculation
ROIList = ['Tumor', 'Liver', 'Lung_L', 'Lung_R']

# 6. Print EUBED and EUD into txt file
createFile = True

# 7. Scale input dose to a max value (float)
maxVoxel = 50

# Main script
calc = EUBEDCalculator(basePath, doseFile, unit, site, maxVoxel)
calc.CalculateBED()
#calc.WriteRTDosebED()
#calc.EUBED(RoiList, createFile)
#calc.EUD(RoiList, createFile)