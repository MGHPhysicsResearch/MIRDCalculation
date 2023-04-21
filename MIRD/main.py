#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 28 18:10:01 2021

@author: alejandrobertolet
"""

from MIRDCalculator import *

#############################################
############ USER PARAMETERS ################
#############################################

# 1. Paths of directories where CT and activity maps are located, respectively
ctpath = '/Users/ai925/Dropbox (Partners HealthCare)/RPT Project/workspace/LiverPatients/Patient1/PRETX/CT'
nmpath = '/Users/ai925/Dropbox (Partners HealthCare)/RPT Project/workspace/LiverPatients/Patient1/PRETX/NM'
# 2. Radionuclide (select from '89Sr', '90Y', '131I', '153Sm', '177Lu', 186Re' and '188Re')
radionuclide = 'Y90'
# 3. Type of tissue considered for the S-values (either'Soft' or 'Bone')
tissue = 'Soft'
# 3. If normalize, calcs refer to 1MBq as total absorption
norm = True
# 4. Unit = 1 refers to mGy/MBq. Select the desired unit (using Gy and/or mCi instead)
unit = 'Gy/mCi'
# 5. If accumulated activity, calcs assume stationary activity until complete decay
accum = True
# 6. Threshold of counts to be considered in activity map (to speed up calculation)
countThreshold = 5000
# 7. Name of the RTDOSE output file
nameDcm = 'MIRDDose.dcm'
#############################################3

# Main script
calc = MIRDCalculator(ctpath, nmpath, radionuclide)
calc.CalculateOnActivityMapGrid(countThreshold, tissue, norm, accum)
calc.DoseInterpolationToCTGrid()
calc.WriteRTDoseCT(nameDcm, unit)
dosepre = calc.doseCTgrid