#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 28 18:10:01 2021

@author: alejandrobertolet
"""

from MIRDCalculator import *

ctpath = '2020-12__Studies/DOE^JANE_ANON60446_CT_2020-12-22_131402_Liver.Scan_AC..ABDOMEN..5.0..B08s_n78__00000/'
nmpath = '2020-12__Studies/DOE^JANE_ANON60446_NM_2020-12-22_131402_Liver.Scan_MAA.BRONCHIAL.EMBOLIZATION.19-218.(Recon.-.AC.)_n80__00000/'                                
countThreshold = 200
calc = MIRDCalculator(ctpath,nmpath,'Y90')
calc.CalculateOnActivityMapGrid(countThreshold)
calc.DoseInterpolationToCTGrid(1)