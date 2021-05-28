#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 28 11:45:58 2021

@author: alejandrobertolet
"""

import DicomPatient as dcmpat

class MIRDCalculator:
    def __init__(self, CTpath, NMpath):
        self.patCT = dcmpat.PatientCT(CTpath)
        self.patActMap = dcmpat.Patient3DActivity(NMpath)
        