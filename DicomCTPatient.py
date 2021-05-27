#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 27 16:09:19 2021

@author: alejandrobertolet
"""

import gdcm
import matplotlib.pylab as plt
from os import listdir

class DicomCTPatient:
    def __init__(self, dicomDirectory):
        self.filesInDir = listdir(dicomDirectory)
        
    def ReadImage(self):
        # Set up scanner
        scanner = gdcm.Scanner()
        # Add relevant tags
        modalityTag = gdcm.Tag(0x08, 0x60)
        acquisitionTag = gdcm.Tag(0x20, 0x12)
        scanner.AddTag(modalityTag)
        scanner.AddTag(acquisitionTag)
        
    