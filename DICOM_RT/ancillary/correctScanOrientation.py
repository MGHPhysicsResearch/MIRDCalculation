#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 4/24/23 12:14 PM

@author: alejandrobertolet
"""

import os
import numpy as np
import pydicom
from pydicom import dcmread

input_folder = '/Users/ai925/workspace/Bailey/Images/DICOM_01/0000/001'
output_folder = '/Users/ai925/workspace/Bailey/Images/DICOM_01/0000/PET_corrected/'

# Ensure the output folder exists
os.makedirs(output_folder, exist_ok=True)

# Read all DICOM files from the input folder into a list
dicom_files = [
    dcmread(os.path.join(input_folder, f))
    for f in os.listdir(input_folder)
]

# Sort the list by slice location (Z)
dicom_files.sort(key=lambda ds: ds.ImagePositionPatient[2])

# Reverse the order of Z slices
dicom_files = dicom_files[::-1]

# Process the DICOM files
for ds in dicom_files:
    # Check if the ImageOrientationPatient is [1, 0, 0, 0, -1, 0]
    if ds.ImageOrientationPatient == [1, 0, 0, 0, -1, 0]:
        # Flip the pixel data vertically (along Y-axis)
        ds.PixelData = ds.pixel_array[:, ::-1].tobytes()

        # Update the ImagePositionPatient (origin) for the Y position
        ds.ImagePositionPatient[1] = ds.ImagePositionPatient[1] - (ds.Rows - 1) * ds.PixelSpacing[1]

        # Update the ImageOrientationPatient to [1, 0, 0, 0, 1, 0]
        ds.ImageOrientationPatient = [1, 0, 0, 0, 1, 0]

        # Patient position
        ds.PatientPosition = 'HFS'
        # Patient orientatoin
        ds.PatientOrientation = ['L', 'A']

    # Save the modified DICOM file to the output folder
    output_path = os.path.join(output_folder, ds.SOPInstanceUID + ".dcm")
    ds.save_as(output_path)

print('Processing completed!')