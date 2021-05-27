#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 27 16:09:19 2021

@author: alejandrobertolet
"""

import pydicom
import matplotlib.pylab as plt
from os import listdir
import numpy as np

class DicomCTPatient:
    def __init__(self, dicomDirectory):
        filesInDir = listdir(dicomDirectory)
        self.dcmFiles = []
        for fname in filesInDir:
            self.dcmFiles.append(pydicom.dcmread(dicomDirectory + '/' + fname))
        self.FilterByModality('CT')
        self.GetSlices()
        print('{} CT slices found'.format(len(self.slices)))
        self.GetImageInfo()
        self.ReadPixelValues()
        self.Rescale()
        self.GetFrameOfReference()
        
    def FilterByModality(self, modality):
        modfiles = []
        for f in self.dcmFiles:
            if hasattr(f, 'Modality') and f.Modality == 'CT':
                modfiles.append(f)
        self.dcmFiles = modfiles
        
    def GetSlices(self):
        self.slices = []
        for f in self.dcmFiles:
            if hasattr(f, 'SliceLocation'):
                self.slices.append(f)
        
    def GetImageInfo(self):
        # Assumes that all slices have same characteristics
        self.pixelSpacing = self.slices[0].PixelSpacing
        self.sliceThickness = self.slices[0].SliceThickness
        self.axAspect = self.pixelSpacing[1] / self.pixelSpacing[0]
        self.sagAspect = self.pixelSpacing[1] / self.sliceThickness
        self.corAspect = self.sliceThickness / self.pixelSpacing[0]
        
    def GetFrameOfReference(self):
        self.forUID = self.slices[0].FrameOfReferenceUID
        self.firstVoxelPosDICOMCoordinates = self.slices[0].ImagePositionPatient
        
    def GetVoxelDICOMPosition(self, ix, iy, iz):
        xpos = self.firstVoxelPosDICOMCoordinates[0] + ix * self.pixelSpacing[0]
        ypos = self.firstVoxelPosDICOMCoordinates[1] + iy * self.pixelSpacing[1]
        zpos = self.firstVoxelPosDICOMCoordinates[2] + iz * self.sliceThickness
        return np.array([xpos, ypos, zpos])
        
    def ReadPixelValues(self):
        imgShape = list(self.slices[0].pixel_array.shape)
        imgShape.append(len(self.slices))
        self.img3D = np.zeros(imgShape)
        
        # Fill 3D array with the images from the files
        for i, s in enumerate(self.slices):
            img2D = s.pixel_array
            self.img3D[:, :, i] = img2D
            
    def Rescale(self):
        self.intercept = self.slices[0].RescaleIntercept
        self.slope = self.slices[0].RescaleSlope
        self.img3D = self.img3D * self.slope + self.intercept
                
    def plotAxialSlice(self, sliceNumber):
        minx = self.firstVoxelPosDICOMCoordinates[0]
        miny = self.firstVoxelPosDICOMCoordinates[1]
        maxx = minx + (self.img3D.shape[0]-1)*self.pixelSpacing[0]
        maxy = miny + (self.img3D.shape[1]-1)*self.pixelSpacing[1]
        p = plt.subplot(1,1,1)
        p.imshow(self.img3D[:,:,sliceNumber], extent=[minx,maxx,miny,maxy], cmap='gray')
        p.set_aspect(self.axAspect)
    