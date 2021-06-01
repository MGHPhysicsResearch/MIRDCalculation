#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 27 16:09:19 2021

@author: alejandrobertolet
"""

import pydicom
import matplotlib.pylab as plt
from datetime import datetime
from os import listdir
import numpy as np

class DicomPatient:
    def __init__(self, dicomDirectory):
        filesInDir = listdir(dicomDirectory)
        self.dcmFiles = []
        for fname in filesInDir:
            self.dcmFiles.append(pydicom.dcmread(dicomDirectory + '/' + fname))
        
    def FilterByModality(self, modality):
        modfiles = []
        for f in self.dcmFiles:
            if hasattr(f, 'Modality') and f.Modality == modality:
                modfiles.append(f)
        self.dcmFiles = modfiles
        
    def GetImageInfo(self):
        # Assumes that all slices have same characteristics
        self.pixelSpacing = self.dcmFiles[0].PixelSpacing
        self.sliceThickness = self.dcmFiles[0].SliceThickness
        self.axAspect = self.pixelSpacing[1] / self.pixelSpacing[0]
        self.sagAspect = self.pixelSpacing[1] / self.sliceThickness
        self.corAspect = self.sliceThickness / self.pixelSpacing[0]
        
    def GetVoxelDICOMPosition(self, ix, iy, iz):
        xpos = self.firstVoxelPosDICOMCoordinates[0] + ix * self.pixelSpacing[0]
        ypos = self.firstVoxelPosDICOMCoordinates[1] + iy * self.pixelSpacing[1]
        zpos = self.firstVoxelPosDICOMCoordinates[2] + iz * self.sliceThickness
        return np.array([xpos, ypos, zpos])
    
    def GetLowerIndexesForDicomPosition(self, position):
        xini = self.firstVoxelPosDICOMCoordinates[0]
        yini = self.firstVoxelPosDICOMCoordinates[1]
        zini = self.firstVoxelPosDICOMCoordinates[2]
        dx = self.pixelSpacing[0]
        dy = self.pixelSpacing[1]
        dz = self.sliceThickness
        ix = np.floor((position[0]-xini)/(dx+1e-6))
        iy = np.floor((position[1]-yini)/(dy+1e-6))
        iz = np.floor((position[2]-zini)/(dz+1e-6))
        return np.array([ix, iy, iz])
    
    def Rescale(self):
        self.intercept = self.dcmFiles[0].RescaleIntercept
        self.slope = self.dcmFiles[0].RescaleSlope
        self.img3D = self.img3D * self.slope + self.intercept

    def plotAxialSlice(self, sliceNumber, colormap='gray'):
        minx = self.firstVoxelPosDICOMCoordinates[0]
        miny = self.firstVoxelPosDICOMCoordinates[1]
        maxx = minx + (self.img3D.shape[0]-1)*self.pixelSpacing[0]
        maxy = miny + (self.img3D.shape[1]-1)*self.pixelSpacing[1]
        p = plt.subplot(1,1,1)
        p.imshow(self.img3D[:,:,sliceNumber], extent=[minx,maxx,miny,maxy], cmap=colormap)
        p.set_aspect(self.axAspect)
        
    def WriteRTDose(self, doseGrid, name):
        try:
            base = self.slices[0].copy()
        except:
            base = self.dcmFiles[0].copy()
        rtdoseSOPClassUID = '1.2.840.10008.5.1.4.1.1.481.2'
        rtdosemodality = 'RTDOSE'
        base.SOPClassUID = rtdoseSOPClassUID
        base.Modality = rtdosemodality
        specificRootUID = '1.2.826.0.1.3680043.9.5872.'
        base.SOPInstanceUID = pydicom.uid.generate_uid(specificRootUID)
        base.SeriesInstanceUID = pydicom.uid.generate_uid(specificRootUID)
        base.Manufacturer = 'MIRDCalculator'
        base.ManufacturerModelName = 'MIRDCalculator v1.0 by abertoletreina@mgh.harvard.edu'
        base.SeriesDescription = 'Dose-MIRDCalculator'
        # Date and time
        now = datetime.now()
        base.StudyDate = now.strftime("%Y%M%d")
        base.SeriesDate = now.strftime("%Y%M%d")
        base.AcquisitionDate = now.strftime("%Y%M%d")
        base.ContentDate = now.strftime("%Y%M%d")
        base.StudyTime = now.strftime("%H%M%S")
        base.SeriesTime = now.strftime("%H%M%S")
        base.AcquisitionTime = now.strftime("%H%M%S")
        base.ContentTime = now.strftime("%H%M%S")
        # Reshape dose grid
        doseGrid = self.reshapeZAxis(doseGrid)
        base.LargestImagePixelValue = int(np.max(doseGrid))
        base.SmallestImagePixelValue = int(np.min(doseGrid))
        base.BitsAllocated = 16
        base.BitsStored = 16
        base.HighBit = 15
        [newGrid, slope] = self.convertInt16(doseGrid)
        del base.RescaleSlope
        del base.RescaleIntercept
        base.DoseGridScaling = slope
        base.DoseSummationType = 'PLAN'
        base.DoseUnits = 'RELATIVE'
        base.ImagePositionPatient = self.firstVoxelPosDICOMCoordinates
        base.NumberOfFrames = newGrid.shape[0]
        base.FrameIncrementPointer = (0x3004, 0x000c)
        frame = []
        for i in range(0, newGrid.shape[0]):
            frame.append(i * self.sliceThickness)
        base.GridFrameOffsetVector = frame
        base.PixelData = newGrid.tobytes()
        base.save_as(name)
        
    def reshapeZAxis(self, grid):
        shape = [grid.shape[2]]
        shape.append(grid.shape[0])
        shape.append(grid.shape[1])
        newgrid = np.zeros(shape)
        for i in range(0, shape[0]):
            img2D = grid[:,:,i]
            newgrid[i,:,:] = img2D
        return newgrid
    
    def convertInt16(self, grid):
        # Determine scaling
        maxAbsScoredValue = np.max(grid)
        minScoredValue = np.min(grid)
        useSigned = minScoredValue < 0
        if useSigned:
            outputScaleFactor = (maxAbsScoredValue - minScoredValue) / 32767
            newGrid = np.zeros(grid.shape, dtype='int16')
        else:
            outputScaleFactor = (maxAbsScoredValue - minScoredValue) / 65535
            newGrid = np.zeros(grid.shape, dtype='uint16')
        for i in range(0, grid.shape[0]):
            for j in range(0, grid.shape[1]):
                for k in range(0, grid.shape[2]):
                    newGrid[i,j,k] = int(grid[i,j,k] / outputScaleFactor)
        return [newGrid, outputScaleFactor]
        
class PatientCT(DicomPatient):
    def __init__(self, dicomDirectory):
        DicomPatient.__init__(self, dicomDirectory)
        self.FilterByModality('CT')
        self.GetSlices()
        print('{} CT slices found'.format(len(self.slices)))
        self.GetImageInfo()
        self.ReadPixelValues()
        self.Rescale()
        self.GetFrameOfReference()
        
    def GetSlices(self):
        self.slices = []
        for f in self.dcmFiles:
            if hasattr(f, 'SliceLocation'):
                self.slices.append(f)
        self.slices = sorted(self.slices, key=lambda s: s.SliceLocation, reverse=True)
    
    def GetFrameOfReference(self):
        self.forUID = self.slices[0].FrameOfReferenceUID
        self.firstVoxelPosDICOMCoordinates = self.slices[0].ImagePositionPatient
        
    def ReadPixelValues(self):
        imgShape = list(self.slices[0].pixel_array.shape)
        imgShape.append(len(self.slices))
        self.img3D = np.zeros(imgShape)
        # Fill 3D array with the images from the files
        for i, s in enumerate(self.slices):
            img2D = s.pixel_array
            self.img3D[:, :, i] = img2D
                            

class Patient3DActivity(DicomPatient):
    def __init__(self, dicomDirectory):
        DicomPatient.__init__(self, dicomDirectory)
        self.FilterByModality('NM')
        print('{} NM slices found'.format(len(self.dcmFiles)))
        self.GetImageInfo()
        self.VoxelSize = self.sliceThickness
        self.ReadPixelValues()
        self.GetFrameOfReference()
        self.totalCounts = np.sum(self.img3D)
        
    def GetFrameOfReference(self):
        self.forUID = self.dcmFiles[0].FrameOfReferenceUID
        self.firstVoxelPosDICOMCoordinates = self.dcmFiles[0].DetectorInformationSequence[0].ImagePositionPatient
        
    def ReadPixelValues(self):
        imgShape = list(self.dcmFiles[0].pixel_array.shape[1:])
        imgShape.append(self.dcmFiles[0].pixel_array.shape[0])
        self.img3D = np.zeros(imgShape)
        for i in range(0, self.dcmFiles[0].pixel_array.shape[0]):
            img2D = self.dcmFiles[0].pixel_array[i,:,:]
            self.img3D[:,:,i] = img2D
        
        