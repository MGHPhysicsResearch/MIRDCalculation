#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 27 16:09:19 2021

@author: alejandrobertolet
"""

from os import listdir

import numpy as np
import pydicom
from rt_utils import RTStructBuilder

import matplotlib.pylab as plt
from datetime import datetime

class DicomPatient:
    def __init__(self, dicomDirectory):
        self.dicomDirectory = dicomDirectory
        filesInDir = listdir(dicomDirectory)
        self.dcmFiles = []
        self.quantitiesOfInterest = []
        for fname in filesInDir:
            if fname[-3:] == 'dcm':
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
        
    def WriteRTDose(self, doseGrid = None, name = None, unit = None):
        if doseGrid == None:
            try:
                for q in self.quantitiesOfInterest:
                    if q.quantity == 'Dose':
                        doseGrid = q.array
                        name = 'RTDose_' + datetime.now().strftime("%m%d%y_%H%M%S") + '.dcm'
                        unit = q.unit
            except:
                print("No dose grid was found.")
                return
        if isinstance(doseGrid, str):
            try:
                for q in self.quantitiesOfInterest:
                    if q.quantity == doseGrid:
                        doseGrid = q.array
                        if name == None:
                            name = 'RTDose_' + doseGrid + '_' + datetime.now().strftime("%m%d%y_%H%M%S") + '.dcm'
                        unit = q.unit
            except:
                print("No " + doseGrid + " grid was found.")
                return
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
        base.ManufacturerModelName = 'RT_DICOM v1.2 by abertoletreina@mgh.harvard.edu'
        base.SeriesDescription = 'Dose-RT_DICOM v1.2'
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
        base.DoseUnits = unit
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
    
    def convertFloat64(self, grid, slope):
        newGrid = np.zeros(grid.shape, dtype='float64')
        newGrid = grid * slope
        return newGrid
        
    def LoadStructures(self, RTStructPath, ROIsList=None):
        '''
        Loads structures from DICOM RTStruct file as 3D arrays. Arrays are stored in a dictionary.
        Function loads all structures if ROIsList is not specified.
        
        Args:
            RTStructPath --> path to input RTStruct file (string)
            ROIsList --> list containing structure names (list of strings)
        '''
        rtstruct = RTStructBuilder.create_from(self.dicomDirectory, RTStructPath)
        if ROIsList is None:
            self.ROINames = rtstruct.get_roi_names()
        else:
            self.ROINames = ROIsList
        structures3DList = []
        excludeROIS = []
        for s in self.ROINames:
            try:
                structures3DList.append(rtstruct.get_roi_mask_by_name(s))
            except:
                excludeROIS.append(s)
                print("Structure " + s + " could not be read.")
        # for i, s in enumerate(structures3DList):
        #     structures3DList[i] = np.swapaxes(s, 0, 2)
        #     structures3DList[i] = np.swapaxes(s, 0, 1)
        self.ROINames = list(set(self.ROINames) - set(excludeROIS))
        self.structures3D = dict(zip(self.ROINames, structures3DList))
        print('CTV shape', self.structures3D['CTV'].shape)
        print('Structures loaded.')
        
    def LoadRTDose(self, RTDosePath, quantity = 'Dose', unit = None, doseScale=1):
        '''
        Loads dose from DICOM RTDose file as 3D array.
        
        Args:
            RTDosePath --> path to input RTDose file (string)
            doseScale --> scale to apply to dose distribution (int / float)
        '''
        ds = pydicom.read_file(RTDosePath)
        dose_arr = ds.pixel_array*doseScale
        dose_arr = np.swapaxes(dose_arr, 0,2)
        dose_arr = np.swapaxes(dose_arr, 0,1)
        slope = ds.DoseGridScaling
        darr = self.convertFloat64(dose_arr, slope)
        qoi = QoIDistribution()
        dx = ds.PixelSpacing[0]
        dy = ds.PixelSpacing[1]
        dz = np.abs(ds.GridFrameOffsetVector[1] - ds.GridFrameOffsetVector[0])
        initPos = ds.ImagePositionPatient
        if darr.shape == self.img3D.shape:
            qoi.array = darr
        else:
            qoi.array = self.DoseInterpolationToCTGrid(darr, dx, dy, dz, initPos)
        qoi.quantity = quantity
        if unit is not None:
            qoi.unit = unit
        else:
            try:
                qoi.unit = ds.DoseUnits
            except:
                qoi.unit = 'arb. unit'
        self.quantitiesOfInterest.append(qoi)
        print(quantity + ' array loaded.')
    
    def DoseInterpolationToCTGrid(self, dosegrid, dx, dy, dz, iniPos, threshold = None):
        shape = self.img3D.shape
        doseCTgrid = np.zeros(shape)
        if threshold == None:
            threshold = 0.01 * np.max(dosegrid)
        minx = int((iniPos[0] - self.firstVoxelPosDICOMCoordinates[0])/(self.pixelSpacing[0]+1e-6))-1
        miny = int((iniPos[1] - self.firstVoxelPosDICOMCoordinates[1])/(self.pixelSpacing[1]+1e-6))-1
        minz = int((iniPos[2] - self.firstVoxelPosDICOMCoordinates[2])/(self.sliceThickness+1e-6))-1
        maxposxCT = self.firstVoxelPosDICOMCoordinates[0] + self.pixelSpacing[0] * shape[0]
        maxposyCT = self.firstVoxelPosDICOMCoordinates[1] + self.pixelSpacing[1] * shape[1]
        maxposzCT = self.firstVoxelPosDICOMCoordinates[2] + self.sliceThickness * shape[2]
        maxposxgrid = iniPos[0] + dx * dosegrid.shape[0]
        maxposygrid = iniPos[1] + dy * dosegrid.shape[1]
        maxposzgrid = iniPos[2] + dz * dosegrid.shape[2]
        maxx = shape[0] - int((maxposxCT - maxposxgrid)/(self.pixelSpacing[0]+1e-6))
        maxy = shape[1] - int((maxposyCT - maxposygrid)/(self.pixelSpacing[1]+1e-6))
        maxz = shape[2] - int((maxposzCT - maxposzgrid)/(self.sliceThickness+1e-6))
        for icx in range(minx, maxx):
            porc = (icx-minx)/(maxx-minx)*100
            print("Interpolating grid... (" + str(round(porc,1))+"%)")
            for icy in range(miny, maxy):
                for icz in range(minz, maxz):
                    position = self.GetVoxelDICOMPosition(icx, icy, icz)
                    iax = int((position[0]-iniPos[0])/(dx+1e-6))
                    iay = int((position[1]-iniPos[1])/(dy+1e-6))
                    iaz = int((position[2]-iniPos[2])/(dz+1e-6))
                    # 8 closest vertices. Weight inversely proportional to the distance to each vertex
                    cumWeight = 0
                    try:
                        if dosegrid[iax, iay, iaz] > threshold:
                            x = iniPos[0] + iax * dx
                            y = iniPos[1] + iay * dy
                            z = iniPos[2] + iaz * dz
                            pos000 = np.array([x,y,z])
                            d000 = self.__distance(position, pos000)
                            if d000 == 0:
                                doseCTgrid[icx, icy, icz] = dosegrid[iax, iay, iaz]
                                break
                            else:
                                doseCTgrid[icx,icy,icz] = doseCTgrid[icx,icy,icz] + dosegrid[iax, iay, iaz] / d000
                                cumWeight = cumWeight + 1/d000
                            if iaz + 1 < dosegrid.shape[2]:
                                pos001 = pos000
                                pos001[2] = pos000[2] + dz
                                d001 = self.__distance(position, pos001)
                                doseCTgrid[icx,icy,icz] = doseCTgrid[icx,icy,icz] + dosegrid[iax, iay, iaz+1] / d001
                                cumWeight = cumWeight + 1/d001
                            if iay + 1 < dosegrid.shape[1]:
                                pos010 = pos000
                                pos010[1] = pos000[1] + dy
                                d010 = self.__distance(position, pos010)
                                doseCTgrid[icx,icy,icz] = doseCTgrid[icx,icy,icz] + dosegrid[iax, iay+1, iaz] / d010
                                cumWeight = cumWeight + 1/d010
                            if iay + 1 < dosegrid.shape[1] and iaz + 1 < dosegrid.shape[2]:
                                pos011 = pos001
                                pos011[1] = pos000[1] + dy
                                d011 = self.__distance(position, pos011)
                                doseCTgrid[icx,icy,icz] = doseCTgrid[icx,icy,icz] + dosegrid[iax, iay+1, iaz+1] / d011
                                cumWeight = cumWeight + 1/d011
                            if iax + 1 < dosegrid.shape[0]:
                                pos100 = pos000
                                pos100[0] = pos000[0] + dx
                                d100 = self.__distance(position, pos100)
                                doseCTgrid[icx,icy,icz] = doseCTgrid[icx,icy,icz] + dosegrid[iax+1, iay, iaz] / d100
                                cumWeight = cumWeight + 1/d100
                            if iax + 1 < dosegrid.shape[0] and iaz + 1 < dosegrid.shape[2]:
                                pos101 = pos001
                                pos101[0] = pos000[0] + dx
                                d101 = self.__distance(position, pos101)
                                doseCTgrid[icx,icy,icz] = doseCTgrid[icx,icy,icz] + dosegrid[iax+1, iay, iaz+1] / d101
                                cumWeight = cumWeight + 1/d101
                            if iax + 1 < dosegrid.shape[0] and iay + 1 < dosegrid.shape[1]:
                                pos110 = pos010
                                pos110[0] = pos000[0] + dx
                                d110 = self.__distance(position, pos110)
                                doseCTgrid[icx,icy,icz] = doseCTgrid[icx,icy,icz] + dosegrid[iax+1, iay+1, iaz] / d110
                                cumWeight = cumWeight + 1/d110
                            if iax + 1 < dosegrid.shape[0] and iay + 1 < dosegrid.shape[1] and iaz + 1 < dosegrid.shape[2]:
                                pos111 = pos011
                                pos111[0] = pos000[0] + dx
                                d111 = self.__distance(position, pos111)
                                doseCTgrid[icx,icy,icz] = doseCTgrid[icx,icy,icz] + dosegrid[iax+1, iay+1, iaz+1] / d111
                                cumWeight = cumWeight + 1/d111
                            doseCTgrid[icx,icy,icz] = doseCTgrid[icx,icy,icz] / cumWeight
                    except:
                        print("Error at: ", iax, iay, iaz)
                        pass
        return doseCTgrid
                
    def __distance(self, pos1, pos2):
        pos1 = np.array(pos1)
        pos2 = np.array(pos2)
        return np.sqrt(np.sum(np.power(pos1-pos2, 2)))
        
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
        print('CT shape', self.img3D.shape)
        
    def GetSlices(self):
        self.slices = []
        for f in self.dcmFiles:
            if hasattr(f, 'SliceLocation'):
                self.slices.append(f)
        self.slices = sorted(self.slices, key=lambda s: s.ImagePositionPatient[2], reverse=False)
    
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
            
    
   
class QoIDistribution:
    def __init__(self, array = None, quantity = None, unit = None):
        self.array = array
        self.quantity = quantity
        self.unit = unit

        