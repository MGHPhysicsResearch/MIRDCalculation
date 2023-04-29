#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 27 16:09:19 2021

@author: alejandrobertolet
"""
import os
from os import listdir
from typing import List, Any

import numpy as np
import pydicom
from rt_utils import RTStructBuilder

import matplotlib.pylab as plt
from datetime import datetime

from scipy import ndimage

from DICOM_RT.StructureManager import Operations, RTStructWriter

class DicomPatient:
    def __init__(self, dicomDirectory):
        self.dicomDirectory = dicomDirectory
        filesInDir = listdir(dicomDirectory)
        self.dcmFiles = []
        self.quantitiesOfInterest = []
        # Assumes all files in dicomDirectory are DICOM files
        for fname in filesInDir:
            self.dcmFiles.append(pydicom.dcmread(dicomDirectory + '/' + fname, force=True))
        self.patientName = self.dcmFiles[0].PatientName
        self.studyDate = self.dcmFiles[0].StudyDate

    def FilterByModality(self, modality):
        modfiles = []
        for f in self.dcmFiles:
            if hasattr(f, 'Modality') and f.Modality == modality:
                modfiles.append(f)
        self.dcmFiles = modfiles
    
    def GetSlices(self):
        self.slices = []
        for f in self.dcmFiles:
            if hasattr(f, 'ImagePositionPatient'):
                self.slices.append(f)
        self.slices = sorted(self.slices, key=lambda s: s.ImagePositionPatient[2], reverse=False)
        
    def GetImageInfo(self):
        # Assumes that all slices have same characteristics
        self.pixelSpacing = self.dcmFiles[0].PixelSpacing
        self.sliceThickness = self.dcmFiles[0].SliceThickness
        self.axAspect = self.pixelSpacing[1] / self.pixelSpacing[0]
        self.sagAspect = self.pixelSpacing[1] / self.sliceThickness
        self.corAspect = self.sliceThickness / self.pixelSpacing[0]
        # Check orientation
        if hasattr(self.dcmFiles[0], 'ImageOrientationPatient'):
            self.orientation = self.dcmFiles[0].ImageOrientationPatient
            if self.orientation[0] == -1:
                self.pixelSpacing[1] = -self.pixelSpacing[1]
            if self.orientation[4] == -1:
                self.pixelSpacing[0] = -self.pixelSpacing[0]
        
    def GetVoxelDICOMPosition(self, ix, iy, iz):
        xpos = self.firstVoxelPosDICOMCoordinates[0] + ix * self.pixelSpacing[1]
        ypos = self.firstVoxelPosDICOMCoordinates[1] + iy * self.pixelSpacing[0]
        zpos = self.firstVoxelPosDICOMCoordinates[2] + iz * self.sliceThickness
        return np.array([xpos, ypos, zpos])
    
    def GetLowerIndexesForDicomPosition(self, position):
        xini = self.firstVoxelPosDICOMCoordinates[0]
        yini = self.firstVoxelPosDICOMCoordinates[1]
        zini = self.firstVoxelPosDICOMCoordinates[2]
        dx = self.pixelSpacing[1]
        dy = self.pixelSpacing[0]
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
        
    def WriteRTDose(self, name = None, doseGrid = None, unit = None, seriesDescription = None):
        if doseGrid is None:
            try:
                for q in self.quantitiesOfInterest:
                    if q.quantity == 'Dose':
                        doseGrid = q.array
                        if name is None:
                            name = 'RTDose_' + datetime.now().strftime("%m%d%y_%H%M%S") + '.dcm'
                        unit = q.unit
            except:
                print("No dose grid was found.")
                return
        if isinstance(doseGrid, str):
            try:
                for q in self.quantitiesOfInterest:
                    if q.quantity == doseGrid:
                        if name is None:
                            name = 'RTDose_' + doseGrid + '_' + datetime.now().strftime("%m%d%y_%H%M%S") + '.dcm'
                        doseGrid = q.array
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
        base.FrameOfReferenceUID = pydicom.uid.generate_uid(specificRootUID)
        base.Manufacturer = 'MIRDCalculator'
        base.ManufacturerModelName = 'DICOM_RT v1.2 by abertoletreina@mgh.harvard.edu'
        if seriesDescription is None:
            base.SeriesDescription = 'Dose-DICOM_RT'
        else:
            base.SeriesDescription = seriesDescription + 'DICOM_RT'
        # Date and time
        now = datetime.now()
        base.StudyDate = now.strftime("%Y%m%d")
        base.SeriesDate = now.strftime("%Y%m%d")
        base.AcquisitionDate = now.strftime("%Y%m%d")
        base.ContentDate = now.strftime("%Y%m%d")
        base.StudyTime = now.strftime("%H%M%S")
        base.SeriesTime = now.strftime("%H%M%S")
        base.AcquisitionTime = now.strftime("%H%M%S")
        base.ContentTime = now.strftime("%H%M%S")
        # Reshape dose grid
        doseGrid = self.reshapeZAxis(doseGrid)
        # Check orientation
        orientation = base.ImageOrientationPatient
        row_orientation, column_orientation = orientation[:3], orientation[3:]
        # Check if the CT slices are flipped along the z-axis
        if np.cross(row_orientation, column_orientation)[-1] < 0:
            doseGrid = np.flip(doseGrid, axis=0)
        base.PixelRepresentation = 0
        base.LargestImagePixelValue = int(np.ceil(np.max(doseGrid)))
        base['LargestImagePixelValue'].VR = 'US'
        base.SmallestImagePixelValue = int(np.min(doseGrid))
        base['SmallestImagePixelValue'].VR = 'US'
        base.BitsAllocated = 16
        base.BitsStored = 16
        base.HighBit = 15
        [newGrid, slope] = self.convertInt16(doseGrid)
        try:
            del base.RescaleSlope
            del base.RescaleIntercept
        except:
            pass
        base.DoseGridScaling = slope
        base.DoseSummationType = 'PLAN'
        base.DoseUnits = unit
        base.NumberOfFrames = newGrid.shape[0]
        base.FrameIncrementPointer = (0x3004, 0x000c)
        # Get the z-coordinates of the CT slices
        ct_slice_positions = [float(slice.ImagePositionPatient[2]) for slice in self.slices]
        firstslice = 0
        if (orientation[0] == 1 and orientation[-2] == -1) or (orientation[0] == -1 and orientation[-2] == 1):
            firstslice = -1
        # Set ImagePositionPatient to the position of the first CT slice
        base.ImagePositionPatient = self.slices[firstslice].ImagePositionPatient
        # Calculate the GridFrameOffsetVector based on the CT slice positions
        frame = [pos - ct_slice_positions[firstslice] for pos in ct_slice_positions]
        # If the dose grid has more slices than the CT, extend the frame with equal spacing
        if newGrid.shape[0] > len(frame):
            for i in range(len(frame), newGrid.shape[0]):
                frame.append(frame[-1] + self.sliceThickness)
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
            ROINames = rtstruct.get_roi_names()
        else:
            ROINames = ROIsList
        structures3DList = []
        self.ROINames = []
        for s in ROINames:
            try:
                structures3DList.append(rtstruct.get_roi_mask_by_name(s))
                self.ROINames.append(s)
            except:
                print("Structure " + s + " could not be read.")
        self.structures3D = dict(zip(self.ROINames, structures3DList))
        print('Structures loaded.')

    def addNewBooleanStructure(self, operation, ROI1, ROI2, newname=None):
        if type(ROI2) is list:
            struct2 = self.structures3D[ROI2[0]]
            for i in range(1, len(ROI2)):
                struct2 = Operations.Union(struct2, self.structures3D[ROI2[i]])
            name2 = 'tumors'
        else:
            struct2 = self.structures3D[ROI2]
            name2 = ROI2
        if operation.lower() == 'subtraction':
            res = Operations.Subtraction(self.structures3D[ROI1], struct2)
            name = ROI1 + '-' + name2
        elif operation.lower() == 'union' or operation.lower() == 'addition':
            res = Operations.Union(self.structures3D[ROI1], struct2)
            name = ROI1 + '+' + name2
        elif operation.lower() == 'intersection':
            res = Operations.Intersection(self.structures3D[ROI1], struct2)
            name = ROI1 + '-int-' + name2
        if newname is not None:
            name = newname
        self.structures3D[name] = res

    def LoadRTDose(self, RTDosePath, quantity = 'Dose', unit=None, desiredUnit='Gy/GBq', nHistories=0):
        '''
        Loads dose from DICOM RTDose file as 3D array.
        
        Args:
            RTDosePath --> path to input RTDose file (string)
            doseScale --> scale to apply to dose distribution (int / float)
        '''
        ds = pydicom.read_file(RTDosePath)
        dose_arr = ds.pixel_array
        dose_arr = np.swapaxes(dose_arr, 0,2)
        dose_arr = np.swapaxes(dose_arr, 0,1)
        slope = ds.DoseGridScaling
        darr = self.convertFloat64(dose_arr, slope)
        qoi = QoIDistribution()
        dx = ds.PixelSpacing[1]
        dy = ds.PixelSpacing[0]
        dz = np.abs(ds.GridFrameOffsetVector[1] - ds.GridFrameOffsetVector[0])
        initPos = np.array(ds.ImagePositionPatient).copy()
        initPos[0] = ds.ImagePositionPatient[1]
        initPos[1] = ds.ImagePositionPatient[0]
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
        print(quantity + ' array loaded with units ' + qoi.unit)
        self._convertDoseUnits(desiredUnit, nHistories)

    def _convertDoseUnits(self, unit, nHistories):
        cumulatedActivityPermCi = 12337446 # MBq s
        GBqInmCi = 1/0.037
        for i, q in enumerate(self.quantitiesOfInterest):
            unitInRTDose = q.unit
            if nHistories > 0 and unitInRTDose == 'arb. unit':
                simulatedActivity = nHistories / 1e6  # MBq
                if unit == 'Gy/GBq':
                    self.quantitiesOfInterest[i].array = cumulatedActivityPermCi/simulatedActivity * GBqInmCi * q.array
                if unit == 'Gy/mCi':
                    self.quantitiesOfInterest[i].array = cumulatedActivityPermCi / simulatedActivity * q.array
                if unit == 'mGy/mCi':
                    self.quantitiesOfInterest[i].array = cumulatedActivityPermCi / simulatedActivity / 1000 * q.array
                print(str(q.quantity) + " units converted from " + unitInRTDose + " to " + unit)
            elif unitInRTDose != unit:
                if unitInRTDose == 'Gy/GBq' and unit == 'Gy/mCi':
                    self.quantitiesOfInterest[i].array = 1/GBqInmCi * q.array
                elif unitInRTDose == "Gy/GBq" and unit == "mGy/mCi":
                    self.quantitiesOfInterest[i].array = 1/GBqInmCi / 1000 * q.array
                elif unitInRTDose == "Gy/mCi" and unit == 'Gy/GBq':
                    self.quantitiesOfInterest[i].array = GBqInmCi * q.array
                elif unitInRTDose == "Gy/mCi" and unit == 'mGy/mCi':
                    self.quantitiesOfInterest[i].array = 1000 * q.array
                elif unitInRTDose == "mGy/mCi" and unit == 'Gy/GBq':
                    self.quantitiesOfInterest[i].array = GBqInmCi * 1000 * q.array
                elif unitInRTDose == "mGy/mCi" and unit == 'Gy/mCi':
                    self.quantitiesOfInterest[i].array = 1000 * q.array
                print(str(q.quantity) + " units converted from " + unitInRTDose + " to " + unit)
            self.quantitiesOfInterest[i].unit = unit


    def DoseInterpolationToCTGrid(self, dosegrid, dx, dy, dz, iniPos, threshold = None):
        shape = self.img3D.shape
        doseCTgrid = np.zeros(shape)
        if threshold == None:
            threshold = 0.01 * np.max(dosegrid)
        minx = int((iniPos[0] - self.firstVoxelPosDICOMCoordinates[0])/(self.pixelSpacing[1]+1e-6))-1
        miny = int((iniPos[1] - self.firstVoxelPosDICOMCoordinates[1])/(self.pixelSpacing[0]+1e-6))-1
        minz = int((iniPos[2] - self.firstVoxelPosDICOMCoordinates[2])/(self.sliceThickness+1e-6))-1
        maxposxCT = self.firstVoxelPosDICOMCoordinates[0] + self.pixelSpacing[1] * shape[0]
        maxposyCT = self.firstVoxelPosDICOMCoordinates[1] + self.pixelSpacing[0] * shape[1]
        maxposzCT = self.firstVoxelPosDICOMCoordinates[2] + self.sliceThickness * shape[2]
        maxposxgrid = iniPos[0] + dx * dosegrid.shape[0]
        maxposygrid = iniPos[1] + dy * dosegrid.shape[1]
        maxposzgrid = iniPos[2] + dz * dosegrid.shape[2]
        maxx = shape[0] - int((maxposxCT - maxposxgrid)/(self.pixelSpacing[1]+1e-6))
        maxy = shape[1] - int((maxposyCT - maxposygrid)/(self.pixelSpacing[0]+1e-6))
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

    def ApplyValueOutsideMask(self, mask, value, mask_dicomPatient=None):
        img_shape = self.img3D.shape

        # Iterate through all voxels in the image
        for ix in range(img_shape[0]):
            for iy in range(img_shape[1]):
                for iz in range(img_shape[2]):
                    if mask_dicomPatient is None:
                        if not mask[ix, iy, iz]:
                            self.img3D[ix, iy, iz] = value
                    else:
                        mask_shape = mask_dicomPatient.img3D.shape
                        # Get the corresponding DICOM position in the PET image
                        dicom_position = self.GetVoxelDICOMPosition(ix, iy, iz)

                        # Get the corresponding indices in the mask
                        mask_indices = mask_dicomPatient.GetLowerIndexesForDicomPosition(dicom_position)
                        mask_indices = mask_indices.astype(int)

                        # Check if the CT indices are within the CT image bounds
                        if (0 <= mask_indices[0] < mask_shape[0] and
                                0 <= mask_indices[1] < mask_shape[1] and
                                0 <= mask_indices[2] < mask_shape[2]):

                            # If the voxel is outside the body mask, set the activity to value
                            if not mask[mask_indices[0], mask_indices[1], mask_indices[2]]:
                                self.img3D[ix, iy, iz] = value

    def RewriteRTStructAsCompatible(self, rtstruct_file):
        # Load the RTSTRUCT file
        rtstruct = pydicom.dcmread(rtstruct_file)
        ct = self.dcmFiles[0]
        # Set the SOP Class UID of the RTSTRUCT to the same as the CT file
        rtstruct.SOPClassUID = '1.2.840.10008.5.1.4.1.1.481.3'
        rtstruct.StudyInstanceUID = ct.StudyInstanceUID
        rtstruct.SeriesInstanceUID = ct.SeriesInstanceUID
        rtstruct.FrameOfReferenceUID = ct.FrameOfReferenceUID
        rtstruct.ReferencedSOPClassUID = ct.SOPClassUID
        rtstruct.ReferencedSOPInstanceUID = ct.SOPInstanceUID
        rtstruct.ReferencedFrameOfReferenceUID = ct.FrameOfReferenceUID
        # Get the corresponding SOP Instance UID for each pair contour/slice
        for i, structure in enumerate(rtstruct.ROIContourSequence):
            for j, contour in enumerate(rtstruct.ROIContourSequence[i].ContourSequence):
                contour_z = float(contour.ContourData[2])
                # Initialize variabels to store the closest slice
                closest_slice = None
                smallest_difference = float('inf')
                # Loop through all slices
                for slice in self.dcmFiles:
                    # Get the slice position
                    slice_z = float(slice.ImagePositionPatient[2])
                    # Calculate the difference between the slice position and the contour position
                    difference = abs(slice_z - contour_z)
                    # If the difference is smaller than the previous one, store the slice
                    if difference < smallest_difference:
                        closest_slice = slice
                        smallest_difference = difference
                # Set the corresponding SOP Instance UID
                rtstruct.ROIContourSequence[i].ContourSequence[j].ReferencedSOPInstanceUID = closest_slice.SOPInstanceUID
                for k, imagecontour in enumerate(rtstruct.ROIContourSequence[i].ContourSequence[j].ContourImageSequence):
                    rtstruct.ROIContourSequence[i].ContourSequence[j].ContourImageSequence[k].ReferencedSOPInstanceUID = closest_slice.SOPInstanceUID

        # Getting ordered list of SOP Instance UIDs
        sop_list = []
        for slice in self.dcmFiles:
            sop_list.append(slice.SOPInstanceUID)
        sop_list.sort()

        for ii, rfofseq in enumerate(rtstruct.ReferencedFrameOfReferenceSequence):
            rtstruct.ReferencedFrameOfReferenceSequence[ii].FrameOfReferenceUID = ct.FrameOfReferenceUID
            for jj, rssq in enumerate(rtstruct.ReferencedFrameOfReferenceSequence[ii].RTReferencedStudySequence):
                for kk, rsseq in enumerate(rtstruct.ReferencedFrameOfReferenceSequence[ii].RTReferencedStudySequence[jj].RTReferencedSeriesSequence):
                    # Get a list of the indexes of the sorted SOP Instance UIDs
                    refsop_list = []
                    for i, contour in enumerate(rtstruct.ReferencedFrameOfReferenceSequence[ii].RTReferencedStudySequence[jj].RTReferencedSeriesSequence[kk].ContourImageSequence):
                        refsop_list.append((contour.ReferencedSOPInstanceUID, i))
                    refsop_list.sort(key=lambda tup: tup[0])
                    for i, contour in enumerate(rtstruct.ReferencedFrameOfReferenceSequence[ii].RTReferencedStudySequence[jj].RTReferencedSeriesSequence[kk].ContourImageSequence):
                        index = refsop_list[i][1]
                        rtstruct.ReferencedFrameOfReferenceSequence[ii].RTReferencedStudySequence[jj].RTReferencedSeriesSequence[kk].ContourImageSequence[index].ReferencedSOPInstanceUID = sop_list[i]
        # Save the updated RTSTRUCT file
        rtstruct.save_as(rtstruct_file)

    def __distance(self, pos1, pos2):
        pos1 = np.array(pos1)
        pos2 = np.array(pos2)
        return np.sqrt(np.sum(np.power(pos1-pos2, 2)))
        
class PatientCT(DicomPatient):
    def __init__(self, dicomDirectory):
        DicomPatient.__init__(self, dicomDirectory)
        self.bodyMask = None
        self.FilterByModality('CT')
        self.GetSlices()
        print('{} CT slices found'.format(len(self.slices)))
        self.GetImageInfo()
        self.ReadPixelValues()
        self.Rescale()
        self.GetFrameOfReference()
    
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

    def CreateBodyMask(self, threshold=-300, write=False):
        # Apply threshold to create a binary mask
        binaryMask = self.img3D > threshold
        # Remove small objects
        structElement = ndimage.generate_binary_structure(3, 2)

        # Process intermediate slices
        for z in range(1, binaryMask.shape[2] - 1):
            binaryMask[:, :, z] = ndimage.binary_opening(binaryMask[:, :, z], structure=structElement[1])
            binaryMask[:, :, z] = ndimage.binary_closing(binaryMask[:, :, z], structure=structElement[1])

        # Fill holes inside the body using binary closing
        closedMask = binaryMask
        # Label connected components
        labeledMask, numComponents = ndimage.label(closedMask)
        # Get the largest connected component
        largestLabel = np.argmax(np.bincount(labeledMask.flatten())[1:]) + 1
        # Create the final body mask
        self.bodyMask = labeledMask == largestLabel

        if write:
            structs = {}
            structs['BODY'] = self.bodyMask
            writer = RTStructWriter(structs, self, 'Body_Contour')
            writer.write(self.dicomDirectory)

class Patient3DActivity(DicomPatient):
    def __init__(self, dicomDirectory):
        DicomPatient.__init__(self, dicomDirectory)
        modalities = list(set([f.Modality for f in self.dcmFiles]))
        if len(modalities) > 1:
            modName = input('Multiple modalities found: {}\nType in and enter modality name...\n'.format(modalities))
        elif len(modalities) == 1:
            modName = modalities[0]
        else:
            raise Exception('Did not find any DICOM files with a given Modality tag!')
        self.FilterByModality(modName)
        print('{} {} slices found'.format(len(self.dcmFiles), modName))
        self.GetImageInfo()
        self.VoxelSize = self.sliceThickness
        if len(self.dcmFiles) == 1:
            self.ReadPixelValues3D()
            self.GetFrameOfReference3D()
        else:
            self.GetSlices()
            # self.slices = self.slices[::-1]
            self.ReadPixelValuesSlices()
            self.Rescale()
            self.GetFrameOfReferenceSlices()
        self.totalCounts = np.sum(self.img3D)
        
    def GetFrameOfReference3D(self):
        self.forUID = self.dcmFiles[0].FrameOfReferenceUID
        self.firstVoxelPosDICOMCoordinates = self.dcmFiles[0].DetectorInformationSequence[0].ImagePositionPatient
    
    def GetFrameOfReferenceSlices(self):
        self.forUID = self.slices[0].FrameOfReferenceUID
        self.firstVoxelPosDICOMCoordinates = self.slices[0].ImagePositionPatient
    
    def ReadPixelValues3D(self):
        imgShape = list(self.dcmFiles[0].pixel_array.shape[1:])
        imgShape.append(self.dcmFiles[0].pixel_array.shape[0])
        self.img3D = np.zeros(imgShape)
        for i in range(0, self.dcmFiles[0].pixel_array.shape[0]):
            img2D = self.dcmFiles[0].pixel_array[i,:,:]
            self.img3D[:,:,i] = img2D

    def ReadPixelValuesSlices(self):
        imgShape = list(self.slices[0].pixel_array.shape)
        imgShape.append(len(self.slices))
        self.img3D = np.zeros(imgShape)
        # Fill 3D array with the images from the files
        for i, s in enumerate(self.slices):
            img2D = s.pixel_array
            self.img3D[:, :, i] = img2D

    def WriteWithNewImg3D(self, path=None):
        if path is None:
            path = self.dicomDirectory
        else:
            if not os.path.exists(path):
                # Create writeable directory
                os.makedirs(path, exist_ok=True, mode=0o777)
            else:
                # Check if directory is writeable and make it writeable if not.
                if not os.access(path, os.W_OK):
                    os.chmod(path, 0o777)
        # Normalize the pixel values in self.img3D
        img_min = self.img3D.min()
        img_max = self.img3D.max()
        normalized_img = (self.img3D - img_min) / (img_max - img_min)

        if hasattr(self, 'slices'):
            slices = self.slices
        else:
            slices = self.dcmFiles
        if len(slices) > 1:
            for i, s in enumerate(slices):
                # Strip filename and get rid of the path
                filename = str(i) + '.dcm'
                filepath = os.path.join(path, filename)
                if os.access(path, os.W_OK):
                    # Create a new DICOM slice and update its metadata
                    new_slice = pydicom.dcmread(s.filename)
                    new_slice.SOPInstanceUID = pydicom.uid.generate_uid()  # Generate a new UID
                    # Conditionally set the TransferSyntaxUID if it exists in the original file
                    if hasattr(s.file_meta, 'TransferSyntaxUID'):
                        new_slice.file_meta.TransferSyntaxUID = s.file_meta.TransferSyntaxUID
                    # Update pixel_array with the correct data type
                    if s.pixel_array.dtype == np.uint16:
                        new_pixel_data = (normalized_img[:, :, i] * np.iinfo(np.uint16).max).astype(np.uint16).tobytes()
                    elif s.pixel_array.dtype == np.uint8:
                        new_pixel_data = (normalized_img[:, :, i] * np.iinfo(np.uint8).max).astype(np.uint8).tobytes()
                    else:
                        raise Exception('Unknown data type: {}'.format(s.pixel_array.dtype))
                    new_slice.PixelData = new_pixel_data
                    new_slice.is_little_endian = s.is_little_endian
                    new_slice.is_implicit_VR = s.is_implicit_VR
                    new_slice.save_as(filepath)
                else:
                    print('Cannot write to file {}'.format(filepath))
        else:
            filename = '0.dcm'
            filepath = os.path.join(path, filename)
            # Normalize the pixel values in self.img3D
            img_min = self.img3D.min()
            img_max = self.img3D.max()
            normalized_img = (self.img3D - img_min) / (img_max - img_min)
            transposed_img = np.transpose(normalized_img, (2, 0, 1))
            if os.access(path, os.W_OK):
                # Create a new DICOM file
                new_file = pydicom.dcmread(slices[0].filename)
                # Conditionally set the TransferSyntaxUID if it exists in the original file
                if hasattr(slices[0].file_meta, 'TransferSyntaxUID'):
                    new_file.file_meta.TransferSyntaxUID = slices[0].file_meta.TransferSyntaxUID
                new_file.update(slices[0])
                new_file.SOPInstanceUID = pydicom.uid.generate_uid()  # Generate a new UID
                # Update pixel_array with the correct data type
                if slices[0].pixel_array.dtype == np.uint16:
                    new_pixel_data = (transposed_img * np.iinfo(np.uint16).max).astype(np.uint16).tobytes()
                elif slices[0].pixel_array.dtype == np.uint8:
                    new_pixel_data = (transposed_img * np.iinfo(np.uint8).max).astype(np.uint8).tobytes()
                else:
                    raise Exception('Unknown data type: {}'.format(slices[0].pixel_array.dtype))
                new_file.PixelData = new_pixel_data
                new_file.Rows, new_file.Columns = transposed_img.shape[:2]
                new_file.is_little_endian = slices[0].is_little_endian
                new_file.is_implicit_VR = slices[0].is_implicit_VR
                new_file.save_as(filepath)

class QoIDistribution:
    def __init__(self, array = None, quantity = None, unit = None):
        self.array = array
        self.quantity = quantity
        self.unit = unit
