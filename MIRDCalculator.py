#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 28 11:45:58 2021

@author: alejandrobertolet
"""

import DicomPatient as dcmpat
import Svalues
import numpy as np

# Units
mCi = 37 #1mCi = 37 MBq
Gy = 1e3 #1 Gy = 1000 mGy

class MIRDCalculator:
    def __init__(self, CTpath, NMpath, radionuclide):
        self.patCT = dcmpat.PatientCT(CTpath)
        self.patActMap = dcmpat.Patient3DActivity(NMpath)
        self.Svalues = Svalues.SValuesData(radionuclide)
        self.accumulate = False
        
    def CalculateOnActivityMapGrid(self, threshold = 0, tissue = 'Soft', normalize = False, accumulate = False):
        shape = self.patActMap.img3D.shape
        self.doseAMGrid = np.zeros(shape)
        maxDistance = self.Svalues.maximumDistanceInVoxels
        for iax in range(0, shape[0]):
            porc = iax/shape[0]*100
            print("Calculating... (" + str(round(porc,1))+"%)")
            for iay in range(0, shape[1]):
                for iaz in range(0, shape[2]):
                    act = self.patActMap.img3D[iax, iay, iaz]
                    if act > threshold:
                        for idx in range(0, maxDistance):
                            for idy in range(0, maxDistance):
                                for idz in range(0, maxDistance):
                                    S = self.Svalues.GetSValue(self.patActMap.VoxelSize, idx, idy, idz, tissue)
                                    if idx == 0 and idy == 0 and idz == 0:
                                        self.doseAMGrid[iax, iay, iaz] = self.doseAMGrid[iax, iay, iaz] + S * act
                                    elif idx == 0 and idy == 0:
                                        if iaz-idz >= 0:
                                             self.doseAMGrid[iax, iay, iaz-idz] = self.doseAMGrid[iax, iay, iaz-idz] + S * act
                                        if iaz+idz < shape[2]:
                                             self.doseAMGrid[iax, iay, iaz+idz] = self.doseAMGrid[iax, iay, iaz+idz] + S * act
                                    elif idx == 0 and idz == 0:
                                        if iay-idy >= 0:
                                             self.doseAMGrid[iax, iay-idy, iaz] = self.doseAMGrid[iax, iay-idy, iaz] + S * act
                                        if iay+idy < shape[1]:
                                             self.doseAMGrid[iax, iay+idy, iaz] = self.doseAMGrid[iax, iay+idy, iaz] + S * act
                                    elif idy == 0 and idz == 0:
                                        if iax-idx >= 0:
                                             self.doseAMGrid[iax-idx, iay, iaz] = self.doseAMGrid[iax-idx, iay, iaz] + S * act
                                        if iax+idx < shape[0]:
                                             self.doseAMGrid[iax+idx, iay, iaz] = self.doseAMGrid[iax+idx, iay, iaz] + S * act
                                    elif idx == 0:
                                        if iay-idy >= 0 and iaz-idz >= 0:
                                            self.doseAMGrid[iax, iay-idy, iaz-idz] = self.doseAMGrid[iax, iay-idy, iaz-idz] + S * act
                                        if iay-idy >= 0 and iaz+idz < shape[2]:
                                            self.doseAMGrid[iax, iay-idy, iaz+idz] = self.doseAMGrid[iax, iay-idy, iaz+idz] + S * act
                                        if iay+idy < shape[1] and iaz-idz >= 0:
                                            self.doseAMGrid[iax, iay+idy, iaz-idz] = self.doseAMGrid[iax, iay+idy, iaz-idz] + S * act
                                        if iay+idy < shape[1] and iaz+idz < shape[2]:
                                            self.doseAMGrid[iax, iay+idy, iaz+idz] = self.doseAMGrid[iax, iay+idy, iaz+idz] + S * act
                                    elif idy == 0:
                                        if iax-idx >= 0 and iaz-idz >= 0:
                                            self.doseAMGrid[iax-idx, iay, iaz-idz] = self.doseAMGrid[iax-idx, iay, iaz-idz] + S * act
                                        if iax-idx >= 0 and iaz+idz < shape[2]:
                                            self.doseAMGrid[iax-idx, iay, iaz+idz] = self.doseAMGrid[iax-idx, iay, iaz+idz] + S * act
                                        if iax+idx < shape[0] and iaz-idz >= 0:
                                            self.doseAMGrid[iax+idx, iay, iaz-idz] = self.doseAMGrid[iax+idx, iay, iaz-idz] + S * act
                                        if iax+idx < shape[0] and iaz+idz < shape[2]:
                                            self.doseAMGrid[iax+idx, iay, iaz+idz] = self.doseAMGrid[iax+idx, iay, iaz+idz] + S * act
                                    elif idz == 0:
                                        if iax-idx >= 0 and iay-idy >= 0:
                                            self.doseAMGrid[iax-idx, iay-idy, iaz] = self.doseAMGrid[iax-idx, iay-idy, iaz] + S * act
                                        if iax-idx >= 0 and iay+idy < shape[1]:
                                            self.doseAMGrid[iax-idx, iay+idy, iaz] = self.doseAMGrid[iax-idx, iay+idy, iaz] + S * act
                                        if iax+idx < shape[0] and iay-idy >= 0:
                                            self.doseAMGrid[iax+idx, iay-idy, iaz] = self.doseAMGrid[iax+idx, iay-idy, iaz] + S * act
                                        if iax+idx < shape[0] and iay+idy < shape[1]:
                                            self.doseAMGrid[iax+idx, iay+idy, iaz] = self.doseAMGrid[iax+idx, iay+idy, iaz] + S * act
                                    else:
                                        if iax-idx >=0 and iay-idy >= 0 and iaz-idz >= 0:
                                            self.doseAMGrid[iax-idx, iay-idy, iaz-idz] = self.doseAMGrid[iax-idx, iay-idy, iaz-idz] + S * act
                                        if iax-idx >=0 and iay-idy >= 0 and iaz+idz < shape[2]:
                                            self.doseAMGrid[iax-idx, iay-idy, iaz+idz] = self.doseAMGrid[iax-idx, iay-idy, iaz+idz] + S * act
                                        if iax-idx >=0 and iay+idy < shape[1] and iaz-idz >= 0:
                                            self.doseAMGrid[iax-idx, iay+idy, iaz-idz] = self.doseAMGrid[iax-idx, iay+idy, iaz-idz] + S * act
                                        if iax-idx >=0 and iay+idy < shape[1] and iaz+idz < shape[2]:
                                            self.doseAMGrid[iax-idx, iay+idy, iaz+idz] = self.doseAMGrid[iax-idx, iay+idy, iaz+idz] + S * act
                                        if iax+idx < shape[0] and iay-idy >= 0 and iaz-idz >= 0:
                                            self.doseAMGrid[iax+idx, iay-idy, iaz-idz] = self.doseAMGrid[iax+idx, iay-idy, iaz-idz] + S * act
                                        if iax+idx < shape[0] and iay-idy >= 0 and iaz+idz < shape[2]:
                                            self.doseAMGrid[iax+idx, iay-idy, iaz+idz] = self.doseAMGrid[iax+idx, iay-idy, iaz+idz] + S * act
                                        if iax+idx < shape[0] and iay+idy < shape[1] and iaz-idz >= 0:
                                            self.doseAMGrid[iax+idx, iay+idy, iaz-idz] = self.doseAMGrid[iax+idx, iay+idy, iaz-idz] + S * act
                                        if iax+idx < shape[0] and iay+idy < shape[1] and iaz+idz < shape[2]:
                                            self.doseAMGrid[iax+idx, iay+idy, iaz+idz] = self.doseAMGrid[iax+idx, iay+idy, iaz+idz] + S * act
        if normalize:
            self.doseAMGrid = self.doseAMGrid / self.patActMap.totalCounts
        if accumulate:
            self.doseAMGrid = self.doseAMGrid / self.Svalues.decayConstant
            self.accumulate = True
                    
    def DoseInterpolationToCTGrid(self, threshold = 0):
        shape = self.patCT.img3D.shape
        self.doseCTgrid = np.zeros(shape)
        for icx in range(0, shape[0]):
            porc = icx/shape[0]*100
            print("Interpolating grid... (" + str(round(porc,1))+"%)")
            for icy in range(0, shape[1]):
                for icz in range(0, shape[2]):
                    position = self.patCT.GetVoxelDICOMPosition(icx, icy, icz)
                    indexes = self.patActMap.GetLowerIndexesForDicomPosition(position)
                    iax = int(indexes[0])
                    iay = int(indexes[1])
                    iaz = int(indexes[2])
                    # 8 closest vertices. Weight inversely proportional to the distance to each vertex
                    cumWeight = 0
                    if iax >= 0 and iay >= 0 and iaz >= 0 and iax < self.patActMap.img3D.shape[0] and iay < self.patActMap.img3D.shape[1] and iaz < self.patActMap.img3D.shape[2]:
                        if self.doseAMGrid[iax, iay, iaz] > threshold:
                            pos000 = self.patActMap.GetVoxelDICOMPosition(iax, iay, iaz)
                            d000 = self.__distance(position, pos000)
                            if d000 == 0:
                                self.doseCTgrid[icx, icy, icz] = self.doseAMGrid[iax, iay, iaz]
                                break
                            else:
                                self.doseCTgrid[icx,icy,icz] = self.doseCTgrid[icx,icy,icz] + self.doseAMGrid[iax, iay, iaz] / d000
                                cumWeight = cumWeight + 1/d000
                            if iaz + 1 < self.patActMap.img3D.shape[2]:
                                pos001 = pos000
                                pos001[2] = pos000[2] + self.patActMap.sliceThickness
                                d001 = self.__distance(position, pos001)
                                self.doseCTgrid[icx,icy,icz] = self.doseCTgrid[icx,icy,icz] + self.doseAMGrid[iax, iay, iaz+1] / d001
                                cumWeight = cumWeight + 1/d001
                            if iay + 1 < self.patActMap.img3D.shape[1]:
                                pos010 = pos000
                                pos010[1] = pos000[1] + self.patActMap.pixelSpacing[1]
                                d010 = self.__distance(position, pos010)
                                self.doseCTgrid[icx,icy,icz] = self.doseCTgrid[icx,icy,icz] + self.doseAMGrid[iax, iay+1, iaz] / d010
                                cumWeight = cumWeight + 1/d010
                            if iay + 1 < self.patActMap.img3D.shape[1] and iaz + 1 < self.patActMap.img3D.shape[2]:
                                pos011 = pos001
                                pos011[1] = pos000[1] + self.patActMap.pixelSpacing[1]
                                d011 = self.__distance(position, pos011)
                                self.doseCTgrid[icx,icy,icz] = self.doseCTgrid[icx,icy,icz] + self.doseAMGrid[iax, iay+1, iaz+1] / d011
                                cumWeight = cumWeight + 1/d011
                            if iax + 1 < self.patActMap.img3D.shape[0]:
                                pos100 = pos000
                                pos100[0] = pos000[0] + self.patActMap.pixelSpacing[0]
                                d100 = self.__distance(position, pos100)
                                self.doseCTgrid[icx,icy,icz] = self.doseCTgrid[icx,icy,icz] + self.doseAMGrid[iax+1, iay, iaz] / d100
                                cumWeight = cumWeight + 1/d100
                            if iax + 1 < self.patActMap.img3D.shape[0] and iaz + 1 < self.patActMap.img3D.shape[2]:
                                pos101 = pos001
                                pos101[0] = pos000[0] + self.patActMap.pixelSpacing[0]
                                d101 = self.__distance(position, pos101)
                                self.doseCTgrid[icx,icy,icz] = self.doseCTgrid[icx,icy,icz] + self.doseAMGrid[iax+1, iay, iaz+1] / d101
                                cumWeight = cumWeight + 1/d101
                            if iax + 1 < self.patActMap.img3D.shape[0] and iay + 1 < self.patActMap.img3D.shape[1]:
                                pos110 = pos010
                                pos110[0] = pos000[0] + self.patActMap.pixelSpacing[0]
                                d110 = self.__distance(position, pos110)
                                self.doseCTgrid[icx,icy,icz] = self.doseCTgrid[icx,icy,icz] + self.doseAMGrid[iax+1, iay+1, iaz] / d110
                                cumWeight = cumWeight + 1/d110
                            if iax + 1 < self.patActMap.img3D.shape[0] and iay + 1 < self.patActMap.img3D.shape[1] and iaz + 1 < self.patActMap.img3D.shape[2]:
                                pos111 = pos011
                                pos111[0] = pos000[0] + self.patActMap.pixelSpacing[0]
                                d111 = self.__distance(position, pos111)
                                self.doseCTgrid[icx,icy,icz] = self.doseCTgrid[icx,icy,icz] + self.doseAMGrid[iax+1, iay+1, iaz+1] / d111
                                cumWeight = cumWeight + 1/d111
                            self.doseCTgrid[icx,icy,icz] = self.doseCTgrid[icx,icy,icz] / cumWeight
                            
    def WriteRTDoseCT(self, name='MIRDDose.dcm', unit = 'mGy/MBq'):
        fn = 1
        if unit == 'mGy/mCi':
            fn = 1/mCi
        elif unit == 'Gy/MBq':
            fn = Gy
        elif unit == 'Gy/mCi':
            fn = Gy/mCi
        if not self.accumulate:
            unit = unit + ' s';
        self.doseCTgrid = self.doseCTgrid / fn
        self.patCT.WriteRTDose(self.doseCTgrid, name, unit)

    def WriteRTDoseAM(self, name='MIRDDose.dcm', unit = 'mGy/MBq'):
        fn = 1
        if unit == 'mGy/mCi':
            fn = 1/mCi
        elif unit == 'Gy/MBq':
            fn = Gy
        elif unit == 'Gy/mCi':
            fn = Gy/mCi
        if not self.accumulate:
            unit = unit + ' s';
        self.doseAMGrid = self.doseAMGrid / fn
        self.patCT.WriteRTDose(self.doseAMGrid, name, unit)
                    
    def __distance(self, pos1, pos2):
        pos1 = np.array(pos1)
        pos2 = np.array(pos2)
        return np.sqrt(np.sum(np.power(pos1-pos2, 2)))