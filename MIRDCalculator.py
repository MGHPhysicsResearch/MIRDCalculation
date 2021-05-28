#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 28 11:45:58 2021

@author: alejandrobertolet
"""

import DicomPatient as dcmpat
import Svalues
import numpy as np

class MIRDCalculator:
    def __init__(self, CTpath, NMpath, radionuclide):
        self.patCT = dcmpat.PatientCT(CTpath)
        self.patActMap = dcmpat.Patient3DActivity(NMpath)
        self.Svalues = Svalues.SValuesData(radionuclide)
        
    def CalculateOnActivityMapGrid(self, threshold = 0):
        shape = self.patActMap.img3D.shape
        self.doseAMGrid = np.zeros(shape)
        maxDistance = self.Svalues.maximumDistanceInVoxels
        for iax in range(0, shape[0]):
            for iay in range(0, shape[1]):
                for iaz in range(0, shape[2]):
                    act = self.patActMap.img3D[iax, iay, iaz]
                    if act > threshold:
                        porc = (iaz+iay*shape[2]+iax*shape[1]*shape[2])/(shape[0]*shape[1]*shape[2])*100
                        print("Calculating... (" + str(round(porc,1))+"%)")
                        for idx in range(0, maxDistance):
                            for idy in range(0, maxDistance):
                                for idz in range(0, maxDistance):
                                    S = self.Svalues.GetSValue(self.patActMap.VoxelSize, idx, idy, idz)
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
                                        