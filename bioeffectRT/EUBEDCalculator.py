#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 1 13:38:00 2022

@authors: mjlindsey, alejandrobertolet
"""

import os
import numpy as np
from DICOM_RT import DicomPatient as dcmpat
from bioData import BioeffectData
from MIRD.Svalues import Radionuclide

class EUBEDCalculator:
    def __init__(self, basepath, dosefile, radionuclide, unit="Gy/GBq", nHistories=0, site=None):
        self.bioeffectData = BioeffectData()
        self.unit = unit
        self.nHistories = nHistories
        rn = Radionuclide(radionuclide)
        self.rnHalfLife = rn.halfLife
        if site is None:
            self.site = 'generic'
        else:
            self.site = site
        ctPath = basepath + "/CT/"
        dosePath = basepath + dosefile
        doseFileFull = os.path.basename(dosefile)
        doseFileSplit = doseFileFull.split('.')[0]
        self.basePath = basepath
        self.doseFileName = doseFileSplit
        self.ctPatient = dcmpat.PatientCT(ctPath)
        self.ctPatient.LoadRTDose(dosePath)
        try:
            structFile = os.listdir(basepath + "/RTSTRUCT/")
            structPath = basepath + "/RTSTRUCT/" + structFile[0]
            self.ctPatient.LoadStructures(structPath)
        except Exception as e1:
            print("ERROR: RTSTRUCT folder was not found. Exception: ", e1)
            try:
                structFile = os.listdir(basepath + "/RTSTRUCT_LUNGSANDLIVER/")
                structPath = basepath + "/RTSTRUCT_LUNGSANDLIVER/" + structFile[0]
                self.ctPatient.LoadStructures(structPath)
                print("RTSTRUCT_LUNGSANDLIVER loaded instead.")
            except:
                pass
        self.ROIs = list(self.ctPatient.structures3D.keys())
        print("ROIs identified: ", self.ROIs)
        self.tumors = []
        for struct in self.ROIs:
            if 'Tumor' in struct:
                self.tumors.append(struct)
        print("Tumor structures identified: ", self.tumors)
        self._convertDoseUnits()
        self.BEDimg3D = np.zeros(self.ctPatient.img3D.shape)

    def CalculateBED(self, doseThreshold=0.01):
        doseArray = self.ctPatient.quantitiesOfInterest[0].array
        threshold = doseThreshold * np.max(doseArray)
        for i in range(doseArray.shape[0]):
            if (i % 20) == 0:
                prog = i/doseArray.shape[0] * 100
                print("Calculating BED... (" + str(round(prog, 1)) + "%)")
            for j in range(doseArray.shape[1]):
                for k in range(doseArray.shape[2]):
                    dose = doseArray[i, j, k]
                    if dose >= threshold:
                        trep = self.bioeffectData.getTRepValue('n', 'generic')
                        alphabeta = self.bioeffectData.getAlphaBetaValue('n', 'generic')
                        voxelBelongsToTumor = False
                        for t in self.tumors:
                            if self.ctPatient.structures3D[t][i, j, k]:
                                voxelBelongsToTumor = True
                                trep = self.bioeffectData.getTRepValue('t', self.site)
                                alphabeta = self.bioeffectData.getAlphaBetaValue('t', self.site)
                                break
                        if not voxelBelongsToTumor:
                            for s in self.ROIs:
                                if self.ctPatient.structures3D[s][i, j, k]:
                                    trep = self.bioeffectData.getTRepValue('n', s)
                                    alphabeta = self.bioeffectData.getAlphaBetaValue('n', s)
                        self.BEDimg3D[i, j, k] = self._BED(dose, trep, alphabeta, self.rnHalfLife)

    def _BED(self, d, Trep, ab, tau):
        return d * (1 + d * Trep / (ab * (Trep + tau)))

    def WriteDICOMRTBED(self):
        name = 'BED_' + self.doseFileName + '.dcm'
        self.ctPatient.WriteRTDose(self.BEDimg3D, self.basePath+name, self.unit)

    def _convertDoseUnits(self):
        cumulatedActivityPermCi = 12337446 # MBq s
        GBqInmCi = 1/0.037
        unitInRTDose = self.ctPatient.quantitiesOfInterest[0].unit
        if self.nHistories > 0 and unitInRTDose == 'arb. unit':
            simulatedActivity = self.nHistories / 1e6  # MBq
            if self.unit == 'Gy/GBq':
                self.ctPatient.quantitiesOfInterest[0].array = cumulatedActivityPermCi/simulatedActivity * GBqInmCi * self.ctPatient.quantitiesOfInterest[0].array
            if self.unit == 'Gy/mCi':
                self.ctPatient.quantitiesOfInterest[0].array = cumulatedActivityPermCi / simulatedActivity * self.ctPatient.quantitiesOfInterest[0].array
            if self.unit == 'mGy/mCi':
                self.ctPatient.quantitiesOfInterest[0].array = cumulatedActivityPermCi / simulatedActivity / 1000 * self.ctPatient.quantitiesOfInterest[0].array
        elif unitInRTDose != self.unit:
            if unitInRTDose == 'Gy/GBq' and self.unit == 'Gy/mCi':
                self.ctPatient.quantitiesOfInterest[0].array = 1/GBqInmCi * self.ctPatient.quantitiesOfInterest[0].array
            elif unitInRTDose == "Gy/GBq" and self.unit == "mGy/mCi":
                self.ctPatient.quantitiesOfInterest[0].array = 1/GBqInmCi / 1000 * self.ctPatient.quantitiesOfInterest[0].array
            elif unitInRTDose == "Gy/mCi" and self.unit == 'Gy/GBq':
                self.ctPatient.quantitiesOfInterest[0].array = GBqInmCi * self.ctPatient.quantitiesOfInterest[0].array
            elif unitInRTDose == "Gy/mCi" and self.unit == 'mGy/mCi':
                self.ctPatient.quantitiesOfInterest[0].array = 1000 * self.ctPatient.quantitiesOfInterest[0].array
            elif unitInRTDose == "mGy/mCi" and self.unit == 'Gy/GBq':
                self.ctPatient.quantitiesOfInterest[0].array = GBqInmCi * 1000 * self.ctPatient.quantitiesOfInterest[0].array
            elif unitInRTDose == "mGy/mCi" and self.unit == 'Gy/mCi':
                self.ctPatient.quantitiesOfInterest[0].array = 1000 * self.ctPatient.quantitiesOfInterest[0].array
        self.ctPatient.quantitiesOfInterest[0].unit = self.unit



