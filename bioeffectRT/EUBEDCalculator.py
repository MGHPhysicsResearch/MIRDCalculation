#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 1 13:38:00 2022

@authors: mjlindsey, alejandrobertolet
"""

import os
import numpy as np
from DICOM_RT import DicomPatient as dcmpat
from DICOM_RT import EvaluationManager as evalman
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
            structFiles = os.listdir(basepath + "/RTSTRUCT/")
            structFile = [f for f in structFiles if 'dcm' in f]
            structPath = basepath + "/RTSTRUCT/" + structFile[0]
            self.ctPatient.LoadStructures(structPath)
        except Exception as e1:
            print("ERROR: RTSTRUCT folder was not found. Exception: ", e1)
            try:
                structFiles = os.listdir(basepath + "/RTSTRUCT_LUNGSANDLIVER/")
                structFile = [f for f in structFiles if 'dcm' in f]
                structPath = basepath + "/RTSTRUCT_LUNGSANDLIVER/" + structFile[0]
                self.ctPatient.LoadStructures(structPath)
                print("RTSTRUCT_LUNGSANDLIVER loaded instead.")
            except:
                pass
        self.ROIs = list(self.ctPatient.structures3D.keys())
        print("ROIs identified: ", self.ROIs)
        self.tumors = []
        for struct in self.ROIs:
            if 'tumor' in struct.lower():
                self.tumors.append(struct)
        print("Tumor structures identified: ", self.tumors)
        self._convertDoseUnits()
        self.EQDXs = []
        self.Xs = []

    def CalculateEQDXs(self, X=[0], doseThreshold=0.01):
        doseArray = self.ctPatient.quantitiesOfInterest[0].array
        threshold = doseThreshold * np.max(doseArray)
        self.Xs = X
        for x in X:
            self.EQDXs.append(np.zeros(self.ctPatient.quantitiesOfInterest[0].array.shape))
        for i in range(doseArray.shape[0]):
            if (i % 20) == 0:
                prog = i/doseArray.shape[0] * 100
                print("Calculating EQDXs... (" + str(round(prog, 1)) + "%)")
            for j in range(doseArray.shape[1]):
                for k in range(doseArray.shape[2]):
                    dose = doseArray[i, j, k]
                    if dose >= threshold:
                        trep = self.bioeffectData.getTRepValue('n', 'generic')
                        alphabeta = self.bioeffectData.getAlphaBetaValue('n', 'generic')
                        voxelBelongsToTumor = False
                        for it, t in enumerate(self.tumors):
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
                        for ix, x in enumerate(X):
                            self.EQDXs[ix][i, j, k] = self._EQDX(x, dose, trep, alphabeta, self.rnHalfLife)
        for i, x in enumerate(self.Xs):
            if x == 0:
                qoi = dcmpat.QoIDistribution(self.EQDXs[i], 'BED', self.unit + '(BED)')
            else:
                qoi = dcmpat.QoIDistribution(self.EQDXs[i], 'EQD_' + str(x), self.unit + '(EQD_' + str(x) + ')')
            self.ctPatient.quantitiesOfInterest.append(qoi)
        self.eval = evalman.EvaluationManager(self.ctPatient)

    def _EQDX(self, X, d, Trep, ab, tau):
        return d * (ab + Trep/(Trep+tau) * d) / (ab + X)

    def _BED(self, d, Trep, ab, tau):
        return self._EQDX(0, d, Trep, ab, tau)

    def WriteDICOMRTEQDXs(self):
        for i, x in enumerate(self.Xs):
            if x == 0:
                description = 'BED_'
            else:
                description = 'EQD_' + str(x)
            name = description + self.doseFileName + '.dcm'
            self.ctPatient.WriteRTDose(self.EQDXs[i], self.basePath+name, self.unit, description)
            print(self.basePath+name, " file saved.")

    def ShowDVHs(self, path=None):
        for x in self.Xs:
            if x == 0:
                quantity = 'BED'
            else:
                quantity = 'EQD_' + str(x)
            self.eval.PlotDVHs(quantity, self.basePath)
            self.eval.printMainResults(quantity, self.basePath)

    def EUEQDX(self, X, struct = None):
        if struct is None:
            sites = self.tumors
        else:
            sites = [struct]
        for site in sites:
            if site in self.tumors:
                alpha = self.bioeffectData.getAlphaValue('t', site)
                #alpha = -10
            else:
                alpha = self.bioeffectData.getAlphaValue('n', site)
            res = []
            for x in X:
                if x == 0:
                    quantity = 'BED'
                else:
                    quantity = 'EQD_' + str(x)
                for q in self.eval.extraQoIs:
                    if q.quantity == quantity:
                        voxels = q.array[self.eval.structureArrays[site]]
                        unit = q.unit
                        break
                voxels = voxels[voxels > 0]
                sum = np.sum(np.exp(-alpha*voxels))
                N = np.count_nonzero(voxels)
                res.append(-1/alpha*np.log(sum/N))
                print('EU' + quantity + ' for ' + site + " = " + str(res[-1]) + " " + unit)
        return res

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



