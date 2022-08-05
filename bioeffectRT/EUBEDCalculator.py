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
        rn = Radionuclide(radionuclide)
        self.rnHalfLife = rn.halfLife
        self.unit = unit
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
        self.ctPatient.LoadRTDose(dosePath, 'Dose', None, unit, nHistories)
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
        self.EQDXs = []
        self.Xs = []

    def CalculateEQDXs(self, X=[0], doseThreshold=0.001):
        doseArray = self.ctPatient.quantitiesOfInterest[0].array
        threshold = doseThreshold * np.max(doseArray)
        self.Xs = X
        alphabetas = np.ones(doseArray.shape) * self.bioeffectData.getAlphaBetaValue('n', 'generic')
        treps = np.ones(doseArray.shape) * self.bioeffectData.getTRepValue('n', 'generic')
        for s in self.ROIs:
            if s not in self.tumors:
                alphabetas[self.ctPatient.structures3D[s]] = self.bioeffectData.getAlphaBetaValue('n', s)
                treps[self.ctPatient.structures3D[s]] = self.bioeffectData.getTRepValue('n', s)
        for t in self.tumors:
            alphabetas[self.ctPatient.structures3D[t]] = self.bioeffectData.getAlphaBetaValue('t', t)
            treps[self.ctPatient.structures3D[t]] = self.bioeffectData.getTRepValue('t', t)
        for x in X:
            EQDX = self._EQDX(x, doseArray, treps, alphabetas, self.rnHalfLife)
            self.EQDXs.append(EQDX)
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
            self.ctPatient.WriteRTDose(self.basePath+name, self.EQDXs[i], self.unit, description)
            print(self.basePath+name, " file saved.")

    def ShowDVHs(self, path=None):
        self.eval.PlotDVHs('Dose', self.basePath)
        self.eval.printMainResults('Dose', self.basePath)
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



