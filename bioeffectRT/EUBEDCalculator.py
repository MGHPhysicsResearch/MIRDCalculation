#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 1 13:38:00 2022

@authors: mjlindsey, alejandrobertolet
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from DICOM_RT import DicomPatient as dcmpat
from DICOM_RT import EvaluationManager as evalman
from bioeffectRT.bioData import BioeffectData
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
            print("ERROR: RTSTRUCT folder was not found. Exception: ", e1)
            try:
                structFiles = os.listdir(basepath + "/RTSTRUCT_LUNGSANDLIVER/")
                structFile = [f for f in structFiles if 'dcm' in f]
                structPath = basepath + "/RTSTRUCT_LUNGSANDLIVER/" + structFile[0]
                self.ctPatient.LoadStructures(structPath)
                print("RTSTRUCT_LUNGSANDLIVER loaded instead.")
            except:
                pass
        try:
            self.ROIs = list(self.ctPatient.structures3D.keys())
            print("ROIs identified: ", self.ROIs)
        except:
            self.ROIs = []
            print("No structures identified.")
        self.tumors = []
        if len(self.ROIs) > 0:
            for struct in self.ROIs:
                if 'tumor' in struct.lower():
                    self.tumors.append(struct)
            print("Tumor structures identified: ", self.tumors)
        self.EQDXs = []
        self.Xs = []

    def CalculateEQDXs(self, X=[0], activityInjected = None, scaleDose = False):
        self.Xs = X
        self.EQDXs = []
        doseArray = self.ctPatient.quantitiesOfInterest[0].array
        if self.ctPatient.quantitiesOfInterest[0].unit == 'Gy/GBq' and activityInjected is not None:
            doseArray = doseArray * activityInjected
            if scaleDose:
                self.ctPatient.quantitiesOfInterest[0].array = doseArray
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
                qoi = dcmpat.QoIDistribution(self.EQDXs[i], 'BED', 'Gy_BED')
            else:
                qoi = dcmpat.QoIDistribution(self.EQDXs[i], 'EQD_' + str(x), 'Gy_EQD' + str(x))
            found = False

            for iq, q in enumerate(self.ctPatient.quantitiesOfInterest):
                if q.quantity == qoi.quantity:
                    self.ctPatient.quantitiesOfInterest[iq] = qoi
                    found = True
                    break
            if not found:
                self.ctPatient.quantitiesOfInterest.append(qoi)
        if len(self.ROIs) > 0:
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

    def WriteScaledDoseFile(self):
        description = "AccumulatedDose"
        self.ctPatient.WriteRTDose(self.basePath+"/AccumulatedDose.dcm", self.ctPatient.quantitiesOfInterest[0].array, "Gy", description)
        print(self.basePath+"/AccumulatedDose.dcm", " file saved.")

    def ShowDVHs(self, Xs = None, path=None):
        self.eval.PlotDVHs('Dose', self.basePath)
        self.eval.printMainResults('Dose', self.basePath)
        if Xs is not None:
            for x in Xs:
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
        if len(res) == 1:
            res = res[0]
        return res

    def GetPredictiveActivityCurves(self, metrics, structures, Xs, activityRange=[0.001, 1]):
        self.Xs = Xs
        for x in self.Xs:
            if x == 0:
                quantity = 'BED'
            else:
                quantity = 'EQD_' + str(x)
        activity = np.linspace(activityRange[0], activityRange[1], 50)
        results = []
        headers = 'Activity'
        for im, m in enumerate(metrics):
            headers += ',' + str(m) + '_' + str(structures[im])
            results.append(np.zeros(activity.shape))
        for ia, a in enumerate(activity):
            self.CalculateEQDXs(self.Xs, a)
            for im, m in enumerate(metrics):
                if 'EUEQDX' in m:
                    results[im][ia] = self.EUEQDX(self.Xs, structures[im])
                if 'MeanDose' in m:
                    results[im][ia] = self.eval.GetMeanDose(structures[im], quantity)
                if m[0] == 'D' and m[1].isnumeric():
                    num = float(m[1:])/100
                    results[im][ia] = self.eval.EvaluateD(num, structures[im], quantity)
                if m[0] == 'V' and m[1].isnumeric():
                    num = float(m[1:])
                    results[im][ia] = 100*self.eval.EvaluateV(num, structures[im], quantity)
        #headers += '\n'
        lines = []
        for ia, a in enumerate(activity):
            line = str(a.round(2))
            for ir, r in enumerate(results):
                line += ',' + str(r[ia].round(2))
            #line += '\n'
            lines.append(line)
        f = open(self.basePath + '/predictiveActivityCurves.csv', 'w+')
        f.write(headers)
        f.write('\n')
        for l in lines:
            f.write(l)
            f.write('\n')
        f.close()

    def PlotPredictiveActivityCurves(self, path = None):
        df = pd.read_csv(self.basePath + '/predictiveActivityCurves.csv')
        fig = plt.figure(dpi=300)
        ax = fig.add_subplot(1, 1, 1)
        colormap = plt.cm.nipy_spectral
        colors = [colormap(i) for i in np.linspace(1, 0, len(df.columns)-1)]
        ax.set_prop_cycle('color', colors)
        x = np.array(df.Activity)
        for i in range(1, len(df.columns)):
            y = np.array(df[df.columns[i]])
            if df.columns[i][0] == 'V':
                ax2 = ax.twinx()
                ax2.plot(x, y, label=df.columns[i])
                ax2.set_ylabel('Volume (%)')
                ax2.legend(loc=4)
                ax2.set_ylim([0, 100])
            else:
                ax.plot(x, y, label=df.columns[i])
        ax.set_xlabel('Activity (GBq)')
        ax.set_ylabel('EQDX (Gy_EQDX)')
        ax.set_xlim([0, np.max(x)])
        ax.set_ylim([0, None])
        ax.legend()
        plt.grid(alpha=0.7, ls='--')
        if path is not None:
            plt.savefig(path + "/predictiveActivityCurve.png")
            print(path + "/predictiveActivityCurve.png saved.")
        plt.show()
        return fig
