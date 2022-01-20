#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 10:33:18 2021
Last modified on Tue Jan 18 14:47 2022

@author: Mislav Bobić and Alejandro Bertolet
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
import DICOM_RT.DicomPatient
from scipy.interpolate import interp1d

class EvaluationManager:
    def __init__(self, dicomPat):
        self.patient = dicomPat
        self.structureArrays = dicomPat.structures3D
        self.ROINames = [i for i in self.structureArrays]
        self.extraQoIs = []
        for q in dicomPat.quantitiesOfInterest:
            if q.quantity == 'Dose':
                self.doseArray = q.array
                self.doseUnit = q.unit
            else:
                self.extraQoIs.append(q)
        self.maxDose =  np.max(self.doseArray)
        self.CalculateDVHs()
        
    def CalculateDVHs(self, numBins=1000):
        self.DVHDataFrame = pd.DataFrame()
        self.maxDose = np.max(self.doseArray)
        if len(self.extraQoIs) > 0:
            self.QoiDVHDataFrames = []
            for q in self.extraQoIs:
                self.QoiDVHDataFrames.append(pd.DataFrame())
        for ROIName in self.ROINames:
            structureDose = self.GetStructureDose(ROIName)
            histRange = (0, round(self.maxDose))
            histBins  = round(self.maxDose) if numBins is None else numBins-1
            differentialDVH, doseValues = np.histogram(structureDose, histBins, histRange)
            cumulativeDVH = np.zeros(len(differentialDVH)+1)
            index = 0
            cumulativeVoxels = 0
            for i in differentialDVH[::-1]:
                np.put(cumulativeDVH, index, cumulativeVoxels)
                cumulativeVoxels += i
                index += 1
            np.put(cumulativeDVH, index, cumulativeVoxels)
            cumulativeDVH = cumulativeDVH[::-1]/cumulativeVoxels
            self.DVHDataFrame[ROIName] = cumulativeDVH
            if len(self.extraQoIs) > 0:
                for iq, q in enumerate(self.extraQoIs):
                    structQ = q.array[self.structureArrays[ROIName]]
                    histR = (0, round(np.max(q.array)))
                    histB = round(np.max(q.array)) if numBins is None else numBins-1
                    diffDVH, qValues = np.histogram(structQ, histB, histR)
                    cumDVH = np.zeros(len(diffDVH)+1)
                    ind = 0
                    cumVoxls = 0
                    for j in diffDVH[::-1]:
                        np.put(cumDVH, ind, cumVoxls)
                        cumVoxls += j
                        ind += 1
                    np.put(cumDVH, ind, cumVoxls)
                    cumDVH = cumDVH[::-1]/cumVoxls
                    self.QoiDVHDataFrames[iq][ROIName] = cumDVH
        self.DVHDataFrame.insert(0, 'Dose', doseValues)
        if len(self.extraQoIs) > 0:
            for i, q in enumerate(self.extraQoIs):
                self.QoiDVHDataFrames[i].insert(0, q.quantity, qValues)
        print('DVHs calculated.')

    def GetStructureDose(self, ROIName):
        return self.doseArray[self.structureArrays[ROIName]]
    
    def GetMeanDose(self, ROIName):
        return np.average(self.GetStructureDose(ROIName))
    
    def GetMaxDose(self, ROIName):
        return np.max(self.GetStructureDose(ROIName))
    
    def GetMinDose(self, ROIName):
        return np.min(self.GetStructureDose(ROIName))

    def EvaluateV(self, dose, ROIName):
        DVHDose = self.DVHDataFrame['Dose']
        DVHVolume = self.DVHDataFrame[ROIName]
        f = interp1d(DVHDose, DVHVolume)
        return float(f(dose))

    def EvaluateD(self, volume, ROIName):
        DVHDose = self.DVHDataFrame['Dose']
        DVHVolume = self.DVHDataFrame[ROIName]
        f = interp1d(DVHVolume, DVHDose)
        return float(f(volume))
    
    def PlotDVHs(self, quantity = "Dose"):
        fig = plt.figure(figsize=(18,12), dpi=300)
        ax = fig.add_subplot(1,1,1)
        xAxis = None
        if quantity == "Dose":
            xAxis = self.DVHDataFrame['Dose']
            DVHDataFrame = self.DVHDataFrame
            unit = self.doseUnit
        else:
            for i, q in enumerate(self.QoiDVHDataFrames):
                if self.extraQoIs[i].quantity == quantity:
                    xAxis = q[quantity]
                    DVHDataFrame = q
                    unit = self.extraQoIs[i].unit
            if xAxis is None:
                print("Quantity " + quantity + " could not be found.")
                return
        for ROIName in self.ROINames:
            yAxis = DVHDataFrame[ROIName]*100
            ax.plot(xAxis, yAxis, label=ROIName, linewidth=1.5)
            plt.xlabel((quantity + '[{}]').format(unit))
            plt.ylabel('Volume [%]')
            plt.legend()
            plt.grid(alpha=0.7, ls='--')
        return fig
        
    def SaveCSV(self, CSVPath):
        self.DVHDataFrame.to_csv(CSVPath)
        
    def LoadCSV(self, CSVPath):
        self.DVHDataFrame = pd.read_csv(CSVPath)
        
    def SaveVoxelByVoxelCSV(self, CSVPath):
        f = open(CSVPath, 'w')
        writer = csv.writer(f)
        for ROIName in self.ROINames:
            writer.writerow([ROIName])
            headers = ['VoxelID', 'Dose (' + self.doseUnit + ')']
            qrois = []
            for q in self.extraQoIs:
                headers.append(q.quantity + " (" + q.unit + " )")
                qrois.append(q.array[self.structureArrays[ROIName]])
            writer.writerow(headers)
            dose = self.doseArray[self.structureArrays[ROIName]]
            indexes = np.where(self.structureArrays[ROIName])
            shape = self.structureArrays[ROIName].shape
            inds = []
            for i, ind in enumerate(indexes[0]):
                inds.append(indexes[0][i]*shape[1]*shape[2]+indexes[1][i]*shape[2] + indexes[2][i])
            for i, ind in enumerate(inds):
                row = []
                row.append(ind)
                row.append(dose[i])
                for qroi in qrois:
                    row.append(qroi[i])
                writer.writerow(row)
        f.close()
