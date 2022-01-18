#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 10:33:18 2021
Last modified on Tue Jan 18 8:42 2022

@author: Mislav Bobić and Alejandro Bertolet
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

class EvaluationManager:
    def __init__(self, structureArrays, qoiArrays):
        self.structureArrays = structureArrays
        self.ROINames = [i for i in structureArrays]
        self.extraQoIs = []
        for q in qoiArrays:
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
                        cumVoxls += i
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
    
    def PlotDVHs(self):
        fig = plt.figure(figsize=(6,4), dpi=300)
        ax = fig.add_subplot(1,1,1)
        xAxis = self.DVHDataFrame['Dose']
        for ROIName in self.ROINames:
            yAxis = self.DVHDataFrame[ROIName]*100
            ax.plot(xAxis, yAxis, label=ROIName, linewidth=1.5)
            plt.xlabel('Dose [{}]'.format(self.doseUnit))
            plt.ylabel('Volume [%]')
            plt.legend()
            plt.grid(alpha=0.7, ls='--')
        return fig
        
    def SaveCSV(self, CSVPath):
        self.DVHDataFrame.to_csv(CSVPath)
        
    def LoadCSV(self, CSVPath):
        self.DVHDataFrame = pd.read_csv(CSVPath)
        
class QoIDistribution:
    def __init__(self, array, quantity, unit):
        self.array = array
        self.quantity = quantity
        self.unit = unit
