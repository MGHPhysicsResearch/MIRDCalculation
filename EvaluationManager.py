#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 10:33:18 2021

@author: Mislav
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

class EvaluationManager:
    def __init__(self, structureArrays, doseArray):
        self.structureArrays = structureArrays
        self.ROINames = [i for i in structureArrays]
        self.doseArray = doseArray
        self.maxDose =  np.max(doseArray)
        self.CalculateDVHs()
        
    def CalculateDVHs(self, numBins=None):
        self.DVHDataFrame = pd.DataFrame()
        self.maxDose = np.max(self.doseArray)
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
        self.DVHDataFrame.insert(0, 'Dose', doseValues)
        print('DVHs calculated.')

    def GetStructureDose(self, ROIName):
        return self.doseArray[self.structureArrays[ROIName]]
    
    def GetMeanDose(self, ROIName):
        return np.average(self.GetStructureDose(ROIName))
    
    def GetMaxDose(self, ROIName):
        return np.max(self.GetStructureDose(ROIName))
    
    def GetMinDose(self, ROIName):
        return np.min(self.GetStructureDose(ROIName))

    def EvaluateV(self, stat, ROIName, prescription):
        DVHDose = self.DVHDataFrame['Dose']
        DVHVolume = self.DVHDataFrame[ROIName]
        f = interp1d(DVHDose, DVHVolume)
        return float(f(stat/100*prescription))

    def EvaluateD(self, stat, ROIName, prescription):
        DVHDose = self.DVHDataFrame['Dose']
        DVHVolume = self.DVHDataFrame[ROIName]
        f = interp1d(DVHVolume, DVHDose)
        return float(f(stat/100))/prescription
    
    def PlotDVHs(self, doseUnit=None):
        fig = plt.figure(figsize=(6,4), dpi=300)
        ax = fig.add_subplot(1,1,1)
        xAxis = self.DVHDataFrame['Dose']
        for ROIName in self.ROINames:
            yAxis = self.DVHDataFrame[ROIName]*100
            ax.plot(xAxis, yAxis, label=ROIName, linewidth=1.5)
            doseUnit = '?' if doseUnit is None else doseUnit
            plt.xlabel('Dose [{}]'.format(doseUnit))
            plt.ylabel('Volume [%]')
            plt.legend()
            plt.grid(alpha=0.7, ls='--')
        return fig
        
    def SaveCSV(self, CSVPath):
        self.DVHDataFrame.to_csv(CSVPath)
        
    def LoadCSV(self, CSVPath):
        self.DVHDataFrame = pd.read_csv(CSVPath)
