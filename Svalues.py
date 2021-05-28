#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 28 11:51:20 2021

@author: alejandrobertolet
"""

from os import listdir
import numpy as np

class SValuesData:
    def __init__(self, radionuclide = '', datapath = 'VoxelSValues'):
        self.datapath = datapath
        self.datafiles = listdir(datapath)
        self.isthereradionuclide = False
        if radionuclide != '':
            self.SetRadionuclide(radionuclide)
        
    def SetRadionuclide(self, radionuclide):
        if not radionuclide[0].isnumeric():
            radionuclide = self.__reverseRadionuclideName(radionuclide)
        self.datasets = []
        self.maximumDistanceInVoxels = 6
        for filename in self.datafiles:
            if filename[0:len(radionuclide)] == radionuclide:
                self.datasets.append(SValueDataset(self.datapath + '/' + filename))
                shape = self.datasets[len(self.datasets)-1].Svalues.shape
                if shape[0] < self.maximumDistanceInVoxels:
                    self.maximumDistanceInVoxels = shape[0]
                if shape[1] < self.maximumDistanceInVoxels:
                    self.maximumDistanceInVoxels = shape[1]
                if shape[2] < self.maximumDistanceInVoxels:
                    self.maximumDistanceInVoxels = shape[2]
        self.datasets.sort(key=lambda x:x.voxelsize)
        self.isthereradionuclide = True

    def GetSValue(self, voxelSize, voxelsX, voxelsY, voxelsZ, tissue = 'Soft'):
        if not self.isthereradionuclide:
            print("No radionuclide was specified or no data was founded.")
            return
        # Get all values for tissue and voxels X, Y and Z specified
        voxelsizes = []
        svalues = []
        for ds in self.datasets:
            if ds.tissue == tissue:
                voxelsizes.append(ds.voxelsize)
                svalues.append(ds.Svalues[voxelsX, voxelsY, voxelsZ])
        self.voxelSizes = np.array(voxelsizes)
        self.Svalues = np.array(svalues)
        return np.interp(voxelSize, self.voxelSizes, self.Svalues)
    
    def __reverseRadionuclideName(self, name):
        number = ''
        alpha = ''
        for i in range(0, len(name)):
            if name[i].isnumeric():
                number = number + name[i]
            if name[i].isalpha():
                alpha = alpha + name[i]
        return number + alpha
    
    
class SValueDataset:
    def __init__(self, filename):
        file = open(filename, 'r', encoding='latin1')
        lines = file.readlines()
        headers = lines[0].split(' - ')
        self.radionuclide = headers[0]
        self.voxelsize = self.__getnumber(headers[1])
        self.tissue = headers[2].split()[0]
        x = []
        y = []
        z = []
        S = []
        for i in range(2, len(lines)):
            row = lines[i].split('\t')
            x.append(int(row[0]))
            y.append(int(row[1]))
            z.append(int(row[2]))
            S.append(row[3].split()[0])
        self.Svalues = np.zeros([max(x)-min(x)+1, max(y)-min(y)+1, max(z)-min(z)+1])
        for i, s in enumerate(S):
            self.Svalues[x[i], y[i], z[i]] = s

    def __getnumber(self, name):
        number = ''
        for i in range(0, len(name)):
            if name[i].isnumeric() or name[i] == '.':
                number = number + name[i]
        return float(number)
            