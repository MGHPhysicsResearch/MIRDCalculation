#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 28 11:51:20 2021

@author: alejandrobertolet
"""

from os import listdir, path
import numpy as np
import matplotlib.pyplot as plt
import csv
import pkg_resources

DATA_PATH = pkg_resources.resource_filename('MIRD', 'VoxelSValues/')

class SValuesData:
    def __init__(self, radionuclide = '', datapath = '', dataTOPASpath = '../TOPASSvalues'):
        self.datapath = datapath
        self.datafiles = []
        if datapath == '':
            self.datapath = DATA_PATH
        if path.isdir(self.datapath):
            self.datafiles = listdir(self.datapath)
        self.TOPASpath = dataTOPASpath
        self.TOPASfiles = []
        if path.isdir(dataTOPASpath):
            self.TOPASfiles = listdir(dataTOPASpath)
        self.isthereradionuclide = False
        self.isthereTOPASdata = False
        if radionuclide != '':
            self.SetRadionuclide(radionuclide)

    def SetRadionuclide(self, radionuclide):
        rn = Radionuclide(radionuclide)
        self.halflife = rn.halfLife
        self.decayConstant = rn.decayConstant
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
        if len(self.datasets) > 0:
            self.datasets.sort(key=lambda x:x.voxelSize)
            self.isthereradionuclide = True

        self.TOPASdatasets = []
        for filename in self.TOPASfiles:
            cumActMBqs = 10
            if filename.find('MBqs') != -1:
                i = filename.find('MBqs')-1
                nStr = ''
                nChar = filename[i]
                while nChar.isnumeric() or nChar == '.':
                    nStr = nStr + nChar
                    i = i-1
                    nChar = filename[i]
                cumActMBqs = float(nStr[::-1])
            if filename[0:len(radionuclide)] == radionuclide:
                self.TOPASdatasets.append(SValueDatasetTOPAS(self.TOPASpath + '/' + filename, cumActMBqs))
        if len(self.TOPASdatasets) > 0:
            self.isthereTOPASdata = True
            self.TOPASdatasets.sort(key = lambda x:x.voxelSize)

    def GetSValue(self, voxelSize, voxelsX, voxelsY, voxelsZ, tissue = 'Soft', source = 'Lanconelli', physics = 'standard'):
        if source == 'Lanconelli':
            if not self.isthereradionuclide:
                print("No radionuclide was specified or no data was founded.")
                return
            # Get all values for tissue and voxels X, Y and Z specified
            voxelsizes = []
            svalues = []
            for ds in self.datasets:
                if ds.tissue == tissue:
                    voxelsizes.append(ds.voxelSize)
                    svalues.append(ds.Svalues[voxelsX, voxelsY, voxelsZ])
            self.voxelSizes = np.array(voxelsizes)
            self.Svalues = np.array(svalues)
        elif source == 'TOPAS':
            if not self.isthereTOPASdata:
                print("No data was found for TOPAS simulations.")
                return
            voxelsizes = []
            svalues = []
            for ds in self.TOPASdatasets:
                if ds.tissue == tissue and ds.physics == physics:
                    voxelsizes.append(ds.voxelSize)
                    svalues.append(ds.Svalues[voxelsX, voxelsY, voxelsZ])
            self.voxelSizes = np.array(voxelsizes)
            self.Svalues = np.array(svalues)
        return np.interp(voxelSize, self.voxelSizes, self.Svalues)

    def GetStdSValue(self, voxelSize, voxelsX, voxelsY, voxelsZ, tissue = 'Soft', physics = 'standard'):
        if not self.isthereTOPASdata:
            print("No data was found for TOPAS simulations.")
            return
        for ds in self.TOPASdatasets:
            if ds.tissue == tissue and round(ds.voxelSize, 3) == round(voxelSize, 3) and ds.physics == physics:
                return ds.StdSvalues[voxelsX, voxelsY, voxelsZ]

    def plot1D(self, voxelSize, tissue = 'Soft', source = 'Lanconelli', physics = 'standard', color = 'black',
               mk = '*', linestyle = '-'):
        x = []
        S = []
        for i in range(0, 6):
            x.append(i * voxelSize)
            S.append(self.GetSValue(voxelSize, 0, 0, i, tissue, source, physics))
        if source == 'Lanconelli':
            plt.errorbar(x, S, marker=mk, label="Voxel size = " + str(voxelSize) + " mm - Lanconelli et al.",
                        c=color, ls=linestyle)
        if source == 'TOPAS':
            errorS = []
            for i in range(0, 6):
                errorS.append(self.GetStdSValue(voxelSize, i, 0, 0, tissue, physics))
            plt.errorbar(x, S, errorS, fmt=mk, c=color, markersize=5,
                         label="Voxel size = " + str(voxelSize) + " mm - TOPAS-" + physics, ls=linestyle)
        plt.xlabel('Distance (mm)')
        plt.ylabel('S value (mGy/(MBq s))')
        plt.xlim([-max(x)*0.05, max(x)*1.05])
        plt.yscale('log')
    
class SValueDataset:
    def __init__(self, filename):
        file = open(filename, 'r', encoding='latin1')
        lines = file.readlines()
        headers = lines[0].split(' - ')
        self.radionuclide = headers[0]
        self.voxelSize = self.__getnumber(headers[1])
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

class SValueDatasetTOPAS:
    def __init__(self, filename, cumActMBqs = 10):
        xd = []
        yd = []
        zd = []
        dose = []
        meanDoseEvt = []
        stdDoseEvt = []
        with open(filename) as csvfile:
            tissue = filename[-8:-4]
            if tissue == 'soft':
                self.tissue = 'Soft'
            if tissue == 'bone':
                self.tissue = 'Bone'
            #Physics
            if filename.find('opt4') == -1:
                self.physics = 'standard'
            else:
                self.physics = 'option4'
            spamreader = csv.reader(csvfile, delimiter=',')
            for row in spamreader:
                if row[0][0] != "#":
                    xd.append(int(row[0]))
                    yd.append(int(row[1]))
                    zd.append(int(row[2]))
                    dose.append(float(row[3]))
                    meanDoseEvt.append(float(row[4]))
                    stdDoseEvt.append(float(row[5]))
                else:
                    if row[0][2] == "X":
                        line = row[0]
                        self.voxelSize = float(line[line.find('of')+3:-3])
                        unit = line[-2:]
                        if unit == 'cm':
                            self.voxelSize = self.voxelSize * 10
                            unit = 'mm'
        self.AddSpecularResults(xd, yd, zd, dose, meanDoseEvt, stdDoseEvt)
        self.ConvertIntoSValue(cumActMBqs)
        
    def AddSpecularResults(self, xd, yd, zd, dose, meanDoseEvt, stdDoseEvt):
        DoseTOPAS = np.zeros([max(xd)-min(xd)+1, max(yd)-min(yd)+1, max(zd)-min(zd)+1])
        StdDTOPAS = np.zeros([max(xd)-min(xd)+1, max(yd)-min(yd)+1, max(zd)-min(zd)+1])
        for i in range(0, len(dose)):
            DoseTOPAS[xd[i], yd[i], zd[i]] = dose[i]
            StdDTOPAS[xd[i], yd[i], zd[i]] = stdDoseEvt[i]
        binsx = np.max(xd) - np.min(xd)
        binsy = np.max(yd) - np.min(yd)
        binsz = np.max(zd) - np.min(zd)
        centerx = np.max(xd) - int(binsx/2)
        centery = np.max(yd) - int(binsy/2)
        centerz = np.max(zd) - int(binsz/2)
        x = []
        y = []
        z = []
        d = []
        sDe = []
        for ix in range(0, int(binsx/2)+1):
            for iy in range(0, int(binsy/2)+1):
                for iz in range(0, int(binsz/2)+1):
                    x.append(ix)
                    y.append(iy)
                    z.append(iz)
                    xp = centerx+ix
                    yp = centery+iy
                    zp = centerz+iz
                    xn = centerx-ix
                    yn = centery-iy
                    zn = centerz-iz
                    d.append((DoseTOPAS[xp, yp, zp] + DoseTOPAS[xn, yn, zn])/2)
                    sDe.append(np.sqrt((StdDTOPAS[xp, yp, zp]**2 + StdDTOPAS[xn, yn, zn]**2)/2))
        self.Dvalues = np.zeros([max(x)-min(x)+1, max(y)-min(y)+1, max(z)-min(z)+1])
        self.StdDvalues = np.zeros([max(x)-min(x)+1, max(y)-min(y)+1, max(z)-min(z)+1])
        for i in range(0, len(d)):
            self.Dvalues[x[i], y[i], z[i]] = d[i]
            self.StdDvalues[x[i], y[i], z[i]] = sDe[i]
            
    def ConvertIntoSValue(self, cumAct = 10, unit='mGy/(MBq s)'):
        if unit == 'mGy/(MBq s)':
            self.Svalues = self.Dvalues / cumAct * 1000
            errorDValues = self.StdDvalues * (cumAct * 1e6) / np.sqrt(cumAct * 1e6)
            self.StdSvalues = errorDValues / cumAct * 1000

class Radionuclide:
    def __init__(self, name):
        if not name[0].isnumeric():
            self.name = self.__reverseRadionuclideName(name)
        else:
            self.name = name

    @property
    def halfLife(self):
        # Units for half life
        hour = 3600
        day = 24*hour
        if self.name == '89Sr':
            self._halflife = 50.5*day
        elif self.name == '90Y':
            self._halflife = 64.2*hour
        elif self.name == '131I':
            self._halflife = 8.02*day
        elif self.name == '153Sm':
            self._halflife = 46.3*hour
        elif self.name == '177Lu':
            self._halflife = 6.647*day
        elif self.name == '186Re':
            self._halflife = 3.8*day
        elif self.name == '188Re':
            self._halflife = 16.98*hour
        return self._halflife

    @property
    def decayConstant(self):
        return np.log(2) / self.halflife

    def __reverseRadionuclideName(self, name):
        number = ''
        alpha = ''
        for i in range(0, len(name)):
            if name[i].isnumeric():
                number = number + name[i]
            if name[i].isalpha():
                alpha = alpha + name[i]
        return number + alpha

def RunTests():
    y90 = SValuesData('Y90')
    vsizes = [2.21, 2.33, 2.4, 3.0, 3.59, 3.9, 4.0, 4.42, 4.8, 5.0, 6.0, 6.8, 9.28]
    vsizes = [2.21, 9.28]
    colors = [[[0,0,1,0.5], [0,0,1,0.75], [0,0,1,1]], [[1,0,0,0.5], [1,0,0,0.75], [1,0,0,1]]]
    lss = [':', '--', '-']
    for i, vs in enumerate(vsizes):
        y90.plot1D(vs, 'Soft', 'Lanconelli', '', colors[i][0], '.', lss[0])
        #y90.plot1D(vs, 'Soft', 'TOPAS', 'standard', colors[i][1], '*', lss[1])
        #y90.plot1D(vs, 'Soft', 'TOPAS', 'option4', colors[i][2], 'd', lss[2])
    plt.legend()
    plt.grid()

#RunTests()
