#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 2 08:42:03 2022

@authors: mjlindsey, alejandrobertolet
"""
import os
import pandas as pd
import pkg_resources
DATA_PATH = pkg_resources.resource_filename('bioeffectRT', 'data/')

class BioeffectData:
    def __init__(self, datapath=''):
        self.datapath = datapath
        if datapath == '':
            self.datapath = DATA_PATH
        if os.path.isdir(self.datapath):
            self.datafiles = os.listdir(self.datapath)
        self.__loadData()

    def __loadData(self):
        self.abData = pd.read_csv(self.datapath + "/alphabeta.csv")
        self.trepData = pd.read_csv(self.datapath + "/trep.csv")
        self.alphaData = pd.read_csv(self.datapath + "/alpha.csv")

    def getAlphaValue(self, type, site):
        stdSite = self.__checkOtherNamesForSite(site)
        dataOfType = self.alphaData[self.alphaData.type == type]
        if len(dataOfType[dataOfType.site == stdSite]) == 0:
            stdSite = 'generic'
        v = dataOfType[dataOfType.site == stdSite].value.squeeze()
        return v

    def getAlphaBetaValue(self, type, site):
        stdSite = self.__checkOtherNamesForSite(site)
        dataOfType = self.abData[self.abData.type == type]
        if len(dataOfType[dataOfType.site == stdSite]) == 0:
            stdSite = 'generic'
        v = dataOfType[dataOfType.site == stdSite].value.squeeze()
        return v

    def getTRepValue(self, type, site):
        stdSite = self.__checkOtherNamesForSite(site)
        dataOfType = self.trepData[self.trepData.type == type]
        if len(dataOfType[dataOfType.site == stdSite]) == 0:
            stdSite = 'generic'
        v = dataOfType[dataOfType.site == stdSite].value.squeeze() * 3600 # From h to s
        return v

    def __checkOtherNamesForSite(self, site):
        ls = site.lower()
        if 'liver' in ls:
            ls = 'liver'
        if 'lung_l' in ls or 'lung_r' in ls or 'right lung' in ls or 'left lung' in ls:
            ls = 'lung'
        return ls