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

    def getAlphaBetaValue(self, type, site):
        stdSite = self.__checkOtherNamesForSite(site)
        if len(self.abData[self.abData.type == type][self.abData.site == stdSite]) == 0:
            stdSite = 'generic'
        v = self.abData[self.abData.type == type][self.abData.site == stdSite].value.squeeze()
        return v

    def getTRepValue(self, type, site):
        stdSite = self.__checkOtherNamesForSite(site)
        if len(self.trepData[self.trepData.type == type][self.trepData.site == stdSite]) == 0:
            stdSite = 'generic'
        v = self.trepData[self.trepData.type == type][self.trepData.site == stdSite].value.squeeze()
        return v

    def __checkOtherNamesForSite(self, site):
        ls = site.lower()
        if 'liver' in ls:
            ls = 'liver'
        if 'lung_l' in ls or 'lung_r' in ls:
            ls = 'lung'
        return ls