#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 4 14:16:34 2022
@author: alejandrobertolet
"""

from reports import *

# Patient directory
basePath = '/Users/ai925/Dropbox (Partners HealthCare)/RPT Project/BronchialSIR/Patient2/'

# Create report
rep = reportDosePerActivity(basePath, 'lung')
rep.OutputReport()