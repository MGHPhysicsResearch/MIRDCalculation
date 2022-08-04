#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 4 9:56:19 2022

@author: alejandrobertolet
"""
import numpy as np

class Operations:
    def __init__(self, struct1, struct2=None):
        self.struct1 = struct1
        self.struct2 = struct2

    @staticmethod
    def Union(struct1, struct2):
        return np.logical_or(struct1, struct2)

    @staticmethod
    def Subtraction(struct1, struct2):
        return np.logical_and(struct1, ~struct2)

    @staticmethod
    def Intersection(struct1, struct2):
        return np.logical_and(struct1, struct2)