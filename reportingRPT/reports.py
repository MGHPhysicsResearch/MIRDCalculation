#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 4 14:04:02 2022
@author: alejandrobertolet
"""

from fpdf import FPDF
from DICOM_RT import DicomPatient
import datetime
import pandas as pd
from bioeffectRT import bioData

class reportDosePerActivity:
    def __init__(self, patientPath, site='liver'):
        self.path = patientPath
        self.site = site
        self._getPatientData()
        self._getBiophysicalParameters()
        self.SetUpPDF()
        self.LayoutContents()

    def _getPatientData(self):
        ctPath = self.path + "/CT/"
        self.patientCT = DicomPatient.PatientCT(ctPath)
        self.dosecsv = pd.read_csv(self.path + "/Dose.csv")
        self.eqd2csv = pd.read_csv(self.path + "/EQD_2.csv")
        self.dosecsv = self._reorderDataFrame(self.dosecsv)
        self.eqd2csv = self._reorderDataFrame(self.eqd2csv)

    def _getBiophysicalParameters(self):
        bioeffect = bioData.BioeffectData()
        self.alphabetaT = bioeffect.getAlphaBetaValue('t', self.site)
        self.trepT = bioeffect.getTRepValue('t', self.site)
        self.rois = []
        self.alphabetaN = []
        self.trepN = []
        for s in self.dosecsv.Structure:
            if 'tumor' not in s.lower():
                self.rois.append(s)
                self.alphabetaN.append(bioeffect.getAlphaBetaValue('n', s))
                self.trepN.append(bioeffect.getTRepValue('n',s))

    def SetUpPDF(self, font='Arial', fontsize=16, fontstyle='B'):
        self.pdf = FPDF()
        self.pdf.add_page()
        self.pdf.set_font(font, fontstyle, fontsize)

    def _reorderDataFrame(self, df):
        tumorsFirst = []
        indexes = []
        for i, s in enumerate(df.Structure):
            if 'tumor' in s.lower():
                tumorsFirst.append(s)
                indexes.append(i)
        for i, s in enumerate(df.Structure):
            if s not in tumorsFirst:
                tumorsFirst.append(s)
                indexes.append(i)
        df = df.loc[indexes]
        return df

    def LayoutContents(self, font='Arial', fontsize=10, fontstyle=''):
        self.pdf.cell(40, 10, 'Dosimetry Report (pretreatment)') # Title
        self.pdf.ln(10) # Line breaks
        self.pdf.set_font(font, fontstyle, fontsize)
        self.pdf.cell(5, 7, 'Patient Name: ' + str(self.patientCT.patientName))
        self.pdf.ln()  # Line breaks
        date = datetime.datetime.strptime(self.patientCT.studyDate, "%Y%m%d")
        self.pdf.cell(5, 7, 'Study date: ' + date.strftime("%Y-%m-%d"))
        self.pdf.ln()
        now = datetime.datetime.now()
        self.pdf.cell(5, 7, 'Dosimetry date: ' + str(now.strftime("%Y-%m-%d")))
        self.pdf.ln()
        self.pdf.set_font(font, 'B', fontsize)
        self.pdf.cell(10, 10, 'Parameters used for EQD_2 calculation')
        self.pdf.set_font(font, fontstyle, fontsize)
        self.pdf.ln()
        self.pdf.cell(10, 5, 'Alpha/Beta (tumor): ' + str(self.alphabetaT) + ' Gy. Trep (tumor): ' + str(self.trepT / 3600) + ' h')
        for i, r in enumerate(self.rois):
            self.pdf.ln()
            self.pdf.cell(10, 5, 'Alpha/Beta (' + r + '): ' + str(self.alphabetaN[i]) + ' Gy. Trep (' + r + '): ' + str(self.trepN[i] / 3600) + ' h')
        self.pdf.ln()
        self.pdf.set_font(font, 'B', fontsize)
        self.pdf.cell(40, 10, 'Dose-Volume Histogram')  # Title
        self.pdf.ln()
        self.pdf.image(self.path + '/Dose.png', 30, None, 120)
        self.pdf.ln()
        self.pdf.set_font(font, 'B', fontsize)
        self.pdf.cell(40, 10, 'EQD2-Volume Histogram')  # Title
        self.pdf.ln()
        self.pdf.image(self.path + '/EQD_2.png', 30, None, 120)
        self.pdf.ln()
        self.pdf.cell(40, 10, 'EQD2-Axial plane')  # Title
        self.pdf.ln()
        self.pdf.image(self.path + '/EQD2_Axial.png', 30, None, 120)
        self.pdf.cell(40, 10, 'Summary - Absorbed Dose (Gy/GBq)')  # Title
        self.pdf.ln()
        dosecsvpdf = self.dosecsv.reset_index()
        numericCols = dosecsvpdf.select_dtypes(include='number').columns
        dosecsvpdf[numericCols] = dosecsvpdf[numericCols].round(2)
        self._outputdftopdf(dosecsvpdf)
        self.pdf.ln()
        self.pdf.set_font(font, 'B', fontsize)
        self.pdf.cell(40, 10, r'Summary - EQD_2 (Gy/GBq)')  # Title
        self.pdf.ln()
        eqd2csvpdf = self.eqd2csv.reset_index()
        numericCols = eqd2csvpdf.select_dtypes(include='number').columns
        eqd2csvpdf[numericCols] = eqd2csvpdf[numericCols].round(2)
        self._outputdftopdf(eqd2csvpdf)


    def OutputReport(self):
        self.pdf.output(self.path + '/' + str(self.patientCT.patientName) + '-DosimetryReport.pdf', 'F')

    def _outputdftopdf(self, df, font='Arial', fontsize=8, fontstyle=''):
        # A cell is a rectangular area, possibly framed, which contains some text
        # Set the width and height of cell
        table_firstcell_width = 40
        table_cell_width = 15
        table_cell_height = 6
        self.pdf.set_font(font, fontstyle, fontsize)
        # Loop over to print column names
        cols = df.columns[1:]
        for i, col in enumerate(cols):
            if i == 0:
                self.pdf.cell(table_firstcell_width, table_cell_height, col, align='C', border=1)
            else:
                self.pdf.cell(table_cell_width, table_cell_height, col, align='C', border=1)
        # Line break
        self.pdf.ln(table_cell_height)
        # Select a font as Arial, regular, 10
        self.pdf.set_font('Arial', '', 10)
        # Loop over to print each data in the table
        for row in df.itertuples():
            for i, col in enumerate(cols):
                value = str(getattr(row, col))
                if i == 0:
                    self.pdf.cell(table_firstcell_width, table_cell_height, value, align='C', border=1)
                else:
                    self.pdf.cell(table_cell_width, table_cell_height, value, align='C', border=1)
            self.pdf.ln(table_cell_height)