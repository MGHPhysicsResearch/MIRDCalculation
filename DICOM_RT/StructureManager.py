#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 4 9:56:19 2022

@author: alejandrobertolet
"""
import cv2
import numpy as np
import pydicom.uid
from pydicom.dataset import Dataset, FileDataset, FileMetaDataset
from pydicom.sequence import Sequence


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


class RTStructWriter:
    def __init__(self, structs, slices, structuresetname):
        self.structs = structs
        self.slices = slices
        self.structuresetname = structuresetname

    def write(self, path):
        # Create a new RTSTRUCT DICOM dataset
        file_meta = FileMetaDataset()
        file_meta.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.481.3'
        file_meta.MediaStorageSOPInstanceUID = pydicom.uid.generate_uid()
        file_meta.ImplementationClassUID = pydicom.uid.PYDICOM_IMPLEMENTATION_UID

        rtstruct = FileDataset(self.path, {}, file_meta=file_meta, preamble=b"\0" * 128)
        rtstruct.is_implicit_VR = True
        rtstruct.is_little_endian = True
        # Set necessary DICOM tags
        rtstruct.SOPClassUID = '1.2.840.10008.5.1.4.1.1.481.3'
        rtstruct.SOPInstanceUID = rtstruct.file_meta.MediaStorageSOPInstanceUID
        rtstruct.StudyInstanceUID = self.slices[0].StudyInstanceUID
        rtstruct.SeriesInstanceUID = pydicom.uid.generate_uid()
        rtstruct.Modality = 'RTSTRUCT'
        rtstruct.Manufacturer = 'MIRDCalculation'
        rtstruct.StructureSetLabel = self.structuresetname
        rtstruct.StructureSetName = self.structuresetname
        rtstruct.PatientID = self.slices[0].PatientID
        rtstruct.PatientName = self.slices[0].PatientName
        rtstruct.PatientBirthDate = self.slices[0].PatientBirthDate
        rtstruct.PatientSex = self.slices[0].PatientSex
        rtstruct.StudyID = self.slices[0].StudyID
        rtstruct.StudyDate = self.slices[0].StudyDate
        rtstruct.StudyTime = self.slices[0].StudyTime
        rtstruct.FrameOfReferenceUID = self.slices[0].FrameOfReferenceUID

        # Create a structure set ROI sequence
        structure_set_roi_sequence = Sequence()
        # Create the ROI contour sequence
        roi_contour_sequence = Sequence()

        # self.struct is a dictionary
        for name, struct in self.structs.items():
            contours = []
            for z in range(struct.shape[2]):
                sliceMask = struct[:, :, z].astype(np.uint8) * 255
                contours_z, _ = cv2.findContours(sliceMask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
                if contours_z:
                    contours.append((z, contours_z))
            # Create a contour sequence
            contour_sequence = Sequence()
            for z, contours_z in contours:
                for contour in contours_z:
                    contour_data = []
                    for point in contour[:, 0, :]:
                        x, y = point
                        slice_position = self.slices[z].ImagePositionPatient
                        pixel_spacing = self.slices[z].PixelSpacing
                        contour_data.extend([x * pixel_spacing[0] + slice_position[0],
                                             y * pixel_spacing[1] + slice_position[1],
                                             slice_position[2]])
                    contour_item = Dataset()
                    contour_item.ContourGeometricType = 'CLOSED_PLANAR'
                    contour_item.NumberOfContourPoints = len(contour)
                    contour_item.ContourData = contour_data
                    contour_item.ContourImageSequence = Sequence([Dataset()])
                    contour_item.ContourImageSequence[0].ReferencedSOPClassUID = self.slices[z].SOPClassUID
                    contour_item.ContourImageSequence[0].ReferencedSOPInstanceUID = self.slices[z].SOPInstanceUID
                    contour_sequence.append(contour_item)
            structure_set_roi_item = Dataset()
            structure_set_roi_item.ROINumber = 1
            structure_set_roi_item.ROIName = name
            structure_set_roi_item.ROIGenerationAlgorithm = ''
            structure_set_roi_sequence.append(structure_set_roi_item)
            roi_contour_item = Dataset()
            roi_contour_item.ReferencedROINumber = 1
            roi_contour_item.ROIDisplayColor = [255, 0, 0]
            roi_contour_item.ContourSequence = contour_sequence
            roi_contour_sequence.append(roi_contour_item)

        # Set the RTStruct sequences
        rtstruct.StructureSetROISequence = structure_set_roi_sequence
        rtstruct.ROIContourSequence = roi_contour_sequence

        # Save the RTStruct file
        rtstruct.save_as(path)
