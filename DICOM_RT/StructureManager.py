#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 4 9:56:19 2022

@author: alejandrobertolet
"""
import cv2
import numpy as np
import pydicom.uid
from pydicom.dataset import Dataset, FileDataset
from pydicom.sequence import Sequence
from datetime import datetime

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
    def __init__(self, structs, patient, structuresetname):
        self.structs = structs
        self.patient = patient
        self.structuresetname = structuresetname

    def write(self, path):
        slices = self.patient.slices
        filename = path + '/RTSTRUCT.dcm'
        # Create a new RTSTRUCT DICOM dataset
        file_meta = Dataset()
        file_meta.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.481.3'
        file_meta.MediaStorageSOPInstanceUID = pydicom.uid.generate_uid()
        file_meta.ImplementationClassUID = pydicom.uid.PYDICOM_IMPLEMENTATION_UID
        file_meta.ImplementationVersionName = "PYDICOM " + pydicom.__version__
        file_meta.TransferSyntaxUID = pydicom.uid.ExplicitVRLittleEndian

        rtstruct = FileDataset(filename, {}, file_meta=file_meta, preamble=b"\0" * 128)
        rtstruct.is_implicit_VR = False
        rtstruct.is_little_endian = True
        # Set necessary DICOM tags
        rtstruct.SOPClassUID = '1.2.840.10008.5.1.4.1.1.481.3'
        rtstruct.SOPInstanceUID = rtstruct.file_meta.MediaStorageSOPInstanceUID
        rtstruct.StudyInstanceUID = slices[0].StudyInstanceUID
        rtstruct.SeriesInstanceUID = pydicom.uid.generate_uid()
        rtstruct.Modality = 'RTSTRUCT'
        rtstruct.Manufacturer = 'MIRDCalculation'
        rtstruct.StructureSetLabel = self.structuresetname
        rtstruct.StructureSetName = self.structuresetname
        rtstruct.PatientID = slices[0].PatientID
        rtstruct.PatientName = slices[0].PatientName
        rtstruct.PatientBirthDate = slices[0].PatientBirthDate
        rtstruct.PatientSex = slices[0].PatientSex
        rtstruct.StudyID = slices[0].StudyID
        rtstruct.StudyDate = slices[0].StudyDate
        rtstruct.StudyTime = slices[0].StudyTime
        rtstruct.FrameOfReferenceUID = slices[0].FrameOfReferenceUID
        rtstruct.ReferencedFrameOfReferenceSequence = Sequence()
        referenced_frame_of_reference_item = Dataset()
        referenced_frame_of_reference_item.FrameOfReferenceUID = slices[0].FrameOfReferenceUID
        rtstruct.ReferencedFrameOfReferenceSequence.append(referenced_frame_of_reference_item)
        rtstruct.SeriesDescription = 'RT Structure Set'
        now = datetime.now()
        rtstruct.StructureSetDate = now.strftime('%Y%m%d')
        rtstruct.StructureSetTime = now.strftime('%H%M%S')

        # Create a structure set ROI sequence
        structure_set_roi_sequence = Sequence()
        # Create the ROI contour sequence
        roi_contour_sequence = Sequence()
        structures = {}
        # self.struct is a dictionary
        iS = 1
        for name, struct in self.structs.items():
            structures[iS] = {}
            structures[iS]['name'] = name
            structures[iS]['color'] = [255,0,0]
            structures[iS]['contours'] = []
            structures[iS]['number'] = iS

            contours = []
            for z in range(struct.shape[2]):
                sliceMask = struct[:, :, z].astype(np.uint8) * 255
                contours_z, _ = cv2.findContours(sliceMask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
                if contours_z:
                    contours.append((z, contours_z))
            for z, contours_z in contours:
                for contour in contours_z:
                    contour_data = []
                    for point in contour[:, 0, :]:
                        x, y = point
                        slice_position = slices[z].ImagePositionPatient
                        pixel_spacing = self.patient.pixelSpacing
                        contour_data.extend([x * pixel_spacing[1] + slice_position[0],
                                             y * pixel_spacing[0] + slice_position[1],
                                             float(slice_position[2])])
                    contour_data = np.array(contour_data).flatten().tolist()
                    structures[iS]['contours'].append({'image_uid': slices[z].SOPInstanceUID, 'points': contour_data})

            structure_set_roi_item = Dataset()
            structure_set_roi_item.ROINumber = 1
            structure_set_roi_item.ROIName = name
            structure_set_roi_item.ROIGenerationAlgorithm = ''
            structure_set_roi_item.ReferencedFrameOfReferenceUID = slices[0].FrameOfReferenceUID
            structure_set_roi_sequence.append(structure_set_roi_item)
            roi_contour_item = Dataset()
            roi_contour_item.ReferencedROINumber = 1
            roi_contour_item.ROIDisplayColor = structures[iS]['color']
            contour_sequence = Sequence()
            for contour_info in structures[iS]['contours']:
                contour = Dataset()
                contour.ContourGeometricType = 'CLOSED_PLANAR'
                contour.NumberOfContourPoints = int(len(contour_info['points']) / 3)
                contour.ContourData = contour_info['points']

                contour_image_sequence = Sequence()
                contour_image_item = Dataset()
                contour_image_item.ReferencedSOPClassUID = '1.2.840.10008.5.1.4.1.1.2'
                contour_image_item.ReferencedSOPInstanceUID = contour_info['image_uid']
                contour_image_sequence.append(contour_image_item)

                contour.ContourImageSequence = contour_image_sequence
                contour_sequence.append(contour)

            roi_contour_item.ContourSequence = contour_sequence
            roi_contour_sequence.append(roi_contour_item)
            iS += 1

        # Set the RTStruct sequences
        rtstruct.StructureSetROISequence = structure_set_roi_sequence
        rtstruct.ROIContourSequence = roi_contour_sequence

        # Save the RTStruct file

        rtstruct.save_as(filename)
