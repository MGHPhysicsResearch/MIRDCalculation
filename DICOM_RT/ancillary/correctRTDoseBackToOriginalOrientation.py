import os
from pydicom import dcmread

input_file = '/Users/ai925/workspace/Bailey/Images/DICOM_01/0000/DoseOnCTGrid.dcm'
output_folder = '/Users/ai925/workspace/Bailey/Images/DICOM_01/0000/dose_corrected/'

# Ensure the output folder exists
os.makedirs(output_folder, exist_ok=True)

# Read the RTDOSE DICOM file
ds = dcmread(input_file)

# Check if the ImageOrientationPatient is [1, 0, 0, 0, 1, 0]
if ds.ImageOrientationPatient == [1, 0, 0, 0, 1, 0]:
    # Flip the pixel data vertically (along the Y-axis)
    ds.PixelData = ds.pixel_array[::-1, :, :].tobytes()

    # Reorder the dose grid along the Z-axis to reverse the order of Z slices
    ds.PixelData = ds.pixel_array[:, :, ::-1].tobytes()

    # Update the ImagePositionPatient (origin) for the Y position
    ds.ImagePositionPatient[1] = ds.ImagePositionPatient[1] + (ds.Rows - 1) * ds.PixelSpacing[1]

    # Update the Z-coordinate of ImagePositionPatient
    ds.ImagePositionPatient[2] = ds.ImagePositionPatient[2] + (ds.NumberOfFrames - 1) * ds.SliceThickness

    # Update the ImageOrientationPatient to [1, 0, 0, 0, -1, 0]
    ds.ImageOrientationPatient = [1, 0, 0, 0, -1, 0]

    # Patient position
    ds.PatientPosition = 'FFS'
    # Patient orientation
    ds.PatientOrientation = ['L', 'A']

# Save the modified DICOM file to the output folder
output_path = os.path.join(output_folder, ds.SOPInstanceUID + ".dcm")
ds.save_as(output_path)

print('Processing completed!')
