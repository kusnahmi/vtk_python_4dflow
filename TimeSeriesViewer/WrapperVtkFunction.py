# MDSC 689.03
# Advanced Medical Image Processing
#
# Week 4 Assignment - Visualization
# Hansuk_Kim
# February 22, 2020
# --------------------------------------------
#
# Wrapper module for commonly used VTK function.
#
# Tested on Windows 10.
# --------------------------------------------
#
import vtk
import numpy as np
from vtk.util import numpy_support

### Convert vtkImageData and numpy Array back and forth
# Convert vtkImageData to numpy Array
def convert_vtkImageData_numpyArray(imageData):
    # Get the point data object vtkPointData from the image data object
    pointData = imageData.GetPointData()
    # Ensure that only one array exists within the 'vtkPointData' object
    assert (pointData.GetNumberOfArrays()==1)
    # Get the vtkArray object from vtkPointData
    arrayData = pointData.GetArray(0)

    # Load dimensions of the image data.
    _extent = imageData.GetExtent()
    # Calculate number of pixels on each axis. This will be used to reshape numpy array.
    ConstPixelDims = [_extent[1]-_extent[0]+1, _extent[3]-_extent[2]+1, _extent[5]-_extent[4]+1]
    print(_extent)

    ### Convert vtkArray to NumPy array
    # Convert the `vtkArray` to a NumPy array
    npImage = numpy_support.vtk_to_numpy(arrayData)
    # Reshape the NumPy array to 3D using 'ConstPixelDims' as a 'shape'
    npImage = npImage.reshape(ConstPixelDims, order='F')

    return npImage

# Convert numpy array to VTK Image Data
def convert_numpyArray_vtkImageData(npImage, referenceVtkImageData=None):
    # Keep shape.
    npImageShape = npImage.shape
    # Rearrange numpy array.
    npImage = np.transpose(npImage, (2,1,0))
    # Convert numpy array to vtkArray
    # Since numpy_to_vtk function only supports 1-D and 2-D dimension array, numpy array should be flattened.
    # vtk_data = numpy_support.numpy_to_vtk(num_array = npImage.ravel(), deep = True, array_type = vtk.VTK_SHORT)
    vtk_data = numpy_support.numpy_to_vtk(num_array = npImage.ravel(), deep = True)
    # Prepare target vtkImageData
    imageData = vtk.vtkImageData()
    # Set dimension of target vtkImageData from the shape of numpy array.
    imageData.SetDimensions(npImageShape)
    # imageData.SetSpacing([1,1,1])	# original image data values should be recovered.
    # imageData.SetOrigin([0,0,0])
    if (referenceVtkImageData == None):     # If reference image data is not given, default spcaing and origin is applied.
        imageData.SetSpacing([1,1,1])
        imageData.SetOrigin([0,0,0])
    else:                                   # Recover spacing and origin from the reference image data.
        imageData.SetSpacing( referenceVtkImageData.GetSpacing() )
        imageData.SetOrigin( referenceVtkImageData.GetOrigin() )

    # Set converted vtkArray object to the source of vtkImageData
    imageData.GetPointData().SetScalars(vtk_data)

    # print(f"GetScalarTypeMin {imageData.GetScalarTypeMin()}")
    # print(f"GetScalarTypeMax {imageData.GetScalarTypeMax()}")
    # print(f"GetScalarSize    {imageData.GetScalarSize()}")
    # print(f"getscalartypemin {imageData.GetScalarTypeMin()}")
    # print(imageData)
    return imageData

### Thresholding and dilation
def FilterThresholding (imageData, thresholdRange):
    # Make Binary threshold image
    threshold = vtk.vtkImageThreshold()
    threshold.SetInputData(imageData)
    threshold.ThresholdBetween(thresholdRange[0], thresholdRange[1])
    threshold.ReplaceInOn()
    threshold.SetInValue(1.0)
    threshold.SetOutValue(0.0)
    threshold.Update()
    return threshold.GetOutput()

def FilterDilateErode (imageData, kernelSize, dilateValue, erodeValue):
    dilateErode = vtk.vtkImageDilateErode3D()
    dilateErode.SetInputData(imageData)
    dilateErode.SetDilateValue(dilateValue)
    dilateErode.SetErodeValue(erodeValue)
    dilateErode.SetKernelSize(kernelSize, kernelSize, kernelSize)
    dilateErode.ReleaseDataFlagOff()
    dilateErode.Update()
    return dilateErode.GetOutput()

def FilterGaussian (imageData, gaussianRadius, gaussianSD):
    gaussian = vtk.vtkImageGaussianSmooth()
    gaussian.SetStandardDeviations(gaussianSD, gaussianSD, gaussianSD)
    gaussian.SetRadiusFactors(gaussianRadius, gaussianRadius, gaussianRadius)
    gaussian.SetInputData(imageData)
    gaussian.Update()
    return gaussian.GetOutput()

def FilterMedian (imageData, kernelSize):
    medianFilter = vtk.vtkImageMedian3D()
    medianFilter.SetInputData(imageData)
    medianFilter.SetKernelSize(kernelSize,kernelSize,kernelSize)
    medianFilter.Update()
    return medianFilter.GetOutput()
