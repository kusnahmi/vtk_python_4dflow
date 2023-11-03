import vtk
import numpy as np
from vtk.util import numpy_support

class vtkUtil:

    # Convert numpy array to VTK Image Data
    @staticmethod
    def convert_numpyArray_vtkImageData(npImage, referenceVtkImageData=None):
        # Keep shape.
        npImageShape = npImage.shape
        # Rearrange numpy array.
        npImage = np.transpose(npImage, (2,1,0))
        # Convert numpy array to vtkArray
        # Since numpy_to_vtk function only supports 1-D and 2-D dimension array, numpy array should be flattened.
        vtk_data = numpy_support.numpy_to_vtk(num_array = npImage.ravel(), deep = True, array_type = vtk.VTK_SHORT)
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

        return imageData
