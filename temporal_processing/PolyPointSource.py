## @package ParticleTracingTest
#  Documentation for PolyPointSource module.
#
#  vtkPolyPointSource is not declared for python.
#  Custom implementation according to https://vtk.org/doc/nightly/html/classvtkPolyPointSource.html
#  Implementation reference: vtkPointSource

import vtk
from vtk.util.vtkAlgorithm import VTKPythonAlgorithmBase


## class PolyPointSource
#  vtkPolyPointSource is not declared for python.
#  Create points from a list of input points
#  This can be used where vtkPlaneSource is being used.
#  Usage:
#     sourcePoints = vtk.vtkPoints()
#     sourcePoints.SetDataType(vtk.VTK_DOUBLE)
#     sourcePoints.Allocate(n)
#     sourcePoints.InsertNextPoint(x,y,z) # repeat n times
#     pointSourceMV = PolyPointSource.PolyPointSource()
#     pointSourceMV.SetPoints(sourcePoints)
#     streamer.SetSourceConnection(pointSourceMV.GetOutputPort())
class PolyPointSource(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self,
                                        nInputPorts=0,
                                        nOutputPorts=1, outputType='vtkPolyData')
        self.__ArrayPoints = None

    # Set a points array (vtkPoints) as source points.
    def SetPoints(self, points):
        self.Modified()
        self.__ArrayPoints = points

    # Return vtkCellArray as vtkPlaneSource does.
    def RequestData(self, request, inInfo, outInfo):
        print("requestData")

        opt = vtk.vtkPolyData.GetData(outInfo)

        numPoints = self.__ArrayPoints.GetNumberOfPoints()

        # newPoints = vtk.vtkPoints()
        # newPoints.SetDataType(vtk.VTK_DOUBLE)
        # newPoints.Allocate(self.__NumberOfPoints)
        #
        newVerts = vtk.vtkCellArray()
        newVerts.Allocate(newVerts.EstimateSize(1,numPoints))
        newVerts.InsertNextCell(numPoints)

        for i in range(numPoints):
            newVerts.InsertCellPoint(i)

        opt.SetPoints(self.__ArrayPoints)
        opt.SetVerts(newVerts)

        return 1
