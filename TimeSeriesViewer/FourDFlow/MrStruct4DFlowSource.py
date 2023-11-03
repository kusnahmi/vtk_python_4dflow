## @package MrStruct4DFlowSource
#  Documentation for MrStruct4DFlowSource module.
#
#  More details.

import vtk
from vtk.util.vtkAlgorithm import VTKPythonAlgorithmBase
from vtk.numpy_interface import dataset_adapter as dsa
import numpy as np
from vtk.util import numpy_support
# import MrStruct
import FourDFlow
import VtkCustom.vtkUtil


## Documentation for a class.
#
#  More details.
class MrStruct4DFlowSource(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self,
                                        nInputPorts=0,
                                        nOutputPorts=1, outputType='vtkImageData')

        self.__Size = (0, 0, 0)

        self.__TimeValue = 0
        self.__TimeSteps = np.linspace(0, 10, 20)
        self.__TimeOffset = 0

        self.__TimeSpaceInverse = False

        self.__FilepathMagnitude = None
        self.__FilepathVelocity = None
        self.__FilepathMask = None
        self.__FlagMagnitudeChanged = False
        self.__FlagVectorChanged = False
        self.__FlagMaskChanged = False
        self.__ArrayScalar = None
        self.__ArrayVector = None
        self.__ArrayMask = None
        self.__Vox = (1, 1, 1)

    def RequestInformation(self, request, inInfo, outInfo):
        print("requestInformation")

        if ((self.__FlagMagnitudeChanged) or (self.__FlagVectorChanged) or (self.__FlagMaskChanged)):
            self.UpdateMrstructData()

        dims = self.__Size
        info = outInfo.GetInformationObject(0)
        info.Set(vtk.vtkStreamingDemandDrivenPipeline.WHOLE_EXTENT(),
                 (0, dims[0] - 1, 0, dims[1] - 1, 0, dims[2] - 1), 6)
        info.Set(vtk.vtkStreamingDemandDrivenPipeline.TIME_STEPS(),
                 self.__TimeSteps, len(self.__TimeSteps))
        info.Set(vtk.vtkStreamingDemandDrivenPipeline.TIME_RANGE(),
                 [self.__TimeSteps[0], self.__TimeSteps[-1]], 2)

        print(f"Time Steps: {self.__TimeSteps} / len: {len(self.__TimeSteps)}")
        print(f"Time Range: {[self.__TimeSteps[0], self.__TimeSteps[-1]]}")
        return 1

    def RequestData(self, request, inInfo, outInfo):
        print("requestData")

        if ((self.__FlagMagnitudeChanged) or (self.__FlagVectorChanged) or (self.__FlagMaskChanged)):
            self.UpdateMrstructData()

        info = outInfo.GetInformationObject(0)
        # The time step requested
        t = info.Get(vtk.vtkStreamingDemandDrivenPipeline.UPDATE_TIME_STEP())

        self.GetImageDataAtTime(t, outInfo)

        return 1

    def UpdateMrstructData(self):
        if (self.__FlagMagnitudeChanged):
            # Set self.__ArrayScalar
            self.SetMrStructScalarFromFile(self.__FilepathMagnitude)
            self.__FlagMagnitudeChanged = False

        if (self.__FlagVectorChanged):
            # Set self.__ArrayVector
            self.__ArrayVector = self.SetMrStructVectorFromFile(self.__FilepathVelocity)
            print(f"Original Magnitude array shape{self.__ArrayScalar.shape}")
            # [x, y, z, echo, t] ==> [x, y, z, t]
            self.__ArrayScalar = np.sqrt(np.sum(self.__ArrayVector ** 2, axis=3))
            print(f"RMS velocity array shape{self.__ArrayScalar.shape}")

            self.__FlagVectorChanged = False

        if (self.__FlagMaskChanged):
            self.__ArrayMask = self.ReadMrStructMask(self.__FilepathMask)
            print(
                f"mag {self.__ArrayScalar.shape} / vel {self.__ArrayVector.shape} / mask {self.__ArrayMask.shape}")
            arrayMaskTrepeated = np.repeat(self.__ArrayMask[:, :, :, np.newaxis], self.__ArrayScalar.shape[3],
                                           axis=3)
            # self.__ArrayMagnitude = self.__ArrayMagnitude * arrayMaskTrepeated
            arrayMaskEchoTrepeated = np.repeat(arrayMaskTrepeated[:, :, :, np.newaxis, :], 3, axis=3)
            # self.__ArrayVelocity = self.__ArrayVelocity * arrayMaskEchoTrepeated
            self.__FlagMaskChanged = False
            print(
                f"mag {self.__ArrayScalar.shape} / vel {self.__ArrayVector.shape} / mask {self.__ArrayMask.shape}")

    def SetMrStructScalarFromFile(self, filepath):
        mrStructScalar = FourDFlow.MrStruct(filepath)
        self.SetMrStructScalar(mrStructScalar)  # [x,y,z,t]

    def SetMrStructScalar(self, mrStructScalar):
        dataArray = mrStructScalar.dataAy
        timeStep = mrStructScalar.user['size_t_step']
        self.__Size = dataArray.shape[0:3]
        # self.__TimeSteps = np.linspace(0, (dataArray.shape[3]-1) * timeStep, dataArray.shape[3])
        self.__TimeSteps = np.linspace(0, (dataArray.shape[3]) * timeStep, dataArray.shape[3] + 1)
        self.__Vox = mrStructScalar.vox[0:3]

        print(f"vox: {self.__Vox}")

        self.__ArrayScalar = mrStructScalar.dataAy  # [x,y,z,t]
        self.__FlagMagnitudeChanged = False

    def SetMrStructVectorFromFile(self, filepath):
        mrStructVector = FourDFlow.MrStruct(filepath)
        self.SetMrStructVector(mrStructVector)  # [x,y,z,echo,t]

    def SetMrStructVector(self, mrStructVector):
        self.__ArrayVector = mrStructVector.dataAy  # [x,y,z,t]
        self.__FlagVectorChanged = False

    def ReadMrStructMask(self, filepath):
        mrStruct_mask = FourDFlow.MrStruct(filepath)
        return mrStruct_mask.dataAy  # [x,y,z]

    def SetFilepathMagnitude(self, filepath):
        if filepath != self.__FilepathMagnitude:
            self.Modified()
            self.__FilepathMagnitude = filepath
            self.__FlagMagnitudeChanged = True

    def SetFilepathVelocity(self, filepath):
        if filepath != self.__FilepathVelocity:
            self.Modified()
            self.__FilepathVelocity = filepath
            self.__FlagVectorChanged = True

    def SetFilepathMask(self, filepath):
        if filepath != self.__FilepathMask:
            self.Modified()
            self.__FilepathMask = filepath
            self.__FlagMaskChanged = True

    def SetSize(self, size):
        if size != self.__Size:
            self.Modified()
            self.__Size = size

    def SetTimeValue(self, time):
        if time != self.__TimeValue:
            self.Modified()
            self.__TimeValue = time

    def SetTimeOffset(self, offset):
        if offset != self.__TimeOffset:
            self.Modified()
            self.__TimeOffset = offset

    def GetTimeSpaceInverse(self):
        return self.__TimeSpaceInverse

    def SetTimeSpaceInverse(self, bInverse):
        self.__TimeSpaceInverse = bInverse

    ## Get the time step between each time frame
    def GetTimeStep(self):
        return self.__TimeSteps[1] - self.__TimeSteps[0]

    ## How many time frames are in source data?
    def GetNumberOfTimeFrames(self):
        # return len(self.__TimeSteps)
        return self.__ArrayScalar.shape[3]

    def GetMaskArray(self):
        return self.__ArrayMask

    def GetMaskImgData(self):
        maskImgData = vtkUtil.convert_numpyArray_vtkImageData(self.__ArrayMask)
        maskImgData.SetSpacing(self.__Vox)
        return maskImgData

    ## Retrieve single time frame data on specified time in vtkImageData
    #
    #  @param time Specified time (in ms)
    #  @param outInfo Output information including vtkImageData object, if this function is called from pipeline. This will be None if it is called standalone.
    def GetImageDataAtTime(self, time, outInfo=None):
        timeStep = self.GetTimeStep()
        if (not self.__TimeSpaceInverse):
            t_index = (np.round(time / timeStep).astype(int) + self.__TimeOffset) % (len(self.__TimeSteps) - 1)
        else:
            t_index = (np.round(-time / timeStep).astype(int) + self.__TimeOffset) % (len(self.__TimeSteps) - 1)

        print(f"UPDATE_TIME_STEP: {time} / index {t_index}")

        return self.GetImageDataAtTimeindex(t_index, outInfo)

    ## Retrieve single time frame data on specified time frame in vtkImageData
    #  @param time Specified frame (in frame number)
    #  @param outInfo Output information including vtkImageData object, if this function is called from pipeline. This will be None if it is called standalone.
    def GetImageDataAtTimeindex(self, timeIndex, outInfo=None):
        # Update if source is modified
        if ((self.__FlagMagnitudeChanged) or (self.__FlagVectorChanged) or (self.__FlagMaskChanged)):
            self.UpdateMrstructData()

        # Prepare output object
        if (outInfo == None):
            output = vtk.vtkImageData()
        else:
            output = vtk.vtkImageData.GetData(outInfo)

        # ## Magnitude
        npDataMagnitude = self.__ArrayScalar[:, :, :, timeIndex]

        # Set output dimension and spacing
        output.SetDimensions(npDataMagnitude.shape)
        # print(f"Data shape: {npDataMagnitude.shape}")
        output.SetSpacing(self.__Vox)

        # Make a VTK array from the numpy array (using pointers)
        vtkDataArrayScalar = dsa.numpyTovtkDataArray(npDataMagnitude.ravel(order='F'))

        # Set magnitude as scalar to output object
        vtkDataArrayScalar.SetName("Magnitude")
        output.GetPointData().SetScalars(vtkDataArrayScalar)

        # ## Velocity
        npDataVelocity = np.zeros((3, self.__Size[0], self.__Size[1], self.__Size[2]), order='F')

        # Inverse direction if backword tracing is selected.
        if (not self.__TimeSpaceInverse):
            npDataVelocity[0, :] = self.__ArrayVector[:, :, :, 0, timeIndex]
            npDataVelocity[1, :] = self.__ArrayVector[:, :, :, 1, timeIndex]
            npDataVelocity[2, :] = self.__ArrayVector[:, :, :, 2, timeIndex]
        else:
            npDataVelocity[0, :] = self.__ArrayVector[:, :, :, 0, timeIndex] * (-1)
            npDataVelocity[1, :] = self.__ArrayVector[:, :, :, 1, timeIndex] * (-1)
            npDataVelocity[2, :] = self.__ArrayVector[:, :, :, 2, timeIndex] * (-1)

        # Make a VTK array from the numpy array (using pointers)
        vtkDataArrayVector = dsa.numpyTovtkDataArray(npDataVelocity.ravel(order='A').reshape(
            self.__Size[0] * self.__Size[1] * self.__Size[2], 3))

        # Set velocity as vector to output object
        vtkDataArrayVector.SetName("VelocityVector")
        output.GetPointData().SetVectors(vtkDataArrayVector)

        print(f"Get ImageData at index {timeIndex}")

        return output
