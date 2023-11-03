import numpy as np
import numpy.ma as ma
import scipy.signal
from vtk.util import numpy_support
import vtk

def calcKineticEnergy(velMrStruct, npMask):
    # Voxel size
    vox = velMrStruct.vox[0:3]
    cellVolume = np.prod(vox) / 1000.0  # [cm3]

    # Make mask and count cells
    numOfCells = np.count_nonzero(npMask)

    # Volume
    volume = numOfCells * cellVolume  # [cm3]

    # KE = 1/2 * mass * velocity^2
    #    = 1/2 * density * volume * velocity^2      (Assumed blood is not compressible)
    #    = 1/2 * density * SIGMA ( unit_volume * velocity^2 )
    #    = 1/2 * density * unit_volume * SIGMA ( velocity^2 )
    # KE indexed = KE / volume
    #            = KE / (unit_volume * num_of_cells)
    #            = 1/2 * density * unit_volume * SIGMA ( velocity^2 ) / (unit_volume * num_of_cells)
    #            = 1/2 * density * SIGMA ( velocity^2 ) / num_of_cells
    # velocity^2 = v_x^2 + v_y^2 + v_z^2
    # density of blood = 1.06 g/cm3
    density = 1.06  # [g/cm3]
    speedSquaredArray = (np.sum(velMrStruct.dataAy ** 2, axis=3))  # [x,y,z,t] [(m/s)^2]
    KE = np.zeros([speedSquaredArray.shape[3], 2])
    for i in range(speedSquaredArray.shape[3]):
        speedSquaredMasked = ma.array(speedSquaredArray[:, :, :, i], mask=np.logical_not(npMask))
        kineticEnergyIndexed = 0.5 * density * speedSquaredMasked.sum() / numOfCells * 1000  # [g/cm3] * [(m/s)^2] * 1000 = [mJ/cm3] * 1000  = [uJ/mL]
        kineticEnergy = kineticEnergyIndexed * volume  # [uJ]
        KE[i, 0] = kineticEnergy
        KE[i, 1] = kineticEnergyIndexed

    return [cellVolume, numOfCells, volume, KE]


def getMaximaMinima(data):
    max, min = _getMaximaMinima(data)
    if ((len(max) > 2) or (len(min) > 2)):
        filtered = scipy.signal.savgol_filter(data, 5, 3, mode='wrap')  # window size 51, polynomial order 3
        max, min = _getMaximaMinima(filtered)

    return max, min


def _getMaximaMinima(data):
    data2 = np.concatenate([data, data])
    min_a = (np.diff(np.sign(np.diff(data2))) > 0).nonzero()[0] + 1
    min_b = min_a[min_a <= len(data)]
    min = np.where(min_b == len(data), 0, min_b)
    min_sorted = min[np.argsort(data[min])]

    max_a = (np.diff(np.sign(np.diff(data2))) < 0).nonzero()[0] + 1
    max_b = max_a[max_a <= len(data)]
    max = np.where(max_b == len(data), 0, max_b)
    max_sorted = max[np.argsort(-data[max])]

    return max_sorted, min_sorted


def calcKineticEnergyStatistics(KEarray, t_pSyst=-1, t_pDias=-1, t_ED=-1, t_ES=-1):
    # parameters
    #     time avearge KE
    #     peak systole
    #     peak diastole
    #     systole: from maximum peak, trace back and forth
    #     diastole: from second peak, trace back and forth
    if ((t_pSyst == -1) or (t_pDias == -1) or (t_ED == -1) or (t_ES == -1)):
        max, min = getMaximaMinima(KEarray)
        if ((len(max) == 2) or (len(min) == 2)):
            t_pSyst = max[0]
            t_pDias = max[1]
            if ((np.max(min) <= max[0]) or (np.min(min) >= max[0])):
                t_ES = np.min(min)
                t_ED = np.max(min)
            else:
                t_ED = np.min(min)
                t_ES = np.max(min)
        else:
            t_pSyst = -1
            t_pDias = -1
            t_ED = -1
            t_ES = -1

    if not ((t_pSyst == -1) or (t_pDias == -1) or (t_ED == -1) or (t_ES == -1)):
        KEavg = np.mean(KEarray)
        KEpeakSystole = KEarray[t_pSyst]
        KEpeakDiastole = KEarray[t_pDias]
        if (t_ED < t_ES):
            KEavgSystole = np.mean(KEarray[t_ED:t_ES])
            KEavgDiastole = np.mean(np.concatenate((KEarray[t_ES:], KEarray[:t_ED])))
        else:
            KEavgDiastole = np.mean(KEarray[t_ES:t_ED])
            KEavgSystole = np.mean(np.concatenate((KEarray[t_ED:], KEarray[:t_ES])))
    else:
        KEavg = 0
        KEpeakSystole = 0
        KEpeakDiastole = 0
        KEavgSystole = 0
        KEavgDiastole = 0

    return t_pSyst, t_pDias, t_ED, t_ES, KEavg, KEpeakSystole, KEpeakDiastole, KEavgSystole, KEavgDiastole


def convertArrayToString(array, separator=' '):
    outString = ''
    if len(array) == 1:
        outString = outString + str(array[0])
    if len(array) > 1:
        outString = outString + str(array[0])
        for item in array[1:]:
            outString = outString + separator + str(item)
    return outString


def buildStreamActor(npVelocityX, npVelocityY, npVelocityZ, npScalar, spacing, planeSource):
    #Utility.buildStreamActor(npVelocityX, npVelocityY, npVelocityZ, npScalar, spacing)
    npShape = npVelocityX.shape
    np_vector = np.transpose(np.vstack((npVelocityX.ravel(order='F'), npVelocityY.ravel(order='F'), npVelocityZ.ravel(order='F'))))
    vtk_vector = numpy_support.numpy_to_vtk(np_vector, deep=True)
    vtk_vector.SetName("vectors")

    vtk_scalar = numpy_support.numpy_to_vtk(npScalar.ravel(order='F'), deep=True)
    vtk_scalar.SetName("scalars")

    ### Make vtkImageData
    imageData = vtk.vtkImageData()
    imageData.SetSpacing(spacing[0], spacing[1], spacing[2])
    imageData.SetExtent(0, npShape[0]-1, 0, npShape[1]-1, 0, npShape[2]-1)

    imageData.GetPointData().SetScalars(vtk_scalar)
    imageData.GetPointData().SetVectors(vtk_vector)

    # Source Plane Definition
    # planeSource = vtk.vtkPlaneSource()
    # planeSource.SetOrigin(0.05, 0.05, 0.05)
    # planeSource.SetPoint1(0.05, 0.85, 0.05)
    # planeSource.SetPoint2(0.05, 0.05, 0.85)
    # planeSource.SetXResolution(5)
    # planeSource.SetYResolution(1)

    # Render Streamline
    rk = vtk.vtkRungeKutta45()
    # Create source for streamtubes
    streamer = vtk.vtkStreamTracer()
    # streamer.SetInputData(gridDataGen)
    streamer.SetInputData(imageData)
    streamer.SetSourceConnection(planeSource.GetOutputPort())
    streamer.SetMaximumPropagation(100)
    streamer.SetIntegrationStepUnit(2)
    streamer.SetMinimumIntegrationStep(0.1)
    streamer.SetMaximumIntegrationStep(1.0)
    streamer.SetInitialIntegrationStep(0.2)
    streamer.SetIntegrationDirection(0)
    streamer.SetIntegrator(rk)
    streamer.SetRotationScale(0.5)
    streamer.SetMaximumError(1.0e-8)
    mapStream = vtk.vtkPolyDataMapper()
    mapStream.SetInputConnection(streamer.GetOutputPort())
    mapStream.SetScalarRange(imageData.GetScalarRange())
    streamActor = vtk.vtkActor()
    streamActor.SetMapper(mapStream)

    return streamActor, imageData


def buildPathlineActor(npVelocityX, npVelocityY, npVelocityZ, npScalar, spacing, planeSource):
    #Utility.buildStreamActor(npVelocityX, npVelocityY, npVelocityZ, npScalar, spacing)
    npShape = npVelocityX.shape
    np_vector = np.transpose(np.vstack((npVelocityX.ravel(order='F'), npVelocityY.ravel(order='F'), npVelocityZ.ravel(order='F'))))
    vtk_vector = numpy_support.numpy_to_vtk(np_vector, deep=True)
    vtk_vector.SetName("vectors")

    np_vectorMag = np.sqrt((npVelocityX**2 + npVelocityY**2 + npVelocityZ**2)/3)
    np_vectorMag = np.clip( np_vectorMag, np.percentile(np_vectorMag, 30), np.percentile(np_vectorMag,70))
    # vtk_scalar = numpy_support.numpy_to_vtk(npScalar.ravel(order='F'), deep=True)
    vtk_scalar = numpy_support.numpy_to_vtk(np_vectorMag.ravel(order='F'), deep=True)
    vtk_scalar.SetName("scalars")

    ### Make vtkImageData
    imageData = vtk.vtkImageData()
    imageData.SetSpacing(spacing[0], spacing[1], spacing[2])
    imageData.SetExtent(0, npShape[0]-1, 0, npShape[1]-1, 0, npShape[2]-1)

    imageData.GetPointData().SetScalars(vtk_scalar)
    imageData.GetPointData().SetVectors(vtk_vector)

    hedgehog = vtk.vtkHedgeHog()
    hedgehog.SetInputData(imageData)
    hedgehog.SetScaleFactor(1)


    arrow = vtk.vtkArrowSource()
    arrow.SetTipResolution(6)
    arrow.SetTipRadius(0.1)
    arrow.SetTipLength(0.35)
    arrow.SetShaftResolution(6)
    arrow.SetShaftRadius(0.03)


    glyph = vtk.vtkGlyph3D()
    glyph.SetInputData(imageData)
    glyph.SetSourceConnection(arrow.GetOutputPort())
    glyph.SetVectorModeToUseVector()
    glyph.SetColorModeToColorByScalar()
    glyph.SetScaleModeToDataScalingOff()
    glyph.OrientOn()
    glyph.SetScaleFactor(1)

    # glyph.SetScaleModeToScaleByScalar()
    # glyph.ScalingOn()


    mapStream = vtk.vtkPolyDataMapper()
    # mapStream.SetInputConnection(hedgehog.GetOutputPort())
    mapStream.SetInputConnection(glyph.GetOutputPort())
    mapStream.SetScalarRange(imageData.GetScalarRange())
    streamActor = vtk.vtkActor()
    streamActor.SetMapper(mapStream)

    return streamActor, imageData
