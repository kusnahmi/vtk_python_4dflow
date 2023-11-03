## @package ParticleTracingTest
#  Particle Tracing Test
#
#  Read 4D Flow magnitude, velocity, and mask
#  Visualize vector field
#  Particle tracing

import vtk
import mrstruct
import vtkVisualization as vtkViz
import numpy as np
import PolyPointSource
import matplotlib.pyplot as plt
from vtk.util import numpy_support


## main function
#
#  Main flow of particle tracing and visualization
def main():
    ###########################################################################
    # Configuration
    ###########################################################################

    ## Config: 4D flow source data
    #filepathMagnitude = "C:/d/project/Case_test/001A00505620161214/3dpc/mrstruct_subject_20161214_user011/mag_struct.mat"
    #filepathVelocity = "C:/d/project/Case_test/001A00505620161214/3dpc/mrstruct_subject_20161214_user011/vel_struct.mat"
    #filepathMask = "C:/d/project/Case_test/001A00505620161214/3dpc/mrstruct_subject_20161214_user011/LA_LV_volume.mat"

    filepathMagnitude = "D:/data_clinical/control_all/done/001B00849420171207/3dpc/mrstruct_subject_20171207_user011/mag_struct.mat"
    filepathVelocity = "D:/data_clinical/control_all/done/001B00849420171207/3dpc/mrstruct_subject_20171207_user011/vel_struct.mat"
    filepathMask = "D:/data_clinical/control_all/done/001B00849420171207/3dpc/mrstruct_subject_20171207_user011/LA_volume_2.mat"

    t = 0

    ## Config: entrance & exit plane coordinate
    pointsPlaneMV = [
        [155.75200746815736, 133.4443205456904, 57.55322438879402],
        [195.11049180338034, 162.2567972045777, 52.920065217945236],
        [158.37555320991834, 137.10930329145776, 102.63175635712525]
    ]
    pointsNew = [
        np.array(pointsPlaneMV[0]),
        np.array(pointsPlaneMV[1]),
        np.array(pointsPlaneMV[2])
    ]
    pointsNew = [
        (pointsNew[0] * 2 + pointsNew[1] * 3) / 5,
        pointsNew[1] + (pointsNew[1] - pointsNew[0]) / 6 + (pointsNew[2] - pointsNew[0]) / 5 * 2,
        (pointsNew[2] * 5 + pointsNew[0]) / 6
    ]

    # pointsNew = [
    #     np.array(pointsPlaneMV[0]) + (np.array(pointsPlaneMV[1]) - np.array(pointsPlaneMV[0]))/2 + (np.array(pointsPlaneMV[2]) - np.array(pointsPlaneMV[0]))/2,
    #     np.array(pointsPlaneMV[1]) + (np.array(pointsPlaneMV[2]) - np.array(pointsPlaneMV[0]))/2,
    #     np.array(pointsPlaneMV[2]) + (np.array(pointsPlaneMV[1]) - np.array(pointsPlaneMV[0]))/2,
    # ]
    pointsPlaneMV = pointsNew
    print(f"pointsPlaneMV {pointsPlaneMV}")
    print(f"pointsNew {pointsNew}")

    pointsPlaneAV = [
        [187.59152548881804, 155.1431218693424, 60.117124070509334],
        [181.91575128426845, 122.62066865580515, 47.517833373914726],
        [167.03837674355384, 148.56220588247717, 86.36324948493984]
    ]
    pointsNew = [
        np.array(pointsPlaneAV[0]),
        np.array(pointsPlaneAV[0]) + (np.array(pointsPlaneAV[1]) - np.array(pointsPlaneAV[0])) * 2 / 3,
        np.array(pointsPlaneAV[2])
    ]
    print(f"pointsPlaneAV {pointsPlaneAV}")
    pointsPlaneAV = pointsNew
    print(f"pointsPlaneAV {pointsPlaneAV}")

    ###########################################################################
    # Prepare 4D flow source
    ###########################################################################
    s = mrstruct.MrStruct4DFlowSource()
    s.SetSize([50, 30, 40])
    s.SetCenter([10, 15, 20])
    s.SetRadius(10)
    s.SetTimeValue(t)
    s.SetTimeOffset(9)
    s.SetFilepathMagnitude(filepathMagnitude)
    s.SetFilepathVelocity(filepathVelocity)
    s.SetFilepathMask(filepathMask)
    s.SetTimeSpaceInverse(False)
    # s.SetTimeSpaceInverse(True)
    s.Update()









    ## test: 4D flow source shape (mask)
    print(s.GetMaskArray().shape)
    ## test: 4D flow source visualization (mask)
    maskPolyData, maskActor = vtkViz.GetMaskActor(s.GetMaskImgData())
    maskActor.GetProperty().SetOpacity(0.2)
    # ## test: 4D flow source vector visualization
    # imgData = s.GetImageDataAtTimeindex(9)
    # vtkViz.ShowVectorData(imgData, None)
    ## test: 4D flow source statistics
    stat = vtk.vtkImageHistogramStatistics()
    stat.SetInputConnection(s.GetOutputPort())
    stat.Update()
    print(f"time: {t} / mean: {stat.GetMean()}")

    ###########################################################################
    # Vector Screenshot for each time frame
    ###########################################################################
    # for t in range (s.GetNumberOfTimeFrames()):
    #     imgData = s.GetImageDataAtTimeindex(t)
    #     showVectorData(imgData, None)

    #     # stat.Update()
    #
    #     renWin.Render()
    #
    #     # screenshot code:
    #     w2if = vtk.vtkWindowToImageFilter()
    #     w2if.SetInput(renWin)
    #     w2if.Update()
    #     writer = vtk.vtkPNGWriter()
    #     writer.SetFileName("animate%d.png" % t)
    #     writer.SetInputData(w2if.GetOutput())
    #     writer.Write()
    #
    #     time.sleep(0.1)

    # s.SetTimeValue(0)
    # s.Update()

    ###########################################################################
    # Plane Definition
    ###########################################################################

    # planeWidget = vtkViz.GetPlaneWidget(s.GetImageDataAtTimeindex(9))
    planeWidgets = [0] * 5
    # for planeWidget in planeWidgets:
    for i in range(5):
        planeWidgets[i] = vtkViz.GetPlaneWidget(s.GetImageDataAtTimeindex(9))

    pointsPlanePulmonaryVein = [
        [[129.44711312832447, 125.039514028837, 62.170267920047195],
         [90.52017962118393, 173.85819366510196, 86.14737895975493],
         [178.19938152531574, 135.5127089670905, 119.99578157167443]],
        [[98.9940059277023, 163.80695868320024, 42.98303579014223],
         [151.46325364619972, 195.27104575744053, 16.71220568994338],
         [101.01525577726633, 204.51517944979338, 95.7754468637734]],
        [[201.32077056755713, 156.13959108883464, 114.29060734179616],
         [141.53143061826114, 123.69914068873558, 108.96800449679179],
         [168.13495096182305, 218.09125486551213, 109.48438249417087]],
        [[105.09152786201763, 220.66160385591442, 91.56071815870541],
         [134.94244264546884, 208.37770722852173, 13.707840033682963],
         [168.53433672036397, 234.91298380632654, 113.6377876793031]],

        pointsPlaneMV
    ]
    pointsPlanePulmonaryVein = np.array(pointsPlanePulmonaryVein)
    for i in range(len(pointsPlanePulmonaryVein)):
        pointsPlane = pointsPlanePulmonaryVein[i]
        center = (pointsPlane[1] + pointsPlane[2]) / 2
        vectorOP1 = pointsPlane[1] - pointsPlane[0]
        vectorOP2 = pointsPlane[2] - pointsPlane[0]
        vectorOP1Normalized = vectorOP1 / np.linalg.norm(vectorOP1)
        vectorOP2Normalized = vectorOP2 / np.linalg.norm(vectorOP2)
        sizeX = 5.0
        sizeY = 5.0
        # sizeX = 20.0
        # sizeY = 20.0
        # sizeX = 100.0
        # sizeY = 100.0
        newOrigin = center - vectorOP1Normalized * sizeX / 2 - vectorOP2Normalized * sizeY / 2
        newPoint1 = center + vectorOP1Normalized * sizeX / 2 - vectorOP2Normalized * sizeY / 2
        newPoint2 = center - vectorOP1Normalized * sizeX / 2 + vectorOP2Normalized * sizeY / 2
        pointsPlanePulmonaryVein[i] = [newOrigin, newPoint1, newPoint2]

    colorWidgetHandle = [
        [1, 0, 0], [0, 1, 0], [0.5, 0.5, 1], [0.5, 0.5, 0.5], [1, 1, 0]
    ]
    for i in range(len(pointsPlanePulmonaryVein)):
        planeWidgets[i].SetOrigin(pointsPlanePulmonaryVein[i][0])
        planeWidgets[i].SetPoint1(pointsPlanePulmonaryVein[i][1])
        planeWidgets[i].SetPoint2(pointsPlanePulmonaryVein[i][2])
        planeWidgets[i].GetHandleProperty().SetColor(colorWidgetHandle[i])
        print(planeWidgets[i].GetOrigin(), planeWidgets[i].GetPoint1(), planeWidgets[i].GetPoint2())

    # Visualize Plane Definition
    # vtkViz.ShowActors([maskActor], widgets=planeWidgets, coordinateMarkerOffset=pointsPlaneMV[0])
    vtkViz.ShowActors([maskActor])
    vtkViz.ShowActors([maskActor], showCoordinateMaker=False)
    imgData = s.GetImageDataAtTimeindex(4)
    vtkViz.ShowVectorData(imgData, maskActor)


    ###########################################################################
    # Define Source Plane
    ###########################################################################
    #  Pulmonary Vein
    planeSourcePVs = [0] * 4
    for i in range(4):
        planeSourcePVs[i] = vtk.vtkPlaneSource()
        planeSourcePVs[i].SetOrigin(pointsPlanePulmonaryVein[i][0])
        planeSourcePVs[i].SetPoint1(pointsPlanePulmonaryVein[i][1])
        planeSourcePVs[i].SetPoint2(pointsPlanePulmonaryVein[i][2])
        planeSourcePVs[i].SetXResolution(20)
        planeSourcePVs[i].SetYResolution(20)
        planeSourcePVs[i].Update()

    #  Mitral Valve
    planeSourceMV = vtk.vtkPlaneSource()
    planeSourceMV.SetOrigin(pointsPlaneMV[0])
    planeSourceMV.SetPoint1(pointsPlaneMV[1])
    planeSourceMV.SetPoint2(pointsPlaneMV[2])
    planeSourceMV.SetXResolution(20)
    planeSourceMV.SetYResolution(20)
    planeSourceMV.Update()

    #  Aortic Valve
    planeSourceAV = vtk.vtkPlaneSource()
    planeSourceAV.SetOrigin(pointsPlaneAV[0])
    planeSourceAV.SetPoint1(pointsPlaneAV[1])
    planeSourceAV.SetPoint2(pointsPlaneAV[2])
    planeSourceAV.SetXResolution(20)
    planeSourceAV.SetYResolution(20)
    planeSourceAV.Update()

    # # Define a set of points as a source
    # # 201125: Source by Points Definition
    # sourcePoints = vtk.vtkPoints()
    # sourcePoints.SetDataType(vtk.VTK_DOUBLE)
    # sourcePoints.Allocate(4)
    # sourcePoints.InsertNextPoint(179.7460175, 139.2701403, 71.03501532)
    # sourcePoints.InsertNextPoint(174.7460175, 134.2701403, 71.03501532)
    # sourcePoints.InsertNextPoint(179.7460175, 144.2701403, 71.03501532)
    # sourcePoints.InsertNextPoint(169.7460175, 149.2701403, 71.03501532)
    # pointSourceMV = PolyPointSource.PolyPointSource()
    # pointSourceMV.SetPoints(sourcePoints)

    ###########################################################################
    # Visual check up of plane definition
    # Image slicing with arbitrary plane: vtkImageReslice
    ###########################################################################
    # vtkImageReslice requires transform matrix
    #     |  vec0_x   vec1_x   vec2_x   center_x  |
    #     |  vec0_y   vec1_y   vec2_y   center_y  |
    #     |  vec0_z   vec1_z   vec2_z   center_z  |
    #     |       0        0        0          1  |

    # planeSourceSets = [[planeSourceMV, 12], [planeSourceAV, 5]]
    planeSourceSets = [[planeSourcePVs[0], 12], [planeSourcePVs[1], 12], [planeSourcePVs[2], 12],
                       [planeSourcePVs[3], 12], [planeSourceMV, 12], [planeSourceAV, 5]]
    # for planeLoop in range(2):
    for planeSourceSet in planeSourceSets:
        planeSource = planeSourceSet[0]
        planeSourceFrame = planeSourceSet[1]

        # Get vtkMatrix4x4 for reslice
        # planeSource = planeSourceAV
        # planeSource = planeSources[planeLoop]
        # obliqueCustom = GetResliceTransformMatrix(planeSourceMV)
        obliqueCustom = GetResliceTransformMatrix(planeSource)

        # Extract a slice in the desired orientation
        # imgData = s.GetImageDataAtTimeindex(5)  # Retrieve image data at a specific time phase from 4D flow source
        # imgData = s.GetImageDataAtTimeindex(12)  # Retrieve image data at a specific time phase from 4D flow source
        imgData = s.GetImageDataAtTimeindex(
            planeSourceFrame)  # Retrieve image data at a specific time phase from 4D flow source
        reslice = vtk.vtkImageReslice()
        reslice.SetInputData(imgData)
        reslice.SetOutputDimensionality(2)
        reslice.SetResliceAxes(obliqueCustom)
        reslice.SetInterpolationModeToLinear()
        reslice.Update()
        resliceOutput = reslice.GetOutput()

        # # vtkImageReslice slices SCALAR image, but it does not return VECTOR point data.
        # vtkDataArrayVector = resliceOutput.GetPointData().GetVectors()
        # print(f"resliceOutput.GetPointData().GetVectors(): {vtkDataArrayVector}")

        #######################################################
        transform = vtk.vtkTransform()
        transform.SetMatrix(obliqueCustom)
        transform.Inverse()
        ptOriginTrans = transform.TransformPoint(planeSourceMV.GetOrigin())
        ptPoint1Trans = transform.TransformPoint(planeSourceMV.GetPoint1())
        ptPoint2Trans = transform.TransformPoint(planeSourceMV.GetPoint2())

        # Create a vtkPoints object and store the points in it
        points = vtk.vtkPoints()
        points.InsertNextPoint(np.array(transform.TransformPoint(planeSource.GetOrigin())) + [0, 0, 1])
        points.InsertNextPoint(np.array(transform.TransformPoint(planeSource.GetPoint1())) + [0, 0, 1])
        points.InsertNextPoint(np.array(
            transform.TransformPoint(np.array(planeSource.GetPoint1()) + np.array(planeSource.GetPoint2()) - np.array(
                planeSource.GetOrigin()))) + [0, 0, 1])
        points.InsertNextPoint(np.array(transform.TransformPoint(planeSource.GetPoint2())) + [0, 0, 1])

        for ptId in range(4):
            print(f"pt {ptId}: {points.GetPoint(ptId)}")

        polyLine = vtk.vtkPolyLine()
        polyLine.GetPointIds().SetNumberOfIds(5)
        polyLine.GetPointIds().SetId(0, 0)
        polyLine.GetPointIds().SetId(1, 1)
        polyLine.GetPointIds().SetId(2, 2)
        polyLine.GetPointIds().SetId(3, 3)
        polyLine.GetPointIds().SetId(4, 0)

        # Create a cell array to store the lines in and add the lines to it
        cells = vtk.vtkCellArray()
        cells.InsertNextCell(polyLine)

        # Create a polydata to store everything in
        polyData = vtk.vtkPolyData()
        polyData.SetPoints(points)
        polyData.SetLines(cells)

        # Setup actor and mapper
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(polyData)

        actorRectangle = vtk.vtkActor()
        actorRectangle.SetMapper(mapper)
        colors = vtk.vtkNamedColors()
        actorRectangle.GetProperty().SetColor(colors.GetColor3d("Tomato"))

        # Text
        textActor = [0] * 4
        textContents = ["o", "1", "+", "2"]
        for i in range(points.GetNumberOfPoints()):
            # textActor[i] = vtk.vtkTextActor()
            # textActor[i].SetInput(textContents[i])
            # pt = points.GetPoint(i)
            # print(f"{pt} // {pt[0:2]}")
            # textActor[i].SetPosition2(pt[0:2])
            # textActor[i].GetTextProperty().SetFontSize(10)
            # textActor[i].GetTextProperty().SetColor ( 1.0, 0.0, 0.0 )

            textActor[i] = vtk.vtkBillboardTextActor3D()
            textActor[i].SetInput(textContents[i])
            textActor[i].SetPosition(points.GetPoint(i))
            textActor[i].GetTextProperty().SetFontSize(10)
            textActor[i].GetTextProperty().SetColor(1.0, 1.0, .4)
            textActor[i].GetTextProperty().SetJustificationToCentered()

        #######################################################

        # Visualize Sliced Image
        actorReslice = vtkViz.GetImageActor(resliceOutput)
        renderer = vtk.vtkRenderer()
        renderer.AddActor(actorReslice)
        renderer.AddActor(actorRectangle)
        for i in range(points.GetNumberOfPoints()):
            renderer.AddActor2D(textActor[i])
        window = vtk.vtkRenderWindow()
        window.SetSize(1440, 1024)
        window.AddRenderer(renderer)
        interactorStyle = vtk.vtkInteractorStyleImage()
        interactor = vtk.vtkRenderWindowInteractor()
        interactor.SetInteractorStyle(interactorStyle)
        window.SetInteractor(interactor)
        window.Render()
        interactor.Start()

    # # test: compare single point transform from the same matrix
    # transform = vtk.vtkTransform()
    # transform.SetMatrix(obliqueCustom)
    # transform.Inverse()
    # ptOriginTrans = transform.TransformPoint(planeSourceMV.GetOrigin())
    # ptPoint1Trans = transform.TransformPoint(planeSourceMV.GetPoint1())
    # ptPoint2Trans = transform.TransformPoint(planeSourceMV.GetPoint2())
    # print(f"ptOriginTrans {ptOriginTrans}, ptPoint1Trans {ptPoint1Trans}, ptPoint2Trans {ptPoint2Trans}")

    ###########################################################################
    # Calculate MV and AV flux
    # Probe vector data of source points with interpolation: vtkProbeFilter
    ###########################################################################
    numFrames = s.GetNumberOfTimeFrames() - 1
    print(f"num frames = {numFrames}")
    totalFlux = np.zeros((2, numFrames))
    for valve in range(2):
        if valve == 0:
            planeSource = planeSourceMV
        else:
            planeSource = planeSourceAV

        ptOrigin = np.array(planeSource.GetOrigin())
        ptPoint1 = np.array(planeSource.GetPoint1())
        ptPoint2 = np.array(planeSource.GetPoint2())
        vectorAxis1 = ptPoint1 - ptOrigin
        vectorAxis2 = ptPoint2 - ptOrigin
        dx = np.linalg.norm(vectorAxis1) / planeSource.GetXResolution()
        dy = np.linalg.norm(vectorAxis2) / planeSource.GetYResolution()
        timeStep = s.GetTimeStep()
        print(f"Axis1 {vectorAxis1}: norm {np.linalg.norm(vectorAxis1)}")
        print(f"Axis2 {vectorAxis2}: norm {np.linalg.norm(vectorAxis2)}")
        print(f"Resolution X {planeSource.GetXResolution()} / Y {planeSource.GetYResolution()}")
        print(f"dx {dx} / dy {dy} / timestep {timeStep}")

        for t in range(numFrames):
            imgData = s.GetImageDataAtTimeindex(t)

            # vtkProbeFilter test
            probe = vtk.vtkProbeFilter()
            probe.SetSourceData(imgData)
            probe.SetInputData(planeSource.GetOutput())
            probe.Update()
            probeOutput = probe.GetOutput()

            dataScalar = probeOutput.GetPointData().GetScalars()
            nTuples = dataScalar.GetNumberOfTuples()
            # nTuples = 20
            # for i in range(nTuples):
            #     val = dataScalar.GetValue(i)
            #     print(f"scalar GetValue({i}): {val}")

            dataVector = probeOutput.GetPointData().GetVectors()
            nTuples = dataVector.GetNumberOfTuples()
            resultList = np.zeros(shape=(nTuples, 3))
            totalFlux[valve, t] = 0
            for i in range(nTuples):
                val = dataVector.GetTuple(i)
                # Find angle between velocity vectors and the plane normal vector
                angle = np.arccos(
                    np.clip(np.dot(val / np.linalg.norm(val), planeSource.GetNormal()), -1.0, 1.0)) * 180 / np.pi
                # print(f"point {i} {probeOutput.GetPoint(i)} ==> vector {val}, angle {angle}, norm {np.linalg.norm(val)}")
                resultList[i] = [i, angle, np.linalg.norm(val)]
                # Flux: dx * dy * (v_z * t)
                # dx, dy, t: constant
                # v_z: project velocity vector on plane normal vector
                totalFlux[valve, t] += np.dot(val, planeSource.GetNormal()) * dx * dy * timeStep * 0.001 # effective velocity * dx * dy * t
                # 1 cubic mm = 0.001 cubic cm = 0.001 mL

            # print(f"result list {resultList}")
            resultListSorted = resultList[np.argsort(resultList[:, 1])]
            # print(f"sorted list {resultListSorted.tolist()}")
            # vtkViz.ShowVectorData(probeOutput, None)
            # print(f"Time frame {t} total flux = {totalFlux}")

    print(f"total flux = {totalFlux}")
    print(f"total MV flux = {np.sum(totalFlux[0,])} // total AV flux = {np.sum(totalFlux[1,])}")

    # Draw graph of MV and AV flux
    timeFrames = range(numFrames)
    plt.plot(timeFrames, totalFlux[0,], color='green', label="MV flux")
    plt.plot(timeFrames, totalFlux[1,], color='orange', label="AV flux")
    plt.legend()
    plt.show()









    ###########################################################################
    # LA Trace Particle in vector field: vtkParticlePathFilter
    # Forward tracing from the pulmonary veins
    ###########################################################################
    # There are three particle tracers for unsteady vector fields: vtkParticleTracer, vtkParticlePathFilter, vtkStreaklineFilter
    # Pathline: vtkParticlePathFilter

    # Calculate sum of flows
    sumFluxWholeFlow = 0
    sumFluxDirectFlow = 0
    sumFluxRetainedInflow = 0
    sumFluxOtherFlow = 0

    # Gather retained inflow endpoints
    pointsRetainedInflowEnd = vtk.vtkPoints()
    arrayRetainedInflowFlux = vtk.vtkDoubleArray()
    arrayRetainedInflowFlux.SetNumberOfComponents(1)

    # Run tracing from each frame to IVR frame.
    offsetIVR = 9
    for i in range(numFrames):
        frameStart = (i + offsetIVR) % numFrames
        print(f"Tracing {i} / {numFrames}: Frame {frameStart}, {-i}")
        # Run particle tracking in 4D flow vector field (s) from emission plane (planeSourceMV) forward (bInverse = False)
        # from 9th frame (9) for a full cycle (0)
        pathOutput = GetParticlePathTrackingOutput2(s, planeSourcePVs, False, frameStart, -i)

        # Add flux calculation to pathline output
        AddFluxCalculation(pathOutput, planeSourceMV, s, frameStart, 1)

        # Separate flow types
        [pathDirectFlow, pathRetainedInflow, pathOtherFlow] = AnalysisSeparateDirectFlow(pathOutput, planeSourceMV,
                                                                                         planeSourceAV)

        # Calculate sum of flows
        sumFluxWholeFlow += CalculateSumCelldataScalars(pathOutput)
        sumFluxDirectFlow += CalculateSumCelldataScalars(pathDirectFlow)
        sumFluxRetainedInflow += CalculateSumCelldataScalars(pathRetainedInflow)
        sumFluxOtherFlow += CalculateSumCelldataScalars(pathOtherFlow)

        # Gather retained inflow endpoints
        AppendEndpoints(pathRetainedInflow, pointsRetainedInflowEnd, arrayRetainedInflowFlux)

        #######################################################################
        ## Visualize pathlines
        #######################################################################

        # # Create actors to visualize entrance and exit planes
        planeEntranceActor = vtkViz.GetPlaneSourceActor(planeSourceMV, "Green")
        planeExitActor = vtkViz.GetPlaneSourceActor(planeSourceAV, "Yellow")
        #
        # # Visualize whole pathlines
        # actorPathLines = vtkViz.GetPolydataTubeActorSpecifiedRadius(pathOutput, 5, "Tomato")
        # vtkViz.ShowActors([actorPathLines, planeEntranceActor, planeExitActor])
        #
        # # Visualize sorted pathlines
        # actorPathDirectFlow = vtkViz.GetPolydataTubeActorSpecifiedRadius(pathDirectFlow, 5, "Green")
        # actorPathRetainedInflow = vtkViz.GetPolydataTubeActorSpecifiedRadius(pathRetainedInflow, 5, "Yellow")
        # actorPathOtherFlow = vtkViz.GetPolydataTubeActorSpecifiedRadius(pathOtherFlow, 5, "Violet")
        # vtkViz.ShowActors([actorPathDirectFlow, actorPathRetainedInflow, actorPathOtherFlow, planeEntranceActor, planeExitActor])

        ## End ################################################################

    # Calculate sum of flows
    print(f"sumFluxWholeFlow      {sumFluxWholeFlow}")
    print(f"sumFluxDirectFlow     {sumFluxDirectFlow}")
    print(f"sumFluxRetainedInflow {sumFluxRetainedInflow}")
    print(f"sumFluxOtherFlow      {sumFluxOtherFlow}")

    # Draw a line graph

    # Visualize endpoints with mask
    actorPointsRetainedInflowEnd = vtkViz.GetPointsActor(pointsRetainedInflowEnd, [0, 255, 0])
    vtkViz.ShowActors([maskActor, actorPointsRetainedInflowEnd], coordinateMarkerOffset=pointsPlaneAV[0])

    insidePoints, outsidePoints = DetermineInsideMask(pointsRetainedInflowEnd, maskPolyData)
    actorPointsRetainedInflowEndInsideMask = vtkViz.GetPointsActor(insidePoints, [0, 255, 0])
    actorPointsRetainedInflowEndOutsideMask = vtkViz.GetPointsActor(outsidePoints, [255, 0, 0])
    vtkViz.ShowActors([maskActor, actorPointsRetainedInflowEndInsideMask, actorPointsRetainedInflowEndOutsideMask],
                      coordinateMarkerOffset=pointsPlaneAV[0])

    ###########################################################################
    ## Re-run path tracking from retained inflow endpoints
    ###########################################################################

    # Re-run path tracking from retained inflow endpoints
    pointSourceRetainedInflow = PolyPointSource.PolyPointSource()
    pointSourceRetainedInflow.SetPoints(pointsRetainedInflowEnd)
    pointSourceRetainedInflow.Update()

    # Run tracking from IVR for a full heart cycle
    pathOutputRerunFromRetainedInflow = GetParticlePathTrackingOutput(s, pointSourceRetainedInflow, False, offsetIVR, 0)
    # Combine flux array with pathlines
    pathOutputRerunFromRetainedInflow.GetCellData().SetScalars(arrayRetainedInflowFlux)
    MatchFluxScalarAfterTracingFromPointSource(pathOutputRerunFromRetainedInflow, arrayRetainedInflowFlux)

    actorPathLines = vtkViz.GetPolydataTubeActorSpecifiedRadius(pathOutputRerunFromRetainedInflow, 5, "Tomato")
    vtkViz.ShowActors([actorPathLines, planeEntranceActor, planeExitActor], coordinateMarkerOffset=pointsPlaneAV[0])

    # Separate flow types
    [pathPartDelayedEjection, pathResidualVolume, pathRerunOtherFlow] = AnalysisSeparateDirectFlow(
        pathOutputRerunFromRetainedInflow, planeSourceMV, planeSourceAV)

    # Calculate sum of flows
    fluxPartDelayedEjection = CalculateSumCelldataScalars(pathPartDelayedEjection)
    fluxResidualVolume = CalculateSumCelldataScalars(pathResidualVolume)
    fluxRerunOtherFlow = CalculateSumCelldataScalars(pathRerunOtherFlow)
    print(
        f"fluxPartDelayedEjection {fluxPartDelayedEjection}, fluxResidualVolume {fluxResidualVolume}, fluxRerunOtherFlow {fluxRerunOtherFlow}")

    actorPartDelayedEjection = vtkViz.GetPolydataTubeActorSpecifiedRadius(pathPartDelayedEjection, 5, "Blue")
    actorResidualVolume = vtkViz.GetPolydataTubeActorSpecifiedRadius(pathResidualVolume, 5, "Red")
    actorRerunOtherFlow = vtkViz.GetPolydataTubeActorSpecifiedRadius(pathRerunOtherFlow, 5, "Violet")
    vtkViz.ShowActors(
        [actorPartDelayedEjection, actorResidualVolume, actorRerunOtherFlow, planeEntranceActor, planeExitActor],
        coordinateMarkerOffset=pointsPlaneAV[0])

    ###########################################################################
    # End
    ###########################################################################











    # ###########################################################################
    # # Trace Particle in vector field: vtkParticlePathFilter
    # # Backward tracing from the exit plane
    # ###########################################################################
    # # There are three particle tracers for unsteady vector fields: vtkParticleTracer, vtkParticlePathFilter, vtkStreaklineFilter
    # # Pathline: vtkParticlePathFilter
    #
    # # Calculate sum of flows
    # sumFluxWholeFlow = 0
    # sumFluxDirectFlow = 0
    # sumFluxDelayedEjection = 0
    # sumFluxOtherFlow = 0
    #
    # # Gather retained inflow endpoints
    # pointsDelayedEjectionEnd = vtk.vtkPoints()
    # arrayDelayedEjectionFlux = vtk.vtkDoubleArray()
    # arrayDelayedEjectionFlux.SetNumberOfComponents(1)
    #
    # # Run tracing from each frame to IVR frame.
    # offsetIVR = 9
    # iTracingReverse = -1
    # for i in range(numFrames):
    #     frameStart = (i * iTracingReverse + offsetIVR) % numFrames
    #     print(f"Tracing {i} / {numFrames}: Frame {frameStart}, {-i}")
    #     # Run particle tracking in 4D flow vector field (s) from exit plane (planeSourceAV) backward (bInverse = True)
    #     # from 9th frame (9) for a full cycle (0)
    #     pathOutput = GetParticlePathTrackingOutput(s, planeSourceAV, True, frameStart, -i)
    #     print("End tracing / Start Flux calculation")
    #
    #     # Add flux calculation to pathline output
    #     AddFluxCalculation(pathOutput, planeSourceAV, s, frameStart, 1)
    #     print("End Flux calculation")
    #
    #     # Separate flow types
    #     [pathDirectFlow, pathDelayedEjection, pathOtherFlow] = AnalysisSeparateDirectFlow(pathOutput, planeSourceAV,
    #                                                                                       planeSourceMV, bInverse=True)
    #
    #     # Calculate sum of flows
    #     sumFluxWholeFlow += CalculateSumCelldataScalars(pathOutput)
    #     sumFluxDirectFlow += CalculateSumCelldataScalars(pathDirectFlow)
    #     sumFluxDelayedEjection += CalculateSumCelldataScalars(pathDelayedEjection)
    #     sumFluxOtherFlow += CalculateSumCelldataScalars(pathOtherFlow)
    #
    #     # Gather retained inflow endpoints
    #     AppendEndpoints(pathDelayedEjection, pointsDelayedEjectionEnd, arrayDelayedEjectionFlux)
    #
    #     #######################################################################
    #     ## Visualize pathlines
    #     #######################################################################
    #
    #     # Create actors to visualize entrance and exit planes
    #     planeEntranceActor = vtkViz.GetPlaneSourceActor(planeSourceMV, "Green")
    #     planeExitActor = vtkViz.GetPlaneSourceActor(planeSourceAV, "Yellow")
    #
    #     # # Visualize whole pathlines
    #     # actorPathLines = vtkViz.GetPolydataTubeActorSpecifiedRadius(pathOutput, 5, "Tomato")
    #     # vtkViz.ShowActors([actorPathLines, planeEntranceActor, planeExitActor], coordinateMarkerOffset=pointsPlaneAV[0])
    #
    #     # Visualize sorted pathlines
    #     actorPathDirectFlow = vtkViz.GetPolydataTubeActorSpecifiedRadius(pathDirectFlow, 5, "Green")
    #     actorPathDelayedEjection = vtkViz.GetPolydataTubeActorSpecifiedRadius(pathDelayedEjection, 5, "Blue")
    #     actorPathOtherFlow = vtkViz.GetPolydataTubeActorSpecifiedRadius(pathOtherFlow, 5, "Violet")
    #     vtkViz.ShowActors(
    #         [actorPathDirectFlow, actorPathDelayedEjection, actorPathOtherFlow, planeEntranceActor, planeExitActor],
    #         coordinateMarkerOffset=pointsPlaneAV[0])
    #
    #     ## End ################################################################
    #
    # # Calculate sum of flows
    # print(f"sumFluxWholeFlow      {sumFluxWholeFlow}")
    # print(f"sumFluxDirectFlow     {sumFluxDirectFlow}")
    # print(f"sumFluxDelayedEjection {sumFluxDelayedEjection}")
    # print(f"sumFluxOtherFlow      {sumFluxOtherFlow}")
    #
    # # Draw a line graph
    #
    # # Visualize endpoints with mask
    # actorPointsDelayedEjectionEnd = vtkViz.GetPointsActor(pointsDelayedEjectionEnd, [0, 255, 0])
    # vtkViz.ShowActors([maskActor, actorPointsDelayedEjectionEnd], coordinateMarkerOffset=pointsPlaneAV[0])
    #
    # insidePoints, outsidePoints = DetermineInsideMask(pointsDelayedEjectionEnd, maskPolyData)
    # actorPointsDelayedEjectionEndInsideMask = vtkViz.GetPointsActor(insidePoints, [0, 255, 0])
    # actorPointsDelayedEjectionEndOutsideMask = vtkViz.GetPointsActor(outsidePoints, [255, 0, 0])
    # vtkViz.ShowActors([maskActor, actorPointsDelayedEjectionEndInsideMask, actorPointsDelayedEjectionEndOutsideMask],
    #                   coordinateMarkerOffset=pointsPlaneAV[0])
    #
    # ###########################################################################
    # ## Re-run path tracking from retained inflow endpoints
    # ###########################################################################
    #
    # # Re-run path tracking from retained inflow endpoints
    # pointSourceDelayedEjection = PolyPointSource.PolyPointSource()
    # pointSourceDelayedEjection.SetPoints(pointsDelayedEjectionEnd)
    # pointSourceDelayedEjection.Update()
    #
    # # Run tracking from IVR for a full heart cycle
    # pathOutputRerunFromDelayedEjection = GetParticlePathTrackingOutput(s, pointSourceDelayedEjection, False, offsetIVR,
    #                                                                    0)
    # # Combine flux array with pathlines
    # pathOutputRerunFromDelayedEjection.GetCellData().SetScalars(arrayDelayedEjectionFlux)
    # MatchFluxScalarAfterTracingFromPointSource(pathOutputRerunFromDelayedEjection, arrayDelayedEjectionFlux)
    #
    # actorPathLines = vtkViz.GetPolydataTubeActorSpecifiedRadius(pathOutputRerunFromDelayedEjection, 5, "Tomato")
    # vtkViz.ShowActors([actorPathLines, planeEntranceActor, planeExitActor], coordinateMarkerOffset=pointsPlaneAV[0])
    #
    # # Separate flow types
    # [pathPartDelayedEjection, pathResidualVolume, pathRerunOtherFlow] = AnalysisSeparateDirectFlow(
    #     pathOutputRerunFromDelayedEjection, planeSourceMV, planeSourceAV, bInverse=True)
    #
    # # Calculate sum of flows
    # fluxPartDelayedEjection = CalculateSumCelldataScalars(pathPartDelayedEjection)
    # fluxResidualVolume = CalculateSumCelldataScalars(pathResidualVolume)
    # fluxRerunOtherFlow = CalculateSumCelldataScalars(pathRerunOtherFlow)
    # print(
    #     f"fluxPartDelayedEjection {fluxPartDelayedEjection}, fluxResidualVolume {fluxResidualVolume}, fluxRerunOtherFlow {fluxRerunOtherFlow}")
    #
    # actorPartDelayedEjection = vtkViz.GetPolydataTubeActorSpecifiedRadius(pathPartDelayedEjection, 5, "Blue")
    # actorResidualVolume = vtkViz.GetPolydataTubeActorSpecifiedRadius(pathResidualVolume, 5, "Red")
    # actorRerunOtherFlow = vtkViz.GetPolydataTubeActorSpecifiedRadius(pathRerunOtherFlow, 5, "Violet")
    # vtkViz.ShowActors(
    #     [actorPartDelayedEjection, actorResidualVolume, actorRerunOtherFlow, planeEntranceActor, planeExitActor],
    #     coordinateMarkerOffset=pointsPlaneAV[0])
    #
    # ###########################################################################
    # # End
    # ###########################################################################

    ###########################################################################
    # Trace Particle in vector field: vtkParticlePathFilter
    # Forward tracing from the entrance plane
    ###########################################################################
    # There are three particle tracers for unsteady vector fields: vtkParticleTracer, vtkParticlePathFilter, vtkStreaklineFilter
    # Pathline: vtkParticlePathFilter

    # Calculate sum of flows
    sumFluxWholeFlow = 0
    sumFluxDirectFlow = 0
    sumFluxRetainedInflow = 0
    sumFluxOtherFlow = 0

    # Gather retained inflow endpoints
    pointsRetainedInflowEnd = vtk.vtkPoints()
    arrayRetainedInflowFlux = vtk.vtkDoubleArray()
    arrayRetainedInflowFlux.SetNumberOfComponents(1)

    # Run tracing from each frame to IVR frame.
    offsetIVR = 9
    for i in range(numFrames):
        frameStart = (i + offsetIVR) % numFrames
        print(f"Tracing {i} / {numFrames}: Frame {frameStart}, {-i}")
        # Run particle tracking in 4D flow vector field (s) from emission plane (planeSourceMV) forward (bInverse = False)
        # from 9th frame (9) for a full cycle (0)
        pathOutput = GetParticlePathTrackingOutput(s, planeSourceMV, False, frameStart, -i)

        # Add flux calculation to pathline output
        AddFluxCalculation(pathOutput, planeSourceMV, s, frameStart, 1)

        # Separate flow types
        [pathDirectFlow, pathRetainedInflow, pathOtherFlow] = AnalysisSeparateDirectFlow(pathOutput, planeSourceMV,
                                                                                         planeSourceAV)

        # Calculate sum of flows
        sumFluxWholeFlow += CalculateSumCelldataScalars(pathOutput)
        sumFluxDirectFlow += CalculateSumCelldataScalars(pathDirectFlow)
        sumFluxRetainedInflow += CalculateSumCelldataScalars(pathRetainedInflow)
        sumFluxOtherFlow += CalculateSumCelldataScalars(pathOtherFlow)

        # Gather retained inflow endpoints
        AppendEndpoints(pathRetainedInflow, pointsRetainedInflowEnd, arrayRetainedInflowFlux)

        #######################################################################
        ## Visualize pathlines
        #######################################################################

        # # Create actors to visualize entrance and exit planes
        planeEntranceActor = vtkViz.GetPlaneSourceActor(planeSourceMV, "Green")
        planeExitActor = vtkViz.GetPlaneSourceActor(planeSourceAV, "Yellow")
        #
        # # Visualize whole pathlines
        # actorPathLines = vtkViz.GetPolydataTubeActorSpecifiedRadius(pathOutput, 5, "Tomato")
        # vtkViz.ShowActors([actorPathLines, planeEntranceActor, planeExitActor])
        #
        # # Visualize sorted pathlines
        # actorPathDirectFlow = vtkViz.GetPolydataTubeActorSpecifiedRadius(pathDirectFlow, 5, "Green")
        # actorPathRetainedInflow = vtkViz.GetPolydataTubeActorSpecifiedRadius(pathRetainedInflow, 5, "Yellow")
        # actorPathOtherFlow = vtkViz.GetPolydataTubeActorSpecifiedRadius(pathOtherFlow, 5, "Violet")
        # vtkViz.ShowActors([actorPathDirectFlow, actorPathRetainedInflow, actorPathOtherFlow, planeEntranceActor, planeExitActor])

        ## End ################################################################

    # Calculate sum of flows
    print(f"sumFluxWholeFlow      {sumFluxWholeFlow}")
    print(f"sumFluxDirectFlow     {sumFluxDirectFlow}")
    print(f"sumFluxRetainedInflow {sumFluxRetainedInflow}")
    print(f"sumFluxOtherFlow      {sumFluxOtherFlow}")

    # Draw a line graph

    # Visualize endpoints with mask
    actorPointsRetainedInflowEnd = vtkViz.GetPointsActor(pointsRetainedInflowEnd, [0, 255, 0])
    vtkViz.ShowActors([maskActor, actorPointsRetainedInflowEnd], coordinateMarkerOffset=pointsPlaneAV[0])

    insidePoints, outsidePoints = DetermineInsideMask(pointsRetainedInflowEnd, maskPolyData)
    actorPointsRetainedInflowEndInsideMask = vtkViz.GetPointsActor(insidePoints, [0, 255, 0])
    actorPointsRetainedInflowEndOutsideMask = vtkViz.GetPointsActor(outsidePoints, [255, 0, 0])
    vtkViz.ShowActors([maskActor, actorPointsRetainedInflowEndInsideMask, actorPointsRetainedInflowEndOutsideMask],
                      coordinateMarkerOffset=pointsPlaneAV[0])

    ###########################################################################
    ## Re-run path tracking from retained inflow endpoints
    ###########################################################################

    # Re-run path tracking from retained inflow endpoints
    pointSourceRetainedInflow = PolyPointSource.PolyPointSource()
    pointSourceRetainedInflow.SetPoints(pointsRetainedInflowEnd)
    pointSourceRetainedInflow.Update()

    # Run tracking from IVR for a full heart cycle
    pathOutputRerunFromRetainedInflow = GetParticlePathTrackingOutput(s, pointSourceRetainedInflow, False, offsetIVR, 0)
    # Combine flux array with pathlines
    pathOutputRerunFromRetainedInflow.GetCellData().SetScalars(arrayRetainedInflowFlux)
    MatchFluxScalarAfterTracingFromPointSource(pathOutputRerunFromRetainedInflow, arrayRetainedInflowFlux)

    actorPathLines = vtkViz.GetPolydataTubeActorSpecifiedRadius(pathOutputRerunFromRetainedInflow, 5, "Tomato")
    vtkViz.ShowActors([actorPathLines, planeEntranceActor, planeExitActor], coordinateMarkerOffset=pointsPlaneAV[0])

    # Separate flow types
    [pathPartDelayedEjection, pathResidualVolume, pathRerunOtherFlow] = AnalysisSeparateDirectFlow(
        pathOutputRerunFromRetainedInflow, planeSourceMV, planeSourceAV)

    # Calculate sum of flows
    fluxPartDelayedEjection = CalculateSumCelldataScalars(pathPartDelayedEjection)
    fluxResidualVolume = CalculateSumCelldataScalars(pathResidualVolume)
    fluxRerunOtherFlow = CalculateSumCelldataScalars(pathRerunOtherFlow)
    print(
        f"fluxPartDelayedEjection {fluxPartDelayedEjection}, fluxResidualVolume {fluxResidualVolume}, fluxRerunOtherFlow {fluxRerunOtherFlow}")

    actorPartDelayedEjection = vtkViz.GetPolydataTubeActorSpecifiedRadius(pathPartDelayedEjection, 5, "Blue")
    actorResidualVolume = vtkViz.GetPolydataTubeActorSpecifiedRadius(pathResidualVolume, 5, "Red")
    actorRerunOtherFlow = vtkViz.GetPolydataTubeActorSpecifiedRadius(pathRerunOtherFlow, 5, "Violet")
    vtkViz.ShowActors(
        [actorPartDelayedEjection, actorResidualVolume, actorRerunOtherFlow, planeEntranceActor, planeExitActor],
        coordinateMarkerOffset=pointsPlaneAV[0])

    ###########################################################################
    # End
    ###########################################################################

    ## Output analysis: When does it cross the plane?
    innerPoints = AnalysisCrossPlane(pathOutput, planeSourceAV)
    # innerPoints = AnalysisCrossPlane(pathOutput, planeSourceMV)

    ## Re-run path tracking from Inner Points
    pointSourceRetainedInflow = PolyPointSource.PolyPointSource()
    pointSourceRetainedInflow.SetPoints(innerPoints)
    pointSourceRetainedInflow.Update()
    print(innerPoints)
    pathfilter.RemoveAllSources()
    pathfilter.AddSourceConnection(pointSourceRetainedInflow.GetOutputPort())
    print(pointSourceRetainedInflow.GetOutputPort())
    # Run Pathline filter
    pathfilter.Update()
    pathOutput = pathfilter.GetOutput()
    print(pathOutput)
    innerPoints = AnalysisCrossPlane(pathOutput, planeSourceAV)

    # Visualization
    insidePointActor = vtkViz.GetPointsActor(insidePoints, [255, 255, 0])
    outsidePointActor = vtkViz.GetPointsActor(outsidePoints, [0, 0, 255])

    vtkViz.ShowVectorData(imgData, pathlineActor, maskActor, insidePointActor, outsidePointActor)
    vtkViz.ShowActors([maskActor, insidePointActor, outsidePointActor])

    ###########################################################################
    # Unused part: StreamLine
    ###########################################################################
    # ## Streamline
    # # for i in [9]:
    # for direction in [vtk.vtkStreamTracer().BOTH, vtk.vtkStreamTracer().FORWARD, vtk.vtkStreamTracer().BACKWARD]:
    #     i = 3
    #     # for i in range(15):
    #     imgData = s.GetImageDataAtTimeindex(i)
    #     # showVectorData(imgData, None)
    #
    #     # Render Streamline
    #     rk = vtk.vtkRungeKutta45()
    #     # Create source for streamtubes
    #     streamer = vtk.vtkStreamTracer()
    #     streamer.SetInputData(imgData)
    #     streamer.SetSourceConnection(planeSourceMV.GetOutputPort())
    #     # streamer.SetSourceConnection(pointSourceMV.GetOutputPort())
    #     streamer.SetMaximumPropagation(50)
    #     streamer.SetIntegrationStepUnit(2)
    #     streamer.SetMinimumIntegrationStep(0.1)
    #     streamer.SetMaximumIntegrationStep(1.0)
    #     streamer.SetInitialIntegrationStep(0.2)
    #     streamer.SetIntegrationDirection(0)
    #     streamer.SetIntegrator(rk)
    #     streamer.SetRotationScale(0.5)
    #     streamer.SetMaximumError(1.0e-8)
    #     # streamer.SetIntegrationDirection(vtk.vtkStreamTracer().BACKWARD)
    #     streamer.SetIntegrationDirection(direction)
    #     mapStream = vtk.vtkPolyDataMapper()
    #     mapStream.SetInputConnection(streamer.GetOutputPort())
    #     mapStream.SetScalarRange(imgData.GetScalarRange())
    #     streamActor = vtk.vtkActor()
    #     streamActor.SetMapper(mapStream)
    #
    #     # vtkViz.ShowVectorData(imgData, streamActor, maskActor)
    #     # vtkViz.ShowActors([streamActor, maskActor])


def MatchFluxScalarAfterTracingFromPointSource(polyData, arrayScalar):
    # Get the number of tuples in Cell data
    nTuples = arrayScalar.GetNumberOfTuples()
    # Get the number of cells in the polydata. This should be same as the number of Scalar tuples.
    idNumCellsInPolyData = polyData.GetNumberOfCells()

    print(
        f"MatchFluxScalarAfterTracingFromPointSource: numCellDataScalarTuples {nTuples} // numCellsInPolyData {idNumCellsInPolyData}")

    arrayScalarNew = vtk.vtkDoubleArray()
    arrayScalarNew.SetNumberOfComponents(1)

    # Iterate each PolyLine
    for i in range(idNumCellsInPolyData):
        # Retrieve each line
        polyLine = polyData.GetCell(i)  # Get a cell (PolyLine)

        # Find flux of the first point
        ptId = polyLine.GetPointIds().GetId(0)
        # print(f"i {i}: first point {ptId}")
        scalar = arrayScalar.GetValue(ptId)  # Get Scalar Value
        # print(f"i {i}: first point {ptId} ==> {scalar}")

        arrayScalarNew.InsertNextTypedTuple([scalar])

    # Get Scalars of Cell data
    polyData.GetCellData().SetScalars(arrayScalarNew)

    return


def AppendEndpoints(polyData, pointsEndpoints, arrayFlux):
    # Get Scalars of Cell data
    arrayScalar = polyData.GetCellData().GetScalars()
    # Get the number of tuples in Cell data
    nTuples = arrayScalar.GetNumberOfTuples()
    # Get the number of cells in the polydata. This should be same as the number of Scalar tuples.
    idNumCellsInPolyData = polyData.GetNumberOfCells()

    if (nTuples != idNumCellsInPolyData):
        print(f"Error: AppendEndpoints: numCellDataScalarTuples {nTuples} != numCellsInPolyData {idNumCellsInPolyData}")
        return

    # Iterate each PolyLine
    for i in range(idNumCellsInPolyData):
        # Retrieve each line and data
        polyLine = polyData.GetCell(i)  # Get a cell (PolyLine)
        scalar = arrayScalar.GetValue(i)  # Get Scalar Value

        # # Find flux of the first point
        # ptId = polyLine.GetPointIds().GetId(0)
        # print(f"i {i}: first point {ptId}")
        #
        # scalar = arrayScalar.GetValue(ptId)    # Get Scalar Value
        # print(f"i {i}: first point {ptId} ==> {scalar}")

        # Find the last point
        numIds = polyLine.GetPointIds().GetNumberOfIds()
        ptId = polyLine.GetPointIds().GetId(numIds - 1)
        pt = polyData.GetPoint(ptId)

        # Append the endpoint and flux
        pointsEndpoints.InsertNextPoint(pt)
        arrayFlux.InsertNextTypedTuple([scalar])

    # print(f"num points {pointsEndpoints.GetNumberOfPoints()} // num flux {arrayFlux.GetNumberOfTuples()}")
    return


def CalculateSumCelldataScalars(polyData):
    # Get Scalars of Cell data
    arrayScalar = polyData.GetCellData().GetScalars()
    npArrayScalar = numpy_support.vtk_to_numpy(arrayScalar)
    # print(f"Array {npArrayScalar}")
    # print(f"Sum {np.sum(npArrayScalar)}")
    return np.sum(npArrayScalar)


## Report location of points into text file
def ReportFinalPoints(insidePoints, outsidePoints):
    f = open(f'report_finalPoints.txt', 'w')

    # print(f"##### Points inside Mask #####")
    print(f"##### Points inside Mask #####", file=f)
    numPoints = insidePoints.GetNumberOfPoints()
    for i in range(numPoints):
        point = insidePoints.GetPoint(i)
        # print(f"{i}: {point}")
        print(f"{i}: {point}", file=f)

    # print(f"##### Points outside Mask #####")
    print(f"##### Points outside Mask #####", file=f)
    numPoints = outsidePoints.GetNumberOfPoints()
    for i in range(numPoints):
        point = outsidePoints.GetPoint(i)
        # print(f"{i}: {point}")
        print(f"{i}: {point}", file=f)

    f.close()


## Extract Final points from Pathline output (vtkPolyData)
#
#  @param polydata vtkPolyData. Pathline output.
def ExtractFinalPointsFromPathlines(polydata):
    # Get the number of points in the polydata
    idNumPointsInFile = polydata.GetNumberOfPoints()
    # numArray = polydata.GetPointData().GetNumberOfArrays()

    # array 3: SimulationTimeStep
    # Get Maximum from simulation time step ==> Find last time index
    arrayTimeStep = polydata.GetPointData().GetAbstractArray(3)
    # print(f"***** array 3: {arrayTimeStep.GetName()} with {idNumPointsInFile} *****")
    maxValue = arrayTimeStep.GetDataTypeValueMin()
    for j in range(idNumPointsInFile):
        value = arrayTimeStep.GetValue(j)
        if (value > maxValue):
            maxValue = value

    # Gather points only on the last time step.
    # Reiterate the time step array, and add points if it is from the last time step.
    # pointIDs = []
    points = vtk.vtkPoints()
    for j in range(idNumPointsInFile):
        value = arrayTimeStep.GetValue(j)
        if (value == maxValue):
            # pointIDs.append(j)
            points.InsertNextPoint(polydata.GetPoint(j))

    return points


## Determine if points are inside or outside the mask
#
#  @param points vtkPoints. Contains points to be judged.
#  @param mask   3D volume mask marked with 1 and 0.
def DetermineInsideMask(points, mask):
    pointsPolydata = vtk.vtkPolyData()
    pointsPolydata.SetPoints(points)

    # Points inside test
    selectEnclosedPoints = vtk.vtkSelectEnclosedPoints()
    selectEnclosedPoints.SetInputData(pointsPolydata)
    selectEnclosedPoints.SetSurfaceData(mask)
    selectEnclosedPoints.Update()

    insidePoints = vtk.vtkPoints()
    outsidePoints = vtk.vtkPoints()

    for i in range(points.GetNumberOfPoints()):
        if (selectEnclosedPoints.IsInside(i)):
            insidePoints.InsertNextPoint(points.GetPoint(i))
        else:
            outsidePoints.InsertNextPoint(points.GetPoint(i))
        # print(f"Point {i}: {selectEnclosedPoints.IsInside(i)}")

    # insideArray = selectEnclosedPoints.GetOutput().GetPointData().GetArray("SelectedPoints")
    #
    # for i in range(insideArray.GetNumberOfTuples()):
    #     print(f"{i}: {insideArray.GetComponent(i, 0)}")

    return insidePoints, outsidePoints


## Analyze pathlines cross the given plane
#
#  This function can be expanded for arbiturary shape of plane, if polygon is given as the input parameter.
#  @param polydata Pathline output
#  @param planeSource vtkPlaneSource. Defines a finite plane.
def AnalysisCrossPlane(polydata, planeSource):
    ## Prepare rectangular polygon from plane source
    # Create a square in the x-y plane.
    thirdPoint = np.array(planeSource.GetPoint1()) + np.array(planeSource.GetPoint2()) - np.array(
        planeSource.GetOrigin())
    points = vtk.vtkPoints()
    points.InsertNextPoint(planeSource.GetOrigin())
    points.InsertNextPoint(planeSource.GetPoint1())
    points.InsertNextPoint(thirdPoint)
    points.InsertNextPoint(planeSource.GetPoint2())

    # Create the polygon with four corner points
    polygon = vtk.vtkPolygon()
    polygon.GetPoints().DeepCopy(points)
    polygon.GetPointIds().SetNumberOfIds(4)  # The 4 corners of the square
    for i in range(4):
        polygon.GetPointIds().SetId(i, i)

    planeNormal = planeSource.GetNormal()

    polyDataInnerLines = vtk.vtkPolyData()
    polyDataOuterLines = vtk.vtkPolyData()
    polyDataCrossLines = vtk.vtkPolyData()
    pointsInnerLines = vtk.vtkPoints()
    pointsOuterLines = vtk.vtkPoints()
    pointsCrossLines = vtk.vtkPoints()
    cellsInnerLines = vtk.vtkCellArray()
    cellsOuterLines = vtk.vtkCellArray()
    cellsCrossLines = vtk.vtkCellArray()

    ## Iterate each pathline and find intersection with plane.
    # Get the number of points in the polydata
    idNumCellsInFile = polydata.GetNumberOfCells()
    print(f"### Number of Cells {idNumCellsInFile}")
    for i in range(idNumCellsInFile):
        cell = polydata.GetCell(i)  # PolyLine
        # print(f"line {i}: {cell.GetNumberOfPoints()} points")

        # Iterate each line segment
        numIds = cell.GetPointIds().GetNumberOfIds()
        pt = None
        prevState = "init"
        prevCrossPtIndex = 0
        for j in range(numIds):
            ptId = cell.GetPointIds().GetId(j)
            pt_prev = pt
            pt = polydata.GetPoint(ptId)
            if (j > 0):
                tolerance = 0.001
                t = vtk.mutable(
                    0.0)  # Parametric coordinate of intersection (0 (corresponding to p1) to 1 (corresponding to p2))
                x = [0.0] * 3  # The coordinate of the intersection
                pcoords = [0.0] * 3
                subId = vtk.mutable(0)
                iD = polygon.IntersectWithLine(pt, pt_prev, tolerance, t, x, pcoords, subId)
                # If there is an intersection between the given plane and the line segment
                if (iD == 1):
                    vectorLineSegment = np.array(pt) - np.array(pt_prev)
                    innerP = np.inner(vectorLineSegment, planeNormal)
                    strInward = "outward" if innerP > 0 else "inward"
                    # print(f"Intersection line {i} segment {j} at {x}, {strInward}")

                    if (strInward == "outward"):
                        if ((prevState == "init") or (prevState == "inward")):
                            # Before cross point: InnerLines
                            polyLine = vtk.vtkPolyLine()
                            polyLine.GetPointIds().SetNumberOfIds(j - prevCrossPtIndex)
                            for k in range(prevCrossPtIndex, j):
                                # polyDataInnerLines.InsertNextPoint(polydata.GetPoint(cell.GetPointIds().GetId(k)))
                                newPtId = pointsInnerLines.InsertNextPoint(
                                    polydata.GetPoint(cell.GetPointIds().GetId(k)))
                                polyLine.GetPointIds().SetId(k - prevCrossPtIndex, newPtId)
                            cellsInnerLines.InsertNextCell(polyLine)

                            # Cross segment: CrossLines
                            polyLine = vtk.vtkPolyLine()
                            polyLine.GetPointIds().SetNumberOfIds(2)
                            newPtId = pointsCrossLines.InsertNextPoint(
                                polydata.GetPoint(cell.GetPointIds().GetId(j - 1)))
                            polyLine.GetPointIds().SetId(0, newPtId)
                            newPtId = pointsCrossLines.InsertNextPoint(polydata.GetPoint(cell.GetPointIds().GetId(j)))
                            polyLine.GetPointIds().SetId(1, newPtId)
                            cellsCrossLines.InsertNextCell(polyLine)

                        else:  # if prevState == "outward"
                            # outward after outward. something wrong. returns volume out of the valve plane.
                            print(f"*** Error: Outward after outward. line {i}, segment {j}")

                        prevState = "outward"
                        prevCrossPtIndex = j

                    else:  # strInward == "inward"
                        if ((prevState == "init") or (prevState == "outward")):
                            # Before cross point: OuterLines
                            polyLine = vtk.vtkPolyLine()
                            polyLine.GetPointIds().SetNumberOfIds(j - prevCrossPtIndex)
                            for k in range(prevCrossPtIndex, j):
                                # polyDataInnerLines.InsertNextPoint(polydata.GetPoint(cell.GetPointIds().GetId(k)))
                                newPtId = pointsOuterLines.InsertNextPoint(
                                    polydata.GetPoint(cell.GetPointIds().GetId(k)))
                                polyLine.GetPointIds().SetId(k - prevCrossPtIndex, newPtId)
                            cellsOuterLines.InsertNextCell(polyLine)

                            # Cross segment: CrossLines
                            polyLine = vtk.vtkPolyLine()
                            polyLine.GetPointIds().SetNumberOfIds(2)
                            newPtId = pointsCrossLines.InsertNextPoint(
                                polydata.GetPoint(cell.GetPointIds().GetId(j - 1)))
                            polyLine.GetPointIds().SetId(0, newPtId)
                            newPtId = pointsCrossLines.InsertNextPoint(polydata.GetPoint(cell.GetPointIds().GetId(j)))
                            polyLine.GetPointIds().SetId(1, newPtId)
                            cellsCrossLines.InsertNextCell(polyLine)

                        else:  # if prevState == "inward"
                            # inward after inward. something wrong. returns volume out of the valve plane.
                            print(f"*** Error: Inward after inward. line {i}, segment {j}")

                        prevState = "inward"
                        prevCrossPtIndex = j

        # After iteration of each line segment, the last part should be added to inner or outer.
        if (prevState == "inward"):
            # After cross point: InnerLines
            polyLine = vtk.vtkPolyLine()
            polyLine.GetPointIds().SetNumberOfIds(numIds - prevCrossPtIndex)
            for k in range(prevCrossPtIndex, numIds):
                newPtId = pointsInnerLines.InsertNextPoint(
                    polydata.GetPoint(cell.GetPointIds().GetId(k)))
                polyLine.GetPointIds().SetId(k - prevCrossPtIndex, newPtId)
            cellsInnerLines.InsertNextCell(polyLine)
        elif (prevState == "outward"):
            # After cross point: OuterLines
            polyLine = vtk.vtkPolyLine()
            polyLine.GetPointIds().SetNumberOfIds(numIds - prevCrossPtIndex)
            for k in range(prevCrossPtIndex, numIds):
                newPtId = pointsOuterLines.InsertNextPoint(
                    polydata.GetPoint(cell.GetPointIds().GetId(k)))
                polyLine.GetPointIds().SetId(k - prevCrossPtIndex, newPtId)
            cellsOuterLines.InsertNextCell(polyLine)
        else:  # if prevState == "init"
            None

    polyDataInnerLines.SetPoints(pointsInnerLines)
    polyDataInnerLines.SetLines(cellsInnerLines)
    polyDataCrossLines.SetPoints(pointsCrossLines)
    polyDataCrossLines.SetLines(cellsCrossLines)
    polyDataOuterLines.SetPoints(pointsOuterLines)
    polyDataOuterLines.SetLines(cellsOuterLines)

    innerLinesActor = vtkViz.GetPolydataActor(polyDataInnerLines, "Tomato")
    crossLinesActor = vtkViz.GetPolydataActor(polyDataCrossLines, "Green")
    outerLinesActor = vtkViz.GetPolydataActor(polyDataOuterLines, "Yellow")
    planeActor = vtkViz.GetPolygonActor(polygon, "Gray")
    vtkViz.ShowActors([innerLinesActor, crossLinesActor, outerLinesActor, planeActor])

    print(pointsInnerLines.GetNumberOfPoints())
    return pointsInnerLines


## Calculate flux of each pathline, then Add them to the PolyData as Scalars of CellData
#
#  @param polydataPathline    Pathline output
#  @param planeSourceEntrance vtkPlaneSource. Defines a finite plane.
#  @param source4DFlow        4D flow source to extract velocity vector from
#  @param tFrameStart         Frame number to extract velocity vector from
#  @param scaleFactorFlux     Flux = dx * dy * (v_z * t). v_z will be calculated from a velocity vector and the normal vector of planeSourceEntrance.
#                             scaleFactorFlux = dx * dy * t
def AddFluxCalculation(polydataPathline, planeSourceEntrance, source4DFlow, tFrameStart, scaleFactorFlux):
    # Collect starting point of each line
    probePoints = vtk.vtkPoints()  # point array to store starting points
    idNumCellsInPolyData = polydataPathline.GetNumberOfCells()  # Get the number of points in the polydata
    for i in range(idNumCellsInPolyData):
        # Get each PolyLine.
        polyLineSingle = polydataPathline.GetCell(i)
        # Find the starting point of the line.
        ptId = polyLineSingle.GetPointIds().GetId(0)
        pt = polydataPathline.GetPoint(ptId)
        # Add the starting point to the point array.
        probePoints.InsertNextPoint(pt)

    # Run probe filter to get velocity vectors of starting points.
    imgData = source4DFlow.GetImageDataAtTimeindex(tFrameStart)
    probe = vtk.vtkProbeFilter()
    probe.SetSourceData(imgData)
    probePolyData = vtk.vtkPolyData()
    probePolyData.SetPoints(probePoints)
    probe.SetInputData(probePolyData)
    probe.Update()
    probeOutput = probe.GetOutput()
    dataVector = probeOutput.GetPointData().GetVectors()

    # Create a vtkDoubleArray container to store flux of each pathline
    scalarsFluxPerPathline = vtk.vtkDoubleArray()
    scalarsFluxPerPathline.SetNumberOfComponents(1)

    # Iterate each velocity vector
    numTuples = dataVector.GetNumberOfTuples()
    for i in range(numTuples):
        vectorSourcePoint = dataVector.GetTuple(i)
        # # Find the angle between velocity vectors and the plane normal vector
        # angle = np.arccos(np.clip(np.dot(vectorSourcePoint / np.linalg.norm(vectorSourcePoint), planeSourceEntrance.GetNormal()), -1.0, 1.0)) * 180 / np.pi
        # Find Flux at the source point
        flux = np.dot(vectorSourcePoint, planeSourceEntrance.GetNormal()) * scaleFactorFlux
        # Store flux in vtkDoubleArray
        scalarsFluxPerPathline.InsertNextTypedTuple([flux])
        # print(f"cell {i} at entrance norm {np.linalg.norm(vectorSourcePoint):.3f} angle {angle:>3.3f} flux {flux*100:>3.3f}")

    # Store flux array with pathlines
    polydataPathline.GetCellData().SetScalars(scalarsFluxPerPathline)


## Analyze pathlines to separate direct flow from other flow
# PolyDataPathLine contains many PathLines.
# This function will examine each PathLine whether it crosses plane source rectangle.
#
#  This function can be expanded for arbiturary shape of plane, if polygon is given as the input parameter.
#  @param polydata Pathline output
#  @param planeSource vtkPlaneSource. Defines a finite plane.
def AnalysisSeparateDirectFlow(polydataPathline, planeSourceEntrance, planeSourceExit, bInverse=False):
    # Prepare rectangular polygon from the exit plane.
    polygonExit = GetPolygonFromPlaneSource(planeSourceExit)

    # Plane Normal vector will be used to determine flow direction.
    if not bInverse:
        vectorExitNormal = planeSourceExit.GetNormal()
    else:
        vectorExitNormal = np.array(planeSourceExit.GetNormal()) * (-1)  # inverse case
        # vectorExitNormal = planeSourceExit.GetNormal()
    # print(f"vectorExitNormal {vectorExitNormal} // Inverse {np.array(planeSourceExit.GetNormal()) * (-1) }")
    ###########################################################################

    # Prepare variables to store direct flow and retained inflow separately.
    polyDataDirectFlow = vtk.vtkPolyData()
    polyDataDirectFlow.SetPoints(polydataPathline.GetPoints())
    polyDataRetainedInflow = vtk.vtkPolyData()
    polyDataRetainedInflow.SetPoints(polydataPathline.GetPoints())
    polyDataOtherFlow = vtk.vtkPolyData()
    polyDataOtherFlow.SetPoints(polydataPathline.GetPoints())

    cellsDirectFlow = vtk.vtkCellArray()
    cellsRetainedInflow = vtk.vtkCellArray()
    cellsOtherFlow = vtk.vtkCellArray()

    # Create a vtkDoubleArray container to store flux of each pathline
    scalarsDirectFlow = vtk.vtkDoubleArray()
    scalarsDirectFlow.SetNumberOfComponents(1)
    scalarsRetainedInflow = vtk.vtkDoubleArray()
    scalarsRetainedInflow.SetNumberOfComponents(1)
    scalarsOtherFlow = vtk.vtkDoubleArray()
    scalarsOtherFlow.SetNumberOfComponents(1)

    # Get the number of points in the polydata
    numCellsInPolyData = polydataPathline.GetNumberOfCells()
    # Get Scalars of Cell data
    arrayScalar = polydataPathline.GetCellData().GetScalars()

    # Iterlate each polyline.
    for i in range(numCellsInPolyData):
        polyLine = polydataPathline.GetCell(i)
        # Retrieve each line and data
        polyLine = polydataPathline.GetCell(i)  # Get a cell (PolyLine)
        scalar = arrayScalar.GetValue(i)  # Get Scalar Value

        # Iterate each line segment
        numIds = polyLine.GetPointIds().GetNumberOfIds()
        pt = None
        stateCross = "init"
        for j in range(numIds):
            ptId = polyLine.GetPointIds().GetId(j)
            pt_prev = pt
            pt = polydataPathline.GetPoint(ptId)

            # Check whether the line segment intersects exit plane.
            #     - If intersected, check inward or outward. Keep the last determination.
            #     - If not intersected, just skip.
            # At last, the status will be three of below.
            #     - "init"    : Retained inflow
            #     - "outward" : Direct flow
            #     - "inward"  : Retained inflow, or error
            if (j > 0):
                tolerance = 0.001
                t = vtk.mutable(
                    0.0)  # Parametric coordinate of intersection (0 (corresponding to p1) to 1 (corresponding to p2))
                x = [0.0] * 3  # The coordinate of the intersection
                pcoords = [0.0] * 3
                subId = vtk.mutable(0)
                # Check whether the line segment intersects exit plane.
                iD = polygonExit.IntersectWithLine(pt, pt_prev, tolerance, t, x, pcoords, subId)
                # If there is an intersection between the given plane and the line segment
                if (iD == 1):
                    vectorLineSegment = np.array(pt) - np.array(pt_prev)
                    innerP = np.inner(vectorLineSegment, vectorExitNormal)
                    stateCross = "outward" if innerP > 0 else "inward"
                    print(
                        f"Intersection line {i} segment {j} at {x}, {stateCross}: {innerP} = Inner ({vectorLineSegment}, {vectorExitNormal})")

        if (stateCross == "inward"):
            cellsOtherFlow.InsertNextCell(polyLine)
            scalarsOtherFlow.InsertNextTypedTuple([scalar])
        elif (stateCross == "outward"):
            cellsDirectFlow.InsertNextCell(polyLine)
            scalarsDirectFlow.InsertNextTypedTuple([scalar])
        else:  # if stateCross == "init"
            cellsRetainedInflow.InsertNextCell(polyLine)
            scalarsRetainedInflow.InsertNextTypedTuple([scalar])

    # Set PolyLines (pathlines) to each polyData
    polyDataDirectFlow.SetLines(cellsDirectFlow)
    polyDataRetainedInflow.SetLines(cellsRetainedInflow)
    polyDataOtherFlow.SetLines(cellsOtherFlow)

    # Store flux scalars array with pathlines
    polyDataDirectFlow.GetCellData().SetScalars(scalarsDirectFlow)
    polyDataRetainedInflow.GetCellData().SetScalars(scalarsRetainedInflow)
    polyDataOtherFlow.GetCellData().SetScalars(scalarsOtherFlow)

    # return pointsInnerLines
    return [polyDataDirectFlow, polyDataRetainedInflow, polyDataOtherFlow]


def GetPolygonFromPlaneSource(planeSource):
    ## Prepare rectangular polygon from the exit plane.
    # Create a square in the x-y plane.
    thirdPoint = np.array(planeSource.GetPoint1()) + np.array(planeSource.GetPoint2()) - np.array(
        planeSource.GetOrigin())
    points = vtk.vtkPoints()
    points.InsertNextPoint(planeSource.GetOrigin())
    points.InsertNextPoint(planeSource.GetPoint1())
    points.InsertNextPoint(thirdPoint)
    points.InsertNextPoint(planeSource.GetPoint2())

    # Create the polygon with four corner points
    polygon = vtk.vtkPolygon()
    polygon.GetPoints().DeepCopy(points)
    polygon.GetPointIds().SetNumberOfIds(4)  # The 4 corners of the square
    for i in range(4):
        polygon.GetPointIds().SetId(i, i)

    return polygon


def GetMiscPointData(polydata):
    # Get the number of points in the polydata
    idNumPointsInFile = polydata.GetNumberOfPoints()

    numArray = polydata.GetPointData().GetNumberOfArrays()
    for i in range(numArray):
        array = polydata.GetPointData().GetAbstractArray(i)
        if (array == None):
            print(f"##### array {i} returns None. #####")
        else:
            f = open(f'array{i}.txt', 'w')
            print(f"##### array {i}: {array.GetName()} with {idNumPointsInFile} #####")
            print(f"##### array {i}: {array.GetName()} with {idNumPointsInFile} #####", file=f)

            for j in range(idNumPointsInFile):
                value = array.GetValue(j)
                # print(f"    {j}: {value}")
                print(f"    {j}: {value}", file=f)

            f.close()


def GetMiscCellData(polydata):
    # Get the number of points in the polydata
    idNumCellsInFile = polydata.GetNumberOfCells()

    celldata = polydata.GetCellData()
    print(f"### Number of Cells {idNumCellsInFile}")
    print(celldata)

    for i in range(idNumCellsInFile):
        cell = polydata.GetCell(i)
        print(f"line {i}: {cell.GetNumberOfPoints()} points")


## GetResliceTransformMatrix
#  vtkImageReslice requires transform matrix.
#  This function generates transform matrix from planeSource and center point.
def GetResliceTransformMatrix(planeSource):
    # vtkImageReslice requires transform matrix
    #     |  vec0_x   vec1_x   vec2_x   center_x  |
    #     |  vec0_y   vec1_y   vec2_y   center_y  |
    #     |  vec0_z   vec1_z   vec2_z   center_z  |
    #     |       0        0        0          1  |

    ptOrigin = np.array(planeSource.GetOrigin())
    ptPoint1 = np.array(planeSource.GetPoint1())
    ptPoint2 = np.array(planeSource.GetPoint2())
    vectorOP1 = ptPoint1 - ptOrigin
    vectorOP2 = ptPoint2 - ptOrigin
    vectorOP1norm = vectorOP1 / np.linalg.norm(vectorOP1)
    vectorOP2norm = vectorOP2 / np.linalg.norm(vectorOP2)
    crossVectors = np.cross(vectorOP1norm, vectorOP2norm)
    # crossVectors = planeSource.GetNormal()    # or method of vtkPlaneSource can be used.
    arrayOblique = np.vstack((np.column_stack((vectorOP1norm, vectorOP2norm, crossVectors, ptOrigin)), [0, 0, 0, 1]))
    listOblique = arrayOblique.flatten().tolist()
    obliqueMatrix = vtk.vtkMatrix4x4()
    obliqueMatrix.DeepCopy(listOblique)
    # print(f"ptOrigin {ptOrigin}, ptPoint1 {ptPoint1}, ptPoint2 {ptPoint2}")
    # print(f"vectorOP1 {vectorOP1}, vectorOP2 {vectorOP2}")
    # print(f"vectorOP1norm {vectorOP1norm}, vectorOP2norm {vectorOP2norm}")
    # print(f"dotVectors {dotVectors}, crossVectors {crossVectors}")
    # print(f"GetNormal {planeSource.GetNormal()}, crossVectors {crossVectors}")
    # print(f"arrayOblique {arrayOblique}")
    # print(f"listOblique {listOblique}")
    # print(f"obliqueMatrix {obliqueMatrix}")

    return obliqueMatrix


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)


## GetParticlePathTrackingOutput: Tracking particle and return pathline output
#  Trace Particle in vector field using vtkParticlePathFilter
#
#  @param source4DFlow     Source 4D flow data in MrStruct (MrStruct4DFlowSource).
#  @param sourceEmission   Source locations to emit particles
#  @param bInverseTracking Forward tracking (False) or Backward tracking (True). This inverts 4D flow vector field direction and time direction.
#  @param frameStartOffset Frame number for tracking to start from.
#  @param frameEndOffset   If this is zero, track for full cycle. If this is negative number, tracking will be ended at specified earlier frame.
def GetParticlePathTrackingOutput(source4DFlow, sourceEmission, bInverseTracking, frameStartOffset, frameEndOffset):
    # There are three particle tracers for unsteady vector fields: vtkParticleTracer, vtkParticlePathFilter, vtkStreaklineFilter
    # Pathline: vtkParticlePathFilter

    # Create pathline filter
    pathfilter = vtk.vtkParticlePathFilter()

    # Setup sources
    source4DFlow.SetTimeSpaceInverse(bInverseTracking)
    source4DFlow.SetTimeOffset(frameStartOffset)
    pathfilter.SetInputConnection(source4DFlow.GetOutputPort())
    pathfilter.AddSourceConnection(sourceEmission.GetOutputPort())  # Forward tracing from the entrance

    # Setup Integrator
    # rk45 = vtk.vtkRungeKutta45()    # 5th order Runge-Kutta method with adaptive stepsize control
    # pathfilter.SetIntegrator(rk45)
    rk2 = vtk.vtkRungeKutta2()  # 2nd order Runge-Kutta method
    pathfilter.SetIntegrator(rk2)

    # Setup misc. settings
    pathfilter.SetRotationScale(0.5)
    terminationTime = source4DFlow.GetTimeStep() * (source4DFlow.GetNumberOfTimeFrames() + frameEndOffset - 1)
    pathfilter.SetTerminationTime(terminationTime)
    print(f"terminationTime {terminationTime}")

    # Run Pathline filter
    pathfilter.Update()
    pathOutput = pathfilter.GetOutput()  # output is vtkPolyData

    # Particles are traced from offset frame for termination time
    # offset frame: s.SetTimeOffset(9)
    # termination time: pathfilter.SetTerminationTime(terminationTime)

    # Output format
    # vtkPolyData
    #   Point Data:
    #     Number Of Arrays: 13
    #     Array 0 name = Magnitude
    #     Array 1 name = VelocityVector
    #     Array 2 name = SimulationTime
    #     Array 3 name = SimulationTimeStep
    #     Array 4 name = ParticleId
    #     Array 5 name = ParticleSourceId
    #     Array 6 name = InjectedPointId
    #     Array 7 name = InjectionStepId
    #     Array 8 name = ErrorCode
    #     Array 9 name = ParticleAge
    #     Array 10 name = Vorticity
    #     Array 11 name = Rotation
    #     Array 12 name = AngularVelocity
    #     Number Of Components: 15
    #     Number Of Tuples: 7490
    #   Number Of Points: 7490
    #   Number Of Lines: 441

    # GetMiscPointData(pathOutput)
    # GetMiscCellData(pathOutput)

    return pathOutput




## GetParticlePathTrackingOutput: Tracking particle and return pathline output
#  Trace Particle in vector field using vtkParticlePathFilter
#
#  @param source4DFlow     Source 4D flow data in MrStruct (MrStruct4DFlowSource).
#  @param sourcesEmission  Array of Source locations to emit particles
#  @param bInverseTracking Forward tracking (False) or Backward tracking (True). This inverts 4D flow vector field direction and time direction.
#  @param frameStartOffset Frame number for tracking to start from.
#  @param frameEndOffset   If this is zero, track for full cycle. If this is negative number, tracking will be ended at specified earlier frame.
def GetParticlePathTrackingOutput2(source4DFlow, sourcesEmission, bInverseTracking, frameStartOffset, frameEndOffset):
    # There are three particle tracers for unsteady vector fields: vtkParticleTracer, vtkParticlePathFilter, vtkStreaklineFilter
    # Pathline: vtkParticlePathFilter

    # Append source emission planes
    appendFilter = vtk.vtkAppendPolyData()
    appendFilter.AddInputConnection(sourcesEmission[0].GetOutputPort())
    appendFilter.AddInputConnection(sourcesEmission[1].GetOutputPort())

    # Remove any duplicate points.
    cleanFilter = vtk.vtkCleanPolyData()
    cleanFilter.SetInputConnection(appendFilter.GetOutputPort())
    cleanFilter.Update()

    return GetParticlePathTrackingOutput(source4DFlow, cleanFilter, bInverseTracking, frameStartOffset, frameEndOffset)
    # # Create pathline filter
    # pathfilter = vtk.vtkParticlePathFilter()
    #
    # # Setup sources
    # source4DFlow.SetTimeSpaceInverse(bInverseTracking)
    # source4DFlow.SetTimeOffset(frameStartOffset)
    # pathfilter.SetInputConnection(source4DFlow.GetOutputPort())
    # pathfilter.AddSourceConnection(cleanFilter.GetOutputPort())  # Forward tracing from the entrance
    #
    # # Setup Integrator
    # # rk45 = vtk.vtkRungeKutta45()    # 5th order Runge-Kutta method with adaptive stepsize control
    # # pathfilter.SetIntegrator(rk45)
    # rk2 = vtk.vtkRungeKutta2()  # 2nd order Runge-Kutta method
    # pathfilter.SetIntegrator(rk2)
    #
    # # Setup misc. settings
    # pathfilter.SetRotationScale(0.5)
    # terminationTime = source4DFlow.GetTimeStep() * (source4DFlow.GetNumberOfTimeFrames() + frameEndOffset - 1)
    # pathfilter.SetTerminationTime(terminationTime)
    # print(f"terminationTime {terminationTime}")
    #
    # # Run Pathline filter
    # pathfilter.Update()
    # pathOutput = pathfilter.GetOutput()  # output is vtkPolyData
    #
    # return pathOutput




###############################################################################
# Unused Functions
###############################################################################

# Calculate the center of the volume
def GetCenterOfImageData(imgData):
    (xMin, xMax, yMin, yMax, zMin, zMax) = imgData.GetExtent()
    (xSpacing, ySpacing, zSpacing) = imgData.GetSpacing()
    (x0, y0, z0) = imgData.GetOrigin()
    center = [x0 + xSpacing * 0.5 * (xMin + xMax),
              y0 + ySpacing * 0.5 * (yMin + yMax),
              z0 + zSpacing * 0.5 * (zMin + zMax)]
    # print(f"center = {center}")
    return center


if __name__ == "__main__":
    main()
