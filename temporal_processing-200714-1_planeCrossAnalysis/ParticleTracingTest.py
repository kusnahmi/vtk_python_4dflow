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


## main function
#
#  Main flow of particle tracing and visualization
def main():
    t = 0

    s = mrstruct.MrStruct4DFlowSource()
    s.SetSize([50, 30, 40])
    s.SetCenter([10, 15, 20])
    s.SetRadius(10)
    s.SetTimeValue(t)
    s.SetTimeOffset(9)
    s.SetFilepathMagnitude(
        "C:/d/project/Case_test/001A00505620161214/3dpc/mrstruct_subject_20161214_user011/mag_struct.mat")
    s.SetFilepathVelocity(
        "C:/d/project/Case_test/001A00505620161214/3dpc/mrstruct_subject_20161214_user011/vel_struct.mat")
    s.SetFilepathMask(
        "C:/d/project/Case_test/001A00505620161214/3dpc/mrstruct_subject_20161214_user011/LA_LV_volume.mat")

    s.SetTimeSpaceInverse(False)
    # s.SetTimeSpaceInverse(True)

    s.Update()
    print(s.GetMaskArray().shape)

    maskPolyData, maskActor = vtkViz.GetMaskActor(s.GetMaskImgData())
    maskActor.GetProperty().SetOpacity(0.2)

    stat = vtk.vtkImageHistogramStatistics()
    stat.SetInputConnection(s.GetOutputPort())
    stat.Update()

    print(f"time: {t} / mean: {stat.GetMean()}")

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

    ## Source Plane Definition
    #  Mitral Valve
    planeSourceMV = vtk.vtkPlaneSource()

    planeSourceMV.SetOrigin((155.75200746815736, 133.4443205456904, 57.55322438879402))
    planeSourceMV.SetPoint1((195.11049180338034, 162.2567972045777, 52.920065217945236))
    planeSourceMV.SetPoint2((158.37555320991834, 137.10930329145776, 102.63175635712525))

    planeSourceMV.SetXResolution(20)
    planeSourceMV.SetYResolution(20)

    #  Aortic Valve
    planeSourceAV = vtk.vtkPlaneSource()

    planeSourceAV.SetOrigin((187.59152548881804, 155.1431218693424, 60.117124070509334))
    planeSourceAV.SetPoint1((181.91575128426845, 122.62066865580515, 47.517833373914726))
    planeSourceAV.SetPoint2((167.03837674355384, 148.56220588247717, 86.36324948493984))

    planeSourceAV.SetXResolution(20)
    planeSourceAV.SetYResolution(20)

    # ## Streamline
    # # for i in [9]:
    # for direction in [vtk.vtkStreamTracer().BOTH, vtk.vtkStreamTracer().FORWARD, vtk.vtkStreamTracer().BACKWARD]:
    #     i = 9
    #     imgData = s.GetImageDataAtTimeindex(i)
    #     # showVectorData(imgData, None)
    #
    #     # Render Streamline
    #     rk = vtk.vtkRungeKutta45()
    #     # Create source for streamtubes
    #     streamer = vtk.vtkStreamTracer()
    #     streamer.SetInputData(imgData)
    #     streamer.SetSourceConnection(planeSourceMV.GetOutputPort())
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
    #     vtkViz.ShowActors([streamActor, maskActor])

    ## Pathline: vtkParticlePathFilter
    # Setup Pathline filter
    rk45 = vtk.vtkRungeKutta45()
    rk2 = vtk.vtkRungeKutta2()
    # Create source for streamtubes
    pathfilter = vtk.vtkParticlePathFilter()
    # pathfilter = vtk.vtkStreaklineFilter()
    pathfilter.SetInputConnection(s.GetOutputPort())

    pathfilter.AddSourceConnection(planeSourceMV.GetOutputPort())
    s.SetTimeSpaceInverse(False)
    # pathfilter.AddSourceConnection(planeSourceAV.GetOutputPort())
    # s.SetTimeSpaceInverse(True)

    # pathfilter.SetIntegrator(rk45)
    pathfilter.SetIntegrator(rk2)
    pathfilter.SetRotationScale(0.5)
    terminationTime = s.GetTimeStep() * s.GetNumberOfTimeFrames()
    pathfilter.SetTerminationTime(1000)

    # Run Pathline filter
    pathfilter.Update()
    pathOutput = pathfilter.GetOutput()
    # print(pathOutput)
    # GetMiscPointData(pathOutput)
    # GetMiscCellData(pathOutput)

    # Visualization: Get pathline actor
    mapStream = vtk.vtkPolyDataMapper()
    mapStream.SetInputConnection(pathfilter.GetOutputPort())
    # imgData = s.GetOutputDataObject(0)
    imgData = s.GetImageDataAtTimeindex(9)
    mapStream.SetScalarRange(imgData.GetScalarRange())
    pathlineActor = vtk.vtkActor()
    pathlineActor.SetMapper(mapStream)

    ## Output analysis: Inside or outside the mask?
    #  Get final points and determine if they are inside or outside the mask.
    finalPoints = ExtractFinalPointsFromPathlines(pathOutput)
    # print("### Render Particle Path Filter Result ###")
    # print(finalPoints)

    insidePoints, outsidePoints = DetermineInsideMask(finalPoints, maskPolyData)
    # ShowPoints(insidePoints, [255,255,0])
    # ShowPoints(outsidePoints, [0,0,255])

    ReportFinalPoints(insidePoints, outsidePoints)

    ## Output analysis: When does it cross the plane?
    AnalysisCrossPlane(pathOutput, planeSourceAV)

    # Visualization
    insidePointActor = vtkViz.GetPointsActor(insidePoints, [255, 255, 0])
    outsidePointActor = vtkViz.GetPointsActor(outsidePoints, [0, 0, 255])

    vtkViz.ShowVectorData(imgData, pathlineActor, maskActor, insidePointActor, outsidePointActor)
    vtkViz.ShowActors([maskActor, insidePointActor, outsidePointActor])


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

    # Create the polygon
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
                    print(f"Intersection line {i} segment {j} at {x}, {strInward}")

                    if (strInward == "outward"):
                        if ((prevState == "init") or (prevState == "inward")):
                            # Before cross point: InnerLines
                            polyLine = vtk.vtkPolyLine()
                            polyLine.GetPointIds().SetNumberOfIds(j - prevCrossPtIndex)
                            for k in range(prevCrossPtIndex, j):
                                # polyDataInnerLines.InsertNextPoint(polydata.GetPoint(cell.GetPointIds().GetId(k)))
                                newPtId = pointsInnerLines.InsertNextPoint(
                                    polydata.GetPoint(cell.GetPointIds().GetId(k)))
                                polyLine.GetPointIds().SetId(k-prevCrossPtIndex, newPtId)
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
                                polyLine.GetPointIds().SetId(k-prevCrossPtIndex, newPtId)
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
                polyLine.GetPointIds().SetId(k-prevCrossPtIndex, newPtId)
            cellsInnerLines.InsertNextCell(polyLine)
        elif (prevState == "outward"):
            # After cross point: OuterLines
            polyLine = vtk.vtkPolyLine()
            polyLine.GetPointIds().SetNumberOfIds(numIds - prevCrossPtIndex)
            for k in range(prevCrossPtIndex, numIds):
                newPtId = pointsOuterLines.InsertNextPoint(
                    polydata.GetPoint(cell.GetPointIds().GetId(k)))
                polyLine.GetPointIds().SetId(k-prevCrossPtIndex, newPtId)
            cellsOuterLines.InsertNextCell(polyLine)
        else:   # if prevState == "init"
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


if __name__ == "__main__":
    main()
