## @package vtkVisualization
#  VTK Visualization helper
#
#  Easy to use functions to make actors from data
import vtk


## Get actor from binary mask data.
#
#  Mask data is 3-dimensional and each element is 1 or 0.
#  @param maskImgData vtkImageData
def GetMaskActor(maskImgData):
    maskSurface = GetIsoSurface(maskImgData, False)

    colors = vtk.vtkNamedColors()
    renderer = vtk.vtkRenderer()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(maskSurface.GetOutput())
    mapper.SetScalarModeToUseCellData()

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetDiffuseColor(colors.GetColor3d("Magenta"))

    return maskSurface.GetOutput(), actor


## Extract surface poly data from the image
#
#  @param Input  imageData vtkImageAlgorithm
#  @param Input  bDecimate If true, simplify the surface after extract a surface.
#  @param Output vtkPolyDataAlgorithm
def GetIsoSurface(imageData, bDecimate):
    # Create the isosurface
    isosurface = vtk.vtkMarchingCubes()
    # isosurface.SetInputConnection(imageAlgorithm.GetOutputPort())
    isosurface.SetInputData(imageData)
    isosurface.SetValue(0, 1)
    isosurface.ComputeNormalsOn()
    isosurface.Update()

    # Use a decimator to reduce the number of triangles, if required
    if (bDecimate):
        decimator = vtk.vtkDecimatePro()
        decimator.SetInputConnection(isosurface.GetOutputPort())
        decimator.SetTargetReduction(0.5)
        decimator.Update()
        return decimator
    else:
        return isosurface


## Create a window and show actors
#
#  Create a window, add coordinate sphere, and add given actors.
#  The view will be adjusted appropriately. Interactor is in trackball mode.
#  @param actors actors to be displayed.
def ShowActors(actors, widgets=None):
    colors = vtk.vtkNamedColors()

    aRenderer = vtk.vtkRenderer()
    aRenderWindow = vtk.vtkRenderWindow()
    aRenderWindow.AddRenderer(aRenderer)
    aRenderWindow.SetSize(640, 480)
    aRenderer.SetBackground(colors.GetColor3d("Gray"))

    anInteractor = vtk.vtkRenderWindowInteractor()
    anInteractor.SetRenderWindow(aRenderWindow)

    # Actors
    sphereActor = GetCoordinateSphereActor()
    for i in range(5):
        aRenderer.AddActor(sphereActor[i])

    for argActor in actors:
        aRenderer.AddActor(argActor)

    # Widgets
    if (widgets != None):
        for argWidget in widgets:
            argWidget.SetInteractor(anInteractor)
            argWidget.EnabledOn()

    axesMarkerWidget = GetOrientationMarkerActor()
    axesMarkerWidget.SetInteractor(anInteractor)
    axesMarkerWidget.EnabledOn()
    axesMarkerWidget.InteractiveOn()

    # Generate an interesting view.

    aRenderer.GetActiveCamera().SetFocalPoint(0, 0, 0)
    aRenderer.GetActiveCamera().SetPosition(1, 0, 0)
    aRenderer.GetActiveCamera().SetViewUp(0, 0, 1)
    aRenderer.ResetCamera()

    aRenderer.GetActiveCamera().Azimuth(60)
    aRenderer.GetActiveCamera().Elevation(30)
    aRenderer.GetActiveCamera().Dolly(1.1)
    aRenderer.ResetCameraClippingRange()

    aRenderWindow.Render()

    # Set Trackball mode
    style = vtk.vtkInteractorStyleSwitch()
    style.SetCurrentStyleToTrackballCamera()
    anInteractor.SetInteractorStyle(style)

    # Interact with the data.
    anInteractor.Start()


## Create a window and show vector data with other actors
#
#  Create a window, add coordinate sphere, and add given actors.
#  Give vector data will be displayed arrows.
#  The view will be adjusted appropriately. Interactor is in trackball mode.
#  @param imageVectorData vector data to be displayed.
#  @param argv other actors to be displayed.
def ShowVectorData(imageVectorData, *argv):
    # hhogActor = GetHedgehogVectorActor(imageVectorData)
    arrowVectorActor = GetArrowVectorActor(imageVectorData)
    outlineActor = GetOutlineActor(imageVectorData)
    planeWidget = GetPlaneWidget(imageVectorData)
    ShowActors([outlineActor, arrowVectorActor, *argv], [planeWidget])
    # showActors([outlineActor, hhogActor, *argv], [planeWidget])

def GetOutlineActor(imageData):
    outline = vtk.vtkOutlineFilter()
    outline.SetInputData(imageData)

    outlineMapper = vtk.vtkPolyDataMapper()
    outlineMapper.SetInputConnection(outline.GetOutputPort())

    outlineActor = vtk.vtkActor()
    outlineActor.SetMapper(outlineMapper)
    outlineActor.GetProperty().SetColor(0, 0, 0)

    return outlineActor

def GetHedgehogVectorActor(imageVectorData):
    hhog = vtk.vtkHedgeHog()
    hhog.SetInputData(imageVectorData)
    # hhog.SetScaleFactor(0.2)
    hhog.SetScaleFactor(10)

    lut = vtk.vtkLookupTable()
    # lut.SetHueRange(.667, 0.0)
    lut.Build()

    hhogMapper = vtk.vtkPolyDataMapper()
    hhogMapper.SetInputConnection(hhog.GetOutputPort())
    # hhogMapper.SetScalarRange(50, 550)
    # valueRange = imageVectorData.GetPointData().GetArray("scalars").GetRange()
    valueRange = imageVectorData.GetScalarRange()
    print(f"scalar range {valueRange}")
    hhogMapper.SetScalarRange(valueRange[0], valueRange[1])
    hhogMapper.SetLookupTable(lut)

    hhogActor = vtk.vtkActor()
    hhogActor.SetMapper(hhogMapper)

    return hhogActor

def GetArrowVectorActor(imageVectorData):
    lut = vtk.vtkLookupTable()
    # lut.SetHueRange(.667, 0.0)
    lut.Build()

    # Source for the glyph filter
    arrow = vtk.vtkArrowSource()
    arrow.SetTipResolution(16)
    arrow.SetTipLength(0.3)
    arrow.SetTipRadius(0.1)

    glyph = vtk.vtkGlyph3D()

    glyph.SetSourceConnection(arrow.GetOutputPort())
    glyph.SetInputData(imageVectorData)
    glyph.SetVectorModeToUseVector()
    # glyph.SetScaleFactor(0.2)
    glyph.SetScaleFactor(10)
    glyph.SetColorModeToColorByScalar()
    glyph.SetScaleModeToScaleByVector()
    glyph.OrientOn()
    glyph.Update()

    glyphMapper = vtk.vtkPolyDataMapper()
    glyphMapper.SetInputConnection(glyph.GetOutputPort())
    # glyphMapper.SetScalarModeToUsePointFieldData()
    # glyphMapper.SelectColorArray('Elevation')
    glyphMapper.SetColorModeToMapScalars()
    glyphMapper.ScalarVisibilityOn()
    glyphMapper.SetScalarRange(imageVectorData.GetScalarRange())
    glyphMapper.SetLookupTable(lut)

    glyphActor = vtk.vtkActor()
    glyphActor.SetMapper(glyphMapper)

    return glyphActor

## Get vtkPlaneWidget with a given Image Data
#
#  @param imageData Source image to display on the plane widget
def GetPlaneWidget(imageData):
    # The plane widget is used probe the dataset.
    planeWidget = vtk.vtkPlaneWidget()
    planeWidget.SetInputData(imageData)
    planeWidget.NormalToXAxisOn()
    planeWidget.SetResolution(20)
    planeWidget.SetRepresentationToOutline()
    planeWidget.PlaceWidget()
    plane = vtk.vtkPolyData()
    planeWidget.GetPolyData(plane)

    # Handle the events.
    # planeWidget.AddObserver("EnableEvent", BeginInteraction)
    # planeWidget.AddObserver("StartInteractionEvent", BeginInteraction)
    planeWidget.AddObserver("InteractionEvent", cbProbeData)

    return planeWidget


## Call back for vtkPlaneWidget
#
#  When plane widget is modified, print out plane coordinates.
def cbProbeData(obj, event):
    print(obj.GetOrigin(), obj.GetPoint1(), obj.GetPoint2())




## Returns an actor with coordinate sphere
#
#  5 spheres: origin, x-axis, y-axis, z-axis, and diagonal point
#  @param location sphere location in float number
#  @param output   an array of sphere actors
def GetCoordinateSphereActor(location=50.0):
    colors = vtk.vtkNamedColors()

    sphereSource = [0] * 5
    sphereMapper = [0] * 5
    sphereActor = [0] * 5

    for i in range(5):
        # Create a sphere
        sphereSource[i] = vtk.vtkSphereSource()
        sphereSource[i].SetRadius(5)
        # Make the surface smooth.
        sphereSource[i].SetPhiResolution(100)
        sphereSource[i].SetThetaResolution(100)

        sphereMapper[i] = vtk.vtkPolyDataMapper()
        sphereMapper[i].SetInputConnection(sphereSource[i].GetOutputPort())

        sphereActor[i] = vtk.vtkActor()
        sphereActor[i].SetMapper(sphereMapper[i])
        sphereActor[i].GetProperty().SetColor(colors.GetColor3d("Cornsilk"))

    sphereSource[0].SetCenter(location, 0.0, 0.0)
    sphereSource[1].SetCenter(0.0, location, 0.0)
    sphereSource[2].SetCenter(0.0, 0.0, location)
    sphereSource[3].SetCenter(location, location, location)
    sphereSource[4].SetCenter(0.0, 0.0, 0.0)
    sphereActor[0].GetProperty().SetColor(colors.GetColor3d("Red"))
    sphereActor[1].GetProperty().SetColor(colors.GetColor3d("Green"))
    sphereActor[2].GetProperty().SetColor(colors.GetColor3d("Blue"))
    sphereActor[4].GetProperty().SetColor(colors.GetColor3d("Black"))

    return sphereActor


## Get actor with given points
#
#  @param points Array points. (vtkPoints)
#  @param color  RGB tuple representing color (ex. [255,0,0] for red)
def GetPointsActor(points, color):
    colors = vtk.vtkUnsignedCharArray()
    colors.SetName("colors")
    colors.SetNumberOfComponents(3)
    # r = [255,0,0]
    for i in range(points.GetNumberOfPoints()):
        # colors.InsertNextTuple(r)
        colors.InsertNextTuple(color)

    # Combine into a polydata
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.GetPointData().SetScalars(colors)

    # Create anything you want here, we will use a cube for the demo.
    cubeSource = vtk.vtkCubeSource()
    sphereSource = vtk.vtkSphereSource()

    glyph3D = vtk.vtkGlyph3D()
    glyph3D.SetColorModeToColorByScalar()
    # glyph3D.SetSourceConnection(cubeSource.GetOutputPort())
    glyph3D.SetSourceConnection(sphereSource.GetOutputPort())
    glyph3D.SetInputData(polydata)
    glyph3D.ScalingOff()
    glyph3D.Update()

    # Create a mapper and actor
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(glyph3D.GetOutputPort())
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    return actor


## Add vtkOrientationMarkerWidget
#
#  Show cube with direction marker on each face.
#  @param iren Interactor which the widget will be shown on.
def AddOrientationMarkerActor(iren):
    axesActor = vtk.vtkAnnotatedCubeActor()
    # axesActor.SetXPlusFaceText('R')
    # axesActor.SetXMinusFaceText('L')
    # axesActor.SetYMinusFaceText('H')
    # axesActor.SetYPlusFaceText('F')
    # axesActor.SetZMinusFaceText('P')
    # axesActor.SetZPlusFaceText('A')
    axesActor.SetXPlusFaceText('F')
    axesActor.SetXMinusFaceText('H')
    axesActor.SetYMinusFaceText('L')
    axesActor.SetYPlusFaceText('R')
    axesActor.SetZMinusFaceText('A')
    axesActor.SetZPlusFaceText('P')
    axesActor.GetTextEdgesProperty().SetColor(1, 1, 0)
    axesActor.GetTextEdgesProperty().SetLineWidth(2)
    axesActor.GetCubeProperty().SetColor(0, 0, 1)
    axes = vtk.vtkOrientationMarkerWidget()
    axes.SetOrientationMarker(axesActor)
    axes.SetInteractor(iren)
    axes.EnabledOn()
    axes.InteractiveOn()
    return axes


## Get vtkOrientationMarkerWidget
#
#  Get orientation marker widget with a direction marker on each face.
def GetOrientationMarkerActor():
    axesActor = vtk.vtkAnnotatedCubeActor()
    # axesActor.SetXPlusFaceText('R')
    # axesActor.SetXMinusFaceText('L')
    # axesActor.SetYMinusFaceText('H')
    # axesActor.SetYPlusFaceText('F')
    # axesActor.SetZMinusFaceText('P')
    # axesActor.SetZPlusFaceText('A')
    axesActor.SetXPlusFaceText('F')
    axesActor.SetXMinusFaceText('H')
    axesActor.SetYMinusFaceText('L')
    axesActor.SetYPlusFaceText('R')
    axesActor.SetZMinusFaceText('A')
    axesActor.SetZPlusFaceText('P')
    axesActor.GetTextEdgesProperty().SetColor(1, 1, 0)
    axesActor.GetTextEdgesProperty().SetLineWidth(2)
    axesActor.GetCubeProperty().SetColor(0, 0, 1)
    axes = vtk.vtkOrientationMarkerWidget()
    axes.SetOrientationMarker(axesActor)
    return axes



## Get actor from PolyData with specified color.
#
#  @param polyData PolyData to draw (ex. PolyLine)
#  @param colorName color name. vtkNamedColors will convert this name to color value.
def GetPolydataActor(polyData, colorName):
    colors = vtk.vtkNamedColors()

    # Setup actor and mapper
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(polyData)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    # print(colors.GetColor3d(colorName))
    # print(colors.GetColor3d(colorName).GetData())
    actor.GetProperty().SetColor(colors.GetColor3d(colorName))

    return actor


## Get actor from polygon with specified color.
#
#  @param polygon polygon to draw
#  @param colorName color name. vtkNamedColors will convert this name to color value.
def GetPolygonActor(polygon, colorName):
    # Add the polygon to a list of polygons
    polygons = vtk.vtkCellArray()
    polygons.InsertNextCell(polygon)

    # Create a PolyData
    polygonPolyData = vtk.vtkPolyData()
    polygonPolyData.SetPoints(polygon.GetPoints())
    polygonPolyData.SetPolys(polygons)

    return GetPolydataActor(polygonPolyData, colorName)
