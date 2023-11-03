import vtk


def main():
    points = vtk.vtkPoints()
    points.InsertPoint(0,1,0,0)
    points.InsertPoint(1,2,0,0)
    points.InsertPoint(2,3,1,0)
    points.InsertPoint(3,4,1,0)
    points.InsertPoint(4,5,0,0)
    points.InsertPoint(5,6,0,0)

    # Fit a spline to the points
    spline = vtk.vtkParametricSpline()
    spline.SetPoints(points)
    functionSource = vtk.vtkParametricFunctionSource()
    functionSource.SetParametricFunction(spline)
    functionSource.SetUResolution(10 * points.GetNumberOfPoints())
    functionSource.Update()

    # Interpolate the scalars
    interpolatedRadius = vtk.vtkTupleInterpolator()
    interpolatedRadius.SetInterpolationTypeToLinear()
    interpolatedRadius.SetNumberOfComponents(1)
    interpolatedRadius.AddTuple(0,[0.2])
    interpolatedRadius.AddTuple(1,[0.2])
    interpolatedRadius.AddTuple(2,[0.2])
    interpolatedRadius.AddTuple(3,[0.1])
    interpolatedRadius.AddTuple(4,[0.1])
    interpolatedRadius.AddTuple(5,[0.1])

    # Generate the radius scalars
    tubeRadius = vtk.vtkDoubleArray()
    n = functionSource.GetOutput().GetNumberOfPoints()
    tubeRadius.SetNumberOfTuples(n)
    tubeRadius.SetName("TubeRadius")
    tMin = interpolatedRadius.GetMinimumT()
    tMax = interpolatedRadius.GetMaximumT()

    r = [0.0]
    for i in range(n):
        t = (tMax - tMin) / (n - 1) * i + tMin
        interpolatedRadius.InterpolateTuple(t, r)
        print(f"interpolated r {r}")
        tubeRadius.SetTuple1(i, r[0])

    # Add the scalars to the polydata
    tubePolyData = functionSource.GetOutput()
    tubePolyData.GetPointData().AddArray(tubeRadius)
    tubePolyData.GetPointData().SetActiveScalars("TubeRadius")

    # Create the tubes
    tuber = vtk.vtkTubeFilter()
    tuber.SetInputData(tubePolyData)
    tuber.SetNumberOfSides(20)
    tuber.SetVaryRadiusToVaryRadiusByAbsoluteScalar()

    #--------------
    # Setup actors and mappers
    lineMapper = vtk.vtkPolyDataMapper()
    lineMapper.SetInputData(tubePolyData)
    lineMapper.SetScalarRange(tubePolyData.GetScalarRange())

    tubeMapper = vtk.vtkPolyDataMapper()
    tubeMapper.SetInputConnection(tuber.GetOutputPort())
    tubeMapper.SetScalarRange(tubePolyData.GetScalarRange())

    lineActor = vtk.vtkActor()
    lineActor.SetMapper(lineMapper)
    lineActor.GetProperty().SetLineWidth(3)

    tubeActor = vtk.vtkActor()
    tubeActor.SetMapper(tubeMapper)
    tubeActor.GetProperty().SetOpacity(.6)

    # Setup render window, renderer, and interactor
    colors = vtk.vtkNamedColors()
    renderer = vtk.vtkRenderer()
    renderer.UseHiddenLineRemovalOn()
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindow.SetSize(640, 480)

    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)
    renderer.AddActor(lineActor)
    renderer.AddActor(tubeActor)
    renderer.SetBackground(colors.GetColor3d("SlateGray"))

    renderWindow.Render()
    renderWindowInteractor.Start()




if __name__ == '__main__':
    main()
