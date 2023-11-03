import vtk


def main():
    points = vtk.vtkPoints()
    points.InsertPoint(0,1,0,0)
    points.InsertPoint(1,2,0,0)
    points.InsertPoint(2,3,1,0)
    points.InsertPoint(3,4,1,0)
    points.InsertPoint(4,5,0,0)
    points.InsertPoint(5,6,0,0)
    points.InsertPoint(6,6,0,1)
    points.InsertPoint(7,6,0,2)

    # Generate the radius scalars
    tubeRadius = vtk.vtkDoubleArray()
    n = points.GetNumberOfPoints()
    tubeRadius.SetNumberOfTuples(n)
    tubeRadius.SetName("TubeRadius")
    # tubeRadius.SetTuple1(0, [0.1])
    # tubeRadius.SetTuple1(1, [0.2])
    # tubeRadius.SetTuple1(2, [0.3])
    # tubeRadius.SetTuple1(3, [0.4])
    # tubeRadius.SetTuple1(4, [0.5])
    # tubeRadius.SetTuple1(5, [0.6])
    tubeRadius.SetTuple1(0, 0.1)
    tubeRadius.SetTuple1(1, 0.2)
    tubeRadius.SetTuple1(2, 0.3)
    tubeRadius.SetTuple1(3, 0.4)
    tubeRadius.SetTuple1(4, 0.5)
    tubeRadius.SetTuple1(5, 0.6)

    # Add the scalars to the polydata
    cellsPolyLines = vtk.vtkCellArray()
    polyLine1 = vtk.vtkPolyLine()
    polyLine2 = vtk.vtkPolyLine()
    polyLine1.GetPointIds().SetNumberOfIds(6)
    polyLine1.GetPointIds().SetId(0,0)
    polyLine1.GetPointIds().SetId(1,1)
    polyLine1.GetPointIds().SetId(2,2)
    polyLine1.GetPointIds().SetId(3,3)
    polyLine1.GetPointIds().SetId(4,4)
    polyLine1.GetPointIds().SetId(5,5)
    polyLine2.GetPointIds().SetNumberOfIds(4)
    polyLine2.GetPointIds().SetId(0,0)
    polyLine2.GetPointIds().SetId(1,2)
    polyLine2.GetPointIds().SetId(2,6)
    polyLine2.GetPointIds().SetId(3,7)
    cellsPolyLines.InsertNextCell(polyLine1)
    cellsPolyLines.InsertNextCell(polyLine2)

    tubePolyData = vtk.vtkPolyData()
    tubePolyData.SetPoints(points)
    tubePolyData.SetLines(cellsPolyLines)

    # tubePolyData = functionSource.GetOutput()
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
