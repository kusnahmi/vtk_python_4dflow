# XML parsing: refer to
#     https://www.geeksforgeeks.org/xml-parsing-python/
#     https://docs.python.org/3/library/xml.etree.elementtree.html

import vtk
import xml.etree.ElementTree as ET
import numpy as np


def main():
    filepath = "D:\\data_clinical\\RASTA\\RASTAX003\\Raw\\Velocity_Export\\study_dwsG700471_2020_10_28_08_38_34\\2020_10_28_13_39_31\\dif001.xml"
    # filepath = "D:\\data_clinical\\RASTA\\RASTAX003\\Raw\\Velocity_Export\\study_dwsG700471_2020_10_28_08_38_34\\2020_10_28_13_39_31\\ModelGroups.xml"

    print(filepath)

    # parse xml file to Volumes
    volumes = parseXML(filepath)
    print(volumes)

    # create a rendering window and renderer
    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)

    # create a renderwindowinteractor
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)

    for volume in volumes:
        # mapper
        mapper = vtk.vtkPolyDataMapper()
        # mapper.SetInputConnection(source.GetOutputPort())
        mapper.SetInputData(volume)

        # actor
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

        # assign actor to the renderer
        ren.AddActor(actor)

    # enable user interface interactor
    iren.Initialize()
    renWin.Render()
    iren.Start()

    for volume in volumes:
        # create a rendering window and renderer
        ren = vtk.vtkRenderer()
        renWin = vtk.vtkRenderWindow()
        renWin.AddRenderer(ren)

        # create a renderwindowinteractor
        iren = vtk.vtkRenderWindowInteractor()
        iren.SetRenderWindow(renWin)

        # mapper
        mapper = vtk.vtkPolyDataMapper()
        # mapper.SetInputConnection(source.GetOutputPort())
        mapper.SetInputData(volume)

        # actor
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

        # assign actor to the renderer
        ren.AddActor(actor)

        # enable user interface interactor
        iren.Initialize()
        renWin.Render()
        iren.Start()


def parseXML(xmlfile):
    # create element tree object
    tree = ET.parse(xmlfile)

    # get root element
    root = tree.getroot()

    # prepare polydata list
    volumes = []

    print("iter Volume")
    for volume in root.iter('Volume'):
        print(volume.attrib)
        vertices = volume.find('Vertices')
        normals = volume.find('Normals')
        polygons = volume.find('Polygons')

        if (vertices == None) or (normals == None) or (polygons == None):
            continue

        # Convert vertices in XML to vtkPoints
        print(vertices.attrib)
        points = vtk.vtkPoints()
        for line in vertices.text.splitlines():
            elementsStr = line.split()
            if (len(elementsStr) != 3):
                continue
            elements = np.asfarray(elementsStr, float)
            points.InsertNextPoint(elements)
        print(f"Vertices Parse Number: {points.GetNumberOfPoints()}")

        # Convert polygons in XML to vtk list of polygons (vtkCellArray)
        polygonsVtk = vtk.vtkCellArray()
        print(polygons.attrib)
        for line in polygons.text.splitlines():
            elementsStr = line.split()
            if (len(elementsStr) != 3):
                continue
            elements = np.asarray(elementsStr, dtype=int)
            elements = elements - 1
            # Add the polygon to a list of polygons
            polygon = vtk.vtkPolygon()
            polygon.GetPointIds().SetNumberOfIds(3)
            polygon.GetPointIds().SetId(0, elements[0])
            polygon.GetPointIds().SetId(1, elements[1])
            polygon.GetPointIds().SetId(2, elements[2])

            # Add the polygon to a list of polygons
            polygonsVtk.InsertNextCell(polygon)
        print(f"Polygons Parse Number: {polygonsVtk.GetNumberOfCells()}")

        # Create a PolyData
        polygonPolyData = vtk.vtkPolyData()
        polygonPolyData.SetPoints(points)
        polygonPolyData.SetPolys(polygonsVtk)

        volumes.append(polygonPolyData)
        # print(f"{points.GetPoint(0)}")
        # print(f"{points.GetPoint(1)}")
        # print(normals.attrib)
        # print(vertices.text)

    print("END: iter Volume")

    return volumes


if __name__ == "__main__":
    main()
