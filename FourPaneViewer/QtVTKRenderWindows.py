# pyuic5 -x QtVTKRenderWindows.ui -o Ui_QtVTKRenderWindows.py

import sys
import vtk
from PyQt5.QtWidgets import QApplication, QMainWindow
from Ui_QtVTKRenderWindows import Ui_QtVTKRenderWindows
import numpy as np
import MrStruct
from vtk.numpy_interface import dataset_adapter as dsa


class QtVTKRenderWindows(QMainWindow, Ui_QtVTKRenderWindows):
    riw = [0] * 3  # ResliceImageViewer
    planeWidget = [0] * 3  # ImagePlaneWidget
    DistanceWidget = [0] * 3
    ResliceMeasurements = 0

    imageReslice = [0] *3 # ImageReslice (for overlay mask)

    def __init__(self, parent=None):

        super().__init__()
        self.setupUi(self)

        # # Read data source
        # reader = vtk.vtkDICOMImageReader()
        # reader.SetDirectoryName(sys.argv[1])
        # reader.Update()
        # imageDims = reader.GetOutput().GetDimensions()
        # imgData = reader.GetOutput()

        # Read MrStruct magnitude
        #filepath = "C:/d/project/Case_test/001A00505620161214/3dpc/mrstruct_subject_20161214_user011/mag_struct.mat"
        filepath = "D:/data_clinical/control_all/done/001B00849420171207/3dpc/mrstruct_subject_20171207_user011/mag_struct.mat"

        mrStruct_mag = MrStruct.MrStruct(filepath)
        data = mrStruct_mag.dataAy[:,:,:,9]
        imgData_Mag = vtk.vtkImageData()
        imgData_Mag.SetDimensions(data.shape)
        print(f"Data shape: {data.shape}")
        imgData_Mag.SetSpacing(mrStruct_mag.vox[0:3])
        # Make a VTK array from the numpy array (using pointers)
        scalarsArray = dsa.numpyTovtkDataArray(data.ravel(order='F'))
        scalarsArray.SetName("magnitude")
        imgData_Mag.GetPointData().SetScalars(scalarsArray)

        imgData = imgData_Mag
        imageDims = imgData_Mag.GetDimensions()

        # Read MrStruct magnitude
        #filepath = "C:/d/project/Case_test/001A00505620161214/3dpc/mrstruct_subject_20161214_user011/LA_LV_volume.mat"
        filepath = "D:/data_clinical/control_all/done/001B00849420171207/3dpc/mrstruct_subject_20171207_user011/LA_volume_2.mat"

        mrStruct_mask = MrStruct.MrStruct(filepath)
        data = mrStruct_mask.dataAy
        imgData_mask = vtk.vtkImageData()
        imgData_mask.SetDimensions(data.shape)
        print(f"Data shape: {data.shape}")
        imgData_mask.SetSpacing(mrStruct_mask.vox[0:3])
        # Make a VTK array from the numpy array (using pointers)
        scalarsArray = dsa.numpyTovtkDataArray(data.ravel(order='F'))
        scalarsArray.SetName("magnitude")
        imgData_mask.GetPointData().SetScalars(scalarsArray)

        # imgData = imgData_mask
        # imageDims = imgData_mask.GetDimensions()





        # Source Data info and stats
        print(f"extent: {imgData.GetExtent()}")

        imgDataStat = vtk.vtkImageHistogramStatistics()
        imgDataStat.SetInputData(imgData)
        imgDataStat.SetAutoRangePercentiles(10, 90)
        imgDataStat.Update()
        print(f"min: {imgDataStat.GetMinimum()}")
        print(f"max: {imgDataStat.GetMaximum()}")
        print(f"mean: {imgDataStat.GetMean()}")
        print(f"median: {imgDataStat.GetMedian()}")
        print(f"autoRange: {imgDataStat.GetAutoRange()}")




        # Setup ResliceImageViewer panes
        for i in range(3):
            self.riw[i] = vtk.vtkResliceImageViewer()

        self.riw[0].SetupInteractor(self.view1)
        self.riw[0].SetRenderWindow(self.view1.GetRenderWindow())

        self.riw[1].SetupInteractor(self.view2)
        self.riw[1].SetRenderWindow(self.view2.GetRenderWindow())

        self.riw[2].SetupInteractor(self.view3)
        self.riw[2].SetRenderWindow(self.view3.GetRenderWindow())


        print(f"ResliceMode: {self.riw[0].GetResliceMode()}")

        #################################################################################
        # Setup ResliceImageViewer panes: view1 - view3
        #################################################################################
        for i in range(3):
            # make them all share the same reslice cursor object.
            self.riw[i].SetResliceCursor(self.riw[0].GetResliceCursor())
            # vtkResliceCursorLineRepresentation *rep
            rep = self.riw[i].GetResliceCursorWidget().GetRepresentation()
            rep.GetResliceCursorActor().GetCursorAlgorithm().SetReslicePlaneNormal(i)

            self.riw[i].SetInputData(imgData)
            self.riw[i].SetColorLevel(np.mean(imgDataStat.GetAutoRange()))
            self.riw[i].SetColorWindow(np.diff(imgDataStat.GetAutoRange()))

            self.riw[i].SetSliceOrientation(i)
            self.riw[i].SetResliceModeToAxisAligned()

            print(f"GetSlice {i}: {self.riw[i].GetSlice()}")
            self.riw[i].SetSlice(int(imgData.GetExtent()[i * 2 + 1] / 2))
            print(f"SetSlice {i}: {int(imgData.GetExtent()[i * 2 + 1] / 2)}")
            self.riw[i].SetResliceModeToAxisAligned()

        for i in range(3):
            self.riw[i].Render()

        # Set background color
        for i in range(3):
            color = [0, 0, 0]
            color[i] = 1
            # Background color on view1 - view3.
            color[0] /= 4.0
            color[1] /= 4.0
            color[2] /= 4.0
            self.riw[i].GetRenderer().SetBackground(color)


        self.resliceMode(1)



        # #################################################################################
        # # Setup ImageReslice reslicers for overlay mask
        # #################################################################################
        # # create lookup table
        # lookupTable = vtk.vtkLookupTable()
        # lookupTable.SetNumberOfTableValues(100)
        # lookupTable.SetRange(0, 185)
        # lookupTable.SetHueRange(0.1,0.9)
        # lookupTable.SetSaturationRange(1,1)
        # lookupTable.SetAlphaRange(1,1)
        # lookupTable.Build()
        #
        # # Setup vtkImageReslice to reslice overlay mask
        # print("#########")
        # for i in range(3):
        #     print(f"### {i} ###")
        #     # Get slice axes
        #     resliceAxes = vtk.vtkMatrix4x4()
        #     resliceAxes.DeepCopy(self.riw[i].GetResliceCursorWidget().GetResliceCursorRepresentation().GetResliceAxes())
        #     # resliceAxes.DeepCopy(self.riw[(1-i)%3].GetResliceCursorWidget().GetResliceCursorRepresentation().GetResliceAxes())
        #     print(resliceAxes)
        #     # vtk.vtkMatrix4x4().Invert(self.riw[i].GetResliceCursorWidget().GetResliceCursorRepresentation().GetResliceAxes(), resliceAxes)
        #     # resliceAxes.Identity()
        #
        #     center = self.riw[0].GetResliceCursor().GetCenter()
        #     print(center)
        #
        #     print(resliceAxes)
        #     resliceAxes.SetElement(0,3,center[0])
        #     resliceAxes.SetElement(1,3,center[1])
        #     resliceAxes.SetElement(2,3,center[2])
        #     # resliceAxes.SetElement(0,3,-center[0])
        #     # resliceAxes.SetElement(1,3,-center[1])
        #     # resliceAxes.SetElement(2,3,-center[2])
        #     print(resliceAxes)
        #     # print(self.riw[i].GetResliceCursor())
        #
        #     self.imageReslice[i] = vtk.vtkImageReslice()
        #     # self.imageReslice[i].SetInputData(imgData_mask)
        #     self.imageReslice[i].SetInputData(imgData_Mag)
        #     scalarRange = self.imageReslice[i].GetInput().GetScalarRange()
        #     self.imageReslice[i].SetBackgroundLevel(scalarRange[0])
        #     self.imageReslice[i].TransformInputSamplingOff()
        #     self.imageReslice[i].AutoCropOutputOn()
        #     self.imageReslice[i].SetResliceAxes(resliceAxes)
        #     self.imageReslice[i].SetOutputSpacing(1,1,1)
        #     self.imageReslice[i].SetOutputOrigin(0,0,0)
        #     self.imageReslice[i].SetOutputExtent(-500,500,-500,500,0,0)
        #
        #     print(f"Scalar Range: {scalarRange}")
        #     self.imageReslice[i].SetOutputDimensionality(2)
        #     # self.imageReslice[i].SetInterpolationModeToLinear()
        #
        #     self.imageReslice[i].Update()
        #
        #     # vtkResliceCursor.vtkPlane.
        #     plane = self.riw[0].GetResliceCursor().GetPlane(i)
        #     print(plane)
        #
        #     # self.riw[i].SetResliceCursor(self.riw[0].GetResliceCursor())
        #     # self.riw[i].GetResliceCursorWidget()
        #
        #     # map colours
        #     mapTransparency2 = vtk.vtkImageMapToColors()
        #     mapTransparency2.SetLookupTable(lookupTable)
        #     mapTransparency2.PassAlphaToOutputOn()
        #     mapTransparency2.SetInputConnection(self.imageReslice[i].GetOutputPort())
        #
        #     # setup image actor for slice
        #     maskActor = vtk.vtkImageActor()
        #     maskActor.GetMapper().SetInputConnection(mapTransparency2.GetOutputPort())
        #
        #     # double pos[3]
        #     # maskActor.GetPosition(pos)
        #     # pos[2] = 55
        #     # maskActor.SetPosition(pos)
        #     pos = maskActor.GetPosition()
        #     pos_list = list(pos)
        #     pos_list[2] = 300
        #     maskActor.SetPosition(pos_list)
        #     maskActor.SetOpacity(0.5)
        #     print(f"maskActor pos {pos}, new_pos {pos_list}")
        #
        #     # Preview Reslice Mask
        #     rep = self.riw[i].GetResliceCursorWidget().GetRepresentation()
        #     # rep.GetResliceCursorActor().GetCursorAlgorithm().SetReslicePlaneNormal(i)
        #     cursorActor = rep.GetResliceCursorActor()
        #
        #
        #     renderer = vtk.vtkRenderer()
        #     renderer.AddActor(maskActor)
        #     renderer.AddActor(cursorActor)
        #     window = vtk.vtkRenderWindow()
        #     window.AddRenderer(renderer)
        #     interactorStyle = vtk.vtkInteractorStyleImage()
        #     interactor = vtk.vtkRenderWindowInteractor()
        #     interactor.SetInteractorStyle(interactorStyle)
        #     window.SetInteractor(interactor)
        #     window.Render()
        #     interactor.Start()
        #
        #
        #
        #
        #
        #     self.riw[i].GetRenderer().AddActor(maskActor)
        #     self.riw[i].GetRenderer().ResetCamera()

        #################################################################################
        # Setup the ImagePlaneWidgets pane: view4
        #################################################################################
        picker = vtk.vtkCellPicker()
        picker.SetTolerance(0.005)

        ipwProp = vtk.vtkProperty()

        ren = vtk.vtkRenderer()
        self.view4.GetRenderWindow().AddRenderer(ren)

        # Set Trackball mode
        style = vtk.vtkInteractorStyleSwitch()
        style.SetCurrentStyleToTrackballCamera()
        self.view4.SetInteractorStyle(style)

        for i in range(3):
            # Setup ImagePlaneWidget on view4
            self.planeWidget[i] = vtk.vtkImagePlaneWidget()
            self.planeWidget[i].SetInteractor(self.view4)
            self.planeWidget[i].SetPicker(picker)
            self.planeWidget[i].RestrictPlaneToVolumeOn()
            color = [0, 0, 0]
            color[i] = 1
            self.planeWidget[i].GetPlaneProperty().SetColor(color)

            self.planeWidget[i].SetTexturePlaneProperty(ipwProp)  # property for the resliced image
            self.planeWidget[i].TextureInterpolateOff()
            self.planeWidget[i].SetResliceInterpolateToLinear()
            # self.planeWidget[i].SetInputConnection(reader.GetOutputPort())
            self.planeWidget[i].SetInputData(imgData)
            # self.planeWidget[i].SetInputData(imgData_mask)
            self.planeWidget[i].SetPlaneOrientation(i)
            self.planeWidget[i].SetSliceIndex(round(imageDims[i] / 2))
            self.planeWidget[i].DisplayTextOn()
            self.planeWidget[i].SetDefaultRenderer(ren)
            self.planeWidget[i].SetWindowLevel(np.diff(imgDataStat.GetAutoRange()), np.mean(imgDataStat.GetAutoRange()))
            # self.planeWidget[i].SetWindowLevel(1,0.5)

            self.planeWidget[i].On()
            self.planeWidget[i].InteractionOn()

        # Added. Initial view4 update.
        # If the reslice plane has modified, update it on the 3D widget
        for i in range(3):
            self.planeWidget[i].UpdatePlacement()
        self.planeWidget[0].GetInteractor().GetRenderWindow().Render()

        #################################################################################
        ### Finalize window
        #################################################################################
        # Setup callbacks
        for i in range(3):
            # ResliceAxesChangedEvent: mouse left button drag on axes
            self.riw[i].GetResliceCursorWidget().AddObserver(vtk.vtkResliceCursorWidget.ResliceAxesChangedEvent,
                                                             self.vtkResliceCursorCallback_ResliceAxesChangedEvent)
            # WindowLevelEvent: mouse left button drag on area other than axes
            self.riw[i].GetResliceCursorWidget().AddObserver(vtk.vtkResliceCursorWidget.WindowLevelEvent,
                                                             self.vtkResliceCursorCallback_WindowLevelEvent)
            # ResliceThicknessChangedEvent: mouse right button drag on axes
            self.riw[i].GetResliceCursorWidget().AddObserver(vtk.vtkResliceCursorWidget.ResliceThicknessChangedEvent,
                                                             self.vtkResliceCursorCallback_ResliceThicknessChangedEvent)
            # ResetCursorEvent: press 'o'
            self.riw[i].GetResliceCursorWidget().AddObserver(vtk.vtkResliceCursorWidget.ResetCursorEvent,
                                                             self.vtkResliceCursorCallback_ResetCursorEvent)
            # vtkCommand::WindowLevelEvent: in NON-oblique mode, mouse left button drag on area other than axes
            self.riw[i].GetInteractorStyle().AddObserver("WindowLevelEvent", self.vtkResliceCursorCallback)

            # Make them all share the same color map.
            self.riw[i].SetLookupTable(self.riw[0].GetLookupTable())
            self.planeWidget[i].GetColorMap().SetLookupTable(self.riw[0].GetLookupTable())
            self.planeWidget[i].SetColorMap(
                self.riw[i].GetResliceCursorWidget().GetResliceCursorRepresentation().GetColorMap())

            self.planeWidget[i].AddObserver("EndInteractionEvent",
                                            self.vtkImagePlaneWidgetCallback_EndInteractionEvent)  # Invoked when StopCursor or StopSliceMotion


            # ResliceAxesChangedEvent: mouse left button drag on axes
            self.riw[i].GetResliceCursorWidget().AddObserver("StartInteractionEvent",
                                                             self.vtkResliceCursorCallback_StartInteractionEvent)



        # Show each panes
        self.view1.show()
        self.view2.show()
        self.view3.show()

        # Set up action signals and slots
        self.resliceModeCheckBox.stateChanged.connect(self.resliceMode)
        self.thickModeCheckBox.stateChanged.connect(self.thickMode)
        self.thickModeCheckBox.setEnabled(0)

        self.radioButton_Max.setChecked(True)
        self.radioButton_Max.pressed.connect(self.SetBlendModeToMaxIP)
        self.radioButton_Min.pressed.connect(self.SetBlendModeToMinIP)
        self.radioButton_Mean.pressed.connect(self.SetBlendModeToMeanIP)
        self.blendModeGroupBox.setEnabled(0)

        self.resetButton.pressed.connect(self.ResetViews)
        self.AddDistance1Button.pressed.connect(self.AddDistanceMeasurementToView1)

        self.btnTest1.pressed.connect(self.cbTest1)
        self.btnTest2.pressed.connect(self.cbTest2)

        self.resliceModeCheckBox.setChecked(True)

        # Show QMainWindow
        self.show()


    #####################################################################################
    ### ImagePlaneWidget Callback
    #####################################################################################

    # Callback: vtkImagePlaneWidgetCallback EndInteractionEvent
    # After pushing or rotating of a ImagePlaneWidget, EndInteractionEvent is invoked.
    # 1. Rotate other ImagePlaneWidget to maintain orthogonality.
    # 2. Apply change to ResliceImageViewer panes.
    def vtkImagePlaneWidgetCallback_EndInteractionEvent(self, caller, ev):
        # Check in which plane this event is called.
        index = -1
        for i in range(3):
            if self.planeWidget[i] == caller:
                index = i
        if i == -1:
            print("(vtkImagePlaneWidgetCallback_LeftButtonReleaseEvent) Invalid caller")
            return

        #################################################################################
        # Response to ROTATE action
        #################################################################################
        # If a ImagePlaneWidget is rotated, rotate one of the other ImagePlaneWidget together
        # to maintain orthogonality.

        # Dot product of normal vectors.
        dot = [0.0] * 3
        baseplane_normal = np.array(self.planeWidget[index].GetPolyDataAlgorithm().GetNormal())
        for i in range(3):
            plane_normal = np.array(self.planeWidget[i].GetPolyDataAlgorithm().GetNormal())
            # dot[i] = np.abs(np.dot(baseplane_normal, plane_normal)) / np.linalg.norm(baseplane_normal) / np.linalg.norm(plane_normal)
            # norm of normal vector is always 1, so no need to divide with norm.
            dot[i] = np.abs(np.dot(baseplane_normal, plane_normal))

        # Sort the result of dot product.
        dot_orderedIndex = np.argsort(dot)

        # The sorted dot array is expected to be [0, 0.xxx, 1]
        # If all planes were orthogonal before manipulation, and only one plane is rotated, one plane is still orthogonal to the plane manipulated. ==> dot product = 0
        # A plane loose orthogonality. ==> dot product = 0.xxx
        # Manipulated plane itself has dot product 1.

        # Update normal vector of non-orthogonal plane.
        # Then update orthogonal plane also - this may not be necessary, but I want to ensure orthogonality.
        for i in [1, 0]:
            k = dot_orderedIndex[i]
            k1 = (k + 1) % 3
            k2 = (k + 2) % 3
            n1 = [0.0] * 3
            n2 = [0.0] * 3
            n1 = self.planeWidget[k1].GetPolyDataAlgorithm().GetNormal()
            n2 = self.planeWidget[k2].GetPolyDataAlgorithm().GetNormal()

            # newNormal = np.cross(self.planeWidget[k1].GetPolyDataAlgorithm().GetNormal(),self.planeWidget[k2].GetPolyDataAlgorithm().GetNormal())
            newNormal = np.cross(n1, n2)
            self.planeWidget[k].GetPolyDataAlgorithm().SetNormal(newNormal)
            self.planeWidget[k].UpdatePlacement()

        #################################################################################
        # Response to PUSH action
        #################################################################################
        # Find new intersection of planes. This will be applied to ResliceImageViewer panes later.

        # Calculate new intersection
        plane0 = self.planeWidget[0].GetPolyDataAlgorithm()
        plane1 = self.planeWidget[1].GetPolyDataAlgorithm()
        plane2 = self.planeWidget[2].GetPolyDataAlgorithm()
        newIntersection = self.CalculateIntersection3Planes(plane0, plane1, plane2)

        #################################################################################
        # Update ResliceImageViewer panes according to ImagePlaneWidget
        #################################################################################
        # Move center of all three ResliceCursor to the new intersection.
        self.riw[0].GetResliceCursor().SetCenter(newIntersection)
        # Rotate ResliceCursorWidget of ResliceImageViewer to match ImagePlaneWidget.
        for i in range(3):
            self.riw[0].GetResliceCursor().GetPlane(i).SetNormal(self.planeWidget[i].GetPolyDataAlgorithm().GetNormal())

        # Render in response to changes.
        self.RenderEverything()

    def CalculateIntersection3Planes(self, plane0, plane1, plane2):
        # Calculate intersect line of plane0 and plane1
        x0 = [0.0] * 3
        x1 = [0.0] * 3
        # vtk.vtkPlane.IntersectWithFinitePlane (plane0.GetNormal(), plane0.GetCenter(), plane1.GetOrigin(), plane1.GetPoint1(), plane1.GetPoint2(), x0, x1)
        self.IntersectWithFinitePlane(plane0.GetNormal(), plane0.GetCenter(), plane1.GetOrigin(), plane1.GetPoint1(),
                                      plane1.GetPoint2(), x0, x1)

        # Calculate intersection point with plane2
        t = vtk.mutable(0)
        x = [0.0] * 3
        vtk.vtkPlane.IntersectWithLine(x0, x1, plane2.GetNormal(), plane2.GetCenter(), t, x)

        return x

    # Redefine IntersectWithFinitePlane because vtk.vtkPlane.IntersectWithFinitePlane returns wrong answer.
    def IntersectWithFinitePlane(self, n, o, pOrigin, px, py, x0, x1):
        # Since we are dealing with convex shapes, if there is an intersection a
        # single line is produced as output. So all this is necessary is to
        # intersect the four bounding lines of the finite line and find the two
        # intersection points.
        numInts = 0
        t = vtk.mutable(0)
        x = x0

        # First line
        if (vtk.vtkPlane.IntersectWithLine(pOrigin, px, n, o, t, x)):
            numInts += 1
            x = x1

        # Second line
        if (vtk.vtkPlane.IntersectWithLine(pOrigin, py, n, o, t, x)):
            numInts += 1
            x = x1
        if (numInts == 2):
            return 1

        # Third line
        # Get third point of rectangle.
        # Original code: xr0[i] = pOrigin[i] + px[i] + py[i]
        # Correct code: xr0[i] = pOrigin[i] + (px[i] - pOrigin[i]) + (py[i] - pOrigin[i]) = px[i] + py[i] - pOrigin[i]
        xr0 = np.array(px) + np.array(py) - np.array(pOrigin)
        if (vtk.vtkPlane.IntersectWithLine(xr0, px, n, o, t, x)):
            numInts += 1
            x = x1
        if (numInts == 2):
            return 1

        # Fourth and last line
        if (vtk.vtkPlane.IntersectWithLine(xr0, py, n, o, t, x)):
            numInts += 1
            x = x1
        if (numInts == 2):
            return 1

        # No intersection has occurred, or a single degenerate point
        return 0


    #####################################################################################
    ### ResliceCursorWidget Callback
    #####################################################################################

    def vtkResliceCursorCallback_ResliceAxesChangedEvent(self, caller, ev):
        self.vtkResliceCursorCallback(caller, "vtkResliceCursorWidget.ResliceAxesChangedEvent")

    def vtkResliceCursorCallback_WindowLevelEvent(self, caller, ev):
        self.vtkResliceCursorCallback(caller, "vtkResliceCursorWidget.WindowLevelEvent")

    def vtkResliceCursorCallback_ResliceThicknessChangedEvent(self, caller, ev):
        self.vtkResliceCursorCallback(caller, "vtkResliceCursorWidget.ResliceThicknessChangedEvent")

    def vtkResliceCursorCallback_ResetCursorEvent(self, caller, ev):
        self.vtkResliceCursorCallback(caller, "vtkResliceCursorWidget.ResetCursorEvent")

    def vtkResliceCursorCallback_StartInteractionEvent(self, caller, ev):
        print("vtkResliceCursorCallback_StartInteractionEvent")
        if (isinstance(caller, vtk.vtkResliceCursorWidget)):
            print("resliceCursorWidget")
            rcw = caller
            rep = rcw.GetRepresentation()
            print(rep.GetManipulationMode())
            rep.SetManipulationMode(vtk.vtkResliceCursorRepresentation.RotateBothAxes)
            print("roateBothAxes")
            print(rep.GetManipulationMode())



    # https:#lorensen.github.io/VTKExamples/site/Python/Interaction/CallBack/
    def vtkResliceCursorCallback(self, caller, ev):
        # Who called callback and the event that triggered it.
        # print(caller.GetClassName(), "Event Id:", ev)

        if (ev == "vtkResliceCursorWidget.WindowLevelEvent") or \
                (ev == "WindowLevelEvent") or \
                (ev == "vtkResliceCursorWidget.ResliceThicknessChangedEvent"):
            # Render everything
            self.RenderEverything()

            return

        # The remained part is for ResliceAxesChangedEvent.
        if (isinstance(caller, vtk.vtkResliceCursorWidget)):
            rcw = caller

            rep = rcw.GetRepresentation()
            # Although the return value is not used, we keep the get calls in case they had side-effects
            cursor = rep.GetResliceCursorActor().GetCursorAlgorithm().GetResliceCursor()
            for i in range(3):
                ps = self.planeWidget[i].GetPolyDataAlgorithm()  ### vtkPlaneSource
                ps.SetOrigin(
                    self.riw[i].GetResliceCursorWidget().GetResliceCursorRepresentation().GetPlaneSource().GetOrigin())
                ps.SetPoint1(
                    self.riw[i].GetResliceCursorWidget().GetResliceCursorRepresentation().GetPlaneSource().GetPoint1())
                ps.SetPoint2(
                    self.riw[i].GetResliceCursorWidget().GetResliceCursorRepresentation().GetPlaneSource().GetPoint2())

                # If the reslice plane has modified, update it on the 3D widget
                self.planeWidget[i].UpdatePlacement()




            # # Setup vtkImageReslice to reslice overlay mask
            # print("#################")
            # for i in range(3):
            #     # Get slice axes
            #     resliceAxes = vtk.vtkMatrix4x4()
            #     resliceAxes.DeepCopy(self.riw[i].GetResliceCursorWidget().GetResliceCursorRepresentation().GetResliceAxes())
            #     print(self.riw[i].GetResliceCursorWidget().GetResliceCursorRepresentation().GetReslice())
            #
            #     resliceAxesTr = vtk.vtkMatrix4x4()
            #     resliceAxesTr.DeepCopy(resliceAxes)
            #     for i in range(3):
            #         for j in range(3):
            #             resliceAxesTr.SetElement(i,j, resliceAxes.GetElement(j,i))
            #
            #     print(f"###{i}###")
            #     # print(self.riw[i].GetResliceCursorWidget().GetResliceCursorRepresentation().GetResliceAxes())
            #     # resliceAxes.SetElement(0,3,198.75)
            #     # resliceAxes.SetElement(1,3,161.25)
            #     # resliceAxes.SetElement(2,3,61.25)
            #     print(resliceAxes)
            #     print(resliceAxesTr)
            #
            #     self.imageReslice[i].SetResliceAxes(resliceAxesTr)




        self.RenderEverything()

    def RenderEverything(self):
        # Render everything
        for i in range(3):
            self.riw[i].GetResliceCursorWidget().Render()
        # self.planeWidget[0].GetInteractor().GetRenderWindow().Render()
        self.view4.GetRenderWindow().Render()


    #####################################################################################
    ### Dialogbox UI Callback
    #####################################################################################

    def resliceMode(self, mode):
        if mode:
            self.thickModeCheckBox.setEnabled(1)
            self.blendModeGroupBox.setEnabled(1)
        else:
            self.thickModeCheckBox.setEnabled(0)
            self.blendModeGroupBox.setEnabled(0)

        print("### resliceMode ###")
        for i in range(3):
            if mode:
                self.riw[i].SetResliceMode(1)
            else:
                self.riw[i].SetResliceMode(0)
            self.riw[i].GetRenderer().ResetCamera()
            print(self.riw[i].GetRenderer().GetActiveCamera())
            self.riw[i].Render()

    def thickMode(self, mode):
        for i in range(3):
            if (mode):
                self.riw[i].SetThickMode(1)
            else:
                self.riw[i].SetThickMode(0)
            self.riw[i].GetRenderer().ResetCamera()
            self.riw[i].Render()

    def SetBlendMode(self, m):
        for i in range(3):
            thickSlabReslice = self.riw[i].GetResliceCursorWidget().GetRepresentation().GetReslice()
            thickSlabReslice.SetBlendMode(m)
            self.riw[i].Render()

    def SetBlendModeToMaxIP(self):
        # VTK_IMAGE_SLAB_MAX == 1
        self.SetBlendMode(1)

    def SetBlendModeToMinIP(self):
        # VTK_IMAGE_SLAB_MIN == 0
        self.SetBlendMode(0)

    def SetBlendModeToMeanIP(self):
        # VTK_IMAGE_SLAB_MEAN == 2
        self.SetBlendMode(2)

    def ResetViews(self):
        # Reset the reslice image views
        for i in range(3):
            self.riw[i].Reset()

        # Also sync the Image plane widget on the 3D top right view with any
        # changes to the reslice cursor.
        for i in range(3):
            ps = self.planeWidget[i].GetPolyDataAlgorithm()  # vtkPlaneSource
            ps.SetNormal(self.riw[0].GetResliceCursor().GetPlane(i).GetNormal())
            ps.SetCenter(self.riw[0].GetResliceCursor().GetPlane(i).GetOrigin())

            # If the reslice plane has modified, update it on the 3D widget
            self.planeWidget[i].UpdatePlacement()

        # Render in response to changes.
        self.RenderEverything()

    def AddDistanceMeasurementToView1(self):
        self.AddDistanceMeasurementToView(0)

    def AddDistanceMeasurementToView(self, i):
        # remove existing widgets.
        if (self.DistanceWidget[i] != 0):
            self.DistanceWidget[i].SetEnabled(0)
            self.DistanceWidget[i] = 0

        # add new widget
        self.DistanceWidget[i] = vtk.vtkDistanceWidget()
        self.DistanceWidget[i].SetInteractor(self.riw[i].GetResliceCursorWidget().GetInteractor())

        # Set a priority higher than our reslice cursor widget
        self.DistanceWidget[i].SetPriority(self.riw[i].GetResliceCursorWidget().GetPriority() + 0.01)

        handleRep = vtk.vtkPointHandleRepresentation2D()
        distanceRep = vtk.vtkDistanceRepresentation2D()
        distanceRep.SetHandleRepresentation(handleRep)
        self.DistanceWidget[i].SetRepresentation(distanceRep)
        distanceRep.InstantiateHandleRepresentation()
        distanceRep.GetPoint1Representation().SetPointPlacer(self.riw[i].GetPointPlacer())
        distanceRep.GetPoint2Representation().SetPointPlacer(self.riw[i].GetPointPlacer())

        # Add the distance to the list of widgets whose visibility is managed based
        # on the reslice plane by the ResliceImageViewerMeasurements class
        self.riw[i].GetMeasurements().AddItem(self.DistanceWidget[i])

        self.DistanceWidget[i].CreateDefaultRepresentation()
        self.DistanceWidget[i].EnabledOn()

    def cbTest1(self):
        # print(self.riw[1].GetRenderer().GetActiveCamera())
        self.cameraPosition = self.riw[1].GetRenderer().GetActiveCamera().GetPosition()
        self.cameraFocalPoint = self.riw[1].GetRenderer().GetActiveCamera().GetFocalPoint()
        self.cameraViewUp = self.riw[1].GetRenderer().GetActiveCamera().GetViewUp()
        self.cameraViewAngle = self.riw[1].GetRenderer().GetActiveCamera().GetViewAngle()
        self.cameraParallelProjection = self.riw[1].GetRenderer().GetActiveCamera().GetParallelProjection()
        self.cameraParallelScale = self.riw[1].GetRenderer().GetActiveCamera().GetParallelScale()


    def cbTest2(self):
        # print(self.riw[1].GetRenderer().GetActiveCamera())
        self.riw[1].GetRenderer().GetActiveCamera().SetPosition(self.cameraPosition)
        self.riw[1].GetRenderer().GetActiveCamera().SetFocalPoint(self.cameraFocalPoint)
        self.riw[1].GetRenderer().GetActiveCamera().SetViewUp(self.cameraViewUp)
        self.riw[1].GetRenderer().GetActiveCamera().SetViewAngle(self.cameraViewAngle)
        self.riw[1].GetRenderer().GetActiveCamera().SetParallelProjection(self.cameraParallelProjection)
        self.riw[1].GetRenderer().GetActiveCamera().SetParallelScale(self.cameraParallelScale)
        self.riw[1].Render()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = QtVTKRenderWindows()
    sys.exit(app.exec())
