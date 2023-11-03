# pyuic5 -x Widget_TimeSeriesView.ui -o Ui_Widget_TimeSeriesView.py

from Ui_Widget_TimeSeriesView import Ui_Widget_TimeSeriesView
from PyQt5.QtWidgets import QWidget, QColorDialog, QFileDialog
import vtk
import WrapperVtkFunction
import numpy as np
from PyQt5 import QtCore
from PyQt5.QtCore import QFileInfo, Qt, QTime, QDate, QTimer

import matplotlib.pyplot as plt

import copy
import Utility
import scipy.signal

import FourDFlow
import VtkCustom


class Widget_TimeSeriesView(QWidget, Ui_Widget_TimeSeriesView):
    sliderColorWindowFactor = 1000.0

    def __init__(self, mrStructScalar, mrStructVector, mrStructMask):
        super().__init__()
        self.setupUi(self)
        self.checkShowLargest.stateChanged.connect(self.on_checkShowLargest_stateChanged)

        self.AdditionalActors = []
        # self.parent = parent

        # Setup self.source
        self.SetSourceData(mrStructScalar, mrStructVector)

        # Setup View
        style = vtk.vtkInteractorStyleSwitch()
        style.SetCurrentStyleToTrackballCamera()
        self.vtkWidgetVolumeRender.SetInteractorStyle(style)

        # Test visualization
        self.volRenWin = self.vtkWidgetVolumeRender.GetRenderWindow()
        self.volRender = vtk.vtkRenderer()
        self.volRenWin.AddRenderer(self.volRender)

        self.UpdateView()

    def SetSourceData(self, mrStructScalar, mrStructVector):
        self.mrStructScalar = mrStructScalar
        self.mrStructVector = mrStructVector

        self.source = FourDFlow.MrStruct4DFlowSource()
        self.source.SetTimeValue(0)
        self.source.SetMrStructScalar(mrStructScalar)
        self.source.SetMrStructVector(mrStructVector)
        self.source.SetTimeSpaceInverse(False)
        self.source.Update()

        self.sliderTimePhase.setRange(0, self.source.GetNumberOfTimeFrames() - 1)
        self.sliderTimePhase.setValue(0)

        self.actorMarchingcube = [None] * self.source.GetNumberOfTimeFrames()
        self.actorRaycast = [None] * self.source.GetNumberOfTimeFrames()
        self.actorHedgehog = [None] * self.source.GetNumberOfTimeFrames()
        self.actorArrow = [None] * self.source.GetNumberOfTimeFrames()

        # self.source.SetTimeOffset(9)
        # self.source.SetFilepathMask()
        # self.source.SetTimeSpaceInverse(True)
        # print(self.source.GetMaskArray().shape)

    def UpdateView(self):
        # print("Update View")
        timePhase = self.sliderTimePhase.value()
        self.source.SetTimeValue(timePhase)
        self.volRender.RemoveAllViewProps()

        imageData = None

        # === Add actors ===
        # Scalar: Marching Cube
        if self.checkBoxShowSurface.checkState() == Qt.Checked:
            if (self.actorMarchingcube[timePhase] == None):
                if imageData == None:
                    imageData = self.source.GetImageDataAtTimeindex(timePhase)
                self.actorMarchingcube[timePhase] = VtkCustom.vtkVisualization.GetMarchingCubeActor(imageData)
            self.volRender.AddActor(self.actorMarchingcube[timePhase])

        # Scalar: Ray Cast
        if self.checkBoxShowRayCast.checkState() == Qt.Checked:
            if (self.actorRaycast[timePhase] == None):
                if imageData == None:
                    imageData = self.source.GetImageDataAtTimeindex(timePhase)
                self.actorRaycast[timePhase] = VtkCustom.vtkVisualization.GetRayCastActor(imageData)
            self.volRender.AddActor(self.actorRaycast[timePhase])

        # Vector: Hedgehog
        if self.checkBoxShowHedgehog.checkState() == Qt.Checked:
            if (self.actorHedgehog[timePhase] == None):
                if imageData == None:
                    imageData = self.source.GetImageDataAtTimeindex(timePhase)
                self.actorHedgehog[timePhase] = VtkCustom.vtkVisualization.GetHedgehogVectorActor(imageData)
            self.volRender.AddActor(self.actorHedgehog[timePhase])

        # Vector: Arrow
        if self.checkBoxShowArrow.checkState() == Qt.Checked:
            if (self.actorArrow[timePhase] == None):
                if imageData == None:
                    imageData = self.source.GetImageDataAtTimeindex(timePhase)
                self.actorArrow[timePhase] = VtkCustom.vtkVisualization.GetArrowVectorActor(imageData)
            self.volRender.AddActor(self.actorArrow[timePhase])

        # Additional Actors
        for actor in self.AdditionalActors:
            self.volRender.AddActor(actor)

        # Render
        self.volRenWin.Render()

    ## callback for Show Surface Checkbox
    @QtCore.pyqtSlot(int)
    def on_checkBoxShowSurface_stateChanged(self, value):
        print("on_checkBoxShowSurface_stateChanged")
        self.UpdateView()
        print("completed")

    ## callback for Show Arrow Checkbox
    @QtCore.pyqtSlot(int)
    def on_checkBoxShowRayCast_stateChanged(self, value):
        print("on_checkBoxShowRayCast_stateChanged")
        self.UpdateView()
        print("completed")

    ## callback for Show Arrow Checkbox
    @QtCore.pyqtSlot(int)
    def on_checkBoxShowHedgehog_stateChanged(self, value):
        print("on_checkBoxShowHedgehog_stateChanged")
        self.UpdateView()
        print("completed")

    ## callback for Show Arrow Checkbox
    @QtCore.pyqtSlot(int)
    def on_checkBoxShowArrow_stateChanged(self, value):
        print("on_checkBoxShowArrow_stateChanged")
        self.UpdateView()
        print("completed")

    # callback for Time Phase Slider
    @QtCore.pyqtSlot(int)
    def on_sliderTimePhase_valueChanged(self, value):
        self.labelTimePhase.setText(f"{value + 1} / {self.source.GetNumberOfTimeFrames()}")
        self.UpdateView()

    ## Callback for Animation Start/Stop Button
    @QtCore.pyqtSlot()
    def on_btnAnimationStartStop_clicked(self):
        if (self.btnAnimationStartStop.text() == "Start"):
            print("Animation Start")
            self.btnAnimationStartStop.setText("Stop")
            self.timerAnimation = QTimer(self)
            self.timerAnimation.timeout.connect(self.AnimationProgress)
            self.timerAnimation.start(200)
        else:
            print("Animation Stop")
            self.btnAnimationStartStop.setText("Start")
            self.timerAnimation.stop()

    def AnimationProgress(self):
        if (self.sliderTimePhase.value() < self.source.GetNumberOfTimeFrames() - 1):
            self.sliderTimePhase.setValue(self.sliderTimePhase.value() + 1)
        else:
            self.sliderTimePhase.setValue(0)
        # self.UpdateView()
        # UpdateView is included in Time Phase Slider value changed.

    def AddMaskActor(self, mrStructMask):
        imageData = VtkCustom.vtkUtil.convert_numpyArray_vtkImageData(mrStructMask.dataAy)
        imageData.SetSpacing(mrStructMask.vox[0:3])
        maskActor = VtkCustom.vtkVisualization.GetMarchingCubeBinaryMaskActor(imageData)
        self.AdditionalActors.append(maskActor)
        self.UpdateView()

    def MaskCutData(self, mrStructMask):
        mrStructScalar = copy.deepcopy(self.mrStructScalar)
        mrStructVector = copy.deepcopy(self.mrStructVector)
        mrStructScalar.dataAy = np.multiply(mrStructScalar.dataAy, mrStructMask.dataAy[:, :, :, np.newaxis])
        mrStructVector.dataAy = np.multiply(mrStructVector.dataAy, mrStructMask.dataAy[:, :, :, np.newaxis, np.newaxis])
        self.SetSourceData(mrStructScalar, mrStructVector)
        self.UpdateView()

    # ====================================================================================================================
    # ====================================================================================================================
    # ====================================================================================================================
    # ====================================================================================================================
    # ====================================================================================================================

    def setupUiPane3D(self, indexPane):
        self.vtkWidgetVolume = self.paneWidgets[indexPane][1]
        self.labelTitlePaneVolume = self.paneWidgets[indexPane][0]
        self.sliderPaneVolume = self.paneWidgets[indexPane][2]
        self.sliderPaneVolume.setTracking(False)

        self.volRenWin = self.vtkWidgetVolume.GetRenderWindow()
        self.volRender = vtk.vtkRenderer()
        self.volRenWin.AddRenderer(self.volRender)

        style = vtk.vtkInteractorStyleSwitch()
        style.SetCurrentStyleToTrackballCamera()
        self.vtkWidgetVolume.SetInteractorStyle(style)

        ### Main Volume
        # Marching Cubes
        self.marchingCubes = vtk.vtkMarchingCubes()
        self.marchingCubes.ComputeScalarsOff()
        self.marchingCubes.ComputeGradientsOff()  # turn off any time-expensive functions
        self.marchingCubes.ComputeNormalsOff()

        self.connectivityFilter = vtk.vtkPolyDataConnectivityFilter()
        self.connectivityFilter.SetInputConnection(
            self.marchingCubes.GetOutputPort())  # feed result to connectivity filter
        # Filter to keep the largest region
        self.connectivityFilter.SetExtractionModeToLargestRegion()  # Set mode of the extraction of connected surfaces

        # Smooth Surface
        self.smoother = vtk.vtkWindowedSincPolyDataFilter()
        # self.smoother.SetInputConnection(self.connectivityFilter.GetOutputPort())
        if self.checkShowLargest.isChecked():
            self.smoother.SetInputConnection(self.connectivityFilter.GetOutputPort())
        else:
            self.smoother.SetInputConnection(self.marchingCubes.GetOutputPort())
        smoothingIterations = 10
        passBand = 0.001
        featureAngle = 60.0
        self.smoother.SetNumberOfIterations(smoothingIterations)
        self.smoother.BoundarySmoothingOff()
        self.smoother.FeatureEdgeSmoothingOff()
        self.smoother.SetFeatureAngle(featureAngle)
        self.smoother.SetPassBand(passBand)
        self.smoother.NonManifoldSmoothingOn()
        self.smoother.NormalizeCoordinatesOn()

        # make mapper
        self.volumeMapper = vtk.vtkPolyDataMapper()
        self.volumeMapper.SetInputConnection(self.smoother.GetOutputPort())
        self.volumeMapper.ScalarVisibilityOff()

        # maker actor
        self.volActor = vtk.vtkActor()
        self.volActor.SetMapper(self.volumeMapper)
        self.volActor.GetProperty().SetDiffuseColor(0.0, 1.0, 0.0)
        self.volActor.GetProperty().SetOpacity(0.4)

        ### Background Volume
        # Marching Cubes
        self.marchingCubesBackground = vtk.vtkMarchingCubes()
        self.marchingCubesBackground.ComputeScalarsOff()
        self.marchingCubesBackground.ComputeGradientsOff()  # turn off any time-expensive functions
        self.marchingCubesBackground.ComputeNormalsOff()

        connectivityFilterBackground = vtk.vtkPolyDataConnectivityFilter()
        connectivityFilterBackground.SetInputConnection(
            self.marchingCubesBackground.GetOutputPort())  # feed result to connectivity filter
        # Filter to keep the largest region
        connectivityFilterBackground.SetExtractionModeToLargestRegion()  # Set mode of the extraction of connected surfaces

        # Smooth Surface
        smootherBackground = vtk.vtkWindowedSincPolyDataFilter()
        smootherBackground.SetInputConnection(connectivityFilterBackground.GetOutputPort())
        smoothingIterations = 10
        passBand = 0.001
        featureAngle = 60.0
        smootherBackground.SetNumberOfIterations(smoothingIterations)
        smootherBackground.BoundarySmoothingOff()
        smootherBackground.FeatureEdgeSmoothingOff()
        smootherBackground.SetFeatureAngle(featureAngle)
        smootherBackground.SetPassBand(passBand)
        smootherBackground.NonManifoldSmoothingOn()
        smootherBackground.NormalizeCoordinatesOn()
        # smootherBackground.Update()

        # make mapper
        self.volumeMapperBackground = vtk.vtkPolyDataMapper()
        self.volumeMapperBackground.SetInputConnection(smootherBackground.GetOutputPort())
        self.volumeMapperBackground.ScalarVisibilityOff()

        # maker actor
        self.volActorBackground = vtk.vtkActor()
        self.volActorBackground.SetMapper(self.volumeMapperBackground)
        self.volActorBackground.GetProperty().SetDiffuseColor(1.0, 0.0, 0.0)
        self.volActorBackground.GetProperty().SetOpacity(0.4)

        ###
        # Add actor and Set camera on the renderer
        # self.volRender.AddActor(self.volActorBackground)
        self.volRender.AddActor(self.volActor)

        self.sliderPaneVolume.setRange(1, 10000)
        self.sliderPaneVolume.setValue(1000)

    # Actually generate contour lines.
    def BeginInteraction(self, obj, event):
        obj.GetPolyData(self.plane)
        self.contourActor.VisibilityOn()

    def ProbeData(self, obj, event):
        obj.GetPolyData(self.plane)

    def setupUiPane3DPlaneWidget(self, indexPane):
        if (hasattr(self, 'planeWidget') == False):
            self.vtkWidgetVolume = self.paneWidgets[indexPane][1]
            self.labelTitlePaneVolume = self.paneWidgets[indexPane][0]
            self.sliderPaneVolume = self.paneWidgets[indexPane][2]
            self.sliderPaneVolume.setTracking(False)
            self.volRenWin = self.vtkWidgetVolume.GetRenderWindow()

            # The plane widget is used probe the dataset.
            self.planeWidget = vtk.vtkPlaneWidget()
            # self.planeWidget.SetInputData(pl3d_output)
            self.planeWidget.SetInputData(self.imageData[0])
            self.planeWidget.NormalToXAxisOn()
            self.planeWidget.SetResolution(20)
            self.planeWidget.SetRepresentationToOutline()
            self.planeWidget.PlaceWidget()
            self.plane = vtk.vtkPolyData()
            self.planeWidget.GetPolyData(self.plane)

            self.planeWidget.SetPlaceFactor(0.2)
            self.planeWidget.PlaceWidget()

            probe = vtk.vtkProbeFilter()
            probe.SetInputData(self.plane)
            # probe.SetSourceData(pl3d_output)
            probe.SetSourceData(self.imageData[0])

            contourMapper = vtk.vtkPolyDataMapper()
            contourMapper.SetInputConnection(probe.GetOutputPort())
            # contourMapper.SetScalarRange(pl3d_output.GetScalarRange())
            contourMapper.SetScalarRange(self.imageData[0].GetScalarRange())
            self.contourActor = vtk.vtkActor()
            self.contourActor.SetMapper(contourMapper)
            self.contourActor.VisibilityOff()

            ren = self.volRenWin.GetRenderers().GetFirstRenderer()

            # Associate the widget with the interactor
            self.planeWidget.SetInteractor(self.vtkWidgetVolume)
            # Handle the events.
            self.planeWidget.AddObserver("EnableEvent", self.BeginInteraction)
            self.planeWidget.AddObserver("StartInteractionEvent", self.BeginInteraction)
            self.planeWidget.AddObserver("InteractionEvent", self.ProbeData)

            ren.AddActor(self.contourActor)

            self.volRenWin.Render()

            self.planeWidget.On()
        elif (self.planeWidget.GetEnabled()):
            self.planeWidget.Off()
        else:
            self.planeWidget.On()

    def setDataArray(self, mrStruct):
        self.mrStruct = copy.deepcopy(mrStruct)
        self.dimension = self.mrStruct.dataAy.ndim
        if (self.dimension == 4):
            self.dataArray = self.mrStruct.dataAy
        elif (self.dimension == 3):
            self.dataArray = np.reshape(self.mrStruct.dataAy, np.append(np.asarray(self.mrStruct.dataAy.shape), 1))
        elif (self.dimension == 5):
            self.dataArray = self.mrStruct.dataAy[:, :, :, 0, :]

        self.sliceNumber = (np.asarray(self.mrStruct.dataAy.shape) / 2).astype(int)

        # Image viewer windowing is working based on integer. If the value is too small, it is not working properly.
        # Normalize data value to certain large value. Here it is arbitrary number 10,000.
        self.dataRange = [self.dataArray.min(), self.dataArray.max()]
        # print(f"range: {self.dataRange}")

        shapeDataArray = self.dataArray.shape
        self.dataExtent = shapeDataArray[0:3]
        self.timePhaseExtent = shapeDataArray[3]
        self.imageData = [0] * self.timePhaseExtent
        for i in range(self.timePhaseExtent):
            np_temp = self.dataArray[:, :, :, i]
            self.imageData[i] = WrapperVtkFunction.convert_numpyArray_vtkImageData(np_temp)
            self.imageData[i].SetSpacing(self.mrStruct.vox[0:3])
            # self.imageData[i].SetOrigin( referenceVtkImageData.GetOrigin() )

        self.dataRange = [np.min(self.dataArray), np.max(self.dataArray)]

        # Set range and value for windowing sliders
        self.sliderColorLevel.setRange(int(self.dataRange[0] * self.sliderColorWindowFactor),
                                       int(self.dataRange[1] * self.sliderColorWindowFactor))
        self.sliderColorWindow.setRange(1, int((self.dataRange[1] - self.dataRange[0]) * self.sliderColorWindowFactor))
        self.sliderColorLevel.setValue(np.mean(self.dataRange) * self.sliderColorWindowFactor)
        self.sliderColorWindow.setValue(np.diff(self.dataRange) * self.sliderColorWindowFactor)

        self.sliderPaneVolume.setRange(int(self.dataRange[0] * self.sliderColorWindowFactor),
                                       int(self.dataRange[1] * self.sliderColorWindowFactor))
        self.sliderPaneVolume.setValue(
            (self.dataRange[0] + np.diff(self.dataRange) * 0.95) * self.sliderColorWindowFactor)

        self.sliderTimePhase.setRange(0, self.timePhaseExtent - 1)
        self.sliderTimePhase.setValue(0)

        # Set current slice to the middle one
        for i in range(3):
            self.drawSlice(i)

        # Set input for volume renderer
        self.marchingCubes.SetInputData(self.imageData[0])  # set thresholded data as input
        self.marchingCubes.SetValue(0, 0.95)  # set isovalue of surface
        self.marchingCubes.Update()
        self.volRenWin.Render()

    def setDataArrayBackground(self, mrStruct):
        self.bgDimension = mrStruct.dataAy.ndim
        if (self.bgDimension == 4):
            self.bgDataArray = mrStruct.dataAy
        elif (self.bgDimension == 3):
            self.bgDataArray = np.reshape(mrStruct.dataAy, np.append(np.asarray(mrStruct.dataAy.shape), 1))
        elif (self.bgDimension == 5):
            self.bgDataArray = mrStruct.dataAy[:, :, :, 0, :]

        # self.bgDataRange = [self.bgDataArray.min(), self.bgDataArray.max()]
        #
        shapeDataArray = self.bgDataArray.shape
        # self.dataExtent = shapeDataArray[0:3]
        self.bgTimePhaseExtent = shapeDataArray[3]
        self.bgImageData = [0] * self.bgTimePhaseExtent
        for i in range(self.bgTimePhaseExtent):
            np_temp = self.bgDataArray[:, :, :, i]
            self.bgImageData[i] = WrapperVtkFunction.convert_numpyArray_vtkImageData(np_temp)
            self.bgImageData[i].SetSpacing(mrStruct.vox[0:3])
            # self.imageData[i].SetOrigin( referenceVtkImageData.GetOrigin() )

        # self.dataRange = [np.min(self.dataArray), np.max(self.dataArray)]
        #
        # self.sliderPaneVolume.setRange(int(self.dataRange[0] * self.sliderColorWindowFactor),
        #                                int(self.dataRange[1] * self.sliderColorWindowFactor))
        # self.sliderPaneVolume.setValue(
        #     (self.dataRange[0] + np.diff(self.dataRange) * 0.95) * self.sliderColorWindowFactor)

        # Set input for volume renderer
        self.volRender.AddActor(self.volActorBackground)
        timePhase = self.sliderTimePhase.value()
        self.marchingCubesBackground.SetInputData(self.bgImageData[timePhase])  # set thresholded data as input
        self.marchingCubesBackground.SetValue(0, 0.95)  # set isovalue of surface
        self.marchingCubesBackground.Update()
        self.volRenWin.Render()

    def drawSlice(self, axis):
        timePhase = self.sliderTimePhase.value()
        imageData = self.imageData[timePhase]

        viewer = self.viewer[axis]
        viewer.SetInputData(imageData)
        viewer.SetSlice(self.sliceNumber[axis])
        viewer.Render()

        self.sliderPane[axis].setRange(0, self.dataExtent[axis])
        self.sliderPane[axis].setValue(self.sliceNumber[axis])

        # renWin = self.vtkWidget[axis].GetRenderWindow()
        # # renWin.AddRenderer(ren)
        # iren = renWin.GetInteractor()
        # iren.Start()
        #
        # # Connect an interactor to the image viewer
        # iren.SetInteractorStyle(vtk.vtkInteractorStyleImage())

    def drawCone(self, vtkWidget):
        colors = vtk.vtkNamedColors()
        bkg = map(lambda x: x / 255.0, [26, 51, 102, 255])
        colors.SetColor("BkgColor", *bkg)

        cylinder = vtk.vtkCylinderSource()
        cylinder.SetResolution(20)

        cylinderMapper = vtk.vtkPolyDataMapper()
        cylinderMapper.SetInputConnection(cylinder.GetOutputPort())

        cylinderActor = vtk.vtkActor()
        cylinderActor.SetMapper(cylinderMapper)
        cylinderActor.GetProperty().SetColor(colors.GetColor3d("Tomato"))
        cylinderActor.RotateX(30.0)
        cylinderActor.RotateY(-45.0)

        ren = vtk.vtkRenderer()
        ren.AddActor(cylinderActor)
        ren.SetBackground(colors.GetColor3d("BkgColor"))
        ren.ResetCamera()
        ren.GetActiveCamera().Zoom(1.5)

        # Set target vtkWidget
        renWin = vtkWidget.GetRenderWindow()
        renWin.AddRenderer(ren)
        # iren = renWin.GetInteractor()
        # iren.Start()
        self.volRenWin.Render()

    # setup slots for slicing sliders
    @QtCore.pyqtSlot(int)
    def on_sliderPane1_valueChanged(self, value):
        for i in range(3):
            if self.sliderPane[i] == self.sliderPane1:
                self.viewer[i].SetSlice(value)
                self.labelTitlePane[i].setText(f"Pane 1: {value}")

    @QtCore.pyqtSlot(int)
    def on_sliderPane2_valueChanged(self, value):
        for i in range(3):
            if self.sliderPane[i] == self.sliderPane2:
                self.viewer[i].SetSlice(value)
                self.labelTitlePane[i].setText(f"Pane 2: {value}")

    @QtCore.pyqtSlot(int)
    def on_sliderPane3_valueChanged(self, value):
        for i in range(3):
            if self.sliderPane[i] == self.sliderPane3:
                self.viewer[i].SetSlice(value)
                self.labelTitlePane[i].setText(f"Pane 3: {value}")

    @QtCore.pyqtSlot(int)
    def on_sliderPane4_valueChanged(self, value):
        # slider: multiplied by sliderColorWindowFactor (1000.0)
        value_converted = value / self.sliderColorWindowFactor
        self.labelTitlePaneVolume.setText(f"Pane 3: {value} ==> {value_converted}")

        self.marchingCubes.SetValue(0, value_converted)  # set isovalue of surface

        if (self.marchingCubes.GetTotalNumberOfInputConnections() > 0):
            self.marchingCubes.Update()
            self.volRenWin.Render()

    # Setup slots for windowing sliders
    @QtCore.pyqtSlot(int)
    def on_sliderColorLevel_valueChanged(self, value):
        for viewer in self.viewer:
            viewer.SetColorLevel(value / self.sliderColorWindowFactor)
            viewer.Render()

    @QtCore.pyqtSlot(int)
    def on_sliderColorWindow_valueChanged(self, value):
        for viewer in self.viewer:
            viewer.SetColorWindow(value / self.sliderColorWindowFactor)
            viewer.Render()

    @QtCore.pyqtSlot()
    def on_btnSetColor3D_clicked(self):
        color = QColorDialog.getColor()
        # rgbcolor = color.getRgb() # = rgba = between 0 and 255
        rgbfcolor = color.getRgbF()  # = rgba / 255.0 = between 0.0 and 1.0
        self.volActor.GetProperty().SetDiffuseColor(rgbfcolor[0:3])
        self.volRenWin.Render()

    @QtCore.pyqtSlot()
    def on_btnLoadBackground_clicked(self):
        print("load background")

        index = -1
        for i in range(self.parent.tableListMrStruct.rowCount()):
            itemCheck = self.parent.tableListMrStruct.item(i, 0)
            if (itemCheck.checkState() == Qt.Checked):
                index = i
                break

        if (index >= 0):
            itemCheck = self.parent.tableListMrStruct.item(index, 0)
            print(f"checked {index} {itemCheck.text()}")
            itemMrStruct = itemCheck.data(Qt.UserRole)
            print(f"mrStruct {itemMrStruct.name}")
            self.setDataArrayBackground(itemMrStruct)
        else:
            # Remove background actor from volume renderer
            self.volRender.RemoveActor(self.volActorBackground)
            self.volRenWin.Render()

    @QtCore.pyqtSlot()
    def on_btnPlaneWidget_clicked(self):
        self.volRenWin = self.vtkWidgetVolume.GetRenderWindow()
        renderers = self.volRenWin.GetRenderers()
        # for renderer in renderers:
        #     self.volRenWin.RemoveRenderer(renderer)

        # self.drawCone(self.vtkWidgetVolume)
        self.setupUiPane3DPlaneWidget(3)
        # self.planeWidget.On()

    @QtCore.pyqtSlot()
    def on_btnCut_clicked(self):
        timePhase = self.sliderTimePhase.value()
        temp_plane = vtk.vtkPlane()
        if (hasattr(self, 'planeWidget') == False):
            return

        if (self.planeWidget.GetEnabled() == False):
            return

        self.planeWidget.GetPlane(temp_plane)
        print(f"GetPlane: {temp_plane}")
        nxa, nya, nza = temp_plane.GetNormal()
        x0a, y0a, z0a = temp_plane.GetOrigin()
        nx = nxa * self.mrStruct.vox[0]
        ny = nya * self.mrStruct.vox[1]
        nz = nza * self.mrStruct.vox[2]
        x0 = x0a / self.mrStruct.vox[0]
        y0 = y0a / self.mrStruct.vox[1]
        z0 = z0a / self.mrStruct.vox[2]

        extent = self.dataArray.shape[0:3]
        mask = np.zeros(extent)

        print("cut")
        ygrid, xgrid, zgrid = np.meshgrid(np.arange(extent[1]), np.arange(extent[0]), np.arange(extent[2]))
        mask = (xgrid - x0) * nx + (ygrid - y0) * ny + (zgrid - z0) * nz > 0

        # temp_array = self.dataArray[:, :, :, timePhase] * mask
        np_currentMask = WrapperVtkFunction.convert_vtkImageData_numpyArray(self.imageData[timePhase])
        temp_array = np_currentMask * mask

        self.imageData[timePhase] = WrapperVtkFunction.convert_numpyArray_vtkImageData(temp_array)
        self.imageData[timePhase].SetSpacing(self.mrStruct.vox[0:3])
        self.marchingCubes.SetInputData(self.imageData[timePhase])  # set thresholded data as input
        self.marchingCubes.Update()
        self.volRenWin.Render()

    @QtCore.pyqtSlot()
    def on_btnCutReset_clicked(self):
        timePhase = self.sliderTimePhase.value()
        temp_array = self.dataArray[:, :, :, timePhase]
        self.imageData[timePhase] = WrapperVtkFunction.convert_numpyArray_vtkImageData(temp_array)
        self.imageData[timePhase].SetSpacing(self.mrStruct.vox[0:3])
        self.marchingCubes.SetInputData(self.imageData[timePhase])  # set thresholded data as input
        self.marchingCubes.Update()
        self.volRenWin.Render()

    @QtCore.pyqtSlot()
    def on_btnExportSegmentation_clicked(self):
        filename = QFileDialog.getSaveFileName(self, 'Export Segmentation', filter="MATLAB Data (*.mat)")
        if filename[0]:
            temp_mrStruct = copy.deepcopy(self.mrStruct)
            timePhase = self.sliderTimePhase.value()
            temp_mrStruct.dataAy = WrapperVtkFunction.convert_vtkImageData_numpyArray(self.imageData[timePhase])
            temp_mrStruct.exportMat(filename[0])
            print(f"Current segmentation is exported to {filename[0]}")

    @QtCore.pyqtSlot()
    def on_btnCalculate_clicked(self):
        print("calculate")
        index = -1
        for i in range(self.parent.tableListMrStruct.rowCount()):
            itemCheck = self.parent.tableListMrStruct.item(i, 0)
            if (itemCheck.text() == "vel"):
                index = i
                break

        if (index >= 0):
            itemVel = self.parent.tableListMrStruct.item(index, 0)
            velMrStruct = itemVel.data(Qt.UserRole)

            # Make mask and count cells
            timePhase = self.sliderTimePhase.value()
            npMask = WrapperVtkFunction.convert_vtkImageData_numpyArray(self.imageData[timePhase])

            # Calculate Kinetic Energy
            unitVolume, numOfCells, volume, KE = Utility.calcKineticEnergy(velMrStruct, npMask)

            print("=" * 60)
            print(f"Volume: {unitVolume:>10.4f} [cm3/cell] * {numOfCells} cells = {volume:>10.4f} [cm3]")
            print(f"timePhase  KE         KEindexed")
            for i in range(KE.shape[0]):
                kineticEnergy = KE[i, 0]
                kineticEnergyIndexed = KE[i, 1]
                print(f"{i:9} {kineticEnergy:>10.4f} {kineticEnergyIndexed:>10.4f}")

            max, min = Utility.getMaximaMinima(KE[:, 1])

            print(f"minima: {min}")
            print(f"maxima: {max}")

            # parameters
            #     time avearge KE
            #     peak systole
            #     peak diastole
            #     systole: from maximum peak, trace back and forth
            #     diastole: from second peak, trace back and forth
            if ((np.max(min) <= max[0]) or (np.min(min) >= max[0])):
                t_ES = np.min(min)
                t_ED = np.max(min)
            else:
                t_ED = np.min(min)
                t_ES = np.max(min)

            KEavg = np.mean(KE[:, 1])
            KEpeakSystole = KE[max[0], 1]
            KEpeakDiastole = KE[max[1], 1]

            if (t_ED < t_ES):
                KEavgSystole = np.mean(KE[t_ED:t_ES, 1])
                KEavgDiastole = np.mean(np.concatenate((KE[t_ES:, 1], KE[:t_ED, 1])))
            else:
                KEavgDiastole = np.mean(KE[t_ES:t_ED, 1])

                print('-' * 60)
                print(KE[:, 1])
                print('-' * 60)
                print(KE[t_ED:, 1])
                print('-' * 60)
                print(KE[:t_ES, 1])
                print('-' * 60)
                print(np.concatenate((KE[t_ED:, 1], KE[:t_ES, 1])))
                print('-' * 60)
                KEavgSystole = np.mean(np.concatenate((KE[t_ED:, 1], KE[:t_ES, 1])))

            strArrayHeader = Utility.convertArrayToString(np.arange(30))
            strArrayValue = Utility.convertArrayToString(KE[:, 1])
            print(f" KEavg     KEpeakSyst KEpeakDias  KEavgSyst  KEavgDias t_Syst t_Dias t_ED t_ES {strArrayHeader}")
            print(
                f"{KEavg:>10.4f} {KEpeakSystole:>10.4f} {KEpeakDiastole:>10.4f} {KEavgSystole:>10.4f} {KEavgDiastole:>10.4f} {max[0]:>6} {max[1]:>6} {t_ED:>4} {t_ES:>4} {strArrayValue}")

            plt.figure()
            plt.plot(KE[:, 1])
            plt.ylabel('Kinetic Energy')
            plt.xlabel('Time Phase')
            plt.show()

            # filtered = np.savitzky_golay(KE[:,1], 5, 3) # window size 51, polynomial order 3
            filtered = scipy.signal.savgol_filter(KE[:, 1], 5, 3, mode='wrap')  # window size 51, polynomial order 3

            plt.figure()
            plt.plot(filtered)
            plt.ylabel('Filtered Kinetic Energy')
            plt.xlabel('Time Phase')
            plt.show()

    @QtCore.pyqtSlot()
    def on_btnStreamline_clicked(self):
        print("Streamline")

        if (hasattr(self, 'streamActor') == True):
            self.volRender.RemoveActor(self.streamActor)
            self.volRenWin.Render()
            del self.streamActor
            return

        ### Check PlaneWidget
        if (hasattr(self, 'planeWidget') == False):
            return

        if (self.planeWidget.GetEnabled() == False):
            return

        ### Check 4D Flow Magnitude and Velocity Data
        indexVel = -1
        indexMag = -1
        for i in range(self.parent.tableListMrStruct.rowCount()):
            itemCheck = self.parent.tableListMrStruct.item(i, 0)
            if (itemCheck.text() == "vel"):
                indexVel = i
            elif (itemCheck.text() == "mag"):
                indexMag = i

        if (indexVel < 0):
            return
        if (indexMag < 0):
            return

        ### Set Plane Source
        origin = self.planeWidget.GetOrigin()
        point1 = self.planeWidget.GetPoint1()
        point2 = self.planeWidget.GetPoint2()

        print(f"Origin {origin}, point1 {point1}, point2 {point2}")

        # Source Plane Definition
        planeSource = vtk.vtkPlaneSource()
        planeSource.SetOrigin(origin)
        planeSource.SetPoint1(point1)
        planeSource.SetPoint2(point2)
        planeSource.SetXResolution(10)
        planeSource.SetYResolution(10)

        ### Set
        flowTimePhase = 0

        itemVel = self.parent.tableListMrStruct.item(indexVel, 0)
        velMrStruct = itemVel.data(Qt.UserRole)

        itemMag = self.parent.tableListMrStruct.item(indexMag, 0)
        magMrStruct = itemMag.data(Qt.UserRole)

        self.streamActor, imageData = Utility.buildStreamActor(velMrStruct.dataAy[:, :, :, 0, flowTimePhase],
                                                               velMrStruct.dataAy[:, :, :, 1, flowTimePhase],
                                                               velMrStruct.dataAy[:, :, :, 2, flowTimePhase],
                                                               magMrStruct.dataAy[:, :, :, flowTimePhase],
                                                               velMrStruct.vox, planeSource)

        # Set input for volume renderer
        self.volRender.AddActor(self.streamActor)
        self.volRenWin.Render()

    @QtCore.pyqtSlot()
    def on_btnPathline_clicked(self):
        print("Pathline")

        if (hasattr(self, 'streamActor') == True):
            self.volRender.RemoveActor(self.streamActor)
            self.volRenWin.Render()
            del self.streamActor
            return

        ### Check PlaneWidget
        if (hasattr(self, 'planeWidget') == False):
            return

        if (self.planeWidget.GetEnabled() == False):
            return

        ### Check 4D Flow Magnitude and Velocity Data
        indexVel = -1
        indexMag = -1
        for i in range(self.parent.tableListMrStruct.rowCount()):
            itemCheck = self.parent.tableListMrStruct.item(i, 0)
            if (itemCheck.text() == "vel"):
                indexVel = i
            elif (itemCheck.text() == "mag"):
                indexMag = i

        if (indexVel < 0):
            return
        if (indexMag < 0):
            return

        ### Set Plane Source
        origin = self.planeWidget.GetOrigin()
        point1 = self.planeWidget.GetPoint1()
        point2 = self.planeWidget.GetPoint2()

        print(f"Origin {origin}, point1 {point1}, point2 {point2}")

        # Source Plane Definition
        planeSource = vtk.vtkPlaneSource()
        planeSource.SetOrigin(origin)
        planeSource.SetPoint1(point1)
        planeSource.SetPoint2(point2)
        planeSource.SetXResolution(10)
        planeSource.SetYResolution(10)

        ### Set
        flowTimePhase = 0

        itemVel = self.parent.tableListMrStruct.item(indexVel, 0)
        velMrStruct = itemVel.data(Qt.UserRole)

        itemMag = self.parent.tableListMrStruct.item(indexMag, 0)
        magMrStruct = itemMag.data(Qt.UserRole)

        self.streamActor, imageData = Utility.buildPathlineActor(velMrStruct.dataAy[:, :, :, 0, flowTimePhase],
                                                                 velMrStruct.dataAy[:, :, :, 1, flowTimePhase],
                                                                 velMrStruct.dataAy[:, :, :, 2, flowTimePhase],
                                                                 magMrStruct.dataAy[:, :, :, flowTimePhase],
                                                                 velMrStruct.vox, planeSource)

        # Set input for volume renderer
        self.volRender.AddActor(self.streamActor)
        self.volRenWin.Render()

    # @QtCore.pyqtSlot()
    def on_checkShowLargest_stateChanged(self, state):
        print("===check changed===")
        # if self.checkShowLargest.isChecked():
        if state == QtCore.Qt.Checked:
            print(f"checked {state}")
            self.smoother.SetInputConnection(self.connectivityFilter.GetOutputPort())
        else:
            print(f"unchecked {state}")
            self.smoother.SetInputConnection(self.marchingCubes.GetOutputPort())  # feed result to connectivity filter
        self.volRenWin.Render()


# Callback class: Reslice Cursor
class ObserverVtkResliceCursorCallback(object):
    imagePlaneWidget = [0] * 3  # vtkImagePlaneWidget
    resliceCursorWidget = [0] * 3  # vtkResliceCursorWidget

    def __init__(self):
        pass

    def __call__(self, caller, ev, *calldata):
        print(f"caller {caller.GetClassName()} ev {ev} calldata {calldata}")
        if (ev == "WindowLevelEvent") or (ev == "ResliceThicknessChangedEvent"):
            # Render everything
            for i in range(0, 3):
                self.resliceCursorWidget[i].Render()
            self.imagePlaneWidget[0].GetInteractor().GetRenderWindow().Render()
            return

        # assume caller is vtkImagePlaneWidget
        if (caller.GetClassName() == "vtkImagePlaneWidget"):
            wl = calldata
            if (caller == self.imagePlaneWidget[0]):
                self.imagePlaneWidget[1].SetWindowLevel(wl[0], wl[1], 1);
                self.imagePlaneWidget[2].SetWindowLevel(wl[0], wl[1], 1);
            elif (caller == self.imagePlaneWidget[1]):
                self.imagePlaneWidget[0].SetWindowLevel(wl[0], wl[1], 1);
                self.imagePlaneWidget[2].SetWindowLevel(wl[0], wl[1], 1);
            elif (caller == self.imagePlaneWidget[2]):
                self.imagePlaneWidget[0].SetWindowLevel(wl[0], wl[1], 1);
                self.imagePlaneWidget[1].SetWindowLevel(wl[0], wl[1], 1);

        # assume caller is vtkResliceCursorWidget
        if (caller.GetClassName() == "vtkResliceCursorWidget"):
            rep = caller.GetRepresentation()  # vtkResliceCursorLineRepresentation
            # Although the return value is not used, we keep the get calls in case they had side-effects
            rep.GetResliceCursorActor().GetCursorAlgorithm().GetResliceCursor()
            for i in range(3):
                ps = self.imagePlaneWidget[i].GetPolyDataAlgorithm()  # vtkPlaneSource
                ps.SetOrigin(self.resliceCursorWidget[i].GetResliceCursorRepresentation().GetPlaneSource().GetOrigin())
                ps.SetPoint1(self.resliceCursorWidget[i].GetResliceCursorRepresentation().GetPlaneSource().GetPoint1())
                ps.SetPoint2(self.resliceCursorWidget[i].GetResliceCursorRepresentation().GetPlaneSource().GetPoint2())

                # If the reslice plane has modified, update it on the 3D widget
                self.imagePlaneWidget[i].UpdatePlacement()

        # Render everything
        for i in range(3):
            self.resliceCursorWidget[i].Render()
        self.imagePlaneWidget[0].GetInteractor().GetRenderWindow().Render()
