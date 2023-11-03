# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'MainWindow_4DFlow.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_MainWindow_4DFlow(object):
    def setupUi(self, MainWindow_4DFlow):
        MainWindow_4DFlow.setObjectName("MainWindow_4DFlow")
        MainWindow_4DFlow.resize(800, 732)
        self.centralwidget = QtWidgets.QWidget(MainWindow_4DFlow)
        self.centralwidget.setObjectName("centralwidget")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.centralwidget)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.label_2 = QtWidgets.QLabel(self.centralwidget)
        self.label_2.setObjectName("label_2")
        self.verticalLayout.addWidget(self.label_2)
        self.btnClearTable = QtWidgets.QPushButton(self.centralwidget)
        self.btnClearTable.setObjectName("btnClearTable")
        self.verticalLayout.addWidget(self.btnClearTable)
        self.tableWidgetMagVel = QtWidgets.QTableWidget(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.tableWidgetMagVel.sizePolicy().hasHeightForWidth())
        self.tableWidgetMagVel.setSizePolicy(sizePolicy)
        self.tableWidgetMagVel.setMinimumSize(QtCore.QSize(320, 0))
        self.tableWidgetMagVel.setObjectName("tableWidgetMagVel")
        self.tableWidgetMagVel.setColumnCount(0)
        self.tableWidgetMagVel.setRowCount(0)
        self.verticalLayout.addWidget(self.tableWidgetMagVel)
        self.btnPC_MRA = QtWidgets.QPushButton(self.centralwidget)
        self.btnPC_MRA.setObjectName("btnPC_MRA")
        self.verticalLayout.addWidget(self.btnPC_MRA)
        self.btnPC_MRA_Time = QtWidgets.QPushButton(self.centralwidget)
        self.btnPC_MRA_Time.setObjectName("btnPC_MRA_Time")
        self.verticalLayout.addWidget(self.btnPC_MRA_Time)
        self.btnVelocityNorm = QtWidgets.QPushButton(self.centralwidget)
        self.btnVelocityNorm.setObjectName("btnVelocityNorm")
        self.verticalLayout.addWidget(self.btnVelocityNorm)
        self.btnVelocitySum = QtWidgets.QPushButton(self.centralwidget)
        self.btnVelocitySum.setObjectName("btnVelocitySum")
        self.verticalLayout.addWidget(self.btnVelocitySum)
        self.btnPeakVelocity = QtWidgets.QPushButton(self.centralwidget)
        self.btnPeakVelocity.setObjectName("btnPeakVelocity")
        self.verticalLayout.addWidget(self.btnPeakVelocity)
        self.btnShowData = QtWidgets.QPushButton(self.centralwidget)
        self.btnShowData.setObjectName("btnShowData")
        self.verticalLayout.addWidget(self.btnShowData)
        self.tableWidgetMask = QtWidgets.QTableWidget(self.centralwidget)
        self.tableWidgetMask.setObjectName("tableWidgetMask")
        self.tableWidgetMask.setColumnCount(0)
        self.tableWidgetMask.setRowCount(0)
        self.verticalLayout.addWidget(self.tableWidgetMask)
        self.btnMaskAddActor = QtWidgets.QPushButton(self.centralwidget)
        self.btnMaskAddActor.setObjectName("btnMaskAddActor")
        self.verticalLayout.addWidget(self.btnMaskAddActor)
        self.btnMaskCutData = QtWidgets.QPushButton(self.centralwidget)
        self.btnMaskCutData.setObjectName("btnMaskCutData")
        self.verticalLayout.addWidget(self.btnMaskCutData)
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem)
        self.horizontalLayout.addLayout(self.verticalLayout)
        self.mdiArea = QtWidgets.QMdiArea(self.centralwidget)
        self.mdiArea.setObjectName("mdiArea")
        self.horizontalLayout.addWidget(self.mdiArea)
        MainWindow_4DFlow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow_4DFlow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 21))
        self.menubar.setObjectName("menubar")
        self.menuFile = QtWidgets.QMenu(self.menubar)
        self.menuFile.setObjectName("menuFile")
        self.menuSegmentation = QtWidgets.QMenu(self.menubar)
        self.menuSegmentation.setEnabled(False)
        self.menuSegmentation.setObjectName("menuSegmentation")
        self.menuCreate = QtWidgets.QMenu(self.menuSegmentation)
        self.menuCreate.setEnabled(False)
        self.menuCreate.setObjectName("menuCreate")
        self.menuEdit_Segmentation = QtWidgets.QMenu(self.menuSegmentation)
        self.menuEdit_Segmentation.setObjectName("menuEdit_Segmentation")
        self.menuAdd = QtWidgets.QMenu(self.menuEdit_Segmentation)
        self.menuAdd.setObjectName("menuAdd")
        self.menuVisualization = QtWidgets.QMenu(self.menubar)
        self.menuVisualization.setEnabled(False)
        self.menuVisualization.setObjectName("menuVisualization")
        self.menuVelocity_Vectors = QtWidgets.QMenu(self.menuVisualization)
        self.menuVelocity_Vectors.setObjectName("menuVelocity_Vectors")
        self.menuStreamlines = QtWidgets.QMenu(self.menuVisualization)
        self.menuStreamlines.setObjectName("menuStreamlines")
        self.menuQuantification = QtWidgets.QMenu(self.menubar)
        self.menuQuantification.setEnabled(False)
        self.menuQuantification.setObjectName("menuQuantification")
        self.menuCreate_Plane = QtWidgets.QMenu(self.menuQuantification)
        self.menuCreate_Plane.setObjectName("menuCreate_Plane")
        self.menuStasis_Maps = QtWidgets.QMenu(self.menuQuantification)
        self.menuStasis_Maps.setObjectName("menuStasis_Maps")
        self.menuWindow = QtWidgets.QMenu(self.menubar)
        self.menuWindow.setObjectName("menuWindow")
        MainWindow_4DFlow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow_4DFlow)
        self.statusbar.setObjectName("statusbar")
        MainWindow_4DFlow.setStatusBar(self.statusbar)
        self.toolBar = QtWidgets.QToolBar(MainWindow_4DFlow)
        self.toolBar.setObjectName("toolBar")
        MainWindow_4DFlow.addToolBar(QtCore.Qt.TopToolBarArea, self.toolBar)
        self.actionLoad_Dicoms = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionLoad_Dicoms.setObjectName("actionLoad_Dicoms")
        self.actionLoad_MrStruct = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionLoad_MrStruct.setObjectName("actionLoad_MrStruct")
        self.actionLoad_Segmentation = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionLoad_Segmentation.setObjectName("actionLoad_Segmentation")
        self.actionSave_Segmentation = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionSave_Segmentation.setEnabled(False)
        self.actionSave_Segmentation.setObjectName("actionSave_Segmentation")
        self.actionSave_Segmentation_with_info = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionSave_Segmentation_with_info.setEnabled(False)
        self.actionSave_Segmentation_with_info.setObjectName("actionSave_Segmentation_with_info")
        self.actionSave_3D_Movie = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionSave_3D_Movie.setEnabled(False)
        self.actionSave_3D_Movie.setObjectName("actionSave_3D_Movie")
        self.actionView_Segmentation = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionView_Segmentation.setObjectName("actionView_Segmentation")
        self.actionMask_Velocity = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionMask_Velocity.setObjectName("actionMask_Velocity")
        self.actionCenterline = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionCenterline.setObjectName("actionCenterline")
        self.actionVelocity_Cloud = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionVelocity_Cloud.setObjectName("actionVelocity_Cloud")
        self.actionVelocity_MIPS = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionVelocity_MIPS.setObjectName("actionVelocity_MIPS")
        self.action3D_Viewer_Movie = QtWidgets.QAction(MainWindow_4DFlow)
        self.action3D_Viewer_Movie.setObjectName("action3D_Viewer_Movie")
        self.actionDiameters_Centerline = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionDiameters_Centerline.setObjectName("actionDiameters_Centerline")
        self.actionDiameter_Projection = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionDiameter_Projection.setObjectName("actionDiameter_Projection")
        self.actionFlow_Plane_s = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionFlow_Plane_s.setObjectName("actionFlow_Plane_s")
        self.actionNormalized_Histogram = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionNormalized_Histogram.setObjectName("actionNormalized_Histogram")
        self.action3D_WSS = QtWidgets.QAction(MainWindow_4DFlow)
        self.action3D_WSS.setObjectName("action3D_WSS")
        self.actionKinetic_Energy = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionKinetic_Energy.setObjectName("actionKinetic_Energy")
        self.actionViscous_Energy_Loss = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionViscous_Energy_Loss.setObjectName("actionViscous_Energy_Loss")
        self.actionVorticity = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionVorticity.setObjectName("actionVorticity")
        self.actionFull_Cycle = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionFull_Cycle.setObjectName("actionFull_Cycle")
        self.actionSystole = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionSystole.setObjectName("actionSystole")
        self.actionDiastole = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionDiastole.setObjectName("actionDiastole")
        self.actionTime_Range = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionTime_Range.setObjectName("actionTime_Range")
        self.actionTime_Point = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionTime_Point.setObjectName("actionTime_Point")
        self.actionUse_PC_MRA = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionUse_PC_MRA.setObjectName("actionUse_PC_MRA")
        self.actionDynamic_Mask = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionDynamic_Mask.setObjectName("actionDynamic_Mask")
        self.actionSeed_Points = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionSeed_Points.setObjectName("actionSeed_Points")
        self.actionFast_Marching = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionFast_Marching.setObjectName("actionFast_Marching")
        self.actionErase_in_slice = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionErase_in_slice.setObjectName("actionErase_in_slice")
        self.actionErase_with_MIPs = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionErase_with_MIPs.setObjectName("actionErase_with_MIPs")
        self.actionErode_Volume = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionErode_Volume.setObjectName("actionErode_Volume")
        self.actionSmooth_Volume = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionSmooth_Volume.setObjectName("actionSmooth_Volume")
        self.actionLevel_Set = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionLevel_Set.setObjectName("actionLevel_Set")
        self.actionContraction = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionContraction.setObjectName("actionContraction")
        self.actionExpansion = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionExpansion.setObjectName("actionExpansion")
        self.actionFree_Hand = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionFree_Hand.setObjectName("actionFree_Hand")
        self.actionVectors_3D = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionVectors_3D.setObjectName("actionVectors_3D")
        self.actionVectors_2D = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionVectors_2D.setObjectName("actionVectors_2D")
        self.actionStreamlines_3D = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionStreamlines_3D.setObjectName("actionStreamlines_3D")
        self.actionStreamlines_2D = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionStreamlines_2D.setObjectName("actionStreamlines_2D")
        self.actionSingle_Plane = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionSingle_Plane.setObjectName("actionSingle_Plane")
        self.actionPlanes_at_Landmarks = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionPlanes_at_Landmarks.setObjectName("actionPlanes_at_Landmarks")
        self.actionStasis_Map_MIP = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionStasis_Map_MIP.setObjectName("actionStasis_Map_MIP")
        self.action3D_Stasis_Map = QtWidgets.QAction(MainWindow_4DFlow)
        self.action3D_Stasis_Map.setObjectName("action3D_Stasis_Map")
        self.actionCascade = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionCascade.setObjectName("actionCascade")
        self.actionTiled = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionTiled.setObjectName("actionTiled")
        self.actionBatch_work_1 = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionBatch_work_1.setObjectName("actionBatch_work_1")
        self.actionBatch_work_2 = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionBatch_work_2.setObjectName("actionBatch_work_2")
        self.actionLoad_Directory = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionLoad_Directory.setObjectName("actionLoad_Directory")
        self.actionClose_All = QtWidgets.QAction(MainWindow_4DFlow)
        self.actionClose_All.setObjectName("actionClose_All")
        self.menuFile.addAction(self.actionLoad_Dicoms)
        self.menuFile.addAction(self.actionLoad_MrStruct)
        self.menuFile.addAction(self.actionLoad_Segmentation)
        self.menuFile.addAction(self.actionLoad_Directory)
        self.menuFile.addAction(self.actionSave_Segmentation)
        self.menuFile.addAction(self.actionSave_Segmentation_with_info)
        self.menuFile.addAction(self.actionSave_3D_Movie)
        self.menuFile.addAction(self.actionBatch_work_1)
        self.menuFile.addAction(self.actionBatch_work_2)
        self.menuCreate.addAction(self.actionFull_Cycle)
        self.menuCreate.addAction(self.actionSystole)
        self.menuCreate.addAction(self.actionDiastole)
        self.menuCreate.addAction(self.actionTime_Range)
        self.menuCreate.addAction(self.actionTime_Point)
        self.menuCreate.addAction(self.actionUse_PC_MRA)
        self.menuCreate.addAction(self.actionDynamic_Mask)
        self.menuCreate.addAction(self.actionSeed_Points)
        self.menuCreate.addAction(self.actionFast_Marching)
        self.menuAdd.addAction(self.actionContraction)
        self.menuAdd.addAction(self.actionExpansion)
        self.menuAdd.addAction(self.actionFree_Hand)
        self.menuEdit_Segmentation.addAction(self.actionErase_in_slice)
        self.menuEdit_Segmentation.addAction(self.actionErase_with_MIPs)
        self.menuEdit_Segmentation.addAction(self.menuAdd.menuAction())
        self.menuEdit_Segmentation.addAction(self.actionErode_Volume)
        self.menuEdit_Segmentation.addAction(self.actionSmooth_Volume)
        self.menuEdit_Segmentation.addAction(self.actionLevel_Set)
        self.menuSegmentation.addAction(self.menuCreate.menuAction())
        self.menuSegmentation.addAction(self.actionView_Segmentation)
        self.menuSegmentation.addAction(self.menuEdit_Segmentation.menuAction())
        self.menuSegmentation.addAction(self.actionMask_Velocity)
        self.menuSegmentation.addAction(self.actionCenterline)
        self.menuVelocity_Vectors.addAction(self.actionVectors_3D)
        self.menuVelocity_Vectors.addAction(self.actionVectors_2D)
        self.menuStreamlines.addAction(self.actionStreamlines_3D)
        self.menuStreamlines.addAction(self.actionStreamlines_2D)
        self.menuVisualization.addAction(self.actionVelocity_Cloud)
        self.menuVisualization.addAction(self.menuVelocity_Vectors.menuAction())
        self.menuVisualization.addAction(self.menuStreamlines.menuAction())
        self.menuVisualization.addAction(self.actionVelocity_MIPS)
        self.menuVisualization.addAction(self.action3D_Viewer_Movie)
        self.menuCreate_Plane.addAction(self.actionSingle_Plane)
        self.menuCreate_Plane.addAction(self.actionPlanes_at_Landmarks)
        self.menuStasis_Maps.addAction(self.actionStasis_Map_MIP)
        self.menuStasis_Maps.addAction(self.action3D_Stasis_Map)
        self.menuQuantification.addAction(self.actionDiameters_Centerline)
        self.menuQuantification.addAction(self.actionDiameter_Projection)
        self.menuQuantification.addAction(self.menuCreate_Plane.menuAction())
        self.menuQuantification.addAction(self.actionFlow_Plane_s)
        self.menuQuantification.addAction(self.actionNormalized_Histogram)
        self.menuQuantification.addAction(self.action3D_WSS)
        self.menuQuantification.addAction(self.actionKinetic_Energy)
        self.menuQuantification.addAction(self.actionViscous_Energy_Loss)
        self.menuQuantification.addAction(self.actionVorticity)
        self.menuQuantification.addAction(self.menuStasis_Maps.menuAction())
        self.menuWindow.addAction(self.actionCascade)
        self.menuWindow.addAction(self.actionTiled)
        self.menuWindow.addAction(self.actionClose_All)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuSegmentation.menuAction())
        self.menubar.addAction(self.menuVisualization.menuAction())
        self.menubar.addAction(self.menuQuantification.menuAction())
        self.menubar.addAction(self.menuWindow.menuAction())
        self.toolBar.addAction(self.actionLoad_Directory)
        self.toolBar.addAction(self.actionLoad_MrStruct)
        self.toolBar.addAction(self.actionLoad_Segmentation)
        self.toolBar.addAction(self.actionBatch_work_1)
        self.toolBar.addAction(self.actionBatch_work_2)

        self.retranslateUi(MainWindow_4DFlow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow_4DFlow)

    def retranslateUi(self, MainWindow_4DFlow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow_4DFlow.setWindowTitle(_translate("MainWindow_4DFlow", "MainWindow"))
        self.label_2.setText(_translate("MainWindow_4DFlow", "List of MrStruct Loaded"))
        self.btnClearTable.setText(_translate("MainWindow_4DFlow", "Clear"))
        self.btnPC_MRA.setText(_translate("MainWindow_4DFlow", "Generate 3D PC-MRA"))
        self.btnPC_MRA_Time.setText(_translate("MainWindow_4DFlow", "Generate Time-Series PC-MRA"))
        self.btnVelocityNorm.setText(_translate("MainWindow_4DFlow", "Generate 4D Velocity Norm from Vector"))
        self.btnVelocitySum.setText(_translate("MainWindow_4DFlow", "Generate Velocity Sum thorugh the cycle"))
        self.btnPeakVelocity.setText(_translate("MainWindow_4DFlow", "Generate Peak Velocity through the cycle"))
        self.btnShowData.setText(_translate("MainWindow_4DFlow", "Show Data"))
        self.btnMaskAddActor.setText(_translate("MainWindow_4DFlow", "Add Actor"))
        self.btnMaskCutData.setText(_translate("MainWindow_4DFlow", "Cut Data"))
        self.menuFile.setTitle(_translate("MainWindow_4DFlow", "File"))
        self.menuSegmentation.setTitle(_translate("MainWindow_4DFlow", "Segmentation"))
        self.menuCreate.setTitle(_translate("MainWindow_4DFlow", "Create"))
        self.menuEdit_Segmentation.setTitle(_translate("MainWindow_4DFlow", "Edit Segmentation"))
        self.menuAdd.setTitle(_translate("MainWindow_4DFlow", "Add"))
        self.menuVisualization.setTitle(_translate("MainWindow_4DFlow", "Visualization"))
        self.menuVelocity_Vectors.setTitle(_translate("MainWindow_4DFlow", "Velocity Vectors"))
        self.menuStreamlines.setTitle(_translate("MainWindow_4DFlow", "Streamlines"))
        self.menuQuantification.setTitle(_translate("MainWindow_4DFlow", "Quantification"))
        self.menuCreate_Plane.setTitle(_translate("MainWindow_4DFlow", "Create Plane"))
        self.menuStasis_Maps.setTitle(_translate("MainWindow_4DFlow", "Stasis Maps"))
        self.menuWindow.setTitle(_translate("MainWindow_4DFlow", "Window"))
        self.toolBar.setWindowTitle(_translate("MainWindow_4DFlow", "toolBar"))
        self.actionLoad_Dicoms.setText(_translate("MainWindow_4DFlow", "Load Dicoms"))
        self.actionLoad_MrStruct.setText(_translate("MainWindow_4DFlow", "Load MrStructs"))
        self.actionLoad_Segmentation.setText(_translate("MainWindow_4DFlow", "Load Segmentation"))
        self.actionSave_Segmentation.setText(_translate("MainWindow_4DFlow", "Save Segmentation"))
        self.actionSave_Segmentation_with_info.setText(_translate("MainWindow_4DFlow", "Save Segmentation with info"))
        self.actionSave_3D_Movie.setText(_translate("MainWindow_4DFlow", "Save 3D Movie"))
        self.actionView_Segmentation.setText(_translate("MainWindow_4DFlow", "View Segmentation"))
        self.actionMask_Velocity.setText(_translate("MainWindow_4DFlow", "Mask Velocity"))
        self.actionCenterline.setText(_translate("MainWindow_4DFlow", "Centerline"))
        self.actionVelocity_Cloud.setText(_translate("MainWindow_4DFlow", "Velocity Cloud"))
        self.actionVelocity_MIPS.setText(_translate("MainWindow_4DFlow", "Velocity MIPS"))
        self.action3D_Viewer_Movie.setText(_translate("MainWindow_4DFlow", "3D Viewer Movie"))
        self.actionDiameters_Centerline.setText(_translate("MainWindow_4DFlow", "Diameters Centerline"))
        self.actionDiameter_Projection.setText(_translate("MainWindow_4DFlow", "Diameter Projection"))
        self.actionFlow_Plane_s.setText(_translate("MainWindow_4DFlow", "Flow Plane(s)"))
        self.actionNormalized_Histogram.setText(_translate("MainWindow_4DFlow", "Normalized Histogram"))
        self.action3D_WSS.setText(_translate("MainWindow_4DFlow", "3D WSS"))
        self.actionKinetic_Energy.setText(_translate("MainWindow_4DFlow", "Kinetic Energy"))
        self.actionViscous_Energy_Loss.setText(_translate("MainWindow_4DFlow", "Viscous Energy Loss"))
        self.actionVorticity.setText(_translate("MainWindow_4DFlow", "Vorticity"))
        self.actionFull_Cycle.setText(_translate("MainWindow_4DFlow", "Full Cycle"))
        self.actionSystole.setText(_translate("MainWindow_4DFlow", "Systole"))
        self.actionDiastole.setText(_translate("MainWindow_4DFlow", "Diastole"))
        self.actionTime_Range.setText(_translate("MainWindow_4DFlow", "Time Range"))
        self.actionTime_Point.setText(_translate("MainWindow_4DFlow", "Time Point"))
        self.actionUse_PC_MRA.setText(_translate("MainWindow_4DFlow", "Use PC-MRA"))
        self.actionDynamic_Mask.setText(_translate("MainWindow_4DFlow", "Dynamic Mask"))
        self.actionSeed_Points.setText(_translate("MainWindow_4DFlow", "Seed Points"))
        self.actionFast_Marching.setText(_translate("MainWindow_4DFlow", "Fast Marching"))
        self.actionErase_in_slice.setText(_translate("MainWindow_4DFlow", "Erase in slice"))
        self.actionErase_with_MIPs.setText(_translate("MainWindow_4DFlow", "Erase with MIPs"))
        self.actionErode_Volume.setText(_translate("MainWindow_4DFlow", "Erode Volume"))
        self.actionSmooth_Volume.setText(_translate("MainWindow_4DFlow", "Smooth Volume"))
        self.actionLevel_Set.setText(_translate("MainWindow_4DFlow", "Level Set"))
        self.actionContraction.setText(_translate("MainWindow_4DFlow", "Contraction"))
        self.actionExpansion.setText(_translate("MainWindow_4DFlow", "Expansion"))
        self.actionFree_Hand.setText(_translate("MainWindow_4DFlow", "Free Hand"))
        self.actionVectors_3D.setText(_translate("MainWindow_4DFlow", "Vectors 3D"))
        self.actionVectors_2D.setText(_translate("MainWindow_4DFlow", "Vectors 2D"))
        self.actionStreamlines_3D.setText(_translate("MainWindow_4DFlow", "Streamlines 3D"))
        self.actionStreamlines_2D.setText(_translate("MainWindow_4DFlow", "Streamlines 2D"))
        self.actionSingle_Plane.setText(_translate("MainWindow_4DFlow", "Single Plane"))
        self.actionPlanes_at_Landmarks.setText(_translate("MainWindow_4DFlow", "Planes at Landmarks"))
        self.actionStasis_Map_MIP.setText(_translate("MainWindow_4DFlow", "Stasis Map MIP"))
        self.action3D_Stasis_Map.setText(_translate("MainWindow_4DFlow", "3D Stasis Map"))
        self.actionCascade.setText(_translate("MainWindow_4DFlow", "Cascade"))
        self.actionTiled.setText(_translate("MainWindow_4DFlow", "Tiled"))
        self.actionBatch_work_1.setText(_translate("MainWindow_4DFlow", "Batch work 1"))
        self.actionBatch_work_2.setText(_translate("MainWindow_4DFlow", "Batch work 2"))
        self.actionLoad_Directory.setText(_translate("MainWindow_4DFlow", "Load Directory"))
        self.actionClose_All.setText(_translate("MainWindow_4DFlow", "Close All"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow_4DFlow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow_4DFlow()
    ui.setupUi(MainWindow_4DFlow)
    MainWindow_4DFlow.show()
    sys.exit(app.exec_())

