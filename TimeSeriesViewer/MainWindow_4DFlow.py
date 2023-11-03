# Compile UI and resource file
# pyuic5 -x MainWindow_4DFlow.ui -o Ui_MainWindow_4DFlow.py
# pyuic5 -x Widget_TimeSeriesView.ui -o Ui_Widget_TimeSeriesView.py

from PyQt5 import QtWidgets, QtCore
from PyQt5.QtWidgets import QFileDialog, QTableWidgetItem, QAbstractItemView, QCheckBox
from PyQt5.QtCore import Qt, QFileInfo, QTime, QDate
from Ui_MainWindow_4DFlow import Ui_MainWindow_4DFlow
from Widget_TimeSeriesView import Widget_TimeSeriesView

import os
import os.path
import copy
import numpy as np
import Utility
import FourDFlow


## MainWindow Class
#
#  Qt Main window of MDI interface.
#  Load source data and generate analysis data.
#  Create sub windows for visualization.
class MainWindow_4DFlow(QtWidgets.QMainWindow, Ui_MainWindow_4DFlow):
    ## Constructor
    #
    #  Setup UI.
    def __init__(self, basePath="."):
        ### Attribute
        # self.basePath = "C:/d/target/Controls/001A00473420160826/3dpc/mrstruct_subject_20160826_user011"
        self.basePath = basePath

        super().__init__()
        self.setupUi(self)

        ### Additional UI Setup
        # Left pane: table of data list
        self.tableWidgetMagVel.setColumnCount(4)
        self.tableWidgetMagVel.setHorizontalHeaderLabels(('Type', 'PatientID', 'Vector', 'Time'))
        self.tableWidgetMagVel.setColumnWidth(0, 60)
        self.tableWidgetMagVel.setColumnWidth(1, 130)
        self.tableWidgetMagVel.setColumnWidth(2, 50)
        self.tableWidgetMagVel.setColumnWidth(3, 50)
        self.tableWidgetMagVel.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tableWidgetMagVel.doubleClicked.connect(self.on_tableWidgetMagVel_doubleclicked)
        # Left pane: table of mask list
        self.tableWidgetMask.setColumnCount(4)
        self.tableWidgetMask.setHorizontalHeaderLabels(('filename', 'PatientID', 'Vector', 'Time'))
        self.tableWidgetMask.setColumnWidth(0, 60)
        self.tableWidgetMask.setColumnWidth(1, 130)
        self.tableWidgetMask.setColumnWidth(2, 50)
        self.tableWidgetMask.setColumnWidth(3, 50)
        self.tableWidgetMask.setSelectionBehavior(QAbstractItemView.SelectRows)
        # menu callback
        self.actionLoad_Dicoms.triggered.connect(self.fileLoadDicoms)
        self.actionLoad_MrStruct.triggered.connect(self.fileLoadMrStruct)
        self.actionLoad_Segmentation.triggered.connect(self.fileLoadSegmentation)
        self.actionCascade.triggered.connect(self.mdiArea.cascadeSubWindows)
        self.actionTiled.triggered.connect(self.mdiArea.tileSubWindows)

        ### Show UI
        self.resize(1200, 800)
        self.show()

        # Batch job.
        self.fileLoadMrStructMagVel(self.basePath)
        self.on_btnPC_MRA_Time_clicked()
        self.fileLoadMrStructMask(os.path.join(self.basePath, "LA_volume.mat"))
        self.on_btnShowData_clicked()

    # MENU CALLBACK: File - Load MrStruct
    def fileLoadMrStruct(self):
        filename = QFileDialog.getOpenFileName(self, 'Open Magnitude File',
                                               os.path.join(self.basePath, "mag_struct.mat"), "MATLAB Data (*.mat)")
        if filename[0]:
            # self.fileLoadMrStructMagVelSpeed(os.path.dirname(filename[0]))
            self.fileLoadMrStructMagVel(os.path.dirname(filename[0]))

    def fileLoadMrStructMagVel(self, basePath):
        filepath_mag = os.path.join(basePath, "mag_struct.mat")
        filepath_vel = os.path.join(basePath, "vel_struct.mat")

        if (os.path.exists(filepath_mag) == False) or (os.path.exists(filepath_vel) == False):
            return

        # Open magnitude and velocity data
        mrStruct_mag = FourDFlow.MrStruct(filepath_mag)
        mrStruct_vel = FourDFlow.MrStruct(filepath_vel)

        # Extract Patient ID
        bLoop = True
        patientID = ""
        path = basePath
        while (bLoop):
            path, tail = os.path.split(path)
            print(tail)
            if (len(tail) == 18):
                patientID = tail
                bLoop = False

        if (patientID != ""):
            print(f"patientID: {patientID}")

        self.addItemOnTableListMr("mag", True, patientID, mrStruct_mag)
        self.addItemOnTableListMr("vel", True, patientID, mrStruct_vel)

    ## Build subwindow from mrStruct and show it.
    def createSubWindowTimeSeriesView(self, mrStructScalar, mrStructVector):
        widget_TimeSeriesView = Widget_TimeSeriesView(mrStructScalar, mrStructVector, self)
        sub = self.mdiArea.addSubWindow(widget_TimeSeriesView)
        sub.setWindowTitle(f"MrStruct: {mrStructScalar.name}")
        sub.show()
        return sub

    ## Add mrStruct as a row in the table
    #
    #  column 0: Row select + type string + data
    #  column 1: name (patient ID)
    #  column 2: if vector or scalar?
    #  column 3: if time-series data?
    def addItemOnTableListMr(self, type, selected, name, mrStruct):
        idx = self.tableWidgetMagVel.rowCount()
        self.tableWidgetMagVel.insertRow(idx)

        #  column 0: Row select + type string + data
        checkBoxItem = QTableWidgetItem();
        if (selected == True):
            checkBoxItem.setCheckState(Qt.Checked);
        else:
            checkBoxItem.setCheckState(Qt.Unchecked);
        checkBoxItem.setText(type)
        checkBoxItem.setData(Qt.UserRole, mrStruct)
        self.tableWidgetMagVel.setItem(idx, 0, checkBoxItem);

        #  column 1: name (patient ID)
        self.tableWidgetMagVel.setItem(idx, 1, QTableWidgetItem(name))

        #  column 2: if vector or scalar?
        checkBoxVector = QCheckBox()
        checkBoxVector.setEnabled(False)
        if (mrStruct.GetIfVector()):
            checkBoxVector.setChecked(Qt.Checked)
            self.tableWidgetMagVel.setCellWidget(idx, 2, checkBoxVector)
        else:
            checkBoxVector.setChecked(Qt.Unchecked)
            self.tableWidgetMagVel.setCellWidget(idx, 2, checkBoxVector)

        #  column 3: if time-series data?
        checkBoxTimeseries = QCheckBox()
        checkBoxTimeseries.setEnabled(False)
        if (mrStruct.GetIfTimeSeries()):
            checkBoxTimeseries.setChecked(Qt.Checked)
            self.tableWidgetMagVel.setCellWidget(idx, 3, checkBoxTimeseries)
        else:
            checkBoxTimeseries.setChecked(Qt.Unchecked)
            self.tableWidgetMagVel.setCellWidget(idx, 3, checkBoxTimeseries)

    ## Add mrStruct as a row in the table
    #
    #  column 0: Row select + filename string + data
    #  column 1: name (patient ID)
    #  column 2: if vector or scalar?
    #  column 3: if time-series data?
    def addItemOnTableListMask(self, filepath, selected, name, mrStruct):
        idx = self.tableWidgetMask.rowCount()
        self.tableWidgetMask.insertRow(idx)

        #  column 0: Row select + type string + data
        checkBoxItem = QTableWidgetItem();
        if (selected == True):
            checkBoxItem.setCheckState(Qt.Checked);
        else:
            checkBoxItem.setCheckState(Qt.Unchecked);

        path, tail = os.path.split(filepath)
        checkBoxItem.setText(tail)
        checkBoxItem.setData(Qt.UserRole, mrStruct)
        self.tableWidgetMask.setItem(idx, 0, checkBoxItem);

        #  column 1: name (patient ID)
        self.tableWidgetMask.setItem(idx, 1, QTableWidgetItem(name))

        #  column 2: if vector or scalar?
        checkBoxVector = QCheckBox()
        checkBoxVector.setEnabled(False)
        if (mrStruct.GetIfVector()):
            checkBoxVector.setChecked(Qt.Checked)
            self.tableWidgetMask.setCellWidget(idx, 2, checkBoxVector)
        else:
            checkBoxVector.setChecked(Qt.Unchecked)
            self.tableWidgetMask.setCellWidget(idx, 2, checkBoxVector)

        #  column 3: if time-series data?
        checkBoxTimeseries = QCheckBox()
        checkBoxTimeseries.setEnabled(False)
        if (mrStruct.GetIfTimeSeries()):
            checkBoxTimeseries.setChecked(Qt.Checked)
            self.tableWidgetMask.setCellWidget(idx, 3, checkBoxTimeseries)
        else:
            checkBoxTimeseries.setChecked(Qt.Unchecked)
            self.tableWidgetMask.setCellWidget(idx, 3, checkBoxTimeseries)

    ## MENU CALLBACK: File - Load Segmentation
    def fileLoadSegmentation(self):
        filename = QFileDialog.getOpenFileName(self, 'Open Magnitude File',
                                               os.path.join(self.basePath, "mag_struct.mat"), "MATLAB Data (*.mat)")
        if filename[0]:
            self.fileLoadMrStructMask(filename[0])

    ## Load Mask in MrStruct. Add to the MrStruct table.
    def fileLoadMrStructMask(self, filepath):
        # Open mask data
        mrSegmentation = FourDFlow.MrStruct(filepath)

        # Extract Patient ID
        bLoop = True
        patientID = ""
        path = filepath
        while (bLoop):
            path, tail = os.path.split(path)
            print(tail)
            if (len(tail) == 18):
                patientID = tail
                bLoop = False

        if (patientID != ""):
            print(f"patientID: {patientID}")

        self.addItemOnTableListMask(filepath, True, patientID, mrSegmentation)

    # ===============================================================================================================
    # ===============================================================================================================
    # ===============================================================================================================
    # ===============================================================================================================
    # ===============================================================================================================

    def on_tableWidgetMagVel_doubleclicked(self):
        for currentQTableWidgetItem in self.tableWidgetMagVel.selectedItems():
            print(currentQTableWidgetItem.row(), currentQTableWidgetItem.column(), currentQTableWidgetItem.text())

    def fileLoadDicoms(self):
        pass

    def fileLoadMrStructMagVelSpeed(self, basePath):
        filepath_mag = os.path.join(basePath, "mag_struct.mat")
        filepath_vel = os.path.join(basePath, "vel_struct.mat")
        filepath_velSpeed = os.path.join(basePath, "vel_struct_speed.mat")

        if (os.path.exists(filepath_mag) == False) or (os.path.exists(filepath_vel) == False):
            return

        # Open magnitude and velocity data
        mrStruct_mag = FourDFlow.MrStruct(filepath_mag)
        mrStruct_vel = FourDFlow.MrStruct(filepath_vel)
        # Open speed data if exist. If not, generate one and save it for next time use.
        if (os.path.exists(filepath_velSpeed)):
            mrStruct_speed = FourDFlow.MrStruct(filepath_velSpeed)
        else:
            mrStruct_vel.Smooth_Velocity()
            mrStruct_speed = mrStruct_vel.generateSpeedArray(mrStruct_mag, "correlation")
            mrStruct_speed.exportMat(filepath_velSpeed)

        # Build subwindow for each data and add to left-pane table.
        self.createSubWindow2DSlice(mrStruct_mag)
        self.addItemOnTableListMr("mag", False, mrStruct_mag.name, mrStruct_mag)

        self.createSubWindow2DSlice(mrStruct_vel)
        self.addItemOnTableListMr("vel", False, mrStruct_vel.name, mrStruct_vel)

        self.createSubWindow2DSlice(mrStruct_speed)
        self.addItemOnTableListMr("speed", True, mrStruct_speed.name, mrStruct_speed)

    ## Build subwindow from mrStruct and show it.
    def createSubWindow2DSlice(self, mrStruct):
        widget_TimeSeriesView = Widget_TimeSeriesView(mrStruct, self)
        sub = self.mdiArea.addSubWindow(widget_TimeSeriesView)
        sub.setWindowTitle(f"MrStruct: {mrStruct.name}")
        sub.show()
        return sub

    @QtCore.pyqtSlot()
    def on_actionLoad_Directory_triggered(self):
        print("load directory")
        basePath = str(QFileDialog.getExistingDirectory(self, "Select Data Directory", self.basePath))
        if (basePath == ""):
            return
        self.fileLoadMrStructMagVelSpeed(basePath)

        # load segmentations
        filename_LA = "LA_volume.mat"
        filename_LV = "LV_volume.mat"
        filename_LA_nPVLAA = "LA_nPVLAA_volume.mat"

        filepath_LA = os.path.join(basePath, filename_LA)
        filepath_LV = os.path.join(basePath, filename_LV)
        filepath_LA_nPVLAA = os.path.join(basePath, filename_LA_nPVLAA)

        if (os.path.exists(filepath_LA)):
            mr_LA = FourDFlow.MrStruct(filepath_LA)
            self.createSubWindow2DSlice(mr_LA)
        if (os.path.exists(filepath_LV)):
            mr_LV = FourDFlow.MrStruct(filepath_LV)
            self.createSubWindow2DSlice(mr_LV)
        if (os.path.exists(filepath_LA_nPVLAA)):
            mr_LA_nPVLAA = FourDFlow.MrStruct(filepath_LA_nPVLAA)
            self.createSubWindow2DSlice(mr_LA_nPVLAA)

        if (os.path.exists(filepath_LA) and os.path.exists(filepath_LA_nPVLAA)):
            mr_calc_PVLAA = copy.deepcopy(mr_LA_nPVLAA)
            mr_calc_PVLAA.name = "calc: PVLAA = LA - nPVLAA"
            mr_calc_PVLAA.dataAy = ((mr_LA.dataAy.astype(np.int8) - mr_LA_nPVLAA.dataAy.astype(np.int8)) > 0) * 1
            self.createSubWindow2DSlice(mr_calc_PVLAA)

    # BUTTON CALLBACK: Left pane - Clear Table
    # remove all rows but keep column header
    @QtCore.pyqtSlot()
    def on_btnClearTable_clicked(self):
        print("Clear Table")
        self.tableWidgetMagVel.clearContents();
        self.tableWidgetMagVel.setRowCount(0);

    # MENU CALLBACK: File - Batch work 1
    # Visit every sub directory and calculate kinetic energy with mask
    @QtCore.pyqtSlot()
    def on_actionBatch_work_1_triggered(self):
        print("Batch work 1")
        basePath = str(QFileDialog.getExistingDirectory(self, "Select Root Directory for Batch Work"))
        if (basePath == ""):
            return

        filename_vel = "vel_struct.mat"
        filename_lvmask = "LV_volume.mat"
        filename_mask = ["LA_volume.mat", "LA_LV_volume.mat", "Ao_volume.mat"]

        # Report Header
        strArrayHeader = Utility.convertArrayToString(np.arange(30))
        print(
            f"root maskname KEavg     KEpeakSyst KEpeakDias  KEavgSyst  KEavgDias t_Syst t_Dias t_ED t_ES numPhase {strArrayHeader}")

        # Iterate all sub directories
        for root, dirs, files in os.walk(basePath):
            if ((filename_vel in files) and (filename_lvmask in files)):
                # Read velocity data and segmentation
                filepath_vel = os.path.join(root, filename_vel)
                filepath_lvmask = os.path.join(root, filename_lvmask)
                mrStruct_vel = FourDFlow.MrStruct(filepath_vel)
                mrStruct_lvMask = FourDFlow.MrStruct(filepath_lvmask)

                # Calculate Kinetic Energy
                unitVolume, numOfCells, volume, KE = Utility.calcKineticEnergy(mrStruct_vel, mrStruct_lvMask.dataAy)
                # Calculate Statistics
                #     time avearge KE
                #     peak systole
                #     peak diastole
                #     systole: from maximum peak, trace back and forth
                #     diastole: from second peak, trace back and forth
                t_pSyst, t_pDias, t_ED, t_ES, KEavg, KEpeakSystole, KEpeakDiastole, KEavgSystole, KEavgDiastole = Utility.calcKineticEnergyStatistics(
                    KE[:, 1])

                # Report LV mask
                strArrayValue = Utility.convertArrayToString(KE[:, 1])
                print(
                    f"{root} {filename_lvmask} {KEavg:>10.4f} {KEpeakSystole:>10.4f} {KEpeakDiastole:>10.4f} {KEavgSystole:>10.4f} {KEavgDiastole:>10.4f} {t_pSyst:>6} {t_pDias:>6} {t_ED:>4} {t_ES:>4} {len(KE[:, 1]):>8} {strArrayValue}")

                # Process other masks with extracted timing parameters
                for maskname in filename_mask:
                    if (maskname in files):
                        # Read Segmentation
                        filepath_mask = os.path.join(root, maskname)
                        mrMask = FourDFlow.MrStruct(filepath_mask)
                        # Calculate Kinetic Energy
                        unitVolume, numOfCells, volume, KE = Utility.calcKineticEnergy(mrStruct_vel, mrMask.dataAy)
                        # Calculate Statistics
                        _, _, _, _, KEavg, KEpeakSystole, KEpeakDiastole, KEavgSystole, KEavgDiastole = Utility.calcKineticEnergyStatistics(
                            KE[:, 1], t_pSyst, t_pDias, t_ED, t_ES)
                        # Report
                        strArrayValue = Utility.convertArrayToString(KE[:, 1])
                        print(
                            f"{root} {maskname} {KEavg:>10.4f} {KEpeakSystole:>10.4f} {KEpeakDiastole:>10.4f} {KEavgSystole:>10.4f} {KEavgDiastole:>10.4f} {t_pSyst:>6} {t_pDias:>6} {t_ED:>4} {t_ES:>4} {len(KE[:, 1]):>8} {strArrayValue}")

        print("=" * 80)

    # MENU CALLBACK: File - Batch work 1
    # Visit every sub directory and calculate kinetic energy with mask
    @QtCore.pyqtSlot()
    def on_actionBatch_work_2_triggered(self):
        print("Batch work 2")
        basePath = str(QFileDialog.getExistingDirectory(self, "Select Root Directory for Batch Work"))
        if (basePath == ""):
            return

        filename_vel = "vel_struct.mat"
        filename_lvmask = "LV_volume.mat"
        filename_mask = ["laa.mat", "LA_volume.mat", "LA_nPVLAA_volume.mat"]

        # Report Header
        strArrayHeader = Utility.convertArrayToString(np.arange(30))
        print(
            f"root maskname cells volume   KEavg     KEpeakSyst KEpeakDias  KEavgSyst  KEavgDias t_Syst t_Dias t_ED t_ES numPhase {strArrayHeader}")

        # Iterate all sub directories
        for root, dirs, files in os.walk(basePath):
            if ((filename_vel in files) and (filename_lvmask in files) and (filename_mask[0] in files)):
                # Read velocity data and segmentation
                filepath_vel = os.path.join(root, filename_vel)
                filepath_lvmask = os.path.join(root, filename_lvmask)
                mrStruct_vel = MrStruct.MrStruct(filepath_vel)
                mrStruct_lvMask = MrStruct.MrStruct(filepath_lvmask)

                # Calculate Kinetic Energy
                unitVolume, numOfCells, volume, KE = Utility.calcKineticEnergy(mrStruct_vel, mrStruct_lvMask.dataAy)
                # Calculate Statistics
                #     time avearge KE
                #     peak systole
                #     peak diastole
                #     systole: from maximum peak, trace back and forth
                #     diastole: from second peak, trace back and forth
                t_pSyst, t_pDias, t_ED, t_ES, KEavg, KEpeakSystole, KEpeakDiastole, KEavgSystole, KEavgDiastole = Utility.calcKineticEnergyStatistics(
                    KE[:, 1])

                # Report LV mask
                strArrayValue = Utility.convertArrayToString(KE[:, 1])
                print(
                    f"{root} {filename_lvmask} {numOfCells:>5} {volume:>8.1f} {KEavg:>10.4f} {KEpeakSystole:>10.4f} {KEpeakDiastole:>10.4f} {KEavgSystole:>10.4f} {KEavgDiastole:>10.4f} {t_pSyst:>6} {t_pDias:>6} {t_ED:>4} {t_ES:>4} {len(KE[:, 1]):>8} {strArrayValue}")

                # Process other masks with extracted timing parameters
                for maskname in filename_mask:
                    if (maskname in files):
                        # Read Segmentation
                        filepath_mask = os.path.join(root, maskname)
                        mrMask = MrStruct.MrStruct(filepath_mask)
                        # Calculate Kinetic Energy
                        unitVolume, numOfCells, volume, KE = Utility.calcKineticEnergy(mrStruct_vel, mrMask.dataAy)
                        # Calculate Statistics
                        _, _, _, _, KEavg, KEpeakSystole, KEpeakDiastole, KEavgSystole, KEavgDiastole = Utility.calcKineticEnergyStatistics(
                            KE[:, 1], t_pSyst, t_pDias, t_ED, t_ES)
                        # Report
                        strArrayValue = Utility.convertArrayToString(KE[:, 1])
                        print(
                            f"{root} {maskname} {numOfCells:>5} {volume:>8.1f} {KEavg:>10.4f} {KEpeakSystole:>10.4f} {KEpeakDiastole:>10.4f} {KEavgSystole:>10.4f} {KEavgDiastole:>10.4f} {t_pSyst:>6} {t_pDias:>6} {t_ED:>4} {t_ES:>4} {len(KE[:, 1]):>8} {strArrayValue}")

        print("=" * 80)

    @QtCore.pyqtSlot()
    def on_actionClose_All_triggered(self):
        print("Close All")

        self.mdiArea.closeAllSubWindows()

        self.tableWidgetMagVel.clearContents();
        self.tableWidgetMagVel.setRowCount(0);

    ## BUTTON CALLBACK: Left pane - Generate 3D PC-MRA
    #
    #  Generate 3D Time-Averaged PC-MRA
    @QtCore.pyqtSlot()
    def on_btnPC_MRA_clicked(self):
        print("Generate 3D Time-Averaged PC-MRA")

    ## BUTTON CALLBACK: Left pane - Generate Time-Series PC-MRA
    #
    #  Generate 4D Time-Series PC-MRA
    @QtCore.pyqtSlot()
    def on_btnPC_MRA_Time_clicked(self):
        gamma = 0.4

        print("Generate 4D Time-Series PC-MRA")
        mrStructScalar, mrStructVector = self.RetrieveSelectedDataFromTable()

        # Copy structure of existing scalar array
        mrStruct4DPcMra = copy.deepcopy(mrStructScalar)

        # Calculate velocity norm
        npNormV = np.sqrt(np.sum(mrStructVector.dataAy ** 2, axis=3))
        npMag = mrStructScalar.dataAy
        npMra = np.multiply(npMag, np.power(npNormV, gamma))

        # Replace the data array
        mrStruct4DPcMra.dataAy = npMra

        # Post on the list
        self.addItemOnTableListMr("normV", True, "NORM", mrStruct4DPcMra)

    ## BUTTON CALLBACK: Left pane - Generate Velocity Norm Scalar Data
    #
    #  Generate 4D Time-Series PC-MRA
    @QtCore.pyqtSlot()
    def on_btnVelocityNorm_clicked(self):
        print("Generate Velocity Norm Scalar Data")
        mrStructScalar, mrStructVector = self.RetrieveSelectedDataFromTable()

        # Copy structure of existing scalar array
        mrStructVelocityNorm = copy.deepcopy(mrStructScalar)

        # Calculate velocity norm
        npNormV = np.sqrt(np.sum(mrStructVector.dataAy ** 2, axis=3))

        # Replace the data array
        mrStructVelocityNorm.dataAy = npNormV

        # Post on the list
        self.addItemOnTableListMr("normV", True, "NORM", mrStructVelocityNorm)

    ## BUTTON CALLBACK: Left pane - Generate Velocity Sum Vector Time-Averaged Data
    #
    #  Generate 4D Time-Series PC-MRA
    @QtCore.pyqtSlot()
    def on_btnVelocitySum_clicked(self):
        print("Generate Velocity Sum Vector Time-Averaged Data")
        mrStructScalar, mrStructVector = self.RetrieveSelectedDataFromTable()

        # Copy structure of existing scalar array
        mrStructVSumNorm = copy.deepcopy(mrStructScalar)
        mrStructVSum = copy.deepcopy(mrStructVector)

        # Calculate velocity sum through the cardiac cycle
        npVSum = np.sum(mrStructVector.dataAy, axis=4)

        # Calculate norm of velocity sum
        npVSumNorm = np.sqrt(np.sum(npVSum ** 2, axis=3))

        # Replace the data array
        mrStructVSumNorm.dataAy = np.repeat(npVSumNorm[:, :, :, np.newaxis], mrStructScalar.dataAy.shape[3], axis=3)
        mrStructVSum.dataAy = np.repeat(npVSum[:, :, :, :, np.newaxis], mrStructVector.dataAy.shape[4], axis=4)

        # Post on the list
        self.addItemOnTableListMr("|VSum|", True, "|VSum|", mrStructVSumNorm)
        self.addItemOnTableListMr("VSum", True, "VSum", mrStructVSum)

    ## BUTTON CALLBACK: Left pane - Generate Peak Velocity Time-Averaged Data
    #
    @QtCore.pyqtSlot()
    def on_btnPeakVelocity_clicked(self):
        print("Generate Peak Veloicty")
        mrStructScalar, mrStructVector = self.RetrieveSelectedDataFromTable()

        # Copy structure of existing scalar array
        mrStructPeakSpeed = copy.deepcopy(mrStructScalar)
        mrStructPeakVelocity = copy.deepcopy(mrStructVector)

        # Calculate speed for each time frame
        npSpeed = np.sqrt(np.sum(mrStructVector.dataAy ** 2, axis=3))
        indexMax = np.expand_dims(np.argmax(npSpeed, axis=3), axis=3)
        npSpeedPeak = np.squeeze(np.take_along_axis(npSpeed, indexMax, axis=3), axis=3)
        print(f"index {indexMax.shape}, array {npSpeedPeak.shape}")

        print(f"index {indexMax.shape}, array {mrStructVector.dataAy[:, :, :, 0, :].shape}")

        x = np.take_along_axis(mrStructVector.dataAy[:, :, :, 0, :], indexMax, axis=3)
        y = np.take_along_axis(mrStructVector.dataAy[:, :, :, 1, :], indexMax, axis=3)
        z = np.take_along_axis(mrStructVector.dataAy[:, :, :, 2, :], indexMax, axis=3)
        npVpeak = np.squeeze(np.stack([x, y, z], axis=3), axis=4)
        print(f"x {x.shape}, npVpeak {npVpeak.shape}")
        # npVPeak = mrStructVector.dataAy[indexMax]

        # Replace the data array
        mrStructPeakSpeed.dataAy = np.repeat(npSpeedPeak[:, :, :, np.newaxis], mrStructScalar.dataAy.shape[3], axis=3)
        mrStructPeakVelocity.dataAy = np.repeat(npVpeak[:, :, :, :, np.newaxis], mrStructVector.dataAy.shape[4], axis=4)
        print(f"speed {mrStructPeakSpeed.dataAy.shape}, velocity {mrStructPeakVelocity.dataAy.shape}")

        # Post on the list
        self.addItemOnTableListMr("PeakS", True, "PeakS", mrStructPeakSpeed)
        self.addItemOnTableListMr("PeakV", True, "PeakV", mrStructPeakVelocity)

    ## BUTTON CALLBACK: Left pane - Show Data
    #
    #  Visualize selected data
    @QtCore.pyqtSlot()
    def on_btnShowData_clicked(self):
        print("Show Data")
        mrStructScalar, mrStructVector = self.RetrieveSelectedDataFromTable()
        self.createSubWindowTimeSeriesView(mrStructScalar, mrStructVector)

    ## Retrieve selected scalar and vector data from MrStruct table.
    def RetrieveSelectedDataFromTable(self):
        mrStructScalar = None
        mrStructVector = None
        for idx in range(self.tableWidgetMagVel.rowCount()):
            item = self.tableWidgetMagVel.item(idx, 0)
            if (item.checkState() == Qt.Checked):
                if (self.tableWidgetMagVel.cellWidget(idx, 2).checkState() == Qt.Checked):
                    mrStructVector = item.data(Qt.UserRole)
                else:
                    mrStructScalar = item.data(Qt.UserRole)

        return mrStructScalar, mrStructVector

    ## Retrieve selected scalar and vector data from MrStruct table.
    def RetrieveSelectedMaskFromTable(self):
        mrStructMask = None
        for idx in range(self.tableWidgetMask.rowCount()):
            item = self.tableWidgetMask.item(idx, 0)
            if (item.checkState() == Qt.Checked):
                mrStructMask = item.data(Qt.UserRole)

        return mrStructMask

    ## BUTTON CALLBACK: Left pane - Mask - Add Actor
    #
    #  Generate Mask Actor and add to the active sub window
    @QtCore.pyqtSlot()
    def on_btnMaskAddActor_clicked(self):
        print("MASK: Add Actor")
        mrStructMask = self.RetrieveSelectedMaskFromTable()

        subwin = self.mdiArea.activeSubWindow()
        subwin.widget().AddMaskActor(mrStructMask)

    ## BUTTON CALLBACK: Left pane - Mask - Cut Data
    #
    #  Cut Data and regenerate visualization
    @QtCore.pyqtSlot()
    def on_btnMaskCutData_clicked(self):
        print("MASK: Cut Data")
        mrStructMask = self.RetrieveSelectedMaskFromTable()

        subwin = self.mdiArea.activeSubWindow()
        subwin.widget().MaskCutData(mrStructMask)
