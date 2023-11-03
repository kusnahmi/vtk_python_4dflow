# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'QtVTKRenderWindows.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_QtVTKRenderWindows(object):
    def setupUi(self, QtVTKRenderWindows):
        QtVTKRenderWindows.setObjectName("QtVTKRenderWindows")
        QtVTKRenderWindows.resize(851, 583)
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(":/Icons/help.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        QtVTKRenderWindows.setWindowIcon(icon)
        QtVTKRenderWindows.setIconSize(QtCore.QSize(22, 22))
        self.centralwidget = QtWidgets.QWidget(QtVTKRenderWindows)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayoutWidget = QtWidgets.QWidget(self.centralwidget)
        self.gridLayoutWidget.setGeometry(QtCore.QRect(10, 30, 591, 531))
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout_2.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.view2 = QVTKRenderWindowInteractor(self.gridLayoutWidget)
        self.view2.setObjectName("view2")
        self.gridLayout_2.addWidget(self.view2, 1, 0, 1, 1)
        self.view4 = QVTKRenderWindowInteractor(self.gridLayoutWidget)
        self.view4.setObjectName("view4")
        self.gridLayout_2.addWidget(self.view4, 0, 1, 1, 1)
        self.view3 = QVTKRenderWindowInteractor(self.gridLayoutWidget)
        self.view3.setObjectName("view3")
        self.gridLayout_2.addWidget(self.view3, 1, 1, 1, 1)
        self.view1 = QVTKRenderWindowInteractor(self.gridLayoutWidget)
        self.view1.setObjectName("view1")
        self.gridLayout_2.addWidget(self.view1, 0, 0, 1, 1)
        self.frame = QtWidgets.QFrame(self.centralwidget)
        self.frame.setGeometry(QtCore.QRect(629, 39, 201, 521))
        self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame.setObjectName("frame")
        self.resliceModeCheckBox = QtWidgets.QCheckBox(self.frame)
        self.resliceModeCheckBox.setGeometry(QtCore.QRect(30, 20, 141, 17))
        self.resliceModeCheckBox.setObjectName("resliceModeCheckBox")
        self.thickModeCheckBox = QtWidgets.QCheckBox(self.frame)
        self.thickModeCheckBox.setGeometry(QtCore.QRect(30, 80, 101, 17))
        self.thickModeCheckBox.setObjectName("thickModeCheckBox")
        self.blendModeGroupBox = QtWidgets.QGroupBox(self.frame)
        self.blendModeGroupBox.setGeometry(QtCore.QRect(10, 120, 181, 91))
        self.blendModeGroupBox.setObjectName("blendModeGroupBox")
        self.radioButton_Min = QtWidgets.QRadioButton(self.blendModeGroupBox)
        self.radioButton_Min.setGeometry(QtCore.QRect(10, 20, 161, 17))
        self.radioButton_Min.setObjectName("radioButton_Min")
        self.radioButton_Max = QtWidgets.QRadioButton(self.blendModeGroupBox)
        self.radioButton_Max.setGeometry(QtCore.QRect(10, 40, 161, 17))
        self.radioButton_Max.setObjectName("radioButton_Max")
        self.radioButton_Mean = QtWidgets.QRadioButton(self.blendModeGroupBox)
        self.radioButton_Mean.setGeometry(QtCore.QRect(10, 60, 111, 17))
        self.radioButton_Mean.setObjectName("radioButton_Mean")
        self.resetButton = QtWidgets.QPushButton(self.frame)
        self.resetButton.setGeometry(QtCore.QRect(10, 220, 51, 21))
        self.resetButton.setObjectName("resetButton")
        self.frame_2 = QtWidgets.QFrame(self.frame)
        self.frame_2.setGeometry(QtCore.QRect(0, 250, 191, 211))
        self.frame_2.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_2.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_2.setObjectName("frame_2")
        self.AddDistance1Button = QtWidgets.QPushButton(self.frame_2)
        self.AddDistance1Button.setGeometry(QtCore.QRect(0, 0, 131, 21))
        self.AddDistance1Button.setObjectName("AddDistance1Button")
        self.btnSyncCursor = QtWidgets.QPushButton(self.frame_2)
        self.btnSyncCursor.setGeometry(QtCore.QRect(0, 40, 131, 21))
        self.btnSyncCursor.setObjectName("btnSyncCursor")
        self.btnTest1 = QtWidgets.QPushButton(self.frame_2)
        self.btnTest1.setGeometry(QtCore.QRect(20, 90, 75, 23))
        self.btnTest1.setObjectName("btnTest1")
        self.btnTest2 = QtWidgets.QPushButton(self.frame_2)
        self.btnTest2.setGeometry(QtCore.QRect(20, 120, 75, 23))
        self.btnTest2.setObjectName("btnTest2")
        self.btnTest3 = QtWidgets.QPushButton(self.frame_2)
        self.btnTest3.setGeometry(QtCore.QRect(20, 150, 75, 23))
        self.btnTest3.setObjectName("btnTest3")
        QtVTKRenderWindows.setCentralWidget(self.centralwidget)
        self.actionOpenFile = QtWidgets.QAction(QtVTKRenderWindows)
        self.actionOpenFile.setEnabled(True)
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap(":/Icons/fileopen.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionOpenFile.setIcon(icon1)
        self.actionOpenFile.setObjectName("actionOpenFile")
        self.actionExit = QtWidgets.QAction(QtVTKRenderWindows)
        icon2 = QtGui.QIcon()
        icon2.addPixmap(QtGui.QPixmap("."), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionExit.setIcon(icon2)
        self.actionExit.setObjectName("actionExit")
        self.actionPrint = QtWidgets.QAction(QtVTKRenderWindows)
        icon3 = QtGui.QIcon()
        icon3.addPixmap(QtGui.QPixmap(":/Icons/print.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionPrint.setIcon(icon3)
        self.actionPrint.setObjectName("actionPrint")
        self.actionHelp = QtWidgets.QAction(QtVTKRenderWindows)
        self.actionHelp.setIcon(icon)
        self.actionHelp.setObjectName("actionHelp")
        self.actionSave = QtWidgets.QAction(QtVTKRenderWindows)
        icon4 = QtGui.QIcon()
        icon4.addPixmap(QtGui.QPixmap(":/Icons/filesave.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.actionSave.setIcon(icon4)
        self.actionSave.setObjectName("actionSave")

        self.retranslateUi(QtVTKRenderWindows)
        QtCore.QMetaObject.connectSlotsByName(QtVTKRenderWindows)

    def retranslateUi(self, QtVTKRenderWindows):
        _translate = QtCore.QCoreApplication.translate
        QtVTKRenderWindows.setWindowTitle(_translate("QtVTKRenderWindows", "QtVTKRenderWindows"))
        self.resliceModeCheckBox.setText(_translate("QtVTKRenderWindows", "Oblique Reslice"))
        self.thickModeCheckBox.setText(_translate("QtVTKRenderWindows", "Thick Slab"))
        self.blendModeGroupBox.setTitle(_translate("QtVTKRenderWindows", "Blend mode"))
        self.radioButton_Min.setText(_translate("QtVTKRenderWindows", "Min Intensity Blend"))
        self.radioButton_Max.setText(_translate("QtVTKRenderWindows", "Max Intensity Blend"))
        self.radioButton_Mean.setText(_translate("QtVTKRenderWindows", "Mean Blend"))
        self.resetButton.setText(_translate("QtVTKRenderWindows", "Reset"))
        self.AddDistance1Button.setText(_translate("QtVTKRenderWindows", "Add Distance On View 1"))
        self.btnSyncCursor.setText(_translate("QtVTKRenderWindows", "Sync Cursor"))
        self.btnTest1.setText(_translate("QtVTKRenderWindows", "Test 1"))
        self.btnTest2.setText(_translate("QtVTKRenderWindows", "Test 2"))
        self.btnTest3.setText(_translate("QtVTKRenderWindows", "Test 3"))
        self.actionOpenFile.setText(_translate("QtVTKRenderWindows", "Open File..."))
        self.actionExit.setText(_translate("QtVTKRenderWindows", "Exit"))
        self.actionPrint.setText(_translate("QtVTKRenderWindows", "Print"))
        self.actionHelp.setText(_translate("QtVTKRenderWindows", "Help"))
        self.actionSave.setText(_translate("QtVTKRenderWindows", "Save"))

from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    QtVTKRenderWindows = QtWidgets.QMainWindow()
    ui = Ui_QtVTKRenderWindows()
    ui.setupUi(QtVTKRenderWindows)
    QtVTKRenderWindows.show()
    sys.exit(app.exec_())

