import sys
from PyQt5.QtWidgets import QApplication
from MainWindow_4DFlow import MainWindow_4DFlow

app = QApplication(sys.argv)

if (len(sys.argv)>=2):
    mainWindow = MainWindow_4DFlow(sys.argv[1])
else:
    mainWindow = MainWindow_4DFlow()
sys.exit(app.exec())
