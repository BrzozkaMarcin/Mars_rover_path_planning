from PyQt5.QtWidgets import *
from PyQt5 import QtCore, QtGui, QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib
import sys

from main import *
from App import Ui_MainWindow
from Generate_Map import *

matplotlib.use('QT5Agg')

TMatrix = None
HMatrix = None
FirstSolution = None
BestSolution = None


class PlotDisplayer(FigureCanvas):
    def __init__(self, plot_file, map=True, on_map=None):
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
        FigureCanvas.__init__(self, self.fig)
        FigureCanvas.setSizePolicy(
            self, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        self.figure.clear()

        ax = self.figure.add_subplot(111)
        if map:
            ax.pcolormesh(plot_file, cmap='jet')
            self.figure.colorbar(ax.pcolormesh(plot_file, cmap='jet'), fraction=0.1, pad=0.03)
            self.figure.axes[1].tick_params(axis="y", labelsize=8)
        else:
            if on_map is not None:
                x = [x[1] for x in on_map]
                y = [x[0] for x in on_map]
                ax.pcolormesh(plot_file, cmap='jet')
                ax.scatter(x, y, s=2, color='magenta')
            else:
                ax.pcolormesh(plot_file, cmap='Paired')
        ax.axis('off')
        ax.grid('on')


class GoalFunctionDisplayer(FigureCanvas):
    def __init__(self, plot1_file, plot2_file):
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
        FigureCanvas.__init__(self, self.fig)
        FigureCanvas.setSizePolicy(
            self, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        plt.title("Funkcja celu")
        ax.plot(plot1_file)
        ax.plot(plot2_file)
        ax.legend(['Najlepsze rozwiązanie', 'Obecne rozwiązanie'])


# self.M1.setStyleSheet("background-color: rgb(255, 255, 255);")
class MainWindow(QMainWindow, Ui_MainWindow):
    def __init__(self, *args, obj=None, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)
        self.setupUi(self)
        self.disableHTN()
        self.H_select.clicked.connect(lambda: self.show_solutions("H"))
        self.T_select.clicked.connect(lambda: self.show_solutions("T"))
        self.N_select.clicked.connect(lambda: self.show_solutions("N"))
        self.wczytaj_mapy.clicked.connect(self.load_maps)
        self.test1.clicked.connect(lambda: self.load_maps("test1"))
        self.test2.clicked.connect(lambda: self.load_maps("test2"))
        self.test3.clicked.connect(lambda: self.load_maps("test3"))
        self.zapisz_mapy.clicked.connect(self.save_maps)
        self.wygeneruj_mapy.clicked.connect(self.gener_maps)
        self.uruchom_algorytm.clicked.connect(self.run_algorithm)
        self.mapy.triggered.connect(self.show_maps)
        self.rozwiazania.triggered.connect(lambda: self.show_solutions("N"))
        self.funkcja_celu.triggered.connect(self.show_goal_function)
        self.M3.vbl = None
        self.GoalFunctionValues = []
        self.BestSolutionValues = []

    def plot_maps(self):
        global HMatrix
        global TMatrix
        if TMatrix is not None and TMatrix is not None:
            self.M1.canvas = PlotDisplayer(HMatrix)
            self.M1.vbl = QtWidgets.QVBoxLayout()
            self.M1.vbl.addWidget(self.M1.canvas)
            self.M1.setLayout(self.M1.vbl)

            self.M2.canvas = PlotDisplayer(TMatrix)
            self.M2.vbl = QtWidgets.QVBoxLayout()
            self.M2.vbl.addWidget(self.M2.canvas)
            self.M2.setLayout(self.M2.vbl)

    def plot_solutions(self, on_map="N"):
        global FirstSolution
        global BestSolution
        global HMatrix
        global TMatrix

        if FirstSolution is not None and BestSolution is not None:
            img = np.zeros(shape=HMatrix.shape)
            for P in FirstSolution.Solution:
                img[P.x, P.y] = 200
            if on_map == "N":
                self.M1.canvas = PlotDisplayer(img, False)
            elif on_map == "H":
                self.M1.canvas = PlotDisplayer(HMatrix, False, FirstSolution.Solution)
            elif on_map == "T":
                self.M1.canvas = PlotDisplayer(TMatrix, False, FirstSolution.Solution)
            self.M1.vbl = QtWidgets.QVBoxLayout()
            self.M1.vbl.addWidget(self.M1.canvas)
            self.M1.setLayout(self.M1.vbl)

            img = np.zeros(shape=HMatrix.shape)
            for P in BestSolution.Solution:
                img[P.x, P.y] = 200
            if on_map == "N":
                self.M2.canvas = PlotDisplayer(img, False)
            elif on_map == "H":
                self.M2.canvas = PlotDisplayer(HMatrix, False, BestSolution.Solution)
            elif on_map == "T":
                self.M2.canvas = PlotDisplayer(TMatrix, False, BestSolution.Solution)
            self.M2.vbl = QtWidgets.QVBoxLayout()
            self.M2.vbl.addWidget(self.M2.canvas)
            self.M2.setLayout(self.M2.vbl)

    def show_maps(self):
        self.enableM1M2()
        self.disableHTN()
        if TMatrix is not None and TMatrix is not None:
            QWidget().setLayout(self.M1.vbl)
            QWidget().setLayout(self.M2.vbl)
            self.plot_maps()
            self.M1_describe.setText("Mapa wysokości")
            self.M2_describe.setText("Mapa terenu")

    def show_solutions(self, on_map="N"):
        self.enableM1M2()
        self.enableHTN()
        if FirstSolution is not None and BestSolution is not None:
            QWidget().setLayout(self.M1.vbl)
            QWidget().setLayout(self.M2.vbl)
            self.plot_solutions(on_map)
            self.M1_describe.setText('FirstSolution, time = {:.3f}'.format(FirstSolution.time))
            self.M2_describe.setText('BestSolution, time = {:.3f}'.format(BestSolution.time))

    def openFileNameDialog(self, title_window):
        options = QFileDialog.Options()
        # options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self, title_window, "", "CSV Files (*.csv)", options=options)
        if fileName:
            return fileName

    def load_maps(self, test=None):
        global HMatrix
        global TMatrix
        # Removes previous layouts
        if TMatrix is not None and TMatrix is not None:
            QWidget().setLayout(self.M1.vbl)
            QWidget().setLayout(self.M2.vbl)
        if test == "test1":
            HMatrix_file = "test_maps/H_test1.csv"
            TMatrix_file = "test_maps/T_test1.csv"
        elif test == "test2":
            HMatrix_file = "test_maps/H_test2.csv"
            TMatrix_file = "test_maps/T_test2.csv"
        elif test == "test3":
            HMatrix_file = "test_maps/H_test3.csv"
            TMatrix_file = "test_maps/T_test3.csv"
        else:
            HMatrix_file = self.openFileNameDialog("Wczytaj mapę wysokości terenu")
            TMatrix_file = self.openFileNameDialog("Wczytaj mapę jakości terenu")
        HMatrix = np.loadtxt(HMatrix_file, dtype="int", delimiter=",")
        TMatrix = np.loadtxt(TMatrix_file, dtype="int", delimiter=",")
        self.plot_maps()
        self.disableHTN()
        self.size_map.setProperty("value", len(TMatrix))
        self.start_x.setMaximum(len(TMatrix) - 1)
        self.start_y.setMaximum(len(TMatrix) - 1)
        self.end_x.setMaximum(len(TMatrix) - 1)
        self.end_y.setMaximum(len(TMatrix) - 1)
        self.M1_describe.setText("Mapa wysokości")
        self.M2_describe.setText("Mapa terenu")

    def saveFileDialog(self, title_window):
        options = QFileDialog.Options()
        # options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getSaveFileName(self, title_window, "", "CSV Files (*.csv)", options=options)
        if fileName:
            return fileName

    def save_maps(self):
        global HMatrix
        global TMatrix
        if TMatrix is not None and TMatrix is not None:
            HMatrix_file = self.saveFileDialog("Zapisz mapę wysokości terenu")
            np.savetxt(HMatrix_file, HMatrix, fmt='%d', delimiter=",")
            TMatrix_file = self.saveFileDialog("Zapisz mapę wysokości terenu")
            np.savetxt(TMatrix_file, TMatrix, fmt='%d', delimiter=",")

    def gener_maps(self):
        global HMatrix
        global TMatrix
        # Removes previous layouts
        if TMatrix is not None and TMatrix is not None:
            QWidget().setLayout(self.M1.vbl)
            QWidget().setLayout(self.M2.vbl)
        HMatrix, TMatrix, _, _ = Generate_Matrices(self.size_map.value())
        self.plot_maps()
        self.start_x.setMaximum(len(TMatrix) - 1)
        self.start_y.setMaximum(len(TMatrix) - 1)
        self.end_x.setMaximum(len(TMatrix) - 1)
        self.end_y.setMaximum(len(TMatrix) - 1)
        self.M1_describe.setText("Mapa wysokości")
        self.M2_describe.setText("Mapa terenu")

    def run_algorithm(self):
        global HMatrix
        global TMatrix
        global FirstSolution
        global BestSolution
        if TMatrix is None or HMatrix is None:
            return
        self.size_map.setProperty("value", len(TMatrix))
        StartPoint = Point(self.start_y.value(), self.start_x.value())
        EndPoint = Point(self.end_y.value(), self.end_x.value())
        if StartPoint == EndPoint:
            self.message("point")
            return
        if not self.check_distance(StartPoint.x, StartPoint.y, EndPoint.x, EndPoint.y):
            self.message("2close")
            return
        MaxEnergy = self.max_energy.value()
        iteration = self.iter.value()
        HCE = self.HCE.value()
        TCE = self.TCE.value()
        HCT = self.HCT.value()
        TCT = self.TCT.value()
        sn1 = self.sn1.value()
        sn2 = self.sn2.value()
        aspiration = self.aspiration.value()
        random_point = self.random_point.value()
        kt1 = self.kt1.value()
        kt2 = self.kt2.value()
        neighbourhood_size = self.neighb_size.value()
        war_zakon = self.war_zakonczenia.value()
        NeighbourMethodList = []
        if self.s1.isChecked(): NeighbourMethodList.append(GenerateNeighbourSolution1)
        if self.s2.isChecked(): NeighbourMethodList.append(GenerateNeighbourSolution2)
        if self.s3.isChecked(): NeighbourMethodList.append(GenerateNeighbourSolution3)
        if self.s4.isChecked(): NeighbourMethodList.append(GenerateNeighbourSolution4)
        TabuProblem = TabuSearch(HMatrix, TMatrix, StartPoint, EndPoint, MaxEnergy, iteration,
                                 HCE, TCE, HCT, TCT, NeighbourMethodList, sn1, sn2,
                                 aspiration, random_point, neighbourhood_size, kt1, kt2, war_zakon)
        FS, BS = TabuProblem.Tabu_Search_Algorithm()
        self.BestSolutionValues = TabuProblem.BestSolutionValues
        self.GoalFunctionValues = TabuProblem.GoalFunctionValues
        if FS is None and BS is None:
            self.message("energy")
        else:
            FirstSolution, BestSolution = FS, BS
            self.show_solutions()

    def show_goal_function(self):
        self.enableM3()
        self.disableHTN()
        if self.M3.vbl is not None:
            QWidget().setLayout(self.M3.vbl)
        self.M3.canvas = GoalFunctionDisplayer(self.BestSolutionValues, self.GoalFunctionValues)
        self.M3.vbl = QtWidgets.QVBoxLayout()
        self.M3.vbl.addWidget(self.M3.canvas)
        self.M3.setLayout(self.M3.vbl)

    def enableM3(self):
        self.M1.close()
        self.M1_describe.close()
        self.M2.close()
        self.M2_describe.close()
        self.M3.show()

    def enableM1M2(self):
        self.M3.close()
        self.M1.show()
        self.M1_describe.show()
        self.M2.show()
        self.M2_describe.show()

    def enableHTN(self):
        self.H_select.show()
        self.H_select.show()
        self.T_select.show()
        self.T_select.show()
        self.N_select.show()
        self.N_select.show()

    def disableHTN(self):
        self.H_select.close()
        self.T_select.close()
        self.N_select.close()

    def message(self, mes):
        dlg = QMessageBox(self)
        dlg.setWindowTitle("Attention!")
        if mes == "energy":
            dlg.setText("Zwiększ energię początkową")
        elif mes == "2close":
            dlg.setText("Punkty są zbyt blisko siebie")
        else:
            dlg.setText("Punkty muszą być różne")
        dlg.setIcon(QMessageBox.Warning)
        button = dlg.exec()
        if button == QMessageBox.Ok:
            print("OK!")

    def check_distance(self, x0, y0, x1, y1):
        return np.sqrt((x0 - x1) ** 2 + (y0 - y1) ** 2) > 25


app = QApplication(sys.argv)
window = MainWindow()
window.show()

app.exec()
