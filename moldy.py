from scipy.stats import gmean
from periodictable import elements
import periodictableGUI
import itertools
import sys
from os.path import expanduser
from PyQt4.QtCore import *
from PyQt4.Qt import QApplication, QWidget, QTableView, QStandardItem, QStandardItemModel, QColor, QPushButton, QFileDialog, QHBoxLayout, QVBoxLayout
import pyqtgraph.opengl as gl
from zmat import ZMError
from utils import *

# create Qt application
qt_app = QApplication(sys.argv)

# create main widget
class MainWidget(QWidget):
    def __init__(self):
        QWidget.__init__(self)

        # define periodic table widget for element selection
        self.periodicTableWidget = periodictableGUI.PeriodicTableDialog()

        # initial molecule Zmatrix (can be empty)
        self.inp = [[]]
        #self.inp = [['H'],
        #['O', 1, 0.9],
        #['O', 2, 1.4, 1, 105.],
        #['H', 3, 0.9, 2, 105., 1, 120.]]

        # define & initialize model that will contain Zmatrix data
        #if self.inp == [[]]:
            #self.inp = [[self.periodicTable()[1]]]
        self.model = QStandardItemModel(len(self.inp), 7, self)
        self.inputField = QTableView(self)
        self.inputField.setModel(self.model)
        self.inputField.setMinimumWidth(325)
        self.model.setHorizontalHeaderLabels(['atom','','bond','','angle','','dihedral'])
        for j, width in enumerate([40, 22, 65, 22, 65, 22, 65]):
            self.inputField.setColumnWidth(j, width)

        # populate the model
        self.populateModel()

        # define GUI buttons and their actions
        self.readZmatButton = QPushButton('Read Zmat from file', self)
        self.readZmatButton.clicked.connect(self.readZmat)
        self.readXYZButton = QPushButton('Read XYZ from file', self)
        self.readXYZButton.clicked.connect(self.readXYZ)
        self.addRowButton = QPushButton('Add row', self)
        self.addRowButton.clicked.connect(self.addRow)
        self.deleteRowButton = QPushButton('Delete row', self)
        self.deleteRowButton.clicked.connect(self.deleteRow)
        #self.periodicTableButton = QPushButton('periodicTable', self)
        #self.periodicTableButton.clicked.connect(self.periodicTable)
        self.writeZmatButton = QPushButton('Write Zmat to file', self)
        self.writeZmatButton.clicked.connect(self.writeZmat)
        self.writeXYZButton = QPushButton('Write XYZ to file', self)
        self.writeXYZButton.clicked.connect(self.writeXYZ)

        # define GL widget that displays the 3D molecule model
        self.window = gl.GLViewWidget()
        self.window.installEventFilter(self)
        self.window.setMinimumSize(500, 500)
        self.updateView()

        # define save file dialog
        self.fileDialog = QFileDialog(self)

        # define application layout
        self.layout = QHBoxLayout(self)
        self.left = QVBoxLayout()
        self.left.addWidget(self.inputField)
        self.leftBot = QHBoxLayout()
        self.leftBot2 = QHBoxLayout()
        self.leftBot3 = QHBoxLayout()
        self.leftBot.addWidget(self.addRowButton)
        self.leftBot.addWidget(self.deleteRowButton)
        #self.leftBot.addWidget(self.periodicTableButton)
        self.leftBot2.addWidget(self.readZmatButton)
        self.leftBot2.addWidget(self.readXYZButton)
        self.leftBot3.addWidget(self.writeZmatButton)
        self.leftBot3.addWidget(self.writeXYZButton)
        self.left.addLayout(self.leftBot)
        self.left.addLayout(self.leftBot2)
        self.left.addLayout(self.leftBot3)
        self.layout.addLayout(self.left)
        self.layout.addWidget(self.window)

        self.adjustSize()
        self.setWindowTitle('Moldy')

        # start monitoring changes in the model
        self.model.dataChanged.connect(self.updateView)
        self.buildList = []

    # run and show the application
    def run(self):
        self.show()
        qt_app.exec_()

    # fill the model with initial data from 'self.inp'
    def populateModel(self):
        self.model.removeRows(0, self.model.rowCount())
        for i, row in enumerate(self.inp):
            for j, cell in enumerate(row):
                item = QStandardItem(str(cell))
                self.model.setItem(i, j, item)
        # some cells should not be editable, they are disabled
        for i in range(min(len(self.inp), 3)):
            for j in range(2*i+1, 7):
                self.model.setItem(i, j, QStandardItem())
                self.model.item(i, j).setBackground(QColor(150,150,150))
                self.model.item(i, j).setFlags(Qt.ItemIsEnabled)

    # add a row to the bottom of the model
    def addRow(self):
        # temporarily stop updating the GL window
        self.model.dataChanged.disconnect(self.updateView)
        row = self.model.rowCount()
        self.model.insertRow(row)
        # some cells should not be editable
        if row < 3:
            for j in range(2*row+1, 7):
                self.model.setItem(row, j, QStandardItem())
                self.model.item(row, j).setBackground(QColor(150,150,150))
                self.model.item(row, j).setFlags(Qt.ItemIsEnabled)
        # restart GL window updating
        self.model.dataChanged.connect(self.updateView)

    # delete the last row of the model
    def deleteRow(self):
        self.model.removeRow(self.model.rowCount()-1)
        self.updateView()

    # show the periodic table widget
    def periodicTable(self):
        self.periodicTableWidget.exec_()
        selection = self.periodicTableWidget.selection()
        return selection

    # import molecule with zmatrix coordinates
    def readZmat(self):
        self.model.dataChanged.disconnect(self.updateView)
        filename = self.fileDialog.getOpenFileName(self, 'Open file', expanduser('~'), '*.zmat;;*.*')
        self.inp = []
        self.populateModel()
        if filename:
            with open(filename, 'r') as f:
                next(f)
                next(f)
                for row in f:
                    self.inp.append(row.split())
                f.close()
        self.populateModel()
        self.model.dataChanged.connect(self.updateView)
        self.updateView()

    # import molecule with xyz coordinates
    def readXYZ(self):
        self.model.dataChanged.disconnect(self.updateView)
        filename = self.fileDialog.getOpenFileName(self, 'Open file', expanduser('~'), '*.xyz;;*.*')
        xyz = []
        elems = []
        self.inp = []
        self.populateModel()
        if filename:
            with open(filename, 'r') as f:
                next(f)
                next(f)
                for row in f:
                    rs = row.split()
                    if len(rs) == 4:
                        elems.append(rs[0])
                        xyz.append([float(f) for f in rs[1:]])
                f.close()
            self.inp = xyz2zmat(xyz, elems)
        self.populateModel()
        self.model.dataChanged.connect(self.updateView)
        self.updateView()

    # export Zmatrix to csv
    def writeZmat(self):
        zm = model2list(self.model)
        filename = self.fileDialog.getSaveFileName(self, 'Save file', expanduser('~')+'/'+getFormula(list(list(zip(*zm))[0]))+'.zmat', '*.zmat;;*.*')
        try:
            filename
        except NameError:
            pass
        else:
            if filename:
                writeOutput(zm, filename)

    # export XYZ coordinates to csv
    def writeXYZ(self):
        xyz = []
        zm = model2list(self.model)
        for i in range(len(v)):
            xyz.append(np.round(v[i], 7).tolist())
            xyz[i][:0] = zm[i][0]
        if len(v) > 0:
            formula = getFormula(list(list(zip(*xyz))[0]))
        else:
            formula = 'moldy_output'
        filename = self.fileDialog.getSaveFileName(self, 'Save file', expanduser('~')+'/'+formula+'.xyz', '*.xyz;;*.*')
        try:
            filename
        except NameError:
            pass
        else:
            if filename:
                writeOutput(xyz, filename)

    # redraw the 3D molecule in GL widget
    def updateView(self):
        global r
        global c
        global v
        global vs
        global elems
        global nelems
        data = model2list(self.model)
        try:
            # create a list with element coordinates
            v = zmat2xyz(data)
        except (AssertionError, IndexError, ZMError):
            pass
        else:
            # clear the screen before redraw
            for item in reversed(self.window.items):
                self.window.removeItem(item)
            # create a second coordinate list 'vs' that is centered in the GL view
            if len(v) > 0:
                shift = np.mean(v, axis=0)
                vs = np.add(v, -shift)
                elems = [ 1 + next((i for i, sublist in enumerate(colors) if row[0] in sublist), -1) for row in data ]
                nelems = len(elems)
                # define molecule radii and colors
                r = []
                c = []
                for i in elems:
                    r.append(elements[i].covalent_radius)
                    c.append(colors[i-1][-1])
                # draw atoms
                for i in range(nelems):
                    addAtom(self.window, i, r, vs, c)
                # draw bonds where appropriate
                combs = itertools.combinations(range(nelems), 2)
                for i in combs:
                    addBond(self.window, i[0], i[1], r, vs, c)
        if len(v) > 1:
            maxDim = float("-inf")
            for dim in v.T:
                span = max(dim)-min(dim)
                if span > maxDim:
                    maxDim = span
        else: maxDim = 2
        self.window.setCameraPosition(distance=maxDim*1.5)

    # detect mouse clicks in GL window and process them
    def eventFilter(self, obj, event):
        if self.layout.indexOf(obj) != -1:
            if event.type() == event.MouseButtonPress:
                allItems = obj.items
                itms = obj.itemsAt((event.pos().x()-2, event.pos().y()-2, 4, 4))
                # check if something was clicked
                if len(itms) > 0:
                    i = allItems.index(itms[0])
                    e = self.model.index(i, 0).data()
                    # check if an atom was clicked
                    if i < nelems:
                        # highlight the clicked atom
                        addAtom(obj, i, r, vs, c, 'highlight')
                        # populate the list of reference atoms
                        # so we know where to put a new one
                        if len(self.buildList) < min(3, nelems):
                            self.buildList.append(i)
                        # add a new atom
                        if len(self.buildList) == min(3, nelems):
                            selection = self.periodicTable()
                            self.addRow()
                            # temporarily disable GL view updating
                            self.model.dataChanged.disconnect(self.updateView)
                            row = self.model.rowCount()-1
                            newSymbol = selection[1]
                            newBond = round(2.1*gmean([ elements[e].covalent_radius for e in [selection[0], elems[self.buildList[0]]] ]), 4)
                            newData = [newSymbol, self.buildList[0]+1, newBond]
                            if len(self.buildList) >= 2:
                                newAngle = 109.4712
                                newData.append(self.buildList[1]+1)
                                newData.append(newAngle)
                                if len(self.buildList) == 3:
                                    newDihedral = 120.
                                    newData.append(self.buildList[2]+1)
                                    newData.append(newDihedral)
                            for j, cell in enumerate(newData):
                                item = QStandardItem(str(cell))
                                self.model.setItem(row, j, item)
                            # clear reference list and enable GL view updating
                            self.buildList = []
                            self.model.dataChanged.connect(self.updateView)
                            self.updateView()
                    # if clicked, un-highlight the last highlighted atom
                    if i >= len(allItems)-1:
                        obj.removeItem(itms[0])
                        self.buildList.pop()
                # if there are no atoms
                if len(allItems) == 0:
                    if self.model.rowCount() == 0:
                        self.addRow()
                    selection = self.periodicTable()
                    item = QStandardItem(selection[1])
                    self.model.setItem(0, 0, item)
        # also do the default click action
        return super(MainWidget, self).eventFilter(obj, event)

# run the application
app = MainWidget()
app.run()
