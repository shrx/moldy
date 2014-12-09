import numpy as np
from math import degrees, radians

from periodictable import *
import periodictableGUI

import csv

import itertools

import sys
 
# SIP allows us to select the API we wish to use
import sip
 
# use the more modern PyQt API (not enabled by default in Python 2.x);
# must precede importing any module that provides the API specified
sip.setapi('QDate', 2)
sip.setapi('QDateTime', 2)
sip.setapi('QString', 2)
sip.setapi('QTextStream', 2)
sip.setapi('QTime', 2)
sip.setapi('QUrl', 2)
sip.setapi('QVariant', 2)
 
# Import all of Qt
from PyQt4.Qt import *

from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph.opengl as gl

from zmat import ZMat

barve = []
with open('barve.csv') as csvfile:
    csvreader = csv.reader(csvfile)
    for row in csvreader:
        barve.append([int(row[0]),row[1],(float(row[2])/255,float(row[3])/255,float(row[4])/255, 1)])
#print(barve)

def addAtom(w, i, opt=''):
    #r2 = r[i]**1.3*.7
    r2 = r[i]*.7
    if opt == 'highlight':
        r2 += .15
    ms = gl.MeshData.sphere(rows=10, cols=20, radius=r2)
    gs = gl.GLMeshItem(meshdata=ms, smooth=True, drawFaces=True, color=(c[i] if opt != 'highlight' else (0, 1, .2, .5)), drawEdges=False, shader='shaded', glOptions=('opaque' if opt != 'highlight' else 'translucent'))
    gs.translate(vs[i][0], vs[i][1], vs[i][2])
    w.addItem(gs)

def addBond(w, i, j):
    l = np.linalg.norm(np.array(vs[i])-np.array(vs[j]))
    if l < (r[i]+r[j])*1.3:
        xyz = np.subtract(vs[j], vs[i])
        xy = xyz[0]**2 + xyz[1]**2
        # for elevation angle defined from Z-axis down
        s1 = degrees(np.arctan2(np.sqrt(xy), xyz[2]))
        # for elevation angle defined from XY-plane up
        s2 = degrees(np.arctan2(xyz[1], xyz[0]))

        if c[i] == c[j]:
            mc = gl.MeshData.cylinder(rows=2, cols=12, radius=[.1, .1], length=l)
            gc = gl.GLMeshItem(meshdata=mc, smooth=True, drawFaces=True, color=c[i], drawEdges=False, shader='shaded')
            gc.rotate(s1, 0, 1, 0)
            gc.rotate(s2, 0, 0, 1)
            gc.translate(vs[i][0], vs[i][1], vs[i][2])
            w.addItem(gc)
        else:
            rf = r[i]/(r[i]+r[j])
            mc1 = gl.MeshData.cylinder(rows=2, cols=12, radius=[.1, .1], length=l*rf)
            gc1 = gl.GLMeshItem(meshdata=mc1, smooth=True, drawFaces=True, color=c[i], drawEdges=False, shader='shaded')
            gc1.rotate(s1, 0, 1, 0)
            gc1.rotate(s2, 0, 0, 1)
            gc1.translate(vs[i][0], vs[i][1], vs[i][2])
            w.addItem(gc1)

            mc2 = gl.MeshData.cylinder(rows=2, cols=12, radius=[.1, .1], length=l*(1-rf))
            gc2 = gl.GLMeshItem(meshdata=mc2, smooth=True, drawFaces=True, color=c[j], drawEdges=False, shader='shaded')
            gc2.rotate(180, 0, 0, 1)
            gc2.rotate(s1-180, 0, 1, 0)
            gc2.rotate(s2, 0, 0, 1)
            gc2.translate(vs[j][0], vs[j][1], vs[j][2])
            w.addItem(gc2)


def flatten(container):
    for i in container:
        if isinstance(i, list) or isinstance(i, tuple):
            for j in flatten(i):
                yield j
        else:
            yield i

def model2list(model):
    result = []
    num_rows, num_cols = model.rowCount(), model.columnCount()
    for row in range(num_rows):
        cols = []
        for col in range(num_cols):
            item = model.index(row, col).data()
            if item == '':
                item = None
            if col == 0 or item == None:
                text = item
            elif col % 2 == 1:
                text = int(item)
            elif col == num_cols:
                text = radians(float(item))
            else:
                text = float(item)
            cols.append(text)
        if not all(x == None for x in cols):
            result.append([x for x in cols if x != None])
    return result

def zmat2xyz(data):
    zslots = [ tuple(np.add(-1, s[1::2])) for s in data ]
    zvals = list(flatten([ [row[2], [ radians(cell) for cell in row[4::2]]] for row in data[1:] ]))
    return ZMat(zslots)(zvals)

# Every Qt application must have one and only one QApplication object;
# it receives the command line arguments passed to the script, as they
# can be used to customize the application's appearance and behavior
qt_app = QApplication(sys.argv)

class Widget(QWidget):
    def __init__(self):
        # Initialize the object as a QWidget
        QWidget.__init__(self)

        self.model = QtGui.QStandardItemModel(min(len(inp), 1), 7, self)
        self.inputField = QTableView(self)
        self.inputField.setModel(self.model)
        #self.inputField.setColumnCount(7)
        #self.inputField.setRowCount(len(inp))
        self.inputField.setMinimumWidth(325)
        self.model.setHorizontalHeaderLabels(['atom','','bond','','angle','','dihedral'])
        self.populateTable()

        self.periodicTableWidget = periodictableGUI.PeriodicTableDialog()

        for j, width in enumerate([40, 22, 65, 22, 65, 22, 65]):
            self.inputField.setColumnWidth(j, width)
        for i in range(min(len(inp), 3)):
            for j in range(2*i+1, 7):
                self.model.setItem(i, j, QStandardItem())
                self.model.item(i, j).setBackground(QtGui.QColor(150,150,150))
                self.model.item(i, j).setFlags(Qt.ItemIsEnabled)

        self.addRowButton = QtGui.QPushButton('Add row', self)
        self.addRowButton.clicked.connect(self.addRow)

        self.deleteRowButton = QtGui.QPushButton('Delete row', self)
        self.deleteRowButton.clicked.connect(self.deleteRow)

        self.periodicTableButton = QtGui.QPushButton('periodicTable', self)
        self.periodicTableButton.clicked.connect(self.periodicTable)

        #self.updateViewButton = QtGui.QPushButton('Update view', self)
        #self.updateViewButton.clicked.connect(self.updateView)

        self.writeZmatButton = QtGui.QPushButton('Write Zmat', self)
        self.writeZmatButton.clicked.connect(self.writeZmat)

        self.writeXYZButton = QtGui.QPushButton('Write XYZ', self)
        self.writeXYZButton.clicked.connect(self.writeXYZ)

        self.window = gl.GLViewWidget()
        self.window.installEventFilter(self)
        self.window.setMinimumSize(500, 500)
        for i in range(nelems):
            addAtom(self.window, i)

        combs = itertools.combinations(range(nelems), 2)
        for i in combs:
            addBond(self.window, i[0], i[1])

        #print(self.window.items)
        #self.window.addItem(zgrid)

        self.layout = QtGui.QHBoxLayout(self)
        self.left = QtGui.QVBoxLayout()
        self.left.addWidget(self.inputField)
        self.leftBot = QtGui.QHBoxLayout()
        self.leftBot2 = QtGui.QHBoxLayout()
        self.leftBot.addWidget(self.addRowButton)
        self.leftBot.addWidget(self.deleteRowButton)
        self.leftBot.addWidget(self.periodicTableButton)
        self.leftBot2.addWidget(self.writeZmatButton)
        self.leftBot2.addWidget(self.writeXYZButton)
        #hbox.addStretch(1)
        self.left.addLayout(self.leftBot)
        self.left.addLayout(self.leftBot2)
        self.layout.addLayout(self.left)
        self.layout.addWidget(self.window)

        self.adjustSize()
        self.setWindowTitle('Moldy')

        self.model.dataChanged.connect(self.updateView)
        self.buildList = []

    def run(self):
        # Show the form
        self.show()
        # Run the Qt application
        qt_app.exec_()

    def populateTable(self):
        for i, row in enumerate(inp):
            for j, cell in enumerate(row):
                item = QStandardItem(str(cell))
                self.model.setItem(i, j, item)

    def addRow(self):
        print ('add row')
        row = self.model.rowCount()
        self.model.insertRow(row)
        if row < 3:
            for j in range(2*row+1, 7):
                self.model.setItem(row, j, QStandardItem())
                self.model.item(row, j).setBackground(QtGui.QColor(150,150,150))
                self.model.item(row, j).setFlags(Qt.ItemIsEnabled)
        self.updateView()

    def deleteRow(self):
        print ('delete row')
        self.model.removeRow(self.model.rowCount()-1)
        self.updateView()

    def periodicTable(self):
        #global selection
        print('periodicTable')
        self.periodicTableWidget.exec_()
        selection = self.periodicTableWidget.selection()
        return selection

    def writeZmat(self):
        print ('write zmat')
        print(model2list(self.model))

    def writeXYZ(self):
        print ('write xyz')
        xyz = []
        modelList = model2list(self.model)
        for i in range(len(v)):
            xyz.append(np.round(v[i], 7).tolist())
            xyz[i][:0] = modelList[i][0]
        print(xyz)

    def updateView(self):
        global r
        global c
        global v
        global vs
        global nelems
        print ('update view')
        data = model2list(self.model)
        print(data)

        v = zmat2xyz(data)
        shift = np.mean(v, axis=0)
        vs = np.add(v, -shift)

        elems = [ 1 + next((i for i, sublist in enumerate(barve) if row[0] in sublist), -1) for row in data ]
        nelems = len(elems)

        r = []
        c = []
        for i in elems:
            r.append(elements[i].covalent_radius)
            c.append(barve[i-1][-1])

        index = []
        for item in self.window.items:
            index.append(self.window.items.index(item))
        for i in reversed(index):
            self.window.removeItem(self.window.items[i])

        for i in range(nelems):
            addAtom(self.window, i)

        combs = itertools.combinations(range(nelems), 2)
        for i in combs:
            addBond(self.window, i[0], i[1])

    def eventFilter(self, obj, event):
        # you could be doing different groups of actions
        # for different types of widgets and either filtering
        # the event or not.
        # Here we just check if its one of the layout widgets
        if self.layout.indexOf(obj) != -1:
            if event.type() == event.MouseButtonPress:
                #print('click')
                allItems = obj.items
                itms = obj.itemsAt((event.pos().x()-2, event.pos().y()-2, 4, 4))
                if len(itms) > 0:
                    #verts = obj[0].opts['meshdata'].vertexes()
                    #verts[0][-1]
                    i = allItems.index(itms[0])
                    e = self.model.index(i, 0).data()
                    print('Pressed button', event.button(), 'at', (event.pos().x(), event.pos().y()), 'item:', e, i+1, 'nelems:', nelems)
                    if i < nelems:
                        addAtom(obj, i, 'highlight')
                        if len(self.buildList) < min(3, nelems):
                            self.buildList.append(i)
                            print(self.buildList)
                        if len(self.buildList) == min(3, nelems):
                            self.model.dataChanged.disconnect(self.updateView)
                            selection = self.periodicTable()
                            self.addRow()
                            row = self.model.rowCount()-1
                            newSymbol = selection[1]
                            newBond = 3*np.mean([elements[e].covalent_radius for e in [selection[0], self.buildList[0]]])
                            newAngle = 120.
                            newDihedral = 180.
                            for j, cell in enumerate([newSymbol, self.buildList[0]+1, newBond, self.buildList[1]+1, newAngle, self.buildList[2]+1, newDihedral]):
                                item = QStandardItem(str(cell))
                                self.model.setItem(row, j, item)
                            self.buildList = []
                            self.model.dataChanged.connect(self.updateView)
                            self.updateView()
                    if i >= len(allItems)-1:
                        obj.removeItem(itms[0])
                        self.buildList.pop()
                        print(self.buildList)
                # if I returned True right here, the event
                # would be filtered and not reach the obj,
                # meaning that I decided to handle it myself

        # regardless, just do the default
        return super(Widget, self).eventFilter(obj, event)

inp = [['H'],
        ['O', 1, 0.9],
        ['O', 2, 1.4, 1, 105.],
        ['H', 3, 0.9, 2, 105., 1, 120.]]

v = zmat2xyz(inp)
shift = np.mean(v, axis=0)
vs = np.add(v, -shift)

elems = [ 1 + next((i for i, sublist in enumerate(barve) if row[0] in sublist), -1) for row in inp ]
nelems = len(elems)

r = []
c = []
for i in elems:
    r.append(elements[i].covalent_radius)
    c.append(barve[i-1][-1])
# Create an instance of the application window and run it
app = Widget()
app.run()
