# coding=utf-8
import pdb
from scipy.stats import gmean
from periodictable import elements
import widgets
import itertools
import sys
from os.path import expanduser
from PyQt4.QtCore import *
from PyQt4.Qt import QApplication, QWidget, QTableView, QStandardItem, QStandardItemModel, QColor, QFileDialog, QHBoxLayout, QVBoxLayout, QStatusBar, QAction, qApp, QMessageBox, QIcon, QMenuBar, QMenu
from pyqtgraph import GraphicsLayoutWidget, mkPen, ErrorBarItem
import pyqtgraph.opengl as gl
from zmat import ZMError
from utils import *
from cclib.parser import ccopen

#debugging
#pdb.set_trace()

# create Qt application
qt_app = QApplication(sys.argv)

# create main widget
class MainWidget(QWidget):
    def __init__(self):
        QWidget.__init__(self)

        # define periodic table widget for element selection
        self.periodicTableWidget = widgets.PeriodicTableDialog()

        # initial molecule Zmatrix (can be empty)
        # self.inp = []
        self.inp = [['H'],
        ['O', 1, 0.9],
        ['O', 2, 1.4, 1, 105.],
        ['H', 3, 0.9, 2, 105., 1, 120.]]

        self.atomList = []
        self.highList = []
        self.labelList = []
        self.fast = False

        # define & initialize ZMatModel that will contain Zmatrix data
        self.ZMatModel = QStandardItemModel(len(self.inp), 7, self)
        self.ZMatTable = QTableView(self)
        self.ZMatTable.setModel(self.ZMatModel)
        self.ZMatTable.setFixedWidth(325)
        #self.ZMatTable.installEventFilter(self)
        #self.ZMatModel.installEventFilter(self)
        self.ZMatModel.setHorizontalHeaderLabels(['atom','','bond','','angle','','dihedral'])
        for j, width in enumerate([40, 22, 65, 22, 65, 22, 65]):
            self.ZMatTable.setColumnWidth(j, width)
        # populate the ZMatModel
        self.populateZMatModel()

        #define Menu bar menus and their actions
        self.menuBar = QMenuBar(self)
        fileMenu = self.menuBar.addMenu('&File')
        editMenu = self.menuBar.addMenu('&Edit')
        viewMenu = self.menuBar.addMenu('&View')
        measureMenu = self.menuBar.addMenu('&Measure')
        helpMenu = self.menuBar.addMenu('&Help')

        readZmatAction = QAction('&Read &ZMat', self)
        readZmatAction.setShortcut('Ctrl+O')
        readZmatAction.setStatusTip('Read Zmat from file')
        readZmatAction.triggered.connect(self.readZmat)
        fileMenu.addAction(readZmatAction)

        readXYZAction = QAction('&Read &XYZ', self)
        readXYZAction.setShortcut('Ctrl+Shift+O')
        readXYZAction.setStatusTip('Read XYZ from file')
        readXYZAction.triggered.connect(self.readXYZ)
        fileMenu.addAction(readXYZAction)

        readGaussianAction = QAction('&Read &Gaussian log', self)
        readGaussianAction.setShortcut('Ctrl+G')
        readGaussianAction.setStatusTip('Read Gaussian log file')
        readGaussianAction.triggered.connect(self.readGaussian)
        fileMenu.addAction(readGaussianAction)

        writeZmatAction = QAction('&Write &ZMat', self)
        writeZmatAction.setShortcut('Ctrl+S')
        writeZmatAction.setStatusTip('Write Zmat to file')
        writeZmatAction.triggered.connect(self.writeZmat)
        fileMenu.addAction(writeZmatAction)

        writeXYZAction = QAction('&Write &XYZ', self)
        writeXYZAction.setShortcut('Ctrl+Shift+S')
        writeXYZAction.setStatusTip('Write XYZ from file')
        writeXYZAction.triggered.connect(self.writeXYZ)
        fileMenu.addAction(writeXYZAction)

        exitAction = QAction('&Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(qApp.quit)
        fileMenu.addAction(exitAction)

        addRowAction = QAction('&Add &row', self)
        addRowAction.setShortcut('Ctrl+R')
        addRowAction.setStatusTip('Add row to ZMatrix')
        addRowAction.triggered.connect(self.addRow)
        editMenu.addAction(addRowAction)

        deleteRowAction = QAction('&Delete &row', self)
        deleteRowAction.setShortcut('Ctrl+Shift+R')
        deleteRowAction.setStatusTip('Delete row from ZMatrix')
        deleteRowAction.triggered.connect(self.deleteRow)
        editMenu.addAction(deleteRowAction)

        addAtomAction = QAction('&Add &atom', self)
        addAtomAction.setShortcut('Ctrl+A')
        addAtomAction.setStatusTip('Add atom to ZMatrix')
        addAtomAction.triggered.connect(self.buildB)
        editMenu.addAction(addAtomAction)

        drawModeMenu = QMenu('Draw mode', self)
        viewMenu.addMenu(drawModeMenu)
        fastDrawAction = QAction('&Fast draw', self)
        fastDrawAction.triggered.connect(self.fastDraw)
        normalDrawAction = QAction('&Normal draw', self)
        normalDrawAction.triggered.connect(self.normalDraw)
        drawModeMenu.addAction(normalDrawAction)
        drawModeMenu.addAction(fastDrawAction)

        clearHighlightsAction = QAction('&Clear selection', self)
        clearHighlightsAction.setShortcut('Ctrl+C')
        clearHighlightsAction.setStatusTip('Clear highlighted atoms')
        clearHighlightsAction.triggered.connect(self.clearHighlights)
        viewMenu.addAction(clearHighlightsAction)

        clearLabelsAction = QAction('&Clear labels', self)
        clearLabelsAction.setShortcut('Ctrl+Alt+C')
        clearLabelsAction.setStatusTip('Clear labels')
        clearLabelsAction.triggered.connect(self.clearLabels)
        viewMenu.addAction(clearLabelsAction)

        clearUpdateViewAction = QAction('&Clear selection and labels', self)
        clearUpdateViewAction.setShortcut('Ctrl+Shift+C')
        clearUpdateViewAction.setStatusTip('Clear highlighted atoms and labels')
        clearUpdateViewAction.triggered.connect(self.clearUpdateView)
        viewMenu.addAction(clearUpdateViewAction)

        self.showGaussAction = QAction('Show &Gaussian geometry optimization', self)
        self.showGaussAction.setShortcut('Ctrl+G')
        self.showGaussAction.setStatusTip('Show Gaussian geometry optimization plots for energy, force and displacement.')
        self.showGaussAction.setEnabled(False)
        self.showGaussAction.triggered.connect(self.showGauss)
        viewMenu.addAction(self.showGaussAction)
        self.showFreqAction = QAction('Show &IR frequency plot', self)
        self.showFreqAction.setShortcut('Ctrl+I')
        self.showFreqAction.setStatusTip('Show Gaussian calculated IR frequency plot.')
        self.showFreqAction.setEnabled(False)
        self.showFreqAction.triggered.connect(self.showFreq)
        viewMenu.addAction(self.showFreqAction)

        measureDistanceAction = QAction('&Measure &distance', self)
        measureDistanceAction.setShortcut('Ctrl+D')
        measureDistanceAction.setStatusTip('Measure distance between two atoms')
        measureDistanceAction.triggered.connect(self.measureDistanceB)
        measureMenu.addAction(measureDistanceAction)

        measureAngleAction = QAction('&Measure &angle', self)
        measureAngleAction.setShortcut('Ctrl+Shift+D')
        measureAngleAction.setStatusTip('Measure angle between three atoms')
        measureAngleAction.triggered.connect(self.measureAngleB)
        measureMenu.addAction(measureAngleAction)

        aboutAction = QAction('&About', self)
        aboutAction.setStatusTip('About this program...')
        aboutAction.triggered.connect(self.about)
        helpMenu.addAction(aboutAction)

        aboutQtAction = QAction('&About Qt', self)
        aboutQtAction.setStatusTip('About Qt...')
        aboutQtAction.triggered.connect(self.aboutQt)
        helpMenu.addAction(aboutQtAction)

        # define GL widget that displays the 3D molecule model
        self.window = widgets.MyGLView()
        self.window.installEventFilter(self)
        self.window.setMinimumSize(500, 500)
        #self.window.setBackgroundColor((50, 0, 10))
        self.updateView()

        self.gaussianPlot = GraphicsLayoutWidget()
        self.gaussianPlot.resize(750, 250)
        self.gaussianPlot.setWindowTitle('Gaussian geometry optimization')
        #self.gaussianPlot.setAspectLocked(True)
        #self.gaussianPlot.addLayout(rowspan=3, colspan=1)

        self.FreqModel = QStandardItemModel(1, 3, self)
        self.freqTable = QTableView(self)
        self.freqTable.setModel(self.FreqModel)
        self.freqTable.setMinimumWidth(240)
        self.freqTable.installEventFilter(self)
        self.FreqModel.installEventFilter(self)
        self.FreqModel.setHorizontalHeaderLabels(['Frequency','IR Intensity','Raman Intensity'])
        for j, width in enumerate([80, 80, 80]):
            self.freqTable.setColumnWidth(j, width)

        self.freqWidget = QWidget()
        self.freqWidget.setWindowTitle('IR frequency plot & table')
        self.freqWidget.resize(800, 400)
        self.freqWidget.layout = QHBoxLayout(self.freqWidget)
        self.freqWidget.layout.setSpacing(1)
        self.freqWidget.layout.setContentsMargins(1, 1, 1, 1)
        self.freqPlot = GraphicsLayoutWidget()
        self.freqWidget.layout.addWidget(self.freqPlot)
        self.freqWidget.layout.addWidget(self.freqTable)
        self.freqTable.clicked.connect(self.freqCellClicked)

        # define other application parts
        self.statusBar = QStatusBar(self)
        self.fileDialog = QFileDialog(self)

        # define application layout
        self.layout = QVBoxLayout(self)
        self.layout.setSpacing(1)
        self.layout.setContentsMargins(1, 1, 1, 1)
        self.layout1 = QHBoxLayout()
        self.layout1.setSpacing(1)
        self.layout1.addWidget(self.ZMatTable)
        self.layout1.addWidget(self.window)
        self.layout.addWidget(self.menuBar)
        self.layout.addLayout(self.layout1)
        self.layout.addWidget(self.statusBar)

        self.adjustSize()
        self.setWindowTitle('Moldy')
        iconPath = 'icon.png'
        icon = QIcon(iconPath)
        icon.addFile(iconPath, QSize(16, 16))
        icon.addFile(iconPath, QSize(24, 24))
        icon.addFile(iconPath, QSize(32, 32))
        icon.addFile(iconPath, QSize(48, 48))
        icon.addFile(iconPath, QSize(256, 256))
        self.setWindowIcon(icon)

        # start monitoring changes in the ZMatModel
        self.ZMatModel.dataChanged.connect(self.clearUpdateView)

    # run and show the application
    def run(self):
        self.show()
        self.ZMatTable.clicked.connect(self.ZMatCellClicked)
        qt_app.instance().aboutToQuit.connect(self.deleteGLwidget)
        qt_app.exec_()

    # fill the ZMatModel with initial data from 'self.inp'
    def populateZMatModel(self):
        self.ZMatModel.removeRows(0, self.ZMatModel.rowCount())
        for i, row in enumerate(self.inp):
            for j, cell in enumerate(row):
                item = QStandardItem(str(cell))
                self.ZMatModel.setItem(i, j, item)
        # some cells should not be editable, they are disabled
        for i in range(min(len(self.inp), 3)):
            for j in range(2*i+1, 7):
                self.ZMatModel.setItem(i, j, QStandardItem())
                self.ZMatModel.item(i, j).setBackground(QColor(150,150,150))
                self.ZMatModel.item(i, j).setFlags(Qt.ItemIsEnabled)
    
    def populateFreqModel(self):
        self.FreqModel.removeRows(0, self.FreqModel.rowCount())
        for i, row in enumerate(zip(self.vibfreqs, self.vibirs, self.vibramans)):
            for j, cell in enumerate(row):
                item = QStandardItem(str(cell))
                self.FreqModel.setItem(i, j, item)

    # add a row to the bottom of the ZMatModel
    def addRow(self):
        # temporarily stop updating the GL window
        self.ZMatModel.dataChanged.disconnect(self.clearUpdateView)
        row = self.ZMatModel.rowCount()
        self.ZMatModel.insertRow(row)
        # some cells should not be editable
        if row < 3:
            for j in range(2*row+1, 7):
                self.ZMatModel.setItem(row, j, QStandardItem())
                self.ZMatModel.item(row, j).setBackground(QColor(150,150,150))
                self.ZMatModel.item(row, j).setFlags(Qt.ItemIsEnabled)
        # restart GL window updating
        self.ZMatModel.dataChanged.connect(self.clearUpdateView)
        self.statusBar.clearMessage()
        self.statusBar.showMessage('Added 1 row.', 3000)

    # delete the last row of the ZMatModel
    def deleteRow(self):
        xyz = [list(vi) for vi in list(v)]
        atoms = [str(elements[e]) for e in elems]
        oldLen = self.ZMatModel.rowCount()
        idxs = sorted(set(idx.row() for idx in self.ZMatTable.selectedIndexes()), reverse=True)
        newLen = oldLen - len(idxs)
        if newLen == oldLen:
            self.ZMatModel.removeRow(self.ZMatModel.rowCount()-1)
        else:
            self.ZMatModel.dataChanged.disconnect(self.clearUpdateView)
            for idx in idxs:
                self.ZMatModel.removeRow(idx)
                if idx < 3:
                    for i in range(idx, min(3, newLen)):
                        for j in range(2*i+1, 7):
                            self.ZMatModel.setItem(i, j, QStandardItem())
                            self.ZMatModel.item(i, j).setBackground(QColor(150,150,150))
                            self.ZMatModel.item(i, j).setFlags(Qt.ItemIsEnabled)
                if len(xyz) > idx:
                    xyz.pop(idx)
                    atoms.pop(idx)
            self.inp = xyz2zmat(xyz, atoms)
            self.populateZMatModel()
            for i in reversed(self.highList):
                self.window.removeItem(i[1])
            self.highList = []
            self.ZMatModel.dataChanged.connect(self.clearUpdateView)
        self.updateView()
        self.statusBar.clearMessage()
        if idxs:
            self.statusBar.showMessage('Deleted row(s): '+str([i+1 for i in idxs]), 3000)
        else:
            self.statusBar.showMessage('Deleted last row.', 3000)

    # show the periodic table widget
    def periodicTable(self):
        self.statusBar.clearMessage()
        self.statusBar.showMessage('Select element from periodic table.')
        self.periodicTableWidget.exec_()
        selection = self.periodicTableWidget.selection()
        return selection

    # import molecule with zmatrix coordinates
    def readZmat(self):
        self.ZMatModel.dataChanged.disconnect(self.clearUpdateView)
        filename = self.fileDialog.getOpenFileName(self, 'Open file', expanduser('~'), '*.zmat;;*.*')
        self.inp = []
        self.populateZMatModel()
        if filename:
            with open(filename, 'r') as f:
                next(f)
                next(f)
                for row in f:
                    self.inp.append(row.split())
                f.close()
            self.populateZMatModel()
        self.ZMatModel.dataChanged.connect(self.clearUpdateView)
        self.updateView()
        self.statusBar.clearMessage()
        self.statusBar.showMessage('Read molecule from '+filename+'.', 5000)
        self.showGaussAction.setEnabled(False)
        self.showFreqAction.setEnabled(False)

    # import molecule with xyz coordinates
    def readXYZ(self):
        self.ZMatModel.dataChanged.disconnect(self.clearUpdateView)
        filename = self.fileDialog.getOpenFileName(self, 'Open file', expanduser('~'), '*.xyz;;*.*')
        xyz = []
        elems = []
        self.inp = []
        self.populateZMatModel()
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
            self.populateZMatModel()
            #print(elems)
        self.ZMatModel.dataChanged.connect(self.clearUpdateView)
        self.updateView()
        self.statusBar.clearMessage()
        self.statusBar.showMessage('Read molecule from '+filename+'.', 5000)
        self.showGaussAction.setEnabled(False)
        self.showFreqAction.setEnabled(False)

    # import Gaussian log file
    def readGaussian(self):
        global vsShifted
        self.ZMatModel.dataChanged.disconnect(self.clearUpdateView)
        filename = self.fileDialog.getOpenFileName(self, 'Open file', expanduser('~'), '*.log;;*.*')
        if filename:
            self.gaussianPlot.clear()
            self.inp = []
            self.populateZMatModel()
            file = ccopen(filename)
            data = file.parse().getattributes()
            self.natom = data['natom']
            self.atomnos = data['atomnos'].tolist()
            self.atomsymbols = [ str(elements[e]) for e in self.atomnos ]
            self.atomcoords = data['atomcoords'].tolist()
            self.scfenergies = data['scfenergies'].tolist()
            self.geovalues = data['geovalues'].T.tolist()
            self.geotargets = data['geotargets'].tolist()
            if 'vibfreqs' in data.keys():
                self.vibfreqs = data['vibfreqs']
                #print(self.vibfreqs)
                self.vibirs = data['vibirs']
                #print(self.vibirs)
                #print(data.keys())
                if 'vibramans' in data.keys():
                    self.vibramans = data['vibramans']
                else:
                    self.vibramans = [''] * len(self.vibirs)
                self.vibdisps = data['vibdisps']
                #print(self.vibdisps)
            self.inp = xyz2zmat(self.atomcoords[0], self.atomsymbols)
            self.populateZMatModel()

            titles = ['SCF Energies', 'RMS & Max Forces', 'RMS & Max Displacements']
            for i in range(3):
                self.gaussianPlot.addPlot(row=1, col=i+1)
                plot = self.gaussianPlot.getItem(1, i+1)
                plot.setTitle(title=titles[i])
                if i == 0:
                    c = ['c']
                    x = [0]
                    y = [self.scfenergies]
                else:
                    c = ['r', 'y']
                    x = [0, 0]
                    y = [self.geovalues[2*i-2], self.geovalues[2*i-1]]
                    targety = [self.geotargets[2*i-2], self.geotargets[2*i-1]]
                plot.clear()
                plot.maxData = plot.plot(y[0], symbol='o', symbolPen=c[0], symbolBrush=c[0], pen=c[0], symbolSize=5, pxMode=True, antialias=True, autoDownsample=False)
                plot.highlight=plot.plot(x, [ yy[0] for yy in y ], symbol='o', symbolPen='w', symbolBrush=None, pen=None, symbolSize=15, pxMode=True, antialias=True, autoDownsample=False)
                plot.maxData.sigPointsClicked.connect(self.gausclicked)
                if i > 0:
                    for j in range(2):
                        plot.addLine(y=np.log10(targety[j]), pen=mkPen((255, 255*j, 0, int(255/2)), width=1))
                    plot.RMSData=plot.plot(y[1], symbol='o', symbolPen=c[1], symbolBrush=c[1], pen=c[1], symbolSize=5, pxMode=True, antialias=True, autoDownsample=False)
                    plot.RMSData.sigPointsClicked.connect(self.gausclicked)
                    plot.setLogMode(y=True)
            self.showGauss()
            self.updateView()
            self.statusBar.clearMessage()
            self.statusBar.showMessage('Read molecule from '+filename+'.', 5000)
            self.ZMatModel.dataChanged.connect(self.clearUpdateView)
            if self.natom:
                self.showGaussAction.setEnabled(True)
            if 'vibfreqs' in data.keys():
                self.showFreqAction.setEnabled(True)

                # populate the FreqModel
                self.populateFreqModel()

                self.freqPlot.clear()
                irPlot = self.freqPlot.addPlot(row=1, col=1)
                irPlot.clear()
                minFreq = np.min(self.vibfreqs)
                maxFreq = np.max(self.vibfreqs)
                maxInt = np.max(self.vibirs)
                x = np.sort(np.concatenate([np.linspace(minFreq-100, maxFreq+100, num=1000), self.vibfreqs]))
                y = x*0
                for f,i in zip(self.vibfreqs, self.vibirs):
                    y += lorentzv(x, f, 2*np.pi, i)
                #xy = np.array([np.concatenate([x, np.array(self.vibfreqs)]), np.concatenate([y, np.array(self.vibirs)])]).T
                #xysort = xy[xy[:,0].argsort()]
                irPlot.maxData = irPlot.plot(x, y, antialias=True)
                markers = ErrorBarItem(x=self.vibfreqs, y=self.vibirs, top=maxInt/30, bottom=None, pen='r')
                irPlot.addItem(markers)
                self.showFreq()
                #self.vibdisps = np.append(self.vibdisps, [np.mean(self.vibdisps, axis=0)], axis=0)
                maxt = 100
                vsShifted = np.array([ [ vs + self.vibdisps[i]*np.sin(t*2*np.pi/maxt)/3 for t in range(maxt) ] for i in range(len(self.vibfreqs)) ])
            else:
                self.showFreqAction.setEnabled(False)
                self.freqWidget.hide()

    def showGauss(self):
        self.gaussianPlot.show()

    def showFreq(self):
        self.freqWidget.show()

    # export Zmatrix to csv
    def writeZmat(self):
        zm = model2list(self.ZMatModel)
        filename = self.fileDialog.getSaveFileName(self, 'Save file', expanduser('~')+'/'+getFormula(list(list(zip(*zm))[0]))+'.zmat', '*.zmat;;*.*')
        try:
            filename
        except NameError:
            pass
        else:
            if filename:
                writeOutput(zm, filename)
                self.statusBar.clearMessage()
                self.statusBar.showMessage('Wrote molecule to '+filename+'.', 5000)

    # export XYZ coordinates to csv
    def writeXYZ(self):
        xyz = []
        zm = model2list(self.ZMatModel)
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
                self.statusBar.clearMessage()
                self.statusBar.showMessage('Wrote molecule to '+filename+'.', 5000)

    # redraw the 3D molecule in GL widget
    def updateView(self):
        global r
        global c
        global v
        global vs
        global elems
        global nelems
        data = model2list(self.ZMatModel)
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
            self.atomList = []
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
                    addAtom(self.window, i, r, vs, c, fast=self.fast)
                    self.atomList.append([i, self.window.items[-1]])
                #print(self.atomList)
                # draw bonds where appropriate
                combs = list(itertools.combinations(range(nelems), 2))
                bonds = []
                for i in combs:
                    bonds.append(addBond(self.window, i[0], i[1], r, vs, c, fast=self.fast))
                if self.fast:
                    bondedAtoms = set(filter((None).__ne__, flatten(bonds)))
                    for i in set(range(nelems)) - bondedAtoms:
                        addUnbonded(self.window, i, vs, c)
                        self.atomList[i][1]=self.window.items[-1]
                    #print(self.atomList)

                for i in self.highList:
                    self.window.addItem(i[1])
                for i in self.labelList:
                    self.window.addItem(i)
        if len(v) > 1:
            maxDim = float('-inf')
            for dim in v.T:
                span = max(dim)-min(dim)
                if span > maxDim:
                    maxDim = span
        else: maxDim = 2
        self.window.setCameraPosition(distance=maxDim*1.5+1)

    global index
    index = 0
    def updateFreq(self):
        global vsShifted, index, r, c
        index += 1
        index = index % len(vsShifted[0])
        #print(index)
        #print(vsShifted[index])
        for item in reversed(self.window.items):
            self.window.removeItem(item)
        for i in range(nelems):
            addAtom(self.window, i, r, vsShifted[self.freqIndex, index], c, fast=self.fast)
            self.atomList.append([i, self.window.items[-1]])
        combs = itertools.combinations(range(nelems), 2)
        bonds = []
        for i in combs:
            bonds.append(addBond(self.window, i[0], i[1], r, vsShifted[self.freqIndex, index], c, fast=self.fast))
        if self.fast:
            bondedAtoms = set(filter((None).__ne__, flatten(bonds)))
            for i in set(range(nelems)) - bondedAtoms:
                addUnbonded(self.window, i, vsShifted[self.freqIndex, index], c)
                self.atomList[i][1]=self.window.items[-1]

    # detect mouse clicks in GL window and process them
    def eventFilter(self, obj, event):
        if obj == self.window:
            if event.type() == event.MouseButtonPress:
                itms = obj.itemsAt((event.pos().x()-2, event.pos().y()-2, 4, 4))
                if len(itms):
                    self.highlight(obj, [itms[0]])
                elif len(self.atomList) == 0:
                    self.build()
        # also do the default click action
        return super(MainWidget, self).eventFilter(obj, event)

    def ZMatCellClicked(self):
        idxs = sorted(set(idx.row() for idx in self.ZMatTable.selectedIndexes()), reverse=True)
        itms = []
        if self.highList:
            highIdx = list(np.array(self.highList).T[0])
        for idx in idxs:
            if self.highList and idx in highIdx:
                itms.append(self.highList[highIdx.index(idx)][1])
            elif len(self.atomList) > idx:
                itms.append(self.atomList[idx][1])
        self.highlight(self.window, itms)

    def freqCellClicked(self):
        global vsShifted
        self.timer = QTimer()
        self.timer.setInterval(30)
        self.timer.timeout.connect(self.updateFreq)
        idxs = [ idx.row() for idx in self.freqTable.selectedIndexes() ]
        if len(idxs) == 1:
            self.freqIndex = idxs[0]
            self.timer.stop()
            self.timer.timeout.connect(self.updateFreq)
            try:
                self.ZMatModel.dataChanged.disconnect(self.clearUpdateView)
            except TypeError:
                pass
            self.timer.start()
        if len(idxs) != 1:
            self.timer.stop()
            self.freqTable.clearSelection()
            self.timer.timeout.disconnect(self.updateFreq)
            self.ZMatModel.dataChanged.connect(self.clearUpdateView)
            self.clearUpdateView()

    def gausclicked(self, item, point):
        itemdata = item.scatter.data
        points = [ row[7] for row in itemdata ]
        idx = points.index(point[0])
        for i in range(3):
            if i == 0:
                x = [idx]
                y = [self.scfenergies[idx]]
            else:
                x = [idx, idx]
                y = [self.geovalues[2*i-2][idx], self.geovalues[2*i-1][idx]]
            plot = self.gaussianPlot.getItem(1, i+1)
            plot.removeItem(plot.highlight)
            plot.highlight=plot.plot(x, y, symbol='o', symbolPen='w', symbolBrush=None, pen=None, symbolSize=15, pxMode=True, antialias=True, autoDownsample=False)
        self.ZMatModel.dataChanged.disconnect(self.clearUpdateView)
        self.inp = []
        self.populateZMatModel()
        self.inp = xyz2zmat(self.atomcoords[min(idx, len(self.atomcoords)-1)], self.atomsymbols)
        self.populateZMatModel()
        self.ZMatModel.dataChanged.connect(self.clearUpdateView)
        self.updateView()

    def highlight(self, obj, itms):
        for itm in itms:
            idx = next((i for i, sublist in enumerate(self.atomList) if itm in sublist), -1)
            #print(idx)
            if idx != -1:
                addAtom(obj, idx, r, vs, c, opt='highlight', fast=self.fast)
                self.highList.append([idx, obj.items[-1]])
                self.ZMatTable.selectRow(idx)
            idx = next((i for i, sublist in enumerate(self.highList) if itm in sublist), -1)
            if idx != -1:
                obj.removeItem(self.highList[idx][1])
                self.highList.pop(idx)
                self.ZMatTable.clearSelection()
        self.statusBar.clearMessage()
        if len(self.highList) > 0:
            idxs = np.asarray(self.highList).T[0]
            selected = []
            for i in idxs:
                selected.append(str(i+1)+str(elements[elems[i]]))
            self.statusBar.showMessage('Selected atoms: '+str(selected), 5000)

    def buildB(self):
        try:
            nelems
        except NameError:
            self.build()
        else:
            if len(self.highList) <= min(nelems, 3):
                diff = min(nelems, 3) - len(self.highList)
                if diff != 0:
                    self.statusBar.clearMessage()
                    self.statusBar.showMessage('Please select '+str(diff)+' more atom(s).')
                else:
                    self.build()
            else:
                self.statusBar.clearMessage()
                self.statusBar.showMessage('Too many atoms selected.')

    def build(self):
        selection = self.periodicTable()
        row = self.ZMatModel.rowCount()
        self.addRow()
        self.ZMatModel.dataChanged.disconnect(self.clearUpdateView)
        newSymbol = selection[1]
        newData = [newSymbol]
        if len(self.highList) >= 1:
            newBond = round(2.1*gmean([ elements[e].covalent_radius for e in [selection[0], elems[self.highList[0][0]]] ]), 4)
            newData.append(self.highList[0][0]+1)
            newData.append(newBond)
            if len(self.highList) >= 2:
                newAngle = 109.4712
                newData.append(self.highList[1][0]+1)
                newData.append(newAngle)
                if len(self.highList) == 3:
                    newDihedral = 120.
                    newData.append(self.highList[2][0]+1)
                    newData.append(newDihedral)
        for j, cell in enumerate(newData):
            item = QStandardItem(str(cell))
            self.ZMatModel.setItem(row, j, item)
        self.highList = []
        self.ZMatModel.dataChanged.connect(self.clearUpdateView)
        self.updateView()

    def measureDistanceB(self):
        sel = len(self.highList)
        if sel <= 2:
            if sel < 2:
                self.statusBar.clearMessage()
                self.statusBar.showMessage('Please select '+str(2-sel)+' more atom(s).')
            else:
                self.measureDistance()
        else:
            self.statusBar.clearMessage()
            self.statusBar.showMessage('Too many atoms selected.')

    def measureDistance(self):
        pts = []
        for pt in self.highList:
            pts.append(vs[pt[0]])
        pts = np.array(pts)
        self.clearHighlights()
        line = gl.GLLinePlotItem(pos=pts, color=(0., 1., 0., 1.), width=3)
        self.window.addItem(line)
        self.labelList.append(line)
        q = pts[1]-pts[0]
        dist = round(np.sqrt(np.dot(q, q)), 4)
        self.window.labelPos.append(np.mean(pts[0:2], axis=0))
        self.window.labelText.append(str(dist))
        self.statusBar.clearMessage()
        self.statusBar.showMessage('Measured distance: '+str(dist)+' A.', 3000)

    def measureAngleB(self):
        sel = len(self.highList)
        if sel <= 3:
            if sel < 3:
                self.statusBar.clearMessage()
                self.statusBar.showMessage('Please select '+str(3-sel)+' more atom(s).')
            else:
                self.measureAngle()
        else:
            self.statusBar.clearMessage()
            self.statusBar.showMessage('Too many atoms selected.')

    def measureAngle(self):
        pts = []
        for pt in self.highList:
            pts.append(vs[pt[0]])
        pts = np.array(pts)
        q = pts[1]-pts[0]
        r = pts[2]-pts[0]
        q_u = q / np.sqrt(np.dot(q, q))
        r_u = r / np.sqrt(np.dot(r, r))
        angle = round(degrees(acos(np.dot(q_u, r_u))), 1)
        srange = np.array([slerp(q, r, t) for t in np.arange(0.0, 13/12, 1/12)])
        self.clearHighlights()
        for i in range(12):
            mesh = gl.MeshData(np.array([[0,0,0],srange[i],srange[i+1]]))
            tri = gl.GLMeshItem(meshdata=mesh, smooth=False, computeNormals=False, color=(0.3, 1., 0.3, 0.5), glOptions=('translucent'))
            tri.translate(pts[0][0], pts[0][1], pts[0][2])
            self.window.addItem(tri)
            self.labelList.append(tri)
        self.window.labelPos.append(slerp(q, r, 0.5)+pts[0])
        self.window.labelText.append(str(angle))
        self.statusBar.clearMessage()
        self.statusBar.showMessage('Measured angle: '+str(angle)+'°', 3000)

    def clearLabels(self):
        self.window.labelPos = []
        self.window.labelText = []
        self.labelList = []
        self.updateView()

    def clearHighlights(self):
        for item in reversed(self.highList):
                self.window.removeItem(item[1])
        self.highList = []
        self.updateView()

    def clearUpdateView(self):
        self.window.labelPos = []
        self.window.labelText = []
        self.labelList = []
        for item in reversed(self.highList):
                self.window.removeItem(item[1])
        self.highList = []
        self.updateView()
        #print(self.highList)

    def fastDraw(self):
        if not self.fast:
            self.fast = True
            self.updateView()

    def normalDraw(self):
        if self.fast:
            self.fast = False
            self.updateView()

    def about(self):
        QMessageBox.about(self, 'About moldy', 'moldy beta 15. 9. 2015')

    def aboutQt(self):
        QMessageBox.aboutQt(self, 'About Qt')

    def deleteGLwidget(self):
        self.window.setParent(None)
        del self.window


# run the application
app = MainWidget()
app.run()
