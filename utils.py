import numpy as np
from prettytable import PrettyTable
import csv
from math import degrees, radians, acos
import pyqtgraph.opengl as gl
from zmat import ZMat
from math import sin, cos
from numpy.linalg import norm

# import element colors
colors = []
with open('colors.csv') as csvfile:
    csvreader = csv.reader(csvfile)
    for row in csvreader:
        colors.append([int(row[0]), row[1], (float(row[2])/255, float(row[3])/255, float(row[4])/255, 1)])

# draw atoms in GL window
def addAtom(w, i, r, vs, c, opt='', fast=False):
    if fast:
        r2 = .1
        if opt == 'highlight':
            r2 += .05
        ms = gl.MeshData.sphere(rows=10, cols=10, radius=r2)
        gs = gl.GLMeshItem(meshdata=ms, smooth=False, drawFaces=True, color=((1, 1, 1, 0) if opt != 'highlight' else (0, 1, .2, .5)), drawEdges=False, glOptions='translucent')
    else:
        r2 = r[i]*.6
        if opt == 'highlight':
            r2 += .05
        ms = gl.MeshData.sphere(rows=10, cols=20, radius=r2)
        gs = gl.GLMeshItem(meshdata=ms, smooth=True, drawFaces=True, color=(c[i] if opt != 'highlight' else (0, 1, .2, .5)), drawEdges=False, shader='shaded', glOptions=('opaque' if opt != 'highlight' else 'translucent'))
    gs.translate(vs[i][0], vs[i][1], vs[i][2])
    w.addItem(gs)

# draw bonds between atoms in GL window;
# max bond length is configurable with maxLen
def addBond(w, i, j, r, vs, c, fast=False):
    l = np.linalg.norm(np.array(vs[i])-np.array(vs[j]))
    maxLen = (r[i]+r[j])*1.2
    if l < maxLen:
        # convert coordinates from cartesian to spherical
        xyz = np.add(vs[j], -vs[i])
        xy = xyz[0]**2 + xyz[1]**2
        s1 = degrees(np.arctan2(np.sqrt(xy), xyz[2]))
        s2 = degrees(np.arctan2(xyz[1], xyz[0]))
        # if atoms are same element, 1 cylinder needed
        if c[i] == c[j]:
            if fast:
                gc = gl.GLLinePlotItem(pos=np.array([vs[i], vs[j]]), color=c[i], width=3)
            else:
                mc = gl.MeshData.cylinder(rows=2, cols=12, radius=[.1, .1], length=l)
                gc = gl.GLMeshItem(meshdata=mc, smooth=True, drawFaces=True, color=c[i], drawEdges=False, shader='shaded')
                gc.rotate(s1, 0, 1, 0)
                gc.rotate(s2, 0, 0, 1)
                gc.translate(vs[i][0], vs[i][1], vs[i][2])
            w.addItem(gc)
        # atoms are different elements; 2 cylinders needed
        else:
            rf = r[i]/(r[i]+r[j])
            if fast:
                gc1 = gl.GLLinePlotItem(pos=np.array([vs[i], np.mean(np.array([vs[i], vs[j]]), axis=0)]), color=c[i], width=3)
            else:
                mc1 = gl.MeshData.cylinder(rows=2, cols=12, radius=[.1, .1], length=l*rf)
                gc1 = gl.GLMeshItem(meshdata=mc1, smooth=True, drawFaces=True, color=c[i], drawEdges=False, shader='shaded')
                gc1.rotate(s1, 0, 1, 0)
                gc1.rotate(s2, 0, 0, 1)
                gc1.translate(vs[i][0], vs[i][1], vs[i][2])
            w.addItem(gc1)
            if fast:
                gc2 = gl.GLLinePlotItem(pos=np.array([vs[j], np.mean(np.array([vs[i], vs[j]]), axis=0)]), color=c[j], width=3)
            else:
                mc2 = gl.MeshData.cylinder(rows=2, cols=12, radius=[.1, .1], length=l*(1-rf))
                gc2 = gl.GLMeshItem(meshdata=mc2, smooth=True, drawFaces=True, color=c[j], drawEdges=False, shader='shaded')
                gc2.rotate(180, 0, 0, 1)
                gc2.rotate(s1-180, 0, 1, 0)
                gc2.rotate(s2, 0, 0, 1)
                gc2.translate(vs[j][0], vs[j][1], vs[j][2])
            w.addItem(gc2)

# flatten a nested list
def flatten(container):
    for i in container:
        if isinstance(i, list) or isinstance(i, tuple):
            for j in flatten(i):
                yield j
        else:
            yield i

# convert Qt model contents to list
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
            result.append([ x for x in cols if x != None ])
    return result

# convert between Zmatrix and XYZ coordinates
def zmat2xyz(data):
    zslots = [ tuple(np.add(s[1::2], -1)) for s in data ]
    zvals = list(flatten([ [row[2], [ radians(cell) for cell in row[4::2] ]] for row in data[1:] ]))
    return ZMat(zslots)(zvals)

# export data to text file arranged in nice columns
def writeOutput(data, filename):
    if len(data) > 0:
        t = PrettyTable()
        maxLen = max([ len(row) for row in data ])
        for row in data:
            l = len(row)
            if l < maxLen:
                row += ['']*(maxLen - l)
            t.add_row(row)
        t.align = 'l'
        t.header = False
        t.border = False
    else:
        t = ''
    with open(filename, 'w') as f:
        f.write(str(len(data))+'\n')
        f.write('\n')
        f.write(str(t))
        f.write('\n\n')
    f.close()

# try to guess the formula from elements; only usable for filename generation
def getFormula(data):
    d = ''
    for i in set(data):
        d += i
        c = data.count(i)
        if c > 1:
            d += str(c)
    return d

#find the point in a 3D array nearest to the specified point
def nearest(node, nodes):
    nodes = np.asarray(nodes)
    dist_2 = np.sum((nodes - node)**2,axis=1)
    idx = np.argmin(dist_2)
    return idx

# convert between XYZ and Zmatrix coordinates
def xyz2zmat(xyz, atoms):
    zmat = []
    A = np.array(xyz)
    for i, coords in enumerate(xyz):
        if i == 0:
            zmat.append(atoms[i])
        elif i == 1:
            nearest1 = A[0]
            q = coords - nearest1
            distance = np.sqrt(np.dot(q, q))
            zmat.append([atoms[i],i,distance])
        elif i == 2:
            n1 = nearest(A[:i], coords)
            nearest1 = np.array(xyz[n1])
            Atemp = np.delete(A[:i], (n1), axis=0)
            n2 = nearest(Atemp, nearest1)
            nearest2 = Atemp[n2]
            i2 = xyz.index(list(Atemp[n2]))
            q = coords - nearest1
            r = nearest1 - nearest2
            # Create unit vectors
            q_u = q / np.sqrt(np.dot(q, q))
            r_u = r / np.sqrt(np.dot(r, r))
            distance = np.sqrt(np.dot(q, q))
            angle = degrees(acos(np.dot(-q_u, r_u)))
            zmat.append([atoms[i],n1+1,distance,i2+1,angle])
        else:
            n1 = nearest(A[:i], coords)
            nearest1 = np.array(xyz[n1])
            Atemp = np.delete(A[:i], (n1), axis=0)
            n2 = nearest(Atemp, nearest1)
            nearest2 = Atemp[n2]
            Atemp2 = np.delete(Atemp, (n2), axis=0)
            n3 = nearest(Atemp2, nearest2)
            nearest3 = Atemp2[n3]
            i2 = xyz.index(list(Atemp[n2]))
            i3 = xyz.index(list(Atemp2[n3]))
            q = coords - nearest1
            r = nearest1 - nearest2
            s = nearest2 - nearest3
            # Create unit vectors
            q_u = q / np.sqrt(np.dot(q, q))
            r_u = r / np.sqrt(np.dot(r, r))
            distance = np.sqrt(np.dot(q, q))
            angle = degrees(acos(np.dot(-q_u, r_u)))
            plane1 = np.cross(q, r)
            plane2 = np.cross(r, s)
            dihedral = degrees(acos(np.dot(plane1, plane2) / (np.sqrt(np.dot(plane1, plane1)) * np.sqrt(np.dot(plane2, plane2)))))
            if np.dot(np.cross(plane1, plane2), r_u) > 0:
                dihedral = -dihedral
            zmat.append([atoms[i],n1+1,distance,i2+1,angle,i3+1,dihedral])
    return zmat

def slerp(q0, r0, t):
    q = np.array(q0)
    r = np.array(r0)
    l = (np.sqrt(np.dot(q, q))+np.sqrt(np.dot(r, r)))*.45
    q_u = q / np.sqrt(np.dot(q, q))
    r_u = r / np.sqrt(np.dot(r, r))
    angle = acos(np.dot(q_u, r_u))
    sa = sin(angle)
    return sin((1.0-t)*angle)/sa*l * q_u + sin(t*angle)/sa*l * r_u