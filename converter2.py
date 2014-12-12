'''
A module for converting between zmatrices and cartesian coordinates
'''

import math as m
import numpy as np

VERBOSE = True


masses = { 'X': 0, 'Ac': 227.028, 'Al': 26.981539, 'Am': 243, 'Sb': 121.757, 'Ar': 39.948, 'As': 74.92159, 'At': 210, 'Ba': 137.327, 'Bk': 247, 'Be': 9.012182, 'Bi': 208.98037, 'Bh': 262, 'B': 10.811, 'Br': 79.904, 'Cd': 112.411, 'Ca': 40.078, 'Cf': 251, 'C': 12.011, 'Ce': 140.115, 'Cs': 132.90543, 'Cl': 35.4527, 'Cr': 51.9961, 'Co': 58.9332, 'Cu': 63.546, 'Cm': 247, 'Db': 262, 'Dy': 162.5, 'Es': 252, 'Er': 167.26, 'Eu': 151.965, 'Fm': 257, 'F': 18.9984032, 'Fr': 223, 'Gd': 157.25, 'Ga': 69.723, 'Ge': 72.61, 'Au': 196.96654, 'Hf': 178.49, 'Hs': 265, 'He': 4.002602, 'Ho': 164.93032, 'H': 1.00794, 'In': 114.82, 'I': 126.90447, 'Ir': 192.22, 'Fe': 55.847, 'Kr': 83.8, 'La': 138.9055, 'Lr': 262, 'Pb': 207.2, 'Li': 6.941, 'Lu': 174.967, 'Mg': 24.305, 'Mn': 54.93805, 'Mt': 266, 'Md': 258, 'Hg': 200.59, 'Mo': 95.94, 'Nd': 144.24, 'Ne': 20.1797, 'Np': 237.048, 'Ni': 58.6934, 'Nb': 92.90638, 'N': 14.00674, 'No': 259, 'Os': 190.2, 'O': 15.9994, 'Pd': 106.42, 'P': 30.973762, 'Pt': 195.08, 'Pu': 244, 'Po': 209, 'K': 39.0983, 'Pr': 140.90765, 'Pm': 145, 'Pa': 231.0359, 'Ra': 226.025, 'Rn': 222, 'Re': 186.207, 'Rh': 102.9055, 'Rb': 85.4678, 'Ru': 101.07, 'Rf': 261, 'Sm': 150.36, 'Sc': 44.95591, 'Sg': 263, 'Se': 78.96, 'Si': 28.0855, 'Ag': 107.8682, 'Na': 22.989768, 'Sr': 87.62, 'S': 32.066, 'Ta': 180.9479, 'Tc': 98, 'Te': 127.6, 'Tb': 158.92534, 'Tl': 204.3833, 'Th': 232.0381, 'Tm': 168.93421, 'Sn': 118.71, 'Ti': 47.88, 'W': 183.85, 'U': 238.0289, 'V': 50.9415, 'Xe': 131.29, 'Yb': 173.04, 'Y': 88.90585, 'Zn': 65.39, 'Zr': 91.224 }
total_mass = 0
cartesian = []
zmatrix = []

def read_zmatrix( input_file='zmatrix.dat'):
      '''Read the input zmatrix file (assumes no errors and no variables)'''
      '''The zmatrix is a list with each element formatted as follows
      [ name, [[ atom1, distance ], [ atom2, angle ], [ atom3, dihedral ]], mass ]
      The first three atoms have blank lists for the undefined coordinates'''
      zmatrix = []
      with open(input_file, 'r') as f:
            f.readline()
            f.readline()
            name = f.readline().strip()
            zmatrix.append([ name, [], masses[name] ])
            name, atom1, distance = f.readline().split()[:3]
            zmatrix.append([ name,
                                             [ [int(atom1) - 1, float(distance)], [], [] ],
                                             masses[name] ])
            name, atom1, distance, atom2, angle = f.readline().split()[:5]
            zmatrix.append([ name,
                                             [[ int(atom1) - 1, float(distance) ],
                                              [int(atom2) - 1, m.radians(float(angle)) ], []],
                                             masses[name] ])
            for line in f.readlines():
                  # Get the components of each line, dropping anything extra
                  name, atom1, distance, atom2, angle, atom3, dihedral = line.split()[:7]
                  # convert to a base 0 indexing system and use radians
                  atom = [ name,
                              [ [int(atom1) - 1, float(distance) ],
                                [int(atom2) - 1, m.radians(float(angle)) ],
                                [int(atom3) - 1, m.radians(float(dihedral)) ] ],
                              masses[name] ]
                  zmatrix.append(atom)

      return zmatrix

def read_cartesian(data):
      '''Read the cartesian coordinates file (assumes no errors)'''
      '''The cartesian coordiantes consist of a list of atoms formatted as follows
      [ name, np.array([ x, y, z ]), mass ]
      '''
      cartesian = []
      for coords in data:
            if len(coords) == 4:
                  name = coords[0]
                  position = []
                  for i in coords[1:4]:
                        position.append(float(i))
                  cartesian.append([name, np.array(position), masses[name]])

      return cartesian

def rotation_matrix(axis, angle):
      '''Euler-Rodrigues formula for rotation matrix'''
      # Normalize the axis
      axis = axis/np.sqrt(np.dot(axis,axis))
      a = np.cos(angle/2)
      b,c,d = -axis*np.sin(angle/2)
      return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                              [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                              [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])

def add_first_three_to_cartesian():
      '''The first three atoms in the zmatrix need to be treated differently'''
      # First atom
      name, coords, mass = zmatrix[0]
      cartesian = [[ name, np.array([0, 0, 0]), mass ]]

      # Second atom
      name, coords, mass = zmatrix[1]
      distance = coords[0][1]
      cartesian.append([ name, np.array([distance, 0, 0]), masses[name] ])

      # Third atom
      name, coords, mass =  zmatrix[2]
      atom1, atom2 = coords[:2]
      atom1, distance = atom1
      atom2, angle = atom2
      q = np.array(cartesian[atom1][1]) # position of atom 1
      r = np.array(cartesian[atom2][1]) # position of atom 2

      # Vector pointing from q to r
      a = r - q

      # Vector of length distance pointing along the x-axis
      d = distance*a/np.sqrt(np.dot(a,a))

      # Rotate d by the angle around the z-axis
      d = np.dot(rotation_matrix([0,0,1], angle), d)

      # Add d to the position of q to get the new coordinates of the atom
      p = q + d
      atom = [ name, p, masses[name] ]
      cartesian.append(atom)

def add_atom_to_cartesian( coords):
      '''Find the cartesian coordinates of the atom'''
      name, coords, mass = coords
      atom1, distance = coords[0]
      atom2, angle = coords[1]
      atom3, dihedral = coords[2]

      q = cartesian[atom1][1] # atom 1
      r = cartesian[atom2][1] # atom 2
      s = cartesian[atom3][1] # atom 3

      # Vector pointing from q to r
      a = r - q
      # Vector pointing from s to r
      b = r - s

      # Vector of length distance pointing from q to r
      d = distance*a/np.sqrt(np.dot(a,a))

      # Vector normal to plane defined by q,r,s
      normal = np.cross(a, b)
      # Rotate d by the angle around the normal to the plane defined by q,r,s
      d = np.dot(rotation_matrix(normal, angle), d)

      # Rotate d around a by the dihedral
      d = np.dot(rotation_matrix(a, dihedral), d)

      # Add d to the position of q to get the new coordinates of the atom
      p = q + d
      atom = [ name, p, mass ]

      cartesian.append(atom)

def zmatrix_to_cartesian():
      '''Convert the zmartix to cartesian coordinates'''
      # Deal with first three line separately
      add_first_three_to_cartesian()

      for i in range(3, len(zmatrix)):
            add_atom_to_cartesian(zmatrix[i])

      remove_dummy_atoms()

      center_cartesian()

      return cartesian

def add_first_three_to_zmatrix(data):
      '''The first three atoms need to be treated differently'''
      # First atom
      name, position, mass = data[0]
      first = [name, [[],[],[]], mass]

      # Second atom
      name, position, mass = data[1]
      atom1 = data[0]
      pos1 = atom1[1]
      q = pos1 - position
      distance = m.sqrt(np.dot(q, q))
      second = [name, [[1,distance],[],[]], mass]

      # Third atom
      name, position, mass = data[2]
      atom1, atom2 = data[:2]
      pos1, pos2 = atom1[1], atom2[1]
      q = pos1 - position
      r = pos2 - pos1
      q_u = q / np.sqrt(np.dot(q, q))
      r_u = r / np.sqrt(np.dot(r, r))
      distance = np.sqrt(np.dot(q, q))
      # Angle between a and b = acos(dot product of the unit vectors)
      angle = m.acos(np.dot(-q_u, r_u))
      third = [name, [ [ 1, distance ], [ 2, np.degrees(angle) ], [] ], mass ]

      return [first, second, third]

def add_atom_to_zmatrix(i, line, data):
      '''Generates an atom for the zmatrix
      (assumes that three previous atoms have been placed in the cartesian coordiantes)'''
      name, position, mass = line
      atom1, atom2, atom3 = data[:3]
      pos1, pos2, pos3 = atom1[1], atom2[1], atom3[1]
      # Create vectors pointing from one atom to the next
      q = pos1 - position
      r = pos2 - pos1
      s = pos3 - pos2
      position_u = position / np.sqrt(np.dot(position, position))
      # Create unit vectors
      q_u = q / np.sqrt(np.dot(q, q))
      r_u = r / np.sqrt(np.dot(r, r))
      s_u = s / np.sqrt(np.dot(s, s))
      distance = np.sqrt(np.dot(q, q))
      # Angle between a and b = acos(dot(a, b) / (|a| |b|))
      angle = m.acos(np.dot(-q_u, r_u))
      angle_123 = m.acos(np.dot(-r_u, s_u))
      # Dihedral angle = acos(dot(normal_vec1, normal_vec2) / (|normal_vec1| |normal_vec2|))
      plane1 = np.cross(q, r)
      plane2 = np.cross(r, s)
      dihedral = m.acos(np.dot(plane1, plane2) / (np.sqrt(np.dot(plane1, plane1)) * np.sqrt(np.dot(plane2, plane2))))
      # Convert to signed dihedral angle
      if np.dot(np.cross(plane1, plane2), r_u) < 0:
            dihedral = -dihedral

      coords = [ [1, distance],[2, np.degrees(angle)],[3, np.degrees(dihedral)] ]
      atom = [ name, coords, mass ]

      return atom

def cartesian_to_zmatrix(data):
      '''Convert the cartesian coordinates to a zmatrix'''
      zmatrix = []
      firstThree = add_first_three_to_zmatrix(data)
      for i in range(3):
            zmatrix.append(firstThree[i])
      for i in range(3, len(data)):
            line = data[i]
            zmatrix.append(add_atom_to_zmatrix(i+1, line, data))

      return zmatrix

def remove_dummy_atoms():
      '''Delete any dummy atoms that may have been placed in the calculated cartesian coordinates'''
      new_cartesian = []
      for line in cartesian:
            if not line[0] == 'X':
                  new_cartesian.append(line)
      cartesian = new_cartesian

def center_cartesian():
      '''Find the center of mass and move it to the origin'''
      total_mass = 0.0
      center_of_mass = np.array([ 0.0, 0.0, 0.0 ])
      for atom in cartesian:
            mass = atom[2]
            total_mass += mass
            center_of_mass += atom[1]*mass
      center_of_mass = center_of_mass / float(total_mass)

      # Translate each atom by the center of mass
      for atom in cartesian:
            atom[1] = atom[1] - center_of_mass

def cartesian_radians_to_degrees():
      for atom in cartesian:
            atom[1][1][1] = np.degrees(atom[1][1][1])
            atom[1][2][1] = np.degrees(atom[1][2][1])

def output_cartesian( output_file='cartesian.dat'):
      '''Output the cartesian coordinates of the file'''
      with open(output_file, 'w') as f:
            f.write(str(len(cartesian)))
            f.write('\n\n')
            for line in cartesian:
                  name, position, mass = line
                  f.write(name + '\t')
                  f.write('\t'.join(str(x) for x in position))
                  f.write('\n')

def print_cartesian():
      '''Print the cartesian coordinates'''
      for line in cartesian:
            print(line[0] + '\t' + '\t'.join(str(x) for x in line[1]))

def output_zmatrix( output_file='zmatrix.dat'):
      '''Output the zmatrix to the file'''
      with open(output_file, 'w') as f:
            f.write('#ZMATRIX\n#\n')
            for line in zmatrix:
                  name, position, mass = line
                  f.write(name)
                  for i in position:
                        for j in range(0, len(i), 2):
                              f.write('\t' + str(i[j]+1) + '\t' + str(i[j+1]))
                  f.write('\n')

def print_zmatrix(data):
      '''Print the zmatrix'''
      for line in data:
            print(line[0] + '\t' + '\t'.join(str(x) for x in line[1]))

def flatten(container):
    for i in container:
        if isinstance(i, list) or isinstance(i, tuple):
            for j in flatten(i):
                yield j
        else:
            yield i
