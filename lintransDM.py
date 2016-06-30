# Program to take input structures from two Gaussian or ORCA files, and return intermediate geometries along a linear transit pathway
# Framework shared with wellfareFF.py

# --------------
# Set up portion
# --------------

import sys
import getopt
import math
import time
from importlib.util import find_spec

# Check for numpy, exit immediately if not available
module_loader = find_spec('numpy')
found = module_loader is not None
if not found:
    print("ERROR: Module numpy is required")
    timestamp("Program aborted")
    sys.exit()
import numpy as np

# Check for scipy, exit immediately if not available
module_loader = find_spec('scipy')
found = module_loader is not None
if not found:
    ProgramError("Module scipy is required")
    ProgramAbort()
import scipy.optimize

# Check for argparse, exit immediately if not available
module_loader = find_spec('argparse')
found = module_loader is not None
if not found:
    print("ERROR: Module argparse is required")
    timestamp("Program aborted")
    sys.exit()
import argparse

def timestamp(s):
    print(s + time.strftime("%Y/%m/%d %X"))

# -----------------------------
# Define command line arguments
# -----------------------------

parser = argparse.ArgumentParser(
    description="linear transit program (accompanying WellFAReFF)",
    epilog="recognised filetypes: g09, orca")
parser.add_argument("-r", "--reactant", metavar='file', help="input file with qc data of the reactant",
                    default="g09-benzene.log")
parser.add_argument("-p", "--product", metavar='file', help="input file with qc data of the product",
                    default="g09-dielsalder-p.log")
parser.add_argument("-o", "--output", metavar='file', help="output file to which geometries for linear transit can be written if desired") # Note may be easier to print to terminal than actually accommodate this option in the program
parser.add_argument("-v", "--verbosity", help="increase output verbosity", type=int, choices=[0, 1, 2, 3], default=2)
parser.add_argument("-b", "--bondcutoff", help="Cutoff value for bond identification through Mayer bond order",
                    type=float, default=0.45)
parser.add_argument("-n", "--numpoints", help="number of structures to be produced on linear transit path", default=10)

args = parser.parse_args()

# -------------------------------------------
# Define accompanying functions to wellfareFF
# -------------------------------------------

# Dictionary to convert atomic symbols to atomic numbers
SymbolToNumber = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9,
    "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17,
    "Ar": 18, "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24,
    "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30, "Ga": 31,
    "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36, "Rb": 37, "Sr": 38,
    "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44, "Rh": 45,
    "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50, "Sb": 51, "Te": 52,
    "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57, "Ce": 58, "Pr": 59,
    "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66,
    "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71, "Hf": 72, "Ta": 73,
    "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80,
    "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85, "Rn": 86, "Fr": 87,
    "Ra": 88, "Ac": 89, "Th": 90, "Pa": 91, "U": 92, "Np": 93, "Pu": 94,
    "Am": 95, "Cm": 96, "Bk": 97, "Cf": 98, "Es": 99, "Fm": 100, "Md": 101,
    "No": 102, "Lr": 103, "Rf": 104, "Db": 105, "Sg": 106, "Bh": 107,
    "Hs": 108, "Mt": 109, "Ds": 110, "Rg": 111, "Cn": 112, "Uut": 113,
    "Fl": 114, "Uup": 115, "Lv": 116, "Uus": 117, "Uuo": 118}

# Invert the above: atomic numbers to atomic symbols
NumberToSymbol = {v: k for k, v in SymbolToNumber.items()}

# Dictionary to convert atomic symbols to atomic masses
SymbolToMass = {
    "H": 1.00794, "He": 4.002602, "Li": 6.941, "Be": 9.012182, "B": 10.811,
    "C": 12.0107, "N": 14.0067, "O": 15.9994, "F": 18.9984032, "Ne": 20.1797,
    "Na": 22.98976928, "Mg": 24.3050, "Al": 26.9815386, "Si": 28.0855,
    "P": 30.973762, "S": 32.065, "Cl": 35.453, "Ar": 39.948, "K": 39.0983,
    "Ca": 40.078, "Sc": 44.955912, "Ti": 47.867, "V": 50.9415, "Cr": 51.9961,
    "Mn": 54.938045, "Fe": 55.845, "Co": 58.933195, "Ni": 58.6934, "Cu": 63.546,
    "Zn": 65.38, "Ga": 69.723, "Ge": 72.64, "As": 74.92160, "Se": 78.96,
    "Br": 79.904, "Kr": 83.798, "Rb": 85.4678, "Sr": 87.62, "Y": 88.90585,
    "Zr": 91.224, "Nb": 92.90638, "Mo": 95.96, "Tc": 98.0, "Ru": 101.07,
    "Rh": 102.90550, "Pd": 106.42, "Ag": 107.8682, "Cd": 112.411, "In": 114.818,
    "Sn": 118.710, "Sb": 121.760, "Te": 127.60, "I": 126.90447, "Xe": 131.293,
    "Cs": 132.9054519, "Ba": 137.327, "La": 138.90547, "Ce": 140.116,
    "Pr": 140.90765, "Nd": 144.242, "Pm": 145.0, "Sm": 150.36, "Eu": 151.964,
    "Gd": 157.25, "Tb": 158.92535, "Dy": 162.500, "Ho": 164.93032, "Er": 167.259,
    "Tm": 168.93421, "Yb": 173.054, "Lu": 174.9668, "Hf": 178.49, "Ta": 180.94788,
    "W": 183.84, "Re": 186.207, "Os": 190.23, "Ir": 192.217, "Pt": 195.084,
    "Au": 196.966569, "Hg": 200.59, "Tl": 204.3833, "Pb": 207.2, "Bi": 208.98040,
    "Po": 209.0, "At": 210.0, "Rn": 222.0, "Fr": 223.0, "Ra": 226.0, "Ac": 227.0,
    "Th": 232.03806, "Pa": 231.03588, "U": 238.02891, "Np": 237.0, "Pu": 244.0,
    "Am": 243.0, "Cm": 247.0, "Bk": 247.0, "Cf": 251.0, "Es": 252.0, "Fm": 257.0,
    "Md": 258.0, "No": 259.0, "Lr": 262.0, "Rf": 267.0, "Db": 268.0, "Sg": 271.0,
    "Bh": 272.0, "Hs": 270.0, "Mt": 276.0, "Ds": 281.0, "Rg": 280.0, "Cn": 285.0,
    "Uut": 284.0, "Uuq": 289.0, "Uup": 288.0, "Uuh": 293.0, "Uuo": 294.0}

# Define dictionary to convert atomic symbols to covalent radii (in Angstrom)
SymbolToRadius = {
    "H": 0.37, "He": 0.32, "Li": 1.34, "Be": 0.90, "B": 0.82, "C": 0.77,
    "N": 0.75, "O": 0.73, "F": 0.71, "Ne": 0.69, "Na": 1.54, "Mg": 1.30,
    "Al": 1.18, "Si": 1.11, "P": 1.06, "S": 1.02, "Cl": 0.99, "Ar": 0.97,
    "K": 1.96, "Ca": 1.74, "Sc": 1.44, "Ti": 1.36, "V": 1.25, "Cr": 1.27,
    "Mn": 1.39, "Fe": 1.25, "Co": 1.26, "Ni": 1.21, "Cu": 1.38, "Zn": 1.31,
    "Ga": 1.26, "Ge": 1.22, "As": 1.19, "Se": 1.16, "Br": 1.14, "Kr": 1.10,
    "Rb": 2.11, "Sr": 1.92, "Y": 1.62, "Zr": 1.48, "Nb": 1.37, "Mo": 1.45,
    "Tc": 1.56, "Ru": 1.26, "Rh": 1.35, "Pd": 1.31, "Ag": 1.53, "Cd": 1.48,
    "In": 1.44, "Sn": 1.41, "Sb": 1.38, "Te": 1.35, "I": 1.33, "Xe": 1.30,
    "Cs": 2.25, "Ba": 1.98, "La": 1.69, "Ce": 1.70, "Pr": 1.70, "Nd": 1.70,
    "Pm": 1.70, "Sm": 1.70, "Eu": 1.70, "Gd": 1.70, "Tb": 1.70, "Dy": 1.70,
    "Ho": 1.70, "Er": 1.70, "Tm": 1.70, "Yb": 1.70, "Lu": 1.60, "Hf": 1.50,
    "Ta": 1.38, "W": 1.46, "Re": 1.59, "Os": 1.28, "Ir": 1.37, "Pt": 1.28,
    "Au": 1.44, "Hg": 1.49, "Tl": 1.48, "Pb": 1.47, "Bi": 1.46, "Po": 1.50,
    "At": 1.50, "Rn": 1.45, "Fr": 1.50, "Ra": 1.50, "Ac": 1.50, "Th": 1.50,
    "Pa": 1.50, "U": 1.50, "Np": 1.50, "Pu": 1.50, "Am": 1.50, "Cm": 1.50,
    "Bk": 1.50, "Cf": 1.50, "Es": 1.50, "Fm": 1.50, "Md": 1.50, "No": 1.50,
    "Lr": 1.50, "Rf": 1.50, "Db": 1.50, "Sg": 1.50, "Bh": 1.50, "Hs": 1.50,
    "Mt": 1.50, "Ds": 1.50, "Rg": 1.50, "Cn": 1.50, "Uut": 1.50, "Uuq": 1.50,
    "Uup": 1.50, "Uuh": 1.50, "Uus": 1.50, "Uuo": 1.50}

def distance(v1=[0.0, 1.0, 0.0], v2=[1.0, 0.0, 0.0]):
    dist = (v1[0] - v2[0]) * (v1[0] - v2[0])
    dist = dist + (v1[1] - v2[1]) * (v1[1] - v2[1])
    dist = dist + (v1[2] - v2[2]) * (v1[2] - v2[2])
    dist = math.sqrt(dist)
    return dist

def ProgramWarning():
    print("\n###################################################################")
    print(" __          __              _             ")
    print(" \ \        / /             (_)            ")
    print("  \ \  /\  / /_ _ _ __ _ __  _ _ __   __ _ ")
    print("   \ \/  \/ / _` | '__| '_ \| | '_ \ / _` |")
    print("    \  /\  / (_| | |  | | | | | | | | (_| |")
    print("     \/  \/ \__,_|_|  |_| |_|_|_| |_|\__, |")
    print("                                      __/ |")
    print("                                     |___/ ")
    timestamp('Warning time/date: ')
    print("###################################################################")
    return

# ------------------------------------------------------------
# Define simplified molecule architecture based off wellfareFF
# ------------------------------------------------------------

class Atom:
    """ An atom with an atomic symbol and cartesian coordinates"""

    def __init__(self, sym, x, y, z, q):
        """ (Atom, int, str, number, number, number) -> NoneType
    
    Create an Atom with (string) symbol sym,
    and (float) cartesian coordinates (x, y, z).
    (int) charge charge and (float) mass mass are set
    automatically according to symbol.
    """

        self.symbol = sym
        self.charge = SymbolToNumber[sym]
        self.mass = SymbolToMass[sym]
        self.coord = [x, y, z]
        self.QMcharge = q  # Extracting q from input file yet to be implemented
        self.basis = []

    def __str__(self):
        """ (Atom) -> str

        Return a string representation of this Atom in this format:

          (SYMBOL, X, Y, Z, [basis])
        """
        s = ""
        for i in self.basis:
            s = s + i.__str__()
            if i != self.basis[-1]:
                s = s + " "

        return '({0}, ({1}, {2}, {3}), ({4})'.format(self.symbol, self.coord[0], self.coord[1], self.coord[2], s)

    def __repr__(self):
        """ (Atom) -> str
    
        Return a string representation of this Atom in this format:"
    
        Atom("SYMBOL", charge, QM charge, mass, X, Y, Z)
        """

        return '("{0}", {1}, {2}, {3}, {4}, {5} {6})'.format(self.symbol, self.charge, self.QMcharge, self.mass,
                                                             self.coord[0], self.coord[1], self.coord[2])
   
    def x(self, x):
        """ (Atom) -> NoneType
    
    Set x coordinate to x
    """
        self.coord[0] = x

    def y(self, y):
        """ (Atom) -> NoneType
    
    Set y coordinate to y
    """
        self.coord[1] = y

    def z(self, z):
        """ (Atom) -> NoneType
    
    Set z coordinate to z
    """
        self.coord[2] = z

    def setq(self, q):
        """ (Atom) -> NoneType
    
    Set QMcharge to q
    """
        self.QMcharge = q

class Molecule:
    """A molecule with a name, charge and a list of atoms"""

    def __init__(self, name, charge=0):
        """ (Molecule, str, int) -> NoneType
    
    Create a Molecule named name with charge charge and no atoms
    (int) Multiplicity mult is automatically set to the lowest
    possible value (1 or 2) and the lists of bonds, angles and
    dihedrals are empty.
    """

        self.name = name
        self.charge = charge
        self.mult = 1
        self.Ee_QM = 0.0  # Initially set to a placeholder value for type only
        self.atoms = []
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.threefolds = []
        self.H_QM = np.zeros((3, 3))  # Array size arbitrary, just a placeholder for type 

    def addAtom(self, a):
        """ (Molecule, Atom) -> NoneType
    
    Add a to my list of Atoms.
    """

        self.atoms.append(a)
        nucchg = 0
        for i in self.atoms:
            nucchg = nucchg + i.charge
        if (nucchg - self.charge) % 2 != 0:
            self.mult = 2
        else:
            self.mult = 1

    def mass(self):
        """ (Molecule) -> number
    
    Return the molar mass as sum of atomic masses
    """

        mass = 0.0
        for atom in self.atoms:
            mass = mass + atom.mass

        return mass

    def numatoms(self):
        """ (Molecule) -> int
    
    Return the number of atoms in the molecule
    """

        return int(len(self.atoms))

    def chatom(self, n, at):
        """ (Molecule) -> NoneType
    
    Change the nth atom of the Molecule
    """

        self.atoms[n] = at

    def movatom(self, n, x, y, z):
        """ (Molecule) -> NoneType
    
    Move the nth atom to position x, y, z
    """

        self.atoms[n].x(x)
        self.atoms[n].y(y)
        self.atoms[n].z(z)

    def setGeometry(self, cartCoordinates):
        """ (Molecule) -> NoneType

    Move all atoms in the molecule to give the geometry defined by the input coordinates
    """
        n = int(len(self.atoms))
        # Print number of atoms
        # print(str(n))
        if len(cartCoordinates) != 3 * n:
            ProgramWarning()
            print("Cannot update geometry of " + str(self.name) + ", supplied coordinates do not match number of atoms")
        else:
            for i in range(n):
                self.movatom(i, cartCoordinates[3*i], cartCoordinates[3 * i+1], cartCoordinates[3 * i+2])

    def atmmass(self, n, m):
        """ (Molecule) -> NoneType
    
        Change the mass of the nth atom to m
        """

        self.atoms[n].mass = m

    def atmatmdist(self, i, j):
        """ (Molecule) -> number
    
    Report the distance between atoms i and j
    """

        distance = (self.atoms[i].coord[0] - self.atoms[j].coord[0]) * (self.atoms[i].coord[0] - self.atoms[j].coord[0])
        distance = distance + (self.atoms[i].coord[1] - self.atoms[j].coord[1]) * (
            self.atoms[i].coord[1] - self.atoms[j].coord[1])
        distance = distance + (self.atoms[i].coord[2] - self.atoms[j].coord[2]) * (
            self.atoms[i].coord[2] - self.atoms[j].coord[2])

        return math.sqrt(distance)

    def bondangle(self, i):
        """ (Molecule) -> number (in radians)

    Report the angle described by three atoms in the bonds list
    """

        # Calculate the distance between each pair of atoms
        angle = self.angles[i]
        d_bond_1 = self.atmatmdist(angle[0], angle[1])
        d_bond_2 = self.atmatmdist(angle[1], angle[2])
        d_non_bond = self.atmatmdist(angle[0], angle[2])

        # Use those distances and the cosine rule to calculate bond angle theta
        numerator = d_bond_1 ** 2 + d_bond_2 ** 2 - d_non_bond ** 2
        denominator = 2 * d_bond_1 * d_bond_2
        argument = numerator / denominator
        theta = np.arccos(argument)

        return theta

    def anybondangle(self, i, j, k):
        """ (Molecule) -> number (in radians)

    Report the angle described by three atoms i, j and k
    """

        # Calculate the distance between each pair of atoms
        d_bond_1 = self.atmatmdist(i, j)
        d_bond_2 = self.atmatmdist(j, k)
        d_non_bond = self.atmatmdist(i, k)

        # Use those distances and the cosine rule to calculate bond angle theta
        numerator = d_bond_1 ** 2 + d_bond_2 ** 2 - d_non_bond ** 2
        denominator = 2 * d_bond_1 * d_bond_2
        argument = numerator / denominator
        theta = np.arccos(argument)

        return theta

    def dihedralangle(self, i):
        """ (Molecule) -> number (in radians)

    Report the dihedral angle described by a set of four atoms in the dihedrals list
    """

        # Calculate the vectors lying along bonds, and their cross products
        dihedral = self.dihedrals[i]
        atom_e1 = self.atoms[dihedral[0]]
        atom_b1 = self.atoms[dihedral[1]]
        atom_b2 = self.atoms[dihedral[2]]
        atom_e2 = self.atoms[dihedral[3]]
        end_1 = [atom_e1.coord[i] - atom_b1.coord[i] for i in range(3)]
        bridge = [atom_b1.coord[i] - atom_b2.coord[i] for i in range(3)]
        end_2 = [atom_b2.coord[i] - atom_e2.coord[i] for i in range(3)]
        vnormal_1 = np.cross(end_1, bridge)
        vnormal_2 = np.cross(bridge, end_2)

        # Construct a set of orthogonal basis vectors to define a frame with vnormal_2 as the x axis
        vcross = np.cross(vnormal_2, bridge)
        norm_vn2 = np.linalg.norm(vnormal_2)
        norm_b = np.linalg.norm(bridge)
        norm_vc = np.linalg.norm(vcross)
        basis_vn2 = [vnormal_2[i] / norm_vn2 for i in range(3)]
        basis_b = [bridge[i] / norm_b for i in range(3)]
        basis_cv = [vcross[i] / norm_vc for i in range(3)]

        # Find the signed angle between vnormal_1 and vnormal_2 in the new frame
        vn1_coord_n2 = np.dot(vnormal_1, basis_vn2)
        vn1_coord_vc = np.dot(vnormal_1, basis_cv)
        psi = math.atan2(vn1_coord_vc, vn1_coord_n2)

        return psi

    def anydihedralangle(self, i, j, k, l):
        """ (Molecule) -> number (in radians)

    Report the dihedral angle described by a set of four atoms i, j, k and l
    """

        # Calculate the vectors lying along bonds, and their cross products
        atom_e1 = self.atoms[i]
        atom_b1 = self.atoms[j]
        atom_b2 = self.atoms[k]
        atom_e2 = self.atoms[l]
        end_1 = [atom_e1.coord[i] - atom_b1.coord[i] for i in range(3)]
        bridge = [atom_b1.coord[i] - atom_b2.coord[i] for i in range(3)]
        end_2 = [atom_b2.coord[i] - atom_e2.coord[i] for i in range(3)]
        vnormal_1 = np.cross(end_1, bridge)
        vnormal_2 = np.cross(bridge, end_2)

        # Construct a set of orthogonal basis vectors to define a frame with vnormal_2 as the x axis
        vcross = np.cross(vnormal_2, bridge)
        norm_vn2 = np.linalg.norm(vnormal_2)
        norm_b = np.linalg.norm(bridge)
        norm_vc = np.linalg.norm(vcross)
        basis_vn2 = [vnormal_2[i] / norm_vn2 for i in range(3)]
        basis_b = [bridge[i] / norm_b for i in range(3)]
        basis_cv = [vcross[i] / norm_vc for i in range(3)]

        # Find the signed angle between vnormal_1 and vnormal_2 in the new frame
        vn1_coord_n2 = np.dot(vnormal_1, basis_vn2)
        vn1_coord_vc = np.dot(vnormal_1, basis_cv)
        psi = math.atan2(vn1_coord_vc, vn1_coord_n2)

        return psi

    def rotateMoleculeArbAxis(self, v1=[0.0, 0.0, 0.0], v2=[0.0, 0.0, 1.0], angle=90.0):
        """ (Molecule) -> none

    Rotate all atoms of a molecule around axis defined by v1 and v2 by angle (in degrees)
    """
        # Debug only
        print("\nCalled rotate function with molecule " + str(self.name))
        print("Using input vectors and angle:")
        print(v1)
        print(v2)
        print(angle)

        # Translate v2 into coordinate origin
        #self.translateMolecule(-v2[0], -v2[1], -v2[2])
        # Debug only: print new geometry
        #print("\nrotating, molecule translated to put v2 in coordinate origin")
        #print(self.xyzString())
        print("Skipping translation of v2")
        u = v1[0] - v2[0]
        v = v1[1] - v2[1]
        w = v1[2] - v2[2]
        L = (u ** 2) + (v ** 2) + (w ** 2)
        pi = np.radians(angle)

        RotMatrix = []
        Rxx = (u ** 2 + (v ** 2 + w ** 2) * np.cos(pi)) / L
        Rxy = ((u * v) * (1 - np.cos(pi)) - w * np.sqrt(L) * np.sin(pi)) / L
        Rxz = ((u * w) * (1 - np.cos(pi)) + v * np.sqrt(L) * np.sin(pi)) / L
        # Rxt = (()*(1-np.cos(pi))+())/L

        Ryx = ((u * v) * (1 - np.cos(pi)) + w * np.sqrt(L) * np.sin(pi)) / L
        Ryy = (v ** 2 + (u ** 2 + w ** 2) * np.cos(pi)) / L
        Ryz = ((v * w) * (1 - np.cos(pi)) - u * np.sqrt(L) * np.sin(pi)) / L
        # Ryt = 0.0

        Rzx = ((u * w) * (1 - np.cos(pi)) - v * np.sqrt(L) * np.sin(pi)) / L
        Rzy = ((v * w) * (1 - np.cos(pi)) + u * np.sqrt(L) * np.sin(pi)) / L
        Rzz = (w ** 2 + (u ** 2 + v ** 2) * np.cos(pi)) / L
        # Rzt = 0.0

        # Rtx = 0.0
        # Rty = 0.0
        # Rtz = 0.0
        # Rtt = 1.0

        RotMatrix.append([Rxx, Rxy, Rxz])
        RotMatrix.append([Ryx, Ryy, Ryz])
        RotMatrix.append([Rzx, Rzy, Rzz])
        # RotMatrix.append([Rtx, Rty, Rtz, Rtt])
        RotMatrix = np.matrix(RotMatrix)

        for i in self.atoms:
            vector = [i.coord[0], i.coord[1], i.coord[2]]
            vector = np.matrix(vector)
            vector = RotMatrix.dot(np.matrix.transpose(vector))
            vector = np.array(vector).flatten().tolist()
            i.coord[0] = vector[0]
            i.coord[1] = vector[1]
            i.coord[2] = vector[2]
        # Debug only
        print("Geometry after rotation before the reverse translation")
        print(self.xyzString())

        # Undo: Translate v2 out of coordinate origin
        #self.translateMolecule(v2[0], v2[1], v2[2])
        # Debug only
        #print("Geometry after translating back out of coordinate origin")
        #print(self.xyzString())
        print("Skipping translate step")

        return

    def rotateMoleculeGivenAxis(self, vrot=[0, 0.5, 0.5], angle=90.0):
        """ (Molecule) -> none

    Rotate all atoms of a molecule around a specified vector through the origin
        """
 
        # Convert angle into radians
        theta = math.radians(angle)

        # Normalise the given vector to be axis of rotation
        norm_vrot = np.linalg.norm(vrot)
        vrotn = [vrot[i] / norm_vrot for i in range(len(vrot))]

        # Separate out the coordinates of the resulting unit vector
        v0 = vrotn[0]
        v1 = vrotn[1]
        v2 = vrotn[2]

        # Construct a rotation matrix for rotating about the given vector by the given angle
        R00 = v0**2 + (1 - v0**2)*np.cos(theta)
        R01 = v0*v1*(1 - np.cos(theta)) - v2*np.sin(theta)
        R02 = v0*v2*(1 - np.cos(theta)) + v1*np.sin(theta)
        R10 = v0*v1*(1 - np.cos(theta)) + v2*np.sin(theta)
        R11 = v1**2 + (1 - v1**2)*np.cos(theta)
        R12 = v1*v2*(1 - np.cos(theta)) - v0*np.sin(theta)
        R20 = v0*v2*(1 - np.cos(theta)) - v1*np.sin(theta)
        R21 = v1*v2*(1 - np.cos(theta)) + v0*np.sin(theta)
        R22 = v2**2 + (1 - v2**2)*np.cos(theta)
        RotMat = [[R00, R01, R02], [R10, R11, R12], [R20, R21, R22]]
        RotMat = np.matrix(RotMat)

        # Multiply the coordinates of each atom in the molecule on the left by the rotation matrix
        for i in range(len(self.atoms)):
            coordsi = [self.atoms[i].coord[0], self.atoms[i].coord[1], self.atoms[i].coord[2]]
            coordsi = np.matrix(coordsi)
            coordsf = RotMat.dot(np.matrix.transpose(coordsi))
            coordsf = np.array(coordsf).flatten().tolist()
            self.atoms[i].x(coordsf[0])
            self.atoms[i].y(coordsf[1])
            self.atoms[i].z(coordsf[2])

    def rotateMolecule(self, axis="z", angle=90.0):
        """ (Molecule) -> none

    Rotate all atoms of a molecule around specified coordinate axis
    """
        # First, convert angle to radians
        angle = np.radians(angle)
        if axis == "x":
            # Setup Matrix for rotation around x axis
            RotMatrix = []
            Rxx = 1.0
            Rxy = 0.0
            Rxz = 0.0
            Ryx = 0.0
            Ryy = np.cos(angle)
            Ryz = -np.sin(angle)
            Rzx = 0.0
            Rzy = np.sin(angle)
            Rzz = np.cos(angle)
            RotMatrix.append([Rxx, Rxy, Rxz])
            RotMatrix.append([Ryx, Ryy, Ryz])
            RotMatrix.append([Rzx, Rzy, Rzz])
            RotMatrix = np.matrix(RotMatrix)
        elif axis == "y":
            # Setup Matrix for rotation around y axis
            RotMatrix = []
            Rxx = np.cos(angle)
            Rxy = 0.0
            Rxz = np.sin(angle)
            Ryx = 0.0
            Ryy = 1.0
            Ryz = 0.0
            Rzx = -np.sin(angle)
            Rzy = 0.0
            Rzz = np.cos(angle)
            RotMatrix.append([Rxx, Rxy, Rxz])
            RotMatrix.append([Ryx, Ryy, Ryz])
            RotMatrix.append([Rzx, Rzy, Rzz])
            RotMatrix = np.matrix(RotMatrix)
        else:
            # Default: Setup Matrix for rotation around z axis
            RotMatrix = []
            Rxx = np.cos(angle)
            Rxy = -np.sin(angle)
            Rxz = 0.0
            Ryx = np.sin(angle)
            Ryy = np.cos(angle)
            Ryz = 0.0
            Rzx = 0.0
            Rzy = 0.0
            Rzz = 1.0
            RotMatrix.append([Rxx, Rxy, Rxz])
            RotMatrix.append([Ryx, Ryy, Ryz])
            RotMatrix.append([Rzx, Rzy, Rzz])
            RotMatrix = np.matrix(RotMatrix)
        # Apply the rotation matrix to all atoms
        for i in self.atoms:
            vector = [i.coord[0], i.coord[1], i.coord[2]]
            vector = np.matrix(vector)
            vector = RotMatrix.dot(np.matrix.transpose(vector))
            vector = np.array(vector).flatten().tolist()
            i.coord[0] = vector[0]
            i.coord[1] = vector[1]
            i.coord[2] = vector[2]

        return

    def translateMolecule(self, x, y, z):
        """ (Molecule) -> none

    Translate all atoms of a molecule along the three axes
    """
        for i in self.atoms:
            i.coord[0] = i.coord[0] + x
            i.coord[1] = i.coord[1] + y
            i.coord[2] = i.coord[2] + z

        return

    def orient(self):
        """ (Molecule) -> NoneType
    
    Translate centre of mass to coordinate origin and
    (re-)Orient the molecule along the principal axes of inertia.
    """

        # The molecular center of mass
        xValue = 0.0
        yValue = 0.0
        zValue = 0.0
        for i in self.atoms:
            xValue = xValue + (i.mass * i.coord[0])
            yValue = yValue + (i.mass * i.coord[1])
            zValue = zValue + (i.mass * i.coord[2])
        xValue = xValue / (self.mass())
        yValue = yValue / (self.mass())
        zValue = zValue / (self.mass())

        # Translate whole molecule into the center of mass reference frame
        for i in self.atoms:
            i.coord[0] = i.coord[0] - xValue
            i.coord[1] = i.coord[1] - yValue
            i.coord[2] = i.coord[2] - zValue

        # Build inertia tensor
        inertiaTensor = []
        Ixx = 0.0
        Ixy = 0.0
        Ixz = 0.0
        Iyx = 0.0
        Iyy = 0.0
        Iyz = 0.0
        Izx = 0.0
        Izy = 0.0
        Izz = 0.0
        for i in self.atoms:
            Ixx = Ixx + (
                i.mass * (
                    (Ang2Bohr(i.coord[1]) * Ang2Bohr(i.coord[1])) + (Ang2Bohr(i.coord[2]) * Ang2Bohr(i.coord[2]))))
            Ixy = Ixy - i.mass * Ang2Bohr(i.coord[0]) * Ang2Bohr(i.coord[1])
            Ixz = Ixz - i.mass * Ang2Bohr(i.coord[0]) * Ang2Bohr(i.coord[2])
            Iyx = Iyx - i.mass * Ang2Bohr(i.coord[1]) * Ang2Bohr(i.coord[0])
            Iyy = Iyy + (
                i.mass * (
                    (Ang2Bohr(i.coord[0]) * Ang2Bohr(i.coord[0])) + (Ang2Bohr(i.coord[2]) * Ang2Bohr(i.coord[2]))))
            Iyz = Iyz - i.mass * Ang2Bohr(i.coord[1]) * Ang2Bohr(i.coord[2])
            Izx = Izx - i.mass * Ang2Bohr(i.coord[2]) * Ang2Bohr(i.coord[0])
            Izy = Izy - i.mass * Ang2Bohr(i.coord[2]) * Ang2Bohr(i.coord[1])
            Izz = Izz + (
                i.mass * (
                    (Ang2Bohr(i.coord[0]) * Ang2Bohr(i.coord[0])) + (Ang2Bohr(i.coord[1]) * Ang2Bohr(i.coord[1]))))
        inertiaTensor.append([Ixx, Ixy, Ixz])
        inertiaTensor.append([Iyx, Iyy, Iyz])
        inertiaTensor.append([Izx, Izy, Izz])
        inertiaTensor = np.matrix(inertiaTensor)

        # Diagonalise inertia tensor
        inertiaMoments, inertialAxes = np.linalg.eig(inertiaTensor)

        # Orthogonalise eigenvectors (only sometimes necessary)...
        inertialAxes, r = np.linalg.qr(inertialAxes)

        # Sort moments from highest to lowest
        idx = inertiaMoments.argsort()[::-1]
        inertiaMoments = inertiaMoments[idx]
        inertialAxes = inertialAxes[:, idx]

        # Transform molecular coordinates into new frame of principal axes of inertia
        for i in self.atoms:
            vector = [i.coord[0], i.coord[1], i.coord[2]]
            vector = np.matrix(vector)
            vector = np.matrix.transpose(inertialAxes).dot(np.matrix.transpose(vector))
            vector = np.array(vector).flatten().tolist()
            i.coord[0] = vector[0]
            i.coord[1] = vector[1]
            i.coord[2] = vector[2]

    def addBond(self, a, b):
        """ (Molecule) -> NoneType
    
    Adds a bond between atoms a and b to the list of bonds
    """
        # Make sure a < b
        if a < b:
            c = a
            d = b
        else:
            c = b
            d = a

        # Check if the bond already exists
        exists = False
        for i in self.bonds:
            if i == [c, d]:
                exists = True

        # Append bond to list if doesn't exist and is plausible
        if exists == False and a >= 0 and b >= 0 and a <= len(self.atoms) and b <= len(self.atoms) and c != d:
            self.bonds.append([c, d])

    def delBond(self, a, b):
        """ (Molecule) -> NoneType
    
    Deletes the bond between atoms a and b from the list of bonds
    """
        # Make sure a < b
        if a < b:
            c = a
            d = b
        else:
            c = b
            d = a

        # Check if the bond actually exists
        exists = False
        for i in self.bonds:
            if i == [c, d]:
                exists = True

        # Remove if it does
        if exists == True:
            self.bonds.remove([c, d])

    def addAngle(self, a, b, c):
        """ (Molecule) -> NoneType
    
    Adds an angle between atoms a, b and c to the list of angles
    """

        # Check if the angle already exists
        exists = False
        for i in self.angles:
            if i == [a, b, c]:
                exists = True

        # Append angle to list if doesn't exist and is plausible
        # (it's a bit unsatisfactory, but I don't think there are
        # better sanity checks)
        if exists == False and a >= 0 and b >= 0 and c >= 0 and a <= len(self.atoms) and b <= len(
                self.atoms) and c <= len(self.atoms) and a != b and a != c and b != c:
            self.angles.append([a, b, c])

    def addDihedral(self, a, b, c, d):
        """ (Molecule) -> NoneType
    
    Adds a dihedral between atoms a, b, c and d to the list of dihedrals
    """

        # Check if the dihedral already exists
        exists = False
        for i in self.dihedrals:
            if i == [a, b, c, d]:
                exists = True

        # Append dihedral to list if doesn't exist and is plausible
        # (it's a bit unsatisfactory, but I don't think there are
        # better sanity checks)
        if exists == False and a >= 0 and b >= 0 and c >= 0 and d >= 0 and a <= len(self.atoms) and b <= len(
                self.atoms) and c <= len(self.atoms) and d <= len(
            self.atoms) and a != b and a != c and a != d and b != c and b != d and c != d:
            self.dihedrals.append([a, b, c, d])

    def cartesianCoordinates(self):
        """ (Molecule) ->

    Returns a list containing the cartesian coordinates.
    """

        coord = []
        for i in self.atoms:
            coord.append(i.coord[0])
            coord.append(i.coord[1])
            coord.append(i.coord[2])

        return coord

    def xyzString(self):
        """ (Molecule) -> str
    
    Returns a string containing an xyz file
    """

        s = str(self.numatoms()) + "\n" + self.name + "\n"
        for i in self.atoms:
            t = "{:<3} {: .8f} {: .8f} {: .8f}\n".format(i.symbol, i.coord[0], i.coord[1], i.coord[2])
            s = s + t

        return s

    def gamessString(self):
        """ (Molecule) -> str
    
      Returns a string containing cartesian coordinates in Gamess format
    """

        s = " $DATA\n" + self.name + "\nC1\n"
        for i in self.atoms:
            t = "{:<3} {:<3d} {: .8f} {: .8f} {: .8f}\n".format(i.symbol, i.charge, i.coord[0], i.coord[1], i.coord[2])
            s = s + t
        s = s + " $END\n"
        return s

    def gaussString(self):
        """ (Molecule) -> str
    
      Returns a string containing cartesian coordinates in Gaussian format
    """

        s = "\n" + self.name + "\n\n" + str(self.charge) + " " + str(self.mult) + "\n"
        for i in self.atoms:
            t = "{:<3} {: .8f} {: .8f} {: .8f}\n".format(i.symbol, i.coord[0], i.coord[1], i.coord[2])
            s = s + t
        s = s + "\n"
        return s

    def gaussStringatX(self, cartCoordinates):
        """ (Molecule) -> str

    Returns a string in Gaussian format with the given cartesian coordinates for the molecule
    """

        s = "\n" + self.name + "\n\n" + str(self.charge) + " " + str(self.mult) + "\n"
        for i in range(len(self.atoms)):
            atm = self.atoms[i]
            t = "{:<3} {: .8f} {: .8f} {: .8f}\n".format(atm.symbol, cartCoordinates[3 * i + 0],
                                                         cartCoordinates[3 * i + 1], cartCoordinates[3 * i + 2])
            s = s + t
        s = s + "\n"
        return s

    def distMatrix(self):
        """ Function to return the upper triangular distance matrix """
 
        distmat = np.zeros((self.numatoms(), self.numatoms()))
        for i in range(len(self.atoms)):
            for j in range(i, len(self.atoms)):
                distmat[i][j] = self.atmatmdist(i, j)
        #distmat = []
        #for i in range(len(self.atoms)):
        #    for j in range(i, len(self.atoms)):
        #        dist = self.atmatmdist(i, j)
        #        distmat.append(dist)

        return distmat

    def prettydistMatrix(self):
        """ Function to print out the calculated distance matrix in a 'nice' format"""

        distmat = np.zeros((self.numatoms(), self.numatoms()))
        for i in range(len(self.atoms)):
            for j in range(i, len(self.atoms)):
                distmat[i][j] = self.atmatmdist(i, j)
        print("Distance matrix without row/ column labelling")
        print(distmat)

        #distmat = []
        #for i in range(len(self.atoms)):
        #    for j in range(i, len(self.atoms)):
        #        dist = self.atmatmdist(i, j)
        #        distmat.append(dist)
        #s = "    "
        #for i in range(len(self.atoms)):
        #    t = "  { }  ".format(i)
        #    s = s + t
        #print(s)
        #s = "    "
        #for i in range(len(self.atoms)):
        #    t = "  {0}  ".format(self.atoms[i].symbol)
        #    s = s + t
        #print(s)

    def prettyDMgiven(self, distmat):

        pass

    def setGeomDM(self, distmat, verbosity=args.verbosity):
        """ Function to set the geometry of the molecule from a given distance matrix """    

        # Note - input is a matrix
        # Working by creating a cartesian coordinate string, then using the standard setGeometry function 
        print("\n\nBeginning routine to set geometry from given distance matrix")

        numats = len(self.atoms)
        print(str(numats) + " atoms to place")

        cartcos = []
        # Place atom index 0 at the origin
        cartcos.append(0.0)
        cartcos.append(0.0)
        cartcos.append(0.0)
        if verbosity >= 2:
            print("Coordinates done for atom 0: " + str([0.0, 0.0, 0.0]) + "\n")

        # Place atom index 1 in the positive x axis
        if 1 < numats:
            at1xco = distmat[0][1]
            cartcos.append(at1xco)
            cartcos.append(0.0)
            cartcos.append(0.0)
            if verbosity >= 2:
                print("Coordinates done for atom 1 " + str([at1xco, 0.0, 0.0]) + "\n")
                print("Distance from atom 0: {0} calculated, {1} from matrix".format(distance([at1xco, 0.0, 0.0], [0.0, 0.0, 0.0]), distmat[0][1]))

        # Set up to position atom 2
        if 2 < numats:
            r01 = distmat[0][1]
            r02 = distmat[0][2]
            r12 = distmat[1][2]
            # Check if this atom also lies on the x axis
            lin012 = False
            if r02 == r01 + r12:
                at2xco = r01 + r12
                cartcos.append(at2xco)
                cartcos.append(0.0)
                cartcos.append(0.0)
                if verbosity >= 2:
                    print("Coordinates done for (linear) atom 2 " + str([at2xco, 0.0, 0.0] + "\n"))
                for j in range(0, 2):
                    if verbosity >= 3:
                        print("Atom {0} to {1} distance: {2} calculated, {3} from matrix".format(2, j, distance([at2xco, 0.0, 0.0], [cartcos[3*j], cartcos[3*j + 1], cartcos[3*j + 2]]), distmat[j][2]))
                    lin012 = True
            elif r01 == r02 + r12:
                at2xco = r01 - r12
                cartcos.append(at2xco)
                cartcos.append(0.0)
                cartcos.append(0.0)
                if verbosity >= 2:
                    print("Coordinates done for (linear) atom 2 " + str([at2xco, 0.0, 0.0]) + "\n")
                for j in range(0, 2):
                    if verbosity >= 3:
                        print("Atom {0} to {1} distance: {2} calculated, {3} from matrix".format(2, j, distance([at2xco, 0.0, 0.0], [cartcos[3*j], cartcos[3*j + 1], cartcos[3*j + 2]]), distmat[j][2]))
                lin012 = True
            elif r12 == r01 + r02:
                at2xco = -r02
                cartcos.append(at2xco)
                cartcos.append(0.0)
                cartcos.append(0.0)
                if verbosity >= 2:
                    print("Coordinates done for (linear) atom 2 " + str([at2xco, 0.0, 0.0]) + "\n")
                for j in range(0, 2):
                    if verbosity >= 3:
                        print("Atom {0} to {1} distance: {2} calculated, {3} from matrix".format(2, j, distance([at2xco, 0.0, 0.0], [cartcos[3*j], cartcos[3*j + 1], cartcos[3*j + 2]]), distmat[j][2]))
                lin012 = True
            # If not linear, determine the angle between atom 2 and the positive x axis, then calculate coordinates to place this atom in the xy plane (positive quadrant)
            if not lin012:
                a102 = math.acos((r02**2 + r01**2 - r12**2)/(2 * r02 * r01))
                at2xco = r02 * math.cos(a102)
                at2yco = r02 * math.sin(a102)
                if at2yco < 0: # NOTE: this step may not be needed
                    at2yco = -at2yco
                cartcos.append(at2xco)
                cartcos.append(at2yco)
                cartcos.append(0.0)
                atxyi = 2
                if verbosity >= 2:
                    print("atxyi assigned")
                atxyco = [at2xco, at2yco, 0.0]
                if verbosity >= 2:
                    print("Coordinates done for atom 2: " + str([at2xco, at2yco, 0.0]) + "\n")
                for j in range(0, 2):
                    rj2 = distmat[j][2]
                    if verbosity >= 3:
                        print("rj2 = " + str(rj2))
                    cartvecj = [cartcos[3*j], cartcos[3*j + 1], cartcos[3*j + 2]]
                    distancej2 = distance([at2xco, at2yco, 0.0], cartvecj)
                    if verbosity >= 3:
                        print("distancej2 = " + str(distancej2))
                    if verbosity >= 3:
                        print("Atom {0} to {1} distance: {2} calculated, {3} from matrix".format(2, j, distancej2, rj2))
           
        # If all atoms so far are linear, repeat until an atom is found not lying on the x axis and place it in the xy plane
        linear = lin012
        index = 3
        if verbosity >= 3:
            print("linear set: " + str(linear))
            print("index: " + str(index))
        while linear and index < numats:
            for i in range(index, numats):
                r0i = distmat[0][i]
                r01 = distmat[0][1]
                r1i = distmat[1][i]
                if r0i == r01 + r1i:
                    atxco = r01 + r1i
                    cartcos.append(atxco)
                    cartcos.append(0.0)
                    cartcos.append(0.0)
                    if verbosity >= 2:
                        print("Coordinates done for (linear) atom {0}: {1} \n".format(index, [atxco, 0.0, 0.0]))
                    for j in range(0, i):
                        if verbosity >= 3:
                            print("Atom {0} to {1} distance: {2} calculated, {3} from matrix".format(i, j, distance([atxco, 0.0, 0.0], [cartcos[3*j], cartcos[3*j + 1], cartcos[3*j + 2]]), distmat[j][i]))
                    index += 1
                elif r01 == r0i + r1i:
                    atxco = r01 - r1i
                    cartcos.append(atxco)
                    cartcos.append(0.0)
                    cartcos.append(0.0)
                    if verbosity >= 2:
                        print("Coordinates done for (linear) atom {0}: {1} \n".format(index, [atxco, 0.0, 0.0]))
                    for j in range(0, i):
                        if verbosity >= 3:
                            print("Atom {0} to {1} distance: {2} calculated, {3} from matrix".format(i, j, distance([atxco, 0.0, 0.0], [cartcos[3*j], cartcos[3*j + 1], cartcos[3*j + 2]]), distmat[j][i]))
                    index += 1
                elif r1i == r01 + r0i:
                    atxco = -r0i
                    cartcos.append(atxco)
                    cartcos.append(0.0)
                    cartcos.append(0.0)
                    if verbosity >= 2:
                        print("Coordinates done for (linear) atom {0}: {1} \n".format(index, [atxco, 0.0, 0.0]))
                    for j in range(0, i):
                        if verbosity >= 3:
                            print("Atom {0} to {1} distance: {2} calculated, {3} from matrix".format(i, j, distance([atxco, 0.0, 0.0], [cartcos[3*j], cartcos[3*j + 1], cartcos[3*j + 2]]), distmat[j][i]))
                    index += 1
                else:
                    if verbosity >= 3:
                        print("Atom {0} not lying on x axis".format(i))
                    linear = False
                    a10i = math.acos((r0i**2 + r01**2 - r1i**2)/(2 * r0i * r01))
                    atxco = r0i * math.cos(a10i)
                    atyco = r0i * math.sin(a10i)
                    if atyco < 0: # NOTE: this step may not be needed
                        atyco = -atyco
                    cartcos.append(atxco)
                    cartcos.append(atyco)
                    cartcos.append(0.0)
                    if verbosity >= 2:
                        print("Coordinates done for atom {0}: {1} \n".format(index, [atxco, atyco, 0.0]))
                    for j in range(0, i):
                        if verbosity >= 3:
                            print("Atom {0} to {1} distance: {2} calculated, {3} from matrix".format(i, j, distance([atxco, atyco, 0.0], [cartcos[3*j], cartcos[3*j + 1], cartcos[3*j + 2]]), distmat[j][i]))
                    atxyi = i
                    if verbosity >= 2:
                        print("atxyi assigned")
                    atxyco = [atxco, atyco, 0.0]
                    index += 1
                    break

        # When there is an atom placed in the xy plane (positive quadrant) start placing the remaining atoms, choozing z positive for the first to fall out of the xy plane
        planar = True
        index2 = index
        while planar and not linear:
            for i in range(index2, len(self.atoms)):
                r0i = distmat[0][i]
                r01 = distmat[0][1]
                r1i = distmat[1][i]
                # First check if atom i lies on the x axis
                if r0i == r01 + r1i:
                    atxco = r01 + r1i
                    cartcos.append(atxco)
                    cartcos.append(0.0)
                    cartcos.append(0.0)
                    if verbosity >= 2:
                        print("Coordinates done for atom {0}: {1} \n".format(index, [atxco, 0.0, 0.0]))
                    for j in range(0, i):
                        if verbosity >= 3:
                            print("Atom {0} to {1} distance: {2} calculated, {3} from matrix".format(i, j, distance([atxco, 0.0, 0.0], [cartcos[3*j], cartcos[3*j + 1], cartcos[3*j + 2]]), distmat[j][i]))
                    index2 += 1
                elif r01 == r0i + r1i:
                    atxco = r01 - r1i
                    cartcos.append(atxco)
                    cartcos.append(0.0)
                    cartcos.append(0.0)
                    if verbosity >= 2:
                        print("Coordinates done for atom {0}: {1} \n".format(index, [atxco, 0.0, 0.0]))
                    for j in range(0, i):
                        if verbosity >= 3:
                            print("Atom {0} to {1} distance: {2} calculated, {3} from matrix".format(i, j, distance([atxco, 0.0, 0.0], [cartcos[3*j], cartcos[3*j + 1], cartcos[3*j + 2]]), distmat[j][i]))
                    index2 += 1
                elif r1i == r01 + r0i:
                    atxco = -r02
                    cartcos.append(atxco)
                    cartcos.append(0.0)
                    cartcos.append(0.0)
                    if verbosity >= 2:
                        print("Coordinates done for atom {0}: {1} \n".format(index, [atxco, 0.0, 0.0]))
                    for j in range(0, i):
                        if verbosity >= 3:
                            print("Atom {0} to {1} distance: {2} calculated, {3} from matrix".format(i, j, distance([atxco, 0.0, 0.0], [cartcos[3*j], cartcos[3*j + 1], cartcos[3*j + 2]]), distmat[j][i]))
                    index2 += 1
                else:
                    # If not assume the point is in the xy plane
                    a10i = math.acos((r0i**2 + r01**2 - r1i**2)/(2 * r0i * r01))
                    atxcogs = r0i * math.cos(a10i)
                    atycogs = r0i * math.sin(a10i)
                    # check whether this guess point in the xy plane is consistent with the point already placed there
                    rxyi = distmat[atxyi][i]
                    distxyi = (atxcogs - atxyco[0]) * (atxcogs - atxyco[0])
                    distxyi = distxyi + (atycogs - atxyco[1]) * (atycogs - atxyco[1])
                    distxyi = math.sqrt(distxyi)
                    if rxyi == distxyi:
                        cartcos.append(atxcogs)
                        cartcos.append(atycogs)
                        cartcos.append(0.0)
                        if verbosity >= 2:
                            print("Coordinates done for (xy plane) atom {0}: {1} \n".format(index, [atxcogs, atycogs, 0.0]))
                        for j in range(0, i):
                            if verbosity >= 3:
                                print("Atom {0} to {1} distance: {2} calculated, {3} from matrix".format(i, j, distance([atxcogs, atycogs, 0.0], [cartcos[3*j], cartcos[3*j + 1], cartcos[3*j + 2]]), distmat[j][i]))
                        index2 += 1
                    elif abs(rxyi - distxyi) < 10**(-10):
                        if verbosity >= 3:
                            print("Distance check for atom {0} vs xy: not equal but within 10dp")
                        cartcos.append(atxcogs)
                        cartcos.append(atycogs)
                        cartcos.append(0.0)
                        index2 += 1
                    else:
                        # if not, check whether the mirror point in the xy plane is consistent
                        atycogs2 = -atycogs
                        distxyi2 = (atxcogs - atxyco[0]) * (atxcogs - atxyco[0])
                        distxyi2 = distxyi2 + (atycogs2 - atxyco[1]) * (atycogs2 - atxyco[1])
                        distxyi2 = math.sqrt(distxyi2)
                        if rxyi == distxyi2:
                            cartcos.append(atxcogs)
                            cartcos.append(atycogs2)
                            cartcos.append(0.0)
                            if verbosity >= 2:
                                print("Coordinates done for (xy plane) atom {0}: {1} \n".format(index, [atxcogs, atycogs2, 0.0]))
                            for j in range(0, i):
                                if verbosity >= 3:
                                    print("Atom {0} to {1} distance: {2} calculated, {3} from matrix".format(i, j, distance([atxcogs, atycogs2, 0.0], [cartcos[3*j], cartcos[3*j + 1], cartcos[3*j + 2]]), distmat[j][i]))
                            index2 += 1
                        elif abs(rxyi - distxyi2) < 10**(-10):
                            if verbosity >= 3:
                                print("Distance check for atom {0} vs xy: not equal but within 10dp")
                            cartcos.append(atxcogs)
                            cartcos.append(atycogs)
                            cartcos.append(0.0)
                            index2 += 1
                        else:
                            planar = False
                            # if not, this point must be placed in the positive z hemisphere. 
                            """
                            if verbosity >= 3:
                                print("\nDetermining coordinates for first out of plane point " + str(i))
                            atxco = atxcogs
                            if verbosity >= 3:
                                print("x coordinate: " + str(atxco))
                            # To calculate the y coordinate, find the projection of this point into the xy plane
                            # Start by finding the line defined by distance from point C
                            
                            # One point on the line is where a line through the point meets the vector to defined xy point at right angles
                            r0xy = distmat[0, atxyi]
                            if verbosity >= 3:
                                print("r0xy " + str(r0xy))
                            axy0i = math.acos((r0i**2 + r0xy**2 - rxyi**2)/(2 * r0i * r0xy))
                            if verbosity >= 3:
                                print("axy0i " + str(axy0i))
                            r0p1 = r0i * math.cos(axy0i)
                            if verbosity >= 3:
                                print("r0p1 " + str(r0p1))
                            p1 = atxyco
                            np1 = np.linalg.norm(p1)
                            p1 = [(r0p1 * p1[i]) / np1 for i in range(3)]
                            if verbosity >= 3:
                                print("p1: " + str(p1))
                            # A second point on the line will lie on the x axis, the third point of the origin, p1, p2 right angled triangle
                            r1xy = distmat[1][atxyi]
                            if verbosity >= 3:
                                print("r1xy " + str(r1xy))
                            a10xy = math.acos((r0xy**2 + r01**2 - r1xy**2)/(2 * r0xy * r01))
                            if verbosity >= 3:
                                print("a10xy " + str(a10xy))
                            r0p2 = r0p1 / math.cos(a10xy)
                            if verbosity >= 3:
                                print("r0p2 " + str(r0p2))
                            p2 = [r0p2, 0.0, 0.0]
                            if verbosity >= 3:
                                print("p2: " + str(p2))
                            # Using these two points, the equation of the line can be found
                            grad = (p2[1] - p1[1]) / (p2[0] - p1[0])
                            intcpt = -grad * p2[0]
                            if verbosity >= 3:
                                print("Line equation y = {0}*x + {1}".format(grad, intcpt))
                            # Now the intersection of this line with that from the x coordinate of point i will give the y coordinate
                            atyco = (grad * atxco) + intcpt
                            if verbosity >= 3:
                                print("y coordinate: " + str(atyco))
                            # From here the z coordinate can be calculated using Pythagoras' theorem
                            atzco = math.sqrt(r0i**2 - atxco**2 - atyco**2) 
                            if atzco < 0:
                                atzco = -atzco
                            if verbosity >= 3:
                                print("z coordinate: " + str(atzco))
                            # check if the distances from this new point to each of the existing points are correct
                            match = True
                            for j in range(0, i):
                                if verbosity >= 3:
                                    print("Atom {0} to {1} distance: {2} calculated, {3} from matrix".format(i, j, distance([atxco, atyco, atzco], [cartcos[3*j], cartcos[3*j + 1], cartcos[3*j + 2]]), distmat[j][i]))
                            for j in range(0, i):
                                distji = distance([atxco, atyco, atzco], [cartcos[3*j], cartcos[3*j + 1], cartcos[3*j + 2]])
                            #    distji = (atxco - cartcos[3*j]) * (atxco - cartcos[3*j])
                            #    distji = distji + ((atyco - cartcos[3*j + 1]) * (atyco - cartcos[3*j + 1]))
                            #    distji = distji + ((atzco - cartcos[3*j + 2]) * (atzco - cartcos[3*j + 2]))
                            #    distji = math.sqrt(distji)
                            #    print("Distance from calculated point {0} to point {1}: {2}".format(i, j, distji))
                                rji = distmat[j][i]
                            #    print("and from distance matrix: " + str(rji))
                                if rji == distji:
                                    if verbosity >= 3:
                                        print("Distances match for atom {0} with atom {1}".format(i, j))
                                else:
                                    if verbosity >= 3:
                                        print("Trying 10dp tolerance for distance {0} to {1} check as not equal".format(j, i))
                                    if abs(rji - distji) < 10**(-10):
                                        if verbosity >= 3:
                                            print("Distances within tolerance")
                                    else:
                                        match = False
                                        if verbosity >= 3:
                                            print("Distance does not match for atom {0} with atom {1}".format(i, j))
                                    break
                            """
                            # Use a method based on system of 3 simultaneous equations for the distance from atom i to each of the atoms 0, 1, and xy to give the 3 cartesian coordinates of atom i
                            # The equations in question are: (1) x**2 + y**2 + z**2 = r0i**2, (2) (x - r01)**2 + y**2 + z**2 = r1i**2, (3) (x - x_xy)**2 + (y - y_xy)**2 + z**2 = rxyi**2
                            # Previously calculated x coordinate might be used as a check
                            # Simplify to solve equations (2) - (1) and (3) - (1):
                            #-2*x*r01 = r1i**2 - r0i**2 - r01**2 
                            #-2*x*x_xy -2*y*y_xy = rxyi**2 -r0i**2 - x_xy**2 - y_xy**2
                            x_xy = atxyco[0]
                            y_xy = atxyco[1]
                            rxyi = distmat[atxyi][i]
                            r0xy = distmat[0][atxyi]
                            
                            # First solve for x from equation (2)
                            atxcolgs = (r1i**2 - r0i**2 -r01**2) / (-2 * r01)
                            if verbosity >= 3:
                                print("x coordinate from simultaneous distance equations: " + str(atxcolgs))
                            # Now substitute x value into (3) for y
                            atycolgs = (rxyi**2 - r0i**2 - x_xy**2 - y_xy**2 + (2 * atxcolgs * x_xy)) / (-2 * y_xy)
                            if verbosity >= 3:
                                print("y coordinate from simultaneous distance equations: " + str(atycolgs))
                            # Next substitute these values for x and y into equation (1) for z:
                            atzcolgsq = r0i**2 - atxcolgs**2 - atycolgs**2
                            atzcolgs = math.sqrt(atzcolgsq)
                            # Now choose that the z coordinate is positive, as this is the first atom to have a nonzero z coordinate
                            if atzcolgs < 0:
                                atzcolgs = -atzcolgs 
                            if verbosity >= 3:
                                print("z coordinate (chosen > 0) from simultaneous distance equations: " + str(atzcolgs))
                            """
                            if abs(atxco - atxcolgs) < 10**(-10):
                                if verbosity >= 3:
                                    print("Values for x coordinate by the two methods are consistent to at least 10dp")
                            if abs(atyco - atycolgs) < 10**(-10):
                                if verbosity >= 3:
                                    print("Values for y coordinate by the two methods are consistent to at least 10dp")
                            if abs(atzco - atzcolgs) < 10**(-10):
                                if verbosity >= 3:
                                    print("Values for z coordinate by the two methods are consistent to at least 10dp")
                            """
                            match = True
                            for j in range(0, i):
                                if verbosity >= 3:
                                    print("Atom {0} to {1} distance: {2} calculated, {3} from matrix".format(i, j, distance([atxcolgs, atycolgs, atzcolgs], [cartcos[3*j], cartcos[3*j + 1], cartcos[3*j + 2]]), distmat[j][i]))
                            for j in range(0, i):
                                distji = distance([atxcolgs, atycolgs, atzcolgs], [cartcos[3*j], cartcos[3*j + 1], cartcos[3*j + 2]])
                                rji = distmat[j][i]
                                if rji == distji:
                                    if verbosity >= 3:
                                        print("Distances match for atom {0} with atom {1}".format(i, j))
                                else:
                                    if verbosity >= 3:
                                        print("Trying 10dp tolerance for distance {0} to {1} check as not equal".format(j, i))
                                    if abs(rji - distji) < 10**(-10):
                                        if verbosity >= 3:
                                            print("Distances within tolerance")
                                    else:
                                        match = False
                                        if verbosity >= 3:
                                            print("Distance does not match for atom {0} with atom {1}".format(i, j))
                                    break
                            if match:
                                if verbosity >= 3:
                                    print("Using coordinate values from simultaneous distance equations to assign (first out of plane) atom: " + str(i))
                                atoopi = i
                                atoopco = [atxcolgs, atycolgs, atzcolgs]
                                cartcos.append(atxcolgs)
                                cartcos.append(atycolgs)
                                cartcos.append(atzcolgs)
                                if verbosity >= 2:
                                    print("Coordinates done for (first out of plane) atom {0}: {1} \n".format(index, [atxcolgs, atycolgs, atzcolgs]))
                                index2 += 1
                                planar = False
                                break
                            else:
                                print("Error attempting to place atom in hemisphere z > 0, for index " + str(i))
                                break
                           
        # Once the first out of plane atom has been placed (atoop), place the remaining atoms with reference to atoms 0 and 1, atxy and atoop
        index3 = index2
        for i in range(index3, len(self.atoms)):
            # First check if atom i lies on the x axis
            r0i = distmat[0][i]
            r01 = distmat[0][1]
            r1i = distmat[1][i]
            if r0i == r01 + r1i:
                atxco = r01 + r1i
                cartcos.append(atxco)
                cartcos.append(0.0)
                cartcos.append(0.0)
                index3 += 1
            elif r01 == r0i + r1i:
                atxco = r01 - r1i
                cartcos.append(atxco)
                cartcos.append(0.0)
                cartcos.append(0.0)
                index3 += 1
            elif r1i == r01 + r0i:
                atxco = -r0i
                cartcos.append(atxco)
                cartcos.append(0.0)
                cartcos.append(0.0)
                index3 += 1
            else:
                # If not, assume that the new atom lies in the xy plane
                a10i = math.acos((r0i**2 + r01**2 - r1i**2)/(2 * r0i * r01))
                atxcogs = r0i * math.cos(a10i)
                atycogs = r0i * math.sin(a10i)
                # check whether this guess point in the xy plane is consistent with the point already placed there
                rxyi = distmat[atxyi, i]
                distxyi = (atxcogs - atxyco[0]) * (atxcogs - atxyco[0])
                distxyi = distxyi + (atycogs - atxyco[1]) * (atycogs - atxyco[1])
                distxyi = math.sqrt(distxyi)
                if rxyi == distxyi:
                    cartcos.append(atxcogs)
                    cartcos.append(atycogs)
                    cartcos.append(0.0)
                    index3 += 1
                elif abs(rxyi - distxyi) < 10**(-10):
                    if verbosity >= 3:
                        print("Distance check for atom {0} vs xy: not equal but within 10dp")
                    cartcos.append(atxcogs)
                    cartcos.append(atycogs)
                    cartcos.append(0.0)
                    index3 += 1
                else:
                    # if not, check whether the mirror point in the xy plane is consistent
                    atycogs2 = -atycogs
                    distxyi2 = (atxcogs - atxyco[0]) * (atxcogs - atxyco[0])
                    distxyi2 = distxyi2 + (atycogs2 - atxyco[1]) * (atycogs2 - atxyco[1])
                    distxyi2 = math.sqrt(distxyi2)
                    if rxyi == distxyi2:
                        cartcos.append(atxcogs)
                        cartcos.append(atycogs2)
                        cartcos.append(0.0)
                        index3 += 1
                    elif abs(rxyi - distxyi2) < 10**(-10):
                        if verbosity >= 3:
                            print("Distance check for atom {0} vs xy: not equal but within 10dp")
                        cartcos.append(atxcogs)
                        cartcos.append(atycogs)
                        cartcos.append(0.0)
                        index3 += 1
                    else:
                        # if not, continue to calculate coordinates for a point out of the xy plane with reference to atoms 0, 1 and xy, and use the out of plane (oop) atom to determine z sign 
                        if verbosity >= 3:
                            print("\nDetermining coordinates for out of plane point " + str(i))
                        """
                        atxco = atxcogs
                        if verbosity >= 3:
                            print("x coordinate: " + str(atxco))
                        # To calculate the y coordinate, find the projection of this point into the xy plane
                        # Start by finding the line defined by distance from point C
                        
                        # One point on the line is where a line through the point meets the vector to defined xy point at right angles
                        r0xy = distmat[0, atxyi]
                        if verbosity >= 3:
                            print("r0xy " + str(r0xy))
                        axy0i = math.acos((r0i**2 + r0xy**2 - rxyi**2)/(2 * r0i * r0xy))
                        if verbosity >= 3:
                            print("axy0i " + str(axy0i))
                        r0p1 = r0i * math.cos(axy0i)
                        if verbosity >= 3:
                            print("r0p1 " + str(r0p1))
                        p1 = atxyco
                        np1 = np.linalg.norm(p1)
                        p1 = [(r0p1 * p1[i]) / np1 for i in range(3)]
                        if verbosity >= 3:
                            print("p1: " + str(p1))
                        # A second point on the line will lie on the x axis, the third point of the origin, p1, p2 right angled triangle
                        r1xy = distmat[1][atxyi]
                        if verbosity >= 3:
                            print("r1xy " + str(r1xy))
                        a10xy = math.acos((r0xy**2 + r01**2 - r1xy**2)/(2 * r0xy * r01))
                        if verbosity >= 3:
                            print("a10xy " + str(a10xy))
                        r0p2 = r0p1 / math.cos(a10xy)
                        if verbosity >= 3:
                            print("r0p2 " + str(r0p2))
                        p2 = [r0p2, 0.0, 0.0]
                        if verbosity >= 3:
                            print("p2: " + str(p2))
                        # Using these two points, the equation of the line can be found
                        grad = (p2[1] - p1[1]) / (p2[0] - p1[0])
                        intcpt = -grad * p2[0]
                        if verbosity >= 3:
                            print("Line equation y = {0}*x + {1}".format(grad, intcpt))
                        # Now the intersection of this line with that from the x coordinate of point i will give the y coordinate
                        atyco = (grad * atxco) + intcpt
                        if verbosity >= 3:
                            print("y coordinate: " + str(atyco))
                        # From here the z coordinate can be calculated using Pythagoras' theorem
                        atzco = math.sqrt(r0i**2 - atxco**2 - atyco**2) 
                        roopi = distmat[atoopi][i]
                        distoopi = distance([atxco, atyco, atzco], atoopco)
                        if roopi != distoopi:
                            if verbosity >= 3:
                                print("Check distances for z coordinate of atom {0} not equal. Trying 10dp tolerance.".format(i))
                            if abs(roopi - distoopi) < 10**(-10):
                                if verbosity >= 3:
                                    print("Distances within tolerance")
                            else:
                                if verbosity >= 3:
                                    print("Switching z sign")
                                atzco = -atzco
                                distoopi = distance([atxco, atyco, atzco], atoopco)
                                if roopi != distoopi:
                                    if verbosity >= 3:
                                        print("Check distances for z coordinate of atom {0} not equal. Trying 10dp tolerance.".format(i))
                                    if abs(roopi - distoopi) < 10**(-10):
                                        if verbosity >= 3:
                                            print("Distances within tolerance")
                                    else:
                                        print("WARNING: Switching z sign could not reconcile distance to out of plane atom")
                        if verbosity >= 3:
                            print("z coordinate: " + str(atzco))
                        # check if the distances from this new point to each of the existing points are correct
                        match = True
                        for j in range(0, i):
                             if verbosity >= 3:
                                 print("Atom {0} to {1} distance: {2} calculated, {3} from matrix".format(i, j, distance([atxco, atyco, atzco], [cartcos[3*j], cartcos[3*j + 1], cartcos[3*j + 2]]), distmat[j][i]))
                        for j in range(0, i):
                            distji = distance([atxco, atyco, atzco], [cartcos[3*j], cartcos[3*j + 1], cartcos[3*j + 2]])
                        #    distji = (atxco - cartcos[3*j]) * (atxco - cartcos[3*j])
                        #    distji = distji + ((atyco - cartcos[3*j + 1]) * (atyco - cartcos[3*j + 1]))
                        #    distji = distji + ((atzco - cartcos[3*j + 2]) * (atzco - cartcos[3*j + 2]))
                        #    distji = math.sqrt(distji)
                        #    print("Distance from calculated point {0} to point {1}: {2}".format(i, j, distji))
                            rji = distmat[j][i]
                        #    print("and from distance matrix: " + str(rji))
                            if rji == distji:
                                if verbosity >= 3:
                                    print("Distances match for atom {0} with atom {1}".format(i, j))
                            else:
                                if verbosity >= 3:
                                    print("Trying 10dp tolerance for distance {0} to {1} check as not equal".format(j, i))
                                if abs(rji - distji) < 10**(-10):
                                    if verbosity >= 3:
                                        print("Distances within tolerance")
                                else:
                                    match = False
                                    if verbosity >= 3:
                                        print("Distance does not match for atom {0} with atom {1}".format(i, j))
                                    break                        
                        """
                        # Use a method based on system of 3 simultaneous equations for the distance from atom i to each of the atoms 0, 1, and xy to give the 3 cartesian coordinates of atom i
                        # The equations in question are: (1) x**2 + y**2 + z**2 = r0i**2, (2) (x - r01)**2 + y**2 + z**2 = r1i**2, (3) (x - x_xy)**2 + (y - y_xy)**2 + z**2 = rxyi**2
                        # Previously calculated x coordinate might be used as a check
                        # Simplify to solve equations (2) - (1) and (3) - (1):
                        #-2*x*r01 = r1i**2 - r0i**2 - r01**2 
                        #-2*x*x_xy -2*y*y_xy = rxyi**2 -r0i**2 - x_xy**2 - y_xy**2
                        x_xy = atxyco[0]
                        y_xy = atxyco[1]
                        rxyi = distmat[atxyi][i]
                        r0xy = distmat[0][atxyi]
                        # First solve for x from equation (2)
                        atxcolgs = (r1i**2 - r0i**2 -r01**2) / (-2 * r01)
                        if verbosity >= 3:
                            print("x coordinate from simultaneous distance equations: " + str(atxcolgs))
                        # Now substitute x value into (3) for y
                        atycolgs = (rxyi**2 - r0i**2 - x_xy**2 - y_xy**2 + (2 * atxcolgs * x_xy)) / (-2 * y_xy)
                        if verbosity >= 3:
                            print("y coordinate from simultaneous distance equations: " + str(atycolgs))
                        # Next substitute these values for x and y into equation (1) for z:
                        # z**2 = rxyi**2 - r0xy - x**2 - y**2 + 2*x*x_xy + 2*y*y_xy
                        #atzcolgsq = rxyi**2 - r0xy - atxcolgs**2 - atycolgs**2 + 2*atxcolgs*x_xy + 2*atycolgs*y_xy
                        atzcolgsq = r0i**2 - atxcolgs**2 - atycolgs**2
                        atzcolgs = math.sqrt(atzcolgsq)
                        # Now check distances against first out of plane atom to ensure correct sign of the z coordinate 
                        roopi = distmat[atoopi][i]
                        distoopi = distance([atxcolgs, atycolgs, atzcolgs], atoopco)
                        needopt = False
                        if roopi != distoopi:
                            if verbosity >= 3:
                                print("Check distances for z coordinate of atom {0} not equal. Trying 10dp tolerance.".format(i))
                            if abs(roopi - distoopi) < 10**(-10):
                                if verbosity >= 3:
                                    print("Distances within tolerance")
                            else:
                                atzcolgs = -atzcolgs
                                distoopi = distance([atxcolgs, atycolgs, atzcolgs], atoopco)
                                if roopi != distoopi:
                                    if verbosity >= 3:
                                        print("Check distances for z coordinate of atom {0} not equal. Trying 10dp tolerance.".format(i))
                                    if abs(roopi - distoopi) < 10**(-10):
                                        if verbosity >= 3:
                                            print("Distances within tolerance")
                                    else:
                                        # First try check again with a loosened tolerance, for each sign of z coordinate
                                        if verbosity >= 3:
                                            print("Check distances for z coordinate of atom {0} not within tolerance. Loosening to 5dp.".format(i))
                                        atzcolgs = -atzcolgs
                                        distoopi = distance([atxcolgs, atycolgs, atzcolgs], atoopco)
                                        if abs(roopi - distoopi) < 10**(-5):
                                            if verbosity >= 3:
                                                print("Distances within 5dp tolerance")
                                        else:
                                            atzcolgs = -atzcolgs
                                            distoopi = distance([atxcolgs, atycolgs, atzcolgs], atoopco)
                                            if abs(roopi - distoopi) < 10**(-5):
                                                if verbosity >= 3:
                                                    print("Distances within 5dp tolerance")
                                            else:
                                                if verbosity >= 2:
                                                    print("Switching z sign could not reconcile distance to first out of plane atom for atom " + str(index3))
                                                    print("Distance from matrix: {0}, Distance from calculated position: {1}".format(roopi, distoopi))
                                                    print("Optimisation will be used to find coordinate values that better match distances")
                                                needopt = True
                                                # Choose the sign of z coordinate which minimises the check distance discrepancy to go on to optimisation
                                                atzcolgs1 = atzcolgs
                                                atzcolgs2 = -atzcolgs
                                                distoopi1 = distance([atxcolgs, atycolgs, atzcolgs1], atoopco)
                                                distoopi2 = distance([atxcolgs, atycolgs, atzcolgs2], atoopco)
                                                diff1 = abs(roopi - distoopi1)
                                                diff2 = abs(roopi - distoopi2)
                                                if diff1 <= diff2:
                                                    atzcolgs = atzcolgs1
                                                elif diff2 < diff1:
                                                    atzcolgs = atzcolgs2
                        if verbosity >= 3 and not needopt:
                            print("z coordinate (checked against oop) from simultaneous distance equations: " + str(atzcolgs))
                        """
                        if abs(atxco - atxcolgs) < 10**(-10):
                            if verbosity >= 3:
                                print("Values for x coordinate by the two methods are consistent to at least 10dp")
                        if abs(atyco - atycolgs) < 10**(-10):
                            if verbosity >= 3:
                                print("Values for y coordinate by the two methods are consistent to at least 10dp")
                        if abs(atzco - atzcolgs) < 10**(-10):
                            if verbosity >= 3:
                                print("Values for z coordinate by the two methods are consistent to at least 10dp")
                        """
                        if needopt:
                            if verbosity >= 3:
                                print("Beginning optimisation of coordinates for atom {0} based on the difference in interatomic distances".format(i))
                            #def distdifffunc(coordsopt, ind, coordslist, distmatrix):
                            #    sumddiffs = 0.0
                            #    for j in range(ind):
                            #        r = distmatrix[j][ind]
                            #        d = distance(coordsopt, [coordslist[3*j], coordslist[3*j + 1], coordslist[3*j + 2]])
                            #        ddiff = (r-d)**2
                            #        sumddiffs += ddiff
                            #    return sumddiffs
                            print("Optimising z coordinate only, with x and y coordinates fixed")
                            def distdifffunc(z, ind, coordslist, distmatrix):
                                sumddiffs = 0.0
                                for j in range(ind):
                                    r = distmatrix[j][ind]
                                    d = distance([atxcolgs, atycolgs, z], [coordslist[3*j], coordslist[3*j + 1], coordslist[3*j + 2]])
                                    ddiff = (r-d)**2
                                    sumddiffs += ddiff
                                return sumddiffs
                            guessco = [atxcolgs, atycolgs, atzcolgs]
                            optout = scipy.optimize.minimize(distdifffunc, atzcolgs, args=(i, cartcos, distmat))        
                            optco = optout.x
                            atzcolgs = optco[0]
                            if verbosity >= 2:
                                #print("Coordinates for atom {0} following optimisation: {1}".format(i, optco))
                                print("Z coordinate for atom {0} following optimisation: {1}".format(i, atzcolgs))
                            #atxcolgs = optco[0]
                            #atycolgs = optco[1]
                            #atzcolgs = optco[2]
                            # NOTE: distances will be checked, but no further optimisation carried out if they do not match
                        match = True
                        for j in range(0, i):
                            if verbosity >= 3:
                                print("Atom {0} to {1} distance: {2} calculated, {3} from matrix".format(i, j, distance([atxcolgs, atycolgs, atzcolgs], [cartcos[3*j], cartcos[3*j + 1], cartcos[3*j + 2]]), distmat[j][i]))
                        for j in range(0, i):
                            distji = distance([atxcolgs, atycolgs, atzcolgs], [cartcos[3*j], cartcos[3*j + 1], cartcos[3*j + 2]])
                            rji = distmat[j][i]
                            if rji == distji:
                                if verbosity >= 3:
                                    print("Distances match for atom {0} with atom {1}".format(i, j))
                            else:
                                if verbosity >= 3:
                                    print("Trying 10dp tolerance for distance {0} to {1} check as not equal".format(j, i))
                                if abs(rji - distji) < 10**(-10):
                                    if verbosity >= 3:
                                        print("Distances within tolerance")
                                else:
                                    match = False
                                    if verbosity >= 2:
                                        print("Distance checks do not match for atom {0} with atom {1}".format(i, j))
                                        print("Distance from matrix: {0}, Distance from calculated position: {1}".format(roopi, distoopi))
                                        if needopt:
                                            print("As they have been optimised, current coordinates will still be used")
                                break
                        # If an optimisation has been done, use the coordinates regardless. Otherwise, use only if all distances are within tolerance
                        if match or needopt:
                            if verbosity >= 2:
                                if needopt:
                                    print("Using optimised coordinate values from simultaneous distance equations to assign (out of plane) atom: " + str(i))
                                else:
                                    print("Using coordinate values from simultaneous distance equations to assign (out of plane) atom: " + str(i))
                            cartcos.append(atxcolgs)
                            cartcos.append(atycolgs)
                            cartcos.append(atzcolgs)
                            if verbosity >= 2:
                                print("Coordinates done for (out of plane) atom {0}: {1} \n".format(i, [atxcolgs, atycolgs, atzcolgs]))
                            index3 += 1
                            planar = False
                        else:
                            print("Error attempting to place atom in hemisphere z > 0, for index " + str(i))
                            break
                        """
                        # if not, calculate the z component using the first point in the xy plane
                        r0xy = distmat[0, atxyi]
                        axy0i = math.acos((r0i**2 + r0xy**2 - rxyi**2)/(2 * r0i * r0xy))
                        atzcogs = r0i * math.sin(axy0i)
                        roopi = distmat[atoopi][i]
                        # check if the distance from xy plane point 1 is correct
                        distxyi3 = (atxcogs - atxyco[0]) * (atxcogs - atxyco[0])
                        distxyi3 = distxyi3 + (atycogs - atxyco[1]) * (atycogs - atxyco[1])
                        distxyi3 = distxyi3 + (atzcogs - atxyco[2]) * (atxcogs - atxyco[2])
                        distxyi3 = math.sqrt(distxyi3)
                        if rxyi == distxyi3:
                            # check if the distance from out of plane point is correct
                            distoopi = (atxcogs - atoopco[0]) * (atxcogs - atoopco[0])
                            distoopi = distxyi3 + (atycogs - atoopco[1]) * (atycogs - atoopco[1])
                            distoopi = distxyi3 + (atzcogs - atoopco[2]) * (atxcogs - atoopco[2])
                            distoopi = math.sqrt(distxyi3)
                            if roopi == distoopi: 
                                cartcos.append(atxcogs)
                                cartcos.append(atycogs)
                                cartcos.append(atzcogs)
                                index3 += 1
                            else:
                                atzcogs2 = -atzcogs
                                # check if the distance from out of plane point is correct
                                distoopi2 = (atxcogs - atoopco[0]) * (atxcogs - atoopco[0])
                                distoopi2 = distxyi3 + (atycogs - atoopco[1]) * (atycogs - atoopco[1])
                                distoopi2 = distxyi3 + (atzcogs2 - atoopco[2]) * (atxcogs2 - atoopco[2])
                                distoopi2 = math.sqrt(distxyi3)
                                if roopi == distoopi2:
                                    cartcos.append(atxcogs)
                                    cartcos.append(atycogs)
                                    cartcos.append(atzcogs2)
                                    index3 += 1
                                else:
                                    print("Unable to consistently place atom with index: " + str(i))
                                    break
                        else:
                            distxyi4 = (atxcogs - atxyco[0]) * (atxcogs - atxyco[0])
                            distxyi4 = distxyi4 + (atycogs2 - atxyco[1]) * (atycogs2 - atxyco[1])
                            distxyi4 = distxyi4 + (atzcogs - atxyco[2]) * (atxcogs - atxyco[2])
                            distxyi4 = math.sqrt(distxyi4)
                            if rxyi == distxyi4:
                                # check if the distance from out of plane point is correct
                                distoopi3 = (atxcogs - atoopco[0]) * (atxcogs - atoopco[0])
                                distoopi3 = distxyi3 + (atycogs2 - atoopco[1]) * (atycogs - atoopco[1])
                                distoopi3 = distxyi3 + (atzcogs - atoopco[2]) * (atxcogs - atoopco[2])
                                distoopi3 = math.sqrt(distxyi3)
                                if roopi == distoopi3: 
                                    cartcos.append(atxcogs)
                                    cartcos.append(atycogs2)
                                    cartcos.append(atzcogs)
                                    index3 += 1
                                else:
                                    atzcogs2 = -atzcogs
                                    # check if the distance from out of plane point is correct
                                    distoopi4 = (atxcogs - atoopco[0]) * (atxcogs - atoopco[0])
                                    distoopi4 = distxyi3 + (atycogs - atoopco[1]) * (atycogs - atoopco[1])
                                    distoopi4 = distxyi3 + (atzcogs2 - atoopco[2]) * (atxcogs2 - atoopco[2])
                                    distoopi4 = math.sqrt(distxyi3)
                                    if roopi == distoopi4:
                                        cartcos.append(atxcogs)
                                        cartcos.append(atycogs2)
                                        cartcos.append(atzcogs2)
                                        index3 += 1
                                    else:
                                        print("Unable to consistently place atom with index: " + str(i))
                                        break 
                        """    
        # NOTE: routine should have all parts in place now, but calculation of x coordinate is very likely wrong.

        # Check coordinate string is complete
        if len(cartcos) == 3 * len(self.atoms):
            print("Cartesian coordinates successfully generated from distance matrix")
        else:
            print("ERROR: missing some coordinates after converting distance matrix")

        # Set geometry using the cartesian coordinate string just generated
        self.setGeometry(cartcos)
        
        # (Temporary) check the distance matrix of the newly set geometry against the input
        resultDM = self.distMatrix()
        diffmat = np.subtract(resultDM, distmat)
        if verbosity >= 3:
            print("Distance matrix of new geometry - input distance matrix:")
        print(diffmat)

# -------------------------------
# Set up geometry reading section
# -------------------------------

def extractCoordinates(filename, molecule, verbosity=0, distfactor=1.3, bondcutoff=0.45):
    if verbosity >= 1:
        print("\nSetting up WellFARe molecule: ", molecule.name)
    f = open(filename, 'r')
    program = "N/A"
    # Determine which QM program we're dealing with
    for line in f:
        if line.find("Entering Gaussian System, Link 0=g09") != -1:
            if verbosity >= 1:
                print("Reading Gaussian output file: ", filename)
            program = "g09"
            break
        elif line.find("* O   R   C   A *") != -1:
            if verbosity >= 1:
                print("Reading Orca output file: ", filename)
            program = "orca"
            break
    f.close()

    # GEOMETRY AND ATOMIC CHARGE READING SECTION
    geom = []
    charges = []
    # Read through Gaussian file, read *last* "Input orientation"
    if program == "g09":
        f = open(filename, 'r')
        for line in f:
            if line.find("Input orientation:") != -1:
                if verbosity >= 2:
                    print("\nInput orientation found, reading coordinates")
                del geom[:]
                for i in range(0, 4):
                    readBuffer = f.__next__()
                while True:
                    readBuffer = f.__next__()
                    if readBuffer.find("-----------") == -1:
                        geom.append(readBuffer)
                        if verbosity >= 3:
                            readBuffer = readBuffer.split()
                            print(" Found atom: {:<3} {: .8f} {: .8f} {: .8f} in current Input orientation".format(
                                NumberToSymbol[int(readBuffer[1])], float(readBuffer[3]), float(readBuffer[4]),
                                float(readBuffer[5])))
                    else:
                        break
        if verbosity >= 1:
            print("\nReading of geometry finished.\nAdding atoms to WellFARe molecule: ", molecule.name)
        for i in geom:
            readBuffer = i.split()
            molecule.addAtom(Atom(NumberToSymbol[int(readBuffer[1])], float(readBuffer[3]), float(readBuffer[4]),
                                  float(readBuffer[5]), 0.1))  # 0.1 a placeholder for QM calculated charge on the atom
            if verbosity >= 2:
                print(" {:<3} {: .8f} {: .8f} {: .8f}".format(NumberToSymbol[int(readBuffer[1])], float(readBuffer[3]),
                                                              float(readBuffer[4]), float(readBuffer[5])))
        f.close()
        # Next read through Gaussian file, read Mulliken charges and assign to their correct atoms
        # This part will need to be tested, may be possible to integrate more into geometry reading!
        f = open(filename, 'r')
        for line in f:
            if line.find("Mulliken charges:") != -1:
                if verbosity >= 2:
                    print("\nMulliken charges found, reading charges")  # May not strictly need this
                del charges[:]
                readBuffer = f.__next__()  # Only one line to skip, but check this works with no for loop
                while True:
                    readBuffer = f.__next__()
                    if readBuffer.find(
                            "Sum of Mulliken charges") == -1:  # Check that finding less than the full line works - if not need to include the value of the sum
                        charges.append(readBuffer)
                        if verbosity >= 3:
                            readBuffer = readBuffer.split()
                            print(" Found atomic charge listing: {:<3} {:<3} {: .8f} in Mulliken charges".format(
                                int(readBuffer[0]), str(readBuffer[1]), float(readBuffer[2])))
                    else:
                        break
        if verbosity >= 1:
            print("\nReading of Mulliken charges finished. \nAdding QM atomic charges to atoms in WellFARe molecule: ",
                  molecule.name)
        for i in charges:
            readBuffer = i.split()
            n = int(readBuffer[0]) - 1
            molecule.atoms[n].setq(float(readBuffer[2]))
            if verbosity >= 2:
                print(" {:<3} ({:3d}) (Charge: {: .3f} e)".format(molecule.atoms[n].symbol, n,
                                                                  molecule.atoms[n].QMcharge))
        f.close()
    # Read through ORCA file, read *last* set of cartesian coordinates
    elif program == "orca":
        f = open(filename, 'r')
        for line in f:
            if line.find("CARTESIAN COORDINATES (ANGSTROEM)") != -1:
                if verbosity >= 2:
                    print("\nCartesian Coordinates found")
                del geom[:]
                readBuffer = f.__next__()
                while True:
                    readBuffer = f.__next__()
                    if readBuffer and readBuffer.strip():
                        geom.append(readBuffer)
                        if verbosity >= 3:
                            readBuffer = readBuffer.split()
                            print(" Found atom: {:<3} {: .8f} {: .8f} {: .8f} in current Cartesian Coordinates".format(
                                readBuffer[0], float(readBuffer[1]), float(readBuffer[2]), float(readBuffer[3])))
                    else:
                        break
        if verbosity >= 1:
            print("\nReading of geometry finished.\nAdding atoms to WellFARe molecule: ", molecule.name)
        for i in geom:
            readBuffer = i.split()
            molecule.addAtom(Atom(readBuffer[0], float(readBuffer[1]), float(readBuffer[2]), float(readBuffer[3]),
                                  0.1))  # 0.1 a placeholder for QM computed carge on the atom
            if verbosity >= 2:
                print(" {:<3} {: .8f} {: .8f} {: .8f}".format(readBuffer[0], float(readBuffer[1]), float(readBuffer[2]),
                                                              float(readBuffer[3])))
        f.close()
        # Read through ORCA file to locate and read Mulliken atomic charges
        f = open(filename, 'r')
        for line in f:
            if line.find("MULLIKEN ATOMIC CHARGES") != -1:
                if verbosity >= 2:
                    print("\nMulliken charges found, reading charges")  # May not strictly need this
                del charges[:]
                readBuffer = f.__next__()  # Only one line to skip, but check this works with no for loop
                while True:
                    readBuffer = f.__next__()
                    if readBuffer.find(
                            "Sum of atomic charges:") == -1:  # Check whether full line needed (as per comment in Gaussian09 section)
                        charges.append(readBuffer)
                        if verbosity >= 3:
                            readBuffer = readBuffer.split()
                            print(" Found atomic charge listing: {:<3} {:<3} {: .8f} in Mulliken charges".format(
                                int(readBuffer[0]), str(readBuffer[1]), float(readBuffer[
                                                                                  3])))  # Assuming that ':' is split into its own list entry. May need to check formatting around whitespace
                    else:
                        break
        if verbosity >= 1:
            print("\nReading of Mulliken charges finished. \nAdding charges to atoms in WellFARe molecule: ",
                  molecule.name)
        for i in charges:
            readBuffer = i.split()
            n = int(readBuffer[0]) - 1
            molecule.atoms[n].setq(readBuffer[
                                       3])  # Again assuming that the charge is the 4th list entry, with ':' having been split on its own.
            if verbosity >= 2:
                print(" {:<3} ({:3d}) (Charge: {: .3f} e)".format(molecule.atoms[n].symbol, n,
                                                                  molecule.atoms[n].QMcharge))
        f.close()

    # BOND ORDER READING SECTION
    bo = []
    bo = np.zeros((molecule.numatoms(), molecule.numatoms()))
    if program == "g09":
        f = open(filename, 'r')
        for line in f:
            if line.find("Atomic Valencies and Mayer Atomic Bond Orders:") != -1:
                if verbosity >= 2:
                    print("\nAtomic Valencies and Mayer Atomic Bond Orders found, reading data")
                bo = np.zeros((molecule.numatoms(), molecule.numatoms()))
                while True:
                    readBuffer = f.__next__()
                    # Check if the whole line is integers only (Header line)
                    if isInt("".join(readBuffer.split())) == True:
                        # And use this information to label the columns
                        columns = readBuffer.split()
                    # If we get to the Löwdin charges, we're done reading
                    elif readBuffer.find("Lowdin Atomic Charges") != -1:
                        break
                    else:
                        row = readBuffer.split()
                        j = 1
                        for i in columns:
                            j = j + 1
                            bo[int(row[0]) - 1][int(i) - 1] = float(row[j])
        f.close()
        if verbosity >= 3:
            print("\nBond Orders:")
            np.set_printoptions(suppress=True)
            np.set_printoptions(formatter={'float': '{: 0.3f}'.format})
            print(bo)
    if program == "orca":
        f = open(filename, 'r')
        for line in f:
            if line.find("Mayer bond orders larger than 0.1") != -1:
                if verbosity >= 2:
                    print("\nMayer bond orders larger than 0.1 found, reading data")
                bo = np.zeros((molecule.numatoms(), molecule.numatoms()))
                while True:
                    readBuffer = f.__next__()
                    # Check if the whole line isn't empty (in that case we're done)
                    if readBuffer and readBuffer.strip():
                        # Break the line into pieces
                        readBuffer = readBuffer[1:].strip()
                        readBuffer = readBuffer.split("B")
                        for i in readBuffer:
                            bondpair1 = int(i[1:4].strip())
                            bondpair2 = int(i[8:11].strip())
                            order = i[-9:].strip()
                            bo[bondpair1][bondpair2] = order
                            bo[bondpair2][bondpair1] = order
                    else:
                        break
        f.close()
        if verbosity >= 3:
            print("\nBond Orders:")
            np.set_printoptions(suppress=True)
            np.set_printoptions(formatter={'float': '{: 0.3f}'.format})
            print(bo)
    # Test if we actually have Mayer Bond orders
    if np.count_nonzero(bo) != 0:
        if verbosity >= 1:
            print("\nAdding bonds to WellFARe molecule: ", molecule.name)
            print("(using bond orders with a cutoff of {: .2f}):".format(bondcutoff))
        for i in range(0, molecule.numatoms()):
            for j in range(i + 1, molecule.numatoms()):
                if bo[i][j] >= bondcutoff:
                    molecule.addBond(i, j)
                    if verbosity >= 2:
                        print(
                            " {:<3} ({:3d}) and {:<3} ({:3d}) (Bond order: {: .3f})".format(molecule.atoms[i].symbol, i,
                                                                                            molecule.atoms[j].symbol, j,
                                                                                            bo[i][j]))
    # Else use 130% of the sum of covalent radii as criterion for a bond (user defined: distfactor)
    else:
        if verbosity >= 1:
            print("\nAdding bonds to WellFARe molecule:", molecule.name)
            print("(using covalent radii scaled by {: .2f}):".format(distfactor))
        for i in range(0, molecule.numatoms()):
            for j in range(i + 1, molecule.numatoms()):
                if molecule.atmatmdist(i, j) <= (
                            SymbolToRadius[molecule.atoms[i].symbol] + SymbolToRadius[
                            molecule.atoms[j].symbol]) * distfactor:
                    molecule.addBond(i, j)
                    if verbosity >= 2:
                        print(
                            " {:<3} ({:3d}) and {:<3} ({:3d}) (Distance: {:.3f} A)".format(molecule.atoms[i].symbol, i,
                                                                                           molecule.atoms[j].symbol, j,
                                                                                           molecule.atmatmdist(i, j)))

    # Now that we know where the bonds are, find angles
    if verbosity >= 2:
        print("\nAdding angles to WellFARe molecule: ", molecule.name)
    for i in range(0, len(molecule.bonds)):
        for j in range(i + 1, len(molecule.bonds)):
            if molecule.bonds[i][0] == molecule.bonds[j][0]:
                molecule.addAngle(molecule.bonds[i][1], molecule.bonds[i][0], molecule.bonds[j][1])
                if verbosity >= 2:
                    print(" {:<3} ({:3d}), {:<3} ({:3d}) and {:<3} ({:3d}) ({:6.2f} deg)".format(
                        molecule.atoms[molecule.bonds[i][1]].symbol, molecule.bonds[i][1],
                        molecule.atoms[molecule.bonds[i][0]].symbol, molecule.bonds[i][0],
                        molecule.atoms[molecule.bonds[j][1]].symbol, molecule.bonds[j][1],
                        math.degrees(molecule.bondangle(len(molecule.angles) - 1))))
            if molecule.bonds[i][0] == molecule.bonds[j][1]:
                molecule.addAngle(molecule.bonds[i][1], molecule.bonds[i][0], molecule.bonds[j][0])
                if verbosity >= 2:
                    print(" {:<3} ({:3d}), {:<3} ({:3d}) and {:<3} ({:3d}) ({:6.2f} deg)".format(
                        molecule.atoms[molecule.bonds[i][1]].symbol, molecule.bonds[i][1],
                        molecule.atoms[molecule.bonds[i][0]].symbol, molecule.bonds[i][0],
                        molecule.atoms[molecule.bonds[j][0]].symbol, molecule.bonds[j][0],
                        math.degrees(molecule.bondangle(len(molecule.angles) - 1))))
            if molecule.bonds[i][1] == molecule.bonds[j][0]:
                molecule.addAngle(molecule.bonds[i][0], molecule.bonds[i][1], molecule.bonds[j][1])
                if verbosity >= 2:
                    print(" {:<3} ({:3d}), {:<3} ({:3d}) and {:<3} ({:3d}) ({:6.2f} deg)".format(
                        molecule.atoms[molecule.bonds[i][0]].symbol, molecule.bonds[i][0],
                        molecule.atoms[molecule.bonds[i][1]].symbol, molecule.bonds[i][1],
                        molecule.atoms[molecule.bonds[j][1]].symbol, molecule.bonds[j][1],
                        math.degrees(molecule.bondangle(len(molecule.angles) - 1))))
            if molecule.bonds[i][1] == molecule.bonds[j][1]:
                molecule.addAngle(molecule.bonds[i][0], molecule.bonds[i][1], molecule.bonds[j][0])
                if verbosity >= 2:
                    print(" {:<3} ({:3d}), {:<3} ({:3d}) and {:<3} ({:3d}) ({:6.2f} deg)".format(
                        molecule.atoms[molecule.bonds[i][0]].symbol, molecule.bonds[i][0],
                        molecule.atoms[molecule.bonds[i][1]].symbol, molecule.bonds[i][1],
                        molecule.atoms[molecule.bonds[j][0]].symbol, molecule.bonds[j][0],
                        math.degrees(molecule.bondangle(len(molecule.angles) - 1))))

    # Same for dihedrals: Use angles to determine where they are
    if verbosity >= 2:
        print("\nAdding dihedrals to WellFARe molecule: ", molecule.name)
    for i in range(0, len(molecule.angles)):
        for j in range(i + 1, len(molecule.angles)):
            if molecule.angles[i][1] == molecule.angles[j][0] and molecule.angles[i][2] == molecule.angles[j][1]:
                molecule.addDihedral(molecule.angles[i][0], molecule.angles[i][1], molecule.angles[i][2],
                                     molecule.angles[j][2])
                if verbosity >= 2:
                    print(" {:<3} ({:3d}), {:<3} ({:3d}), {:<3} ({:3d}) and {:<3} ({:3d}) ({: 7.2f} deg)".format(
                        molecule.atoms[molecule.angles[i][0]].symbol, molecule.angles[i][0],
                        molecule.atoms[molecule.angles[i][1]].symbol, molecule.angles[i][1],
                        molecule.atoms[molecule.angles[i][2]].symbol, molecule.angles[i][2],
                        molecule.atoms[molecule.angles[j][2]].symbol, molecule.angles[j][2],
                        math.degrees(molecule.dihedralangle(len(molecule.dihedrals) - 1))))
            if molecule.angles[i][1] == molecule.angles[j][2] and molecule.angles[i][2] == molecule.angles[j][1]:
                molecule.addDihedral(molecule.angles[i][0], molecule.angles[i][1], molecule.angles[i][2],
                                     molecule.angles[j][0])
                if verbosity >= 2:
                    print(" {:<3} ({:3d}), {:<3} ({:3d}), {:<3} ({:3d}) and {:<3} ({:3d}) ({: 7.2f} deg)".format(
                        molecule.atoms[molecule.angles[i][0]].symbol, molecule.angles[i][0],
                        molecule.atoms[molecule.angles[i][1]].symbol, molecule.angles[i][1],
                        molecule.atoms[molecule.angles[i][2]].symbol, molecule.angles[i][2],
                        molecule.atoms[molecule.angles[j][0]].symbol, molecule.angles[j][0],
                        math.degrees(molecule.dihedralangle(len(molecule.dihedrals) - 1))))
            if molecule.angles[i][1] == molecule.angles[j][0] and molecule.angles[i][0] == molecule.angles[j][1]:
                molecule.addDihedral(molecule.angles[i][2], molecule.angles[j][0], molecule.angles[j][1],
                                     molecule.angles[j][2])
                if verbosity >= 2:
                    print(" {:<3} ({:3d}), {:<3} ({:3d}), {:<3} ({:3d}) and {:<3} ({:3d}) ({: 7.2f} deg)".format(
                        molecule.atoms[molecule.angles[i][2]].symbol, molecule.angles[i][2],
                        molecule.atoms[molecule.angles[j][0]].symbol, molecule.angles[j][0],
                        molecule.atoms[molecule.angles[j][1]].symbol, molecule.angles[j][1],
                        molecule.atoms[molecule.angles[j][2]].symbol, molecule.angles[j][2],
                        math.degrees(molecule.dihedralangle(len(molecule.dihedrals) - 1))))
            if molecule.angles[i][1] == molecule.angles[j][2] and molecule.angles[i][0] == molecule.angles[j][1]:
                molecule.addDihedral(molecule.angles[i][2], molecule.angles[j][2], molecule.angles[j][1],
                                     molecule.angles[j][0])
                if verbosity >= 2:
                    print(" {:<3} ({:3d}), {:<3} ({:3d}), {:<3} ({:3d}) and {:<3} ({:3d}) ({: 7.2f} deg)".format(
                        molecule.atoms[molecule.angles[i][2]].symbol, molecule.angles[i][2],
                        molecule.atoms[molecule.angles[j][2]].symbol, molecule.angles[j][2],
                        molecule.atoms[molecule.angles[j][1]].symbol, molecule.angles[j][1],
                        molecule.atoms[molecule.angles[j][0]].symbol, molecule.angles[j][0],
                        math.degrees(molecule.dihedralangle(len(molecule.dihedrals) - 1))))


# -----------------------------------
# Define functions for linear transit
# -----------------------------------

def alignstructures(molecule_R, molecule_P, verbosity=args.verbosity):
    # First check - make sure number of atoms in each molecule is equal
    if not molecule_R.numatoms() == molecule_P.numatoms():
        print("\nERROR: input structures have different numbers of atoms")
        timestamp("Program aborted at: ")
        sys.exit()
    else:
        if verbosity >= 1:
             print("Preparing to align structures")
    # Simplest check for common atom labels - compare symbols
    for i in range(len(molecule_R.atoms)):
        if molecule_R.atoms[i].symbol != molecule_P.atoms[i].symbol:
             print("ERROR: Atom symbols do not match for index " + str(i))
             timestamp("Program aborted at: ")
             sys.exit()
    # Consider adding a further check to determine whether labelling is consistent (eg are all hydrogen atoms tracked correctly)

    # Search lists for bonds
    commonbonds = []
    combdindices_R = []
    combdindices_P = []
    # NOTE: allowing for differences in ordering is more complex, achieved by keeping separate index lists
    for i in range(len(molecule_R.bonds)):
        for j in range(len(molecule_P.bonds)):
            if molecule_R.bonds[i][0] == molecule_P.bonds[j][0] and molecule_R.bonds[i][1] == molecule_P.bonds[j][1]:
                commonbonds.append(molecule_R.bonds[i])
                combdindices_R.append(i)
                combdindices_P.append(j)
            elif molecule_R.bonds[i][0] == molecule_P.bonds[j][1] and molecule_R.bonds[i][1] == molecule_P.bonds[j][0]:
                commonbonds.append(molecule_R.bonds[i])
                combdindices_R.append(i)
                combdindices_P.append(j)
    if verbosity >= 2:
        print("\nList of bonds in common:")
        print(commonbonds)
    # Debug only, print the two atom lists side by side for comparison
    print("Common bond indices, R and P respectively:")
    print(combdindices_R)
    print(combdindices_P)
    print("\nCommon bonds as found in each molecule:")
    print("Index    R bond    P bond    Common bond")
    for i in range(len(commonbonds)):
        print("{:^8} {:^9} {:^9} {:^11}".format(str(i), str(molecule_R.bonds[combdindices_R[i]]), str(molecule_P.bonds[combdindices_P[i]]), str(commonbonds[i])))

    # Carry out an equivalent search for bond angles
    commonangles = []
    comangindices_R = []
    comangindices_P = []
    # NOTE: allowing for differences in ordering is more complex, achieved by keeping separate index lists
    for i in range(len(molecule_R.angles)):
        for j in range(len(molecule_P.angles)):
            if molecule_R.angles[i][0] == molecule_P.angles[j][0] and molecule_R.angles[i][1] == molecule_P.angles[j][1] and molecule_R.angles[i][2] == molecule_P.angles[j][2]:
                commonangles.append(molecule_R.angles[i])
                comangindices_R.append(i)
                comangindices_P.append(j)
            elif molecule_R.angles[i][0] == molecule_P.angles[j][2] and molecule_R.angles[i][1] == molecule_P.angles[j][1] and molecule_R.angles[i][2] == molecule_P.angles[j][0]:
                commonangles.append(molecule_R.angles[i])
                comangindices_R.append(i)
                comangindices_P.append(j)
    if verbosity >= 2:
        print("\nList of angles in common:")
        print(commonangles)
    # Debug only, print the two atom lists side by side for comparison
    print("Common angle indices, R and P respectively:")
    print(comangindices_R)
    print(comangindices_P)
    print("\nCommon angles as found in each molecule:")
    print("Index    R angle      P angle      Common angle")
    for i in range(len(commonbonds)):
        print("{:^8} {:^11} {:^11} {:^13}".format(str(i), str(molecule_R.angles[comangindices_R[i]]), str(molecule_P.angles[comangindices_P[i]]), str(commonangles[i])))

    # Translate to put central atom of angle at the origin (simple vector addition)
    originatom_R = molecule_R.atoms[molecule_R.angles[comangindices_R[0]][1]]
    originatom_P = molecule_P.atoms[molecule_P.angles[comangindices_P[0]][1]]
    # Debug only, print coordinates of the atom used for translation
    print("Identifying atom to translate to the origin, coordinates for R and P:")
    print(originatom_R.coord)
    print(originatom_P.coord)
    shiftx_R = originatom_R.coord[0] * -1
    shifty_R = originatom_R.coord[1] * -1
    shiftz_R = originatom_R.coord[2] * -1
    shiftx_P = originatom_P.coord[0] * -1
    shifty_P = originatom_P.coord[1] * -1
    shiftz_P = originatom_P.coord[2] * -1
    molecule_R.translateMolecule(shiftx_R, shifty_R, shiftz_R)
    molecule_P.translateMolecule(shiftx_P, shifty_P, shiftz_P)
    # Debug only, print coordinates of the focus atom in each molecule after translation
    print("Translated coordinates of chosen atom, in R and P molecules:")
    print(molecule_R.atoms[molecule_R.angles[comangindices_R[0]][1]].coord)
    print(molecule_P.atoms[molecule_P.angles[comangindices_P[0]][1]].coord)
    if verbosity >= 1:
        print("\nCoordinates after translating to common origin:") 
        print(molecule_R.xyzString())
        print(molecule_P.xyzString())

    # Rotate one end atom into the positive x axis (slightly more complicated)
    # Set up vectors lying along the positive x axis and the bond of interest
    if molecule_R.angles[comangindices_R[0]][0] == molecule_P.angles[comangindices_P[0]][0]:
        end1_R = molecule_R.atoms[molecule_R.angles[comangindices_R[0]][0]].coord
        end1_P = molecule_P.atoms[molecule_P.angles[comangindices_P[0]][0]].coord
    elif molecule_R.angles[comangindices_R[0]][0] == molecule_P.angles[comangindices_P[0]][2]:
        end1_R = molecule_R.atoms[molecule_R.angles[comangindices_R[0]][0]].coord
        end1_P = molecule_P.atoms[molecule_P.angles[comangindices_P[0]][2]].coord
    xpos = [1.0, 0.0, 0.0]
    # Calculate the relevant norms and dot products to use in finding angles
    norm_e1R = np.linalg.norm(end1_R)
    norm_e1P = np.linalg.norm(end1_P)
    norm_xp = np.linalg.norm(xpos)
    ctheta_1R = np.dot(end1_R, xpos)
    ctheta_1P = np.dot(end1_P, xpos) 
    # Calculate the angle between each of the bond-vectors and the positive x axis
    theta_1R = math.acos(ctheta_1R / (norm_e1R * norm_xp))
    theta_1P = math.acos(ctheta_1P / (norm_e1P * norm_xp))
    # Debug only - print the calculated angles
    print("Calculated angles to x axis, in degrees, for reactant and product")
    print(math.degrees(theta_1R))
    print(math.degrees(theta_1P))
    # Calculate cross products to give the axes for rotation
    cross_R = np.cross(xpos, end1_R)
    cross_P = np.cross(xpos, end1_P)
    # Carry out the rotation using a built in function
    #molecule_R.rotateMoleculeArbAxis(xpos, end1_R, angle=math.degrees(theta_1R))
    #molecule_P.rotateMoleculeArbAxis(xpos, end1_P, angle=math.degrees(theta_1P))
    molecule_R.rotateMoleculeGivenAxis(-cross_R, angle=math.degrees(theta_1R))
    molecule_P.rotateMoleculeGivenAxis(-cross_P, angle=math.degrees(theta_1P))
    if verbosity >= 1:
        print("\nCoordinates after rotating common bond into x axis")
        print(molecule_R.xyzString())
        print(molecule_P.xyzString())

    # Rotate about the x axis to bring the other end atom into the xy plane 
    if molecule_R.angles[comangindices_R[0]][2] == molecule_P.angles[comangindices_P[0]][2]:
        end2_R = molecule_R.atoms[molecule_R.angles[comangindices_R[0]][2]].coord
        end2_P = molecule_P.atoms[molecule_P.angles[comangindices_P[0]][2]].coord
    elif molecule_R.angles[comangindices_R[0]][2] == molecule_P.angles[comangindices_P[0]][0]:
        end2_R = molecule_R.atoms[molecule_R.angles[comangindices_R[0]][2]].coord
        end2_P = molecule_P.atoms[molecule_P.angles[comangindices_P[0]][0]].coord
    # Calculate the angle to rotate based on the projection into the yz plane
    yco_2R = end2_R[1]
    zco_2R = end2_R[2]
    yco_2P = end2_P[1]
    zco_2P = end2_P[2]
    theta_2R = math.atan(zco_2R / yco_2R)
    theta_2P = math.atan(zco_2P / yco_2P)
    # Debug only - print the calculated angles
    print("Calculated angles to rotate into xy plane, in degrees, for reactant and product")
    print(math.degrees(theta_2R))
    print(math.degrees(theta_2P))
    # Carry out the rotation using a built in function
    molecule_R.rotateMolecule(axis="x", angle=math.degrees(-theta_2R))
    molecule_P.rotateMolecule(axis="x", angle=math.degrees(-theta_2P))
    if verbosity >= 1:
        print("\nCoordinates after rotating second common bond into xy plane")
        print(molecule_R.xyzString())
        print(molecule_P.xyzString())

def interpolatepoints(molecule_R, molecule_P, n, verbosity=args.verbosity):
    if not molecule_R.numatoms() == molecule_P.numatoms():
        print("\nERROR: input structures have different numbers of atoms")
        print("Reactant atom count: " + str(molecule_R.numatoms))
        print("Product atom count: " + str(molecule_P.numatoms))
        timestamp("Program aborted at: ")
        sys.exit()
    else:
        if verbosity >= 1:
             print("\nBeginning interpolation")
             print("Using " + str(n) + " linear transit points")
    reactant_distMat = molecule_R.distMatrix()
    product_distMat = molecule_P.distMatrix()
    if verbosity >= 1:
        print("\nDistance matrix of reactant:")
        print(reactant_distMat)
        print("\nDistance matrix of product:")
        print(product_distMat)
    if verbosity >= 3:
        # As a temporary check, subtract reactant distance matrix from product matrix, sum up the absolute value of entries
        diffrpdistmat = np.subtract(product_distMat, reactant_distMat)
        sumabsrpdiffs = 0.0
        for i in range(molecule_R.numatoms()):
            for j in range(molecule_R.numatoms()):
                sumabsrpdiffs += abs(diffrpdistmat[i][j])
        print("\nSum of absolute differences in product and reactant distance matrices: " + str(sumabsrpdiffs))
    step = 1/(n - 1)
    mix = 0.0
    pointnum = 0
    while pointnum < n:
        mixed_distMat = np.zeros((molecule_R.numatoms(), molecule_R.numatoms()))
        for i in range(molecule_R.numatoms()):
            for j in range(i, molecule_R.numatoms()):
                mixed_distMat[i][j] = (((1 - mix) * reactant_distMat[i][j]) + (mix * product_distMat[i][j]))
        if verbosity >= 3:
            print("\nMixed distance matrix at point: " + str(pointnum))
            print(mixed_distMat)
        # As a temporary check, subtract reactant distance matrix from mixed matrix, sum up the absolute value of entries
        diffdistmat = np.subtract(mixed_distMat, reactant_distMat)
        sumabsdiffs = 0.0
        for i in range(molecule_R.numatoms()):
            for j in range(molecule_R.numatoms()):
                sumabsdiffs += abs(diffdistmat[i][j])
        if verbosity >= 2:
            print("\nSum of absolute differences in distance matrix entries at point {0}: reactant vs mix: {1}".format(pointnum, sumabsdiffs))
        molecule_R.setGeomDM(mixed_distMat)
        print("\nCoordinates at linear transit point " + str(pointnum) + " ")
        #print(mix) # Debug only 
        print(molecule_R.gaussString())
        pointnum += 1
        mix += step
        #print(mix) # Debug only
        print("Done")
    
# --------------------------------------
# Main working part of the program below
# --------------------------------------

print("###################################")
print("            lintransit.py          ")
timestamp("Program start time")
print("###################################")

reactant_mol = Molecule("Reactant", 0)
extractCoordinates(args.reactant, reactant_mol, verbosity=args.verbosity, bondcutoff=args.bondcutoff)

product_mol = Molecule("Product", 0)
extractCoordinates(args.product, product_mol, verbosity=args.verbosity, bondcutoff=args.bondcutoff)

print("Geometry of reactant and product molecules at start of program")
print(reactant_mol.xyzString())
print(product_mol.xyzString())

print("\nDistance matrix for reactant")
reactant_mol.prettydistMatrix()

print("\nDistance matrix for product")
product_mol.prettydistMatrix()

print("Testing geometry from distance matrix - feeding reactant distance matrix back to reactant")
rDistM = reactant_mol.distMatrix()
reactant_mol.setGeomDM(rDistM)

print("\nTesting geometry from distance matrix - feeding product distance matrix back to product")
pDistM = product_mol.distMatrix()
product_mol.setGeomDM(pDistM)

print("\nGeometry of reactant after setting from reactant distance matrix:")
print(reactant_mol.xyzString())

print("\nGeometry of product after setting from product distance matrix:")
print(product_mol.xyzString())

#alignstructures(reactant_mol, product_mol)

interpolatepoints(reactant_mol, product_mol, int(args.numpoints))

print("Program ends")
print("-----------------------------------------")
