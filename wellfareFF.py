import sys
import getopt
import math
import time
from importlib.util import find_spec


def timestamp(s):
    print(s + time.strftime("%Y/%m/%d %X"))


# ASCII FONTS from: http://patorjk.com/software/taag/
# Font = "Big"
def ProgramHeader():
    print("###################################################################")
    print("Wellington Fast Assessment of Reactions using Force Fields")
    print(" __          __  _ _ ______      _____      ______ ______ ")
    print(" \ \        / / | | |  ____/\   |  __ \    |  ____|  ____|")
    print("  \ \  /\  / /__| | | |__ /  \  | |__) |___| |__  | |__   ")
    print("   \ \/  \/ / _ \ | |  __/ /\ \ |  _  // _ \  __| |  __|  ")
    print("    \  /\  /  __/ | | | / ____ \| | \ \  __/ |    | |     ")
    print("     \/  \/ \___|_|_|_|/_/    \_\_|  \_\___|_|    |_|     ")
    print("                                              Version 0.01")
    print("      WellFAReFF Copyright (C) 2015 Matthias Lein         ")
    print("    This program comes with ABSOLUTELY NO WARRANTY        ")
    print("     This is free software, and you are welcome to        ")
    print("       redistribute it under certain conditions.          ")
    timestamp('Program started at: ')
    print("###################################################################\n")


def ProgramFooter():
    print("\n###################################################################")
    print("  _____                                                      _     ")
    print(" |  __ \                                                    | |    ")
    print(" | |__) | __ ___   __ _ _ __ __ _ _ __ ___     ___ _ __   __| |___ ")
    print(" |  ___/ '__/ _ \ / _` | '__/ _` | '_ ` _ \   / _ \ '_ \ / _` / __|")
    print(" | |   | | | (_) | (_| | | | (_| | | | | | | |  __/ | | | (_| \__ \ ")
    print(" |_|   |_|  \___/ \__, |_|  \__,_|_| |_| |_|  \___|_| |_|\__,_|___/")
    print("                   __/ |                                           ")
    print("                  |___/                                            ")
    timestamp('Program terminated at: ')
    print("###################################################################")


def ProgramAbort():
    print("\n###################################################################")
    print("  _____                     _                _           _ ")
    print(" |  __ \                   | |              | |         | |")
    print(" | |__) |   _ _ __     __ _| |__   ___  _ __| |_ ___  __| |")
    print(" |  _  / | | | '_ \   / _` | '_ \ / _ \| '__| __/ _ \/ _` |")
    print(" | | \ \ |_| | | | | | (_| | |_) | (_) | |  | ||  __/ (_| |")
    print(" |_|  \_\__,_|_| |_|  \__,_|_.__/ \___/|_|   \__\___|\__,_|")
    timestamp('Program aborted at: ')
    print("###################################################################")
    sys.exit()
    return


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


def ProgramError():
    print("\n###################################################################")
    print("  ______                     ")
    print(" |  ____|                    ")
    print(" | |__   _ __ _ __ ___  _ __ ")
    print(" |  __| | '__| '__/ _ \| '__|")
    print(" | |____| |  | | | (_) | |   ")
    print(" |______|_|  |_|  \___/|_|   ")
    timestamp('Error time/date: ')
    print("###################################################################")
    return


# Check for numpy, exit immediately if not available
module_loader = find_spec('numpy')
found = module_loader is not None
if not found:
    ProgramError("Module numpy is required")
    ProgramAbort()
import numpy

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
    ProgramError("Module argparse is required")
    ProgramAbort()
import argparse


#############################################################################################################
# This section is for the definition of *all* constants and conversion factors
#############################################################################################################

# Conversion of mass in atomic mass units (AMU) to
# atomic units (electron masses)
def AMU2au(amu):
    return amu * 1822.88839


# Same in reverse
def au2AMU(au):
    return au / 1822.88839


# Conversion of length in Angstroms to  to
# atomic units (Bohrs)
def Ang2Bohr(ang):
    return ang * 1.889725989


# Same in reverse
def Bohr2Ang(bohr):
    return bohr / 1.889725989


# Conversion of energy in Joules to atomic units (Hartrees)
def J2au(J):
    return J / (4.35974394 * (10 ** -18))


# Same in reverse
def au2J(au):
    return au * (4.35974393 * (10 ** -18))


# Conversion of energy in kcal/mol to atomic units (Hartrees)
def kcal_mol2au(kcm):
    return kcm / 627.503


# Same in reverse
def au2kcal_mol(au):
    return au * 627.503


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

# Define dictionary to convert atomic symbols to van der Waals radii (in Angstrom)
SymbolToVdWRadius = {
    "H": 1.10, "He": 1.40, "Li": 1.82, "Be": 1.53, "B": 1.92, "C": 1.70,
    "N": 1.55, "O": 1.52, "F": 1.47, "Ne": 1.54, "Na": 2.27, "Mg": 1.73,
    "Al": 1.84, "Si": 2.10, "P": 1.80, "S": 1.80, "Cl": 1.75, "Ar": 1.88,
    "K": 2.75, "Ca": 2.31, "Sc": 2.15, "Ti": 2.11, "V": 2.07, "Cr": 2.06,
    "Mn": 2.05, "Fe": 2.04, "Co": 2.00, "Ni": 1.97, "Cu": 1.96, "Zn": 2.01,
    "Ga": 1.87, "Ge": 2.11, "As": 1.85, "Se": 1.90, "Br": 1.85, "Kr": 2.02,
    "Rb": 3.03, "Sr": 2.49, "Y": 2.32, "Zr": 2.23, "Nb": 2.18, "Mo": 2.17,
    "Tc": 2.16, "Ru": 2.13, "Rh": 2.10, "Pd": 2.10, "Ag": 2.11, "Cd": 2.18,
    "In": 1.93, "Sn": 2.17, "Sb": 2.06, "Te": 2.06, "I": 1.98, "Xe": 2.16,
    "Cs": 3.43, "Ba": 2.68, "La": 2.43, "Ce": 2.42, "Pr": 2.40, "Nd": 2.39,
    "Pm": 2.38, "Sm": 2.36, "Eu": 2.35, "Gd": 2.34, "Tb": 2.33, "Dy": 2.31,
    "Ho": 2.30, "Er": 2.29, "Tm": 2.27, "Yb": 2.26, "Lu": 2.24, "Hf": 2.23,
    "Ta": 2.22, "W": 2.18, "Re": 2.16, "Os": 2.16, "Ir": 2.13, "Pt": 2.13,
    "Au": 2.14, "Hg": 2.23, "Tl": 1.96, "Pb": 2.02, "Bi": 2.07, "Po": 1.97,
    "At": 2.02, "Rn": 2.20, "Fr": 3.48, "Ra": 2.83, "Ac": 2.47, "Th": 2.45,
    "Pa": 2.43, "U": 2.41, "Np": 2.39, "Pu": 2.43, "Am": 2.44, "Cm": 2.45,
    "Bk": 2.44, "Cf": 2.45, "Es": 2.45, "Fm": 2.45, "Md": 2.46, "No": 2.46,
    "Lr": 2.46, "Rf": "?", "Db": "?", "Sg": "?", "Bh": "?", "Hs": "?",
    "Mt": "?", "Ds": "?", "Rg": "?", "Cn": "?", "Uut": "?", "Uuq": "?",
    "Uup": "?", "Uuh": "?", "Uus": "?", "Uuo": "?"}

# Define dictionary to convert atomic symbols to (Pauling) electronegativity
SymbolToEN = {
    "H": 2.20, "He": 0.00, "Li": 0.98, "Be": 1.57, "B": 2.04, "C": 2.55,
    "N": 3.04, "O": 3.44, "F": 3.98, "Ne": 0.00, "Na": 0.93, "Mg": 1.31,
    "Al": 1.61, "Si": 1.90, "P": 2.19, "S": 2.58, "Cl": 3.16, "Ar": 0.00,
    "K": 0.82, "Ca": 1.00, "Sc": 1.36, "Ti": 1.54, "V": 1.63, "Cr": 1.66,
    "Mn": 1.55, "Fe": 1.83, "Co": 1.88, "Ni": 1.91, "Cu": 1.90, "Zn": 1.65,
    "Ga": 1.81, "Ge": 2.01, "As": 2.18, "Se": 2.55, "Br": 2.96, "Kr": 3.00,
    "Rb": 0.82, "Sr": 0.95, "Y": 1.22, "Zr": 1.33, "Nb": 1.60, "Mo": 2.16,
    "Tc": 1.90, "Ru": 2.00, "Rh": 2.28, "Pd": 2.20, "Ag": 1.93, "Cd": 1.69,
    "In": 1.78, "Sn": 1.96, "Sb": 2.05, "Te": 2.10, "I": 2.66, "Xe": 2.60,
    "Cs": 0.79, "Ba": 0.89, "La": 1.10, "Ce": 1.12, "Pr": 1.13, "Nd": 1.14,
    "Pm": 1.13, "Sm": 1.17, "Eu": 1.20, "Gd": 1.20, "Tb": 1.10, "Dy": 1.22,
    "Ho": 1.23, "Er": 1.24, "Tm": 1.25, "Yb": 1.10, "Lu": 1.27, "Hf": 1.30,
    "Ta": 1.50, "W": 2.36, "Re": 1.90, "Os": 2.20, "Ir": 2.20, "Pt": 2.28,
    "Au": 2.54, "Hg": 2.00, "Tl": 1.62, "Pb": 1.87, "Bi": 2.02, "Po": 2.00,
    "At": 2.20, "Rn": 2.20, "Fr": 0.70, "Ra": 0.90, "Ac": 1.10, "Th": 1.30,
    "Pa": 1.50, "U": 1.38, "Np": 1.36, "Pu": 1.28, "Am": 1.13, "Cm": 1.28,
    "Bk": 1.30, "Cf": 1.30, "Es": 1.30, "Fm": 1.30, "Md": 1.30, "No": 1.30,
    "Lr": 1.30, "Rf": 1.30, "Db": 1.30, "Sg": 1.30, "Bh": 1.30, "Hs": 1.30,
    "Mt": 1.30, "Ds": 1.30, "Rg": 1.30, "Cn": 1.30, "Uut": 1.30, "Uuq": 1.30,
    "Uup": 1.30, "Uuh": 1.30, "Uus": 1.30, "Uuo": 1.30}

# Define dictionary to convert atomic symbods to valence electron number
SymbolToValenceE = {
    "H": 1, "He": 2, "Li": 1, "Be": 2, "B": 3, "C": 4,
    "N": 5, "O": 6, "F": 7, "Ne": 8, "Na": 1, "Mg": 2,
    "Al": 3, "Si": 4, "P": 5, "S": 6, "Cl": 7, "Ar": 8,
    "K": 1, "Ca": 2, "Sc": 3, "Ti": 4, "V": 5, "Cr": 6,
    "Mn": 7, "Fe": 8, "Co": 9, "Ni": 10, "Cu": 11, "Zn": 12}
# Note that this will need to be completed later

# ---------------------------------------------------
# Define global empirical parameters for force field
# ---------------------------------------------------

# Define a dictionary for the element specific parameter k_a
k_a = {
    "H": 1.755, "He": 1.755, "B": 2.287, "C": 2.463, "N": 2.559, "O": 2.579,
    "F": 2.465, "Ne": 2.465, "Al": 2.508, "Si": 2.684, "P": 2.780, "S": 2.800,
    "Cl": 2.686, "Li": 2.20, "Na": 2.20, "K": 2.20, "Rb": 2.20, "Cs": 2.20,
    "Fr": 2.20, "Be": 2.80, "Mg": 2.80, "Ca": 2.80, "Sr": 2.80, "Ba": 2.80,
    "Ra": 2.80, "Ar": 2.75, "Sc": 2.95, "Ti": 2.95, "V": 2.95, "Cr": 2.95,
    "Mn": 2.95, "Fe": 2.95, "Co": 2.95, "Ni": 2.95, "Cu": 2.95, "Zn": 2.95,
    "Ga": 2.95, "Ge": 2.95, "As": 2.95, "Se": 2.95, "Br": 2.95, "Kr": 2.95,
    "Y": 3.15, "Zr": 3.15, "Nb": 3.15, "Mo": 3.15, "Tc": 3.15, "Ru": 3.15,
    "Rh": 3.15, "Pd": 3.15, "Ag": 3.15, "Cd": 3.15, "In": 3.15, "Sn": 3.15,
    "Sb": 3.15, "Te": 3.15, "I": 3.15, "Xe": 3.15, "La": 3.80, "Ce": 3.80,
    "Pr": 3.80, "Nd": 3.80, "Pm": 3.80, "Sm": 3.80, "Eu": 3.80, "Gd": 3.80,
    "Tb": 3.80, "Dy": 3.80, "Ho": 3.80, "Er": 3.80, "Tm": 3.80, "Yb": 3.80,
    "Lu": 3.80, "Hf": 3.80, "Ta": 3.80, "W": 3.80, "Re": 3.80, "Os": 3.80,
    "Ir": 3.80, "Pt": 3.80, "Au": 3.80, "Hg": 3.80, "Tl": 3.80, "Pb": 3.80,
    "Bi": 3.80, "Po": 3.80, "At": 3.80, "Rn": 3.80}

# Define dictionary for the element specific parameter k_z
k_z = {
    "H": 3.00, "He": 2.35, "Li": 1.70, "Be": 5.50, "B": 0.95, "C": 0.95,
    "N": 0.95, "O": 0.95, "F": 0.95, "Ne": 0.95, "Na": 2.50, "Mg": 3.00,
    "Al": 0.75, "Si": 0.75, "P": 0.75, "S": 0.75, "Cl": 0.75, "Ar": 0.75,
    "K": 3.00, "Ca": 3.00, "Sc": 0.65, "Ti": 0.65, "V": 0.65, "Cr": 0.65,
    "Mn": 0.65, "Fe": 0.65, "Co": 0.65, "Ni": 0.65, "Cu": 0.65, "Zn": 0.65,
    "Ga": 0.65, "Ge": 0.65, "As": 0.65, "Se": 0.65, "Br": 0.65, "Kr": 0.65,
    "Rb": 3.00, "Sr": 3.00, "Cs": 0.60, "Ba": 0.60, "La": 0.60, "Ce": 0.60,
    "Pr": 0.60, "Nd": 0.60, "Pm": 0.60, "Sm": 0.60, "Eu": 0.60, "Gd": 0.60,
    "Tb": 0.60, "Dy": 0.60, "Ho": 0.60, "Er": 0.60, "Tm": 0.60, "Yb": 0.60,
    "Lu": 0.60, "Hf": 0.60, "Ta": 0.60, "W": 0.60, "Re": 0.60, "Os": 0.60,
    "Ir": 0.60, "Pt": 0.60, "Au": 0.60, "Hg": 0.60, "Tl": 0.60, "Pb": 0.60,
    "Bi": 0.60, "Po": 0.60, "At": 0.60, "Rn": 0.60, "Fr": 0.60, "Ra": 0.60,
    "Ac": 0.60, "Th": 0.60, "Pa": 0.60, "U": 0.60, "Np": 0.60, "Pu": 0.60,
    "Am": 0.60, "Cm": 0.60, "Bk": 0.60, "Cf": 0.60, "Es": 0.60, "Fm": 0.60,
    "Md": 0.60, "No": 0.60, "Lr": 0.60, "Rf": 0.60, "Db": 0.60, "Sg": 0.60,
    "Bh": 0.60, "Hs": 0.60, "Mt": 0.60, "Ds": 0.60, "Rg": 0.60, "Cn": 0.60,
    "Uut": 0.60, "Uuq": 0.60, "Uup": 0.60, "Uuh": 0.60, "Uus": 0.60, "Uuo": 0.60}
# Note no values specified for row 5 elements outside s block - Grimme gives rows 1, 2, 3 and 4, then Z>54

# Define individual global parameters
k_EN = -0.164
k_a2 = 0.221
k_a13 = 2.81
k_b13 = 0.53
k_13r = 0.7
k_damping = 0.11
k_overlap = 0.5
a1 = 0.45
a2 = 4.0
s8 = 2.7
E_ES_14 = 0.85
E_disp_rep = 0.5
beta_rep = 16.5
k_q = 1.15

# Define dictionary for the element specific hydrogen bonding parameter
k_hbnd = {"N": 0.8,
          "O": 0.3,
          "F": 0.1,
          "P": 2.0,
          "S": 2.0,
          "Cl": 2.0,
          "Se": 2.0,
          "Br": 2.0}

# Define dictionary for the element specific parameter k_X
k_X = {"Cl": 0.3,
       "Br": 0.6,
       "I": 0.8,
       "At": 1.0}

# Define dictionaries for the hydrogen and halogen bonding parameters k_q1 and k_q2
k_q1 = {"hbond": 10,
        "xbond": -6.5}
k_q2 = {"hbond": 5,
        "xbond": 1}

# Define a dictionary for the DFT-D3 C6 coefficients
# (Values taken from www.thch.uni-bonn.de/tc/downloads/DFT-D3/data/refmol.txt)
C6 = {
    "H": 7.5916, "He": 1.5583, "Li": 1163.4454, "Be": 257.4863, "B": 107.1777,
    "C": 49.1130, "N": 25.2685, "O": 15.5059, "F": 9.6916, "Ne": 6.2896,
    "Na": 1608.0286, "Mg": 683.3758, "Al": 540.5406, "Si": 317.8574, "P": 191.6887,
    "S": 134.0066, "Cl": 92.3460, "Ar": 64.6462, "K": 4983.5009, "Ca": 2352.6862,
    "Sc": 1702.6213, "Ti": 1361.9185, "V": 1116.0984, "Cr": 690.7425, "Mn": 802.7484,
    "Fe": 109.5041, "Co": 532.7794, "Ni": 574.7436, "Cu": 337.1808, "Zn": 340.5213,
    "Ga": 483.7516, "Ge": 363.5474, "As": 262.9498, "Se": 213.6738, "Br": 167.1297,
    "Kr": 130.4017, "Rb": 6138.7755, "Sr": 3381.3672, "Y": 2365.8925, "Zr": 1822.7181,
    "Nb": 1475.2500, "Mo": 845.8972, "Tc": 1067.0169, "Ru": 239.0071, "Pd": 608.5041,
    "Ag": 426.7450, "Cd": 468.1900, "In": 757.7397, "Sn": 627.5677, "Sb": 492.9379,
    "Te": 425.5355, "I": 351.9667, "Xe": 290.2223, "Cs": 9330.7294, "Ba": 5726.9887,
    "La": 3990.6172, "Ce": 688.0353, "Pr": 4342.2386, "Nd": 3924.4211, "Pm": 3710.9375,
    "Sm": 3522.0508, "Eu": 3358.3122, "Gd": 1891.6719, "Tb": 2851.6677, "Dy": 2617.3310,
    "Ho": 2664.1668, "Er": 2545.1713, "Tm": 2437.4539, "Yb": 2390.1227, "Lu": 1597.4796,
    "Hf": 1441.2394, "Ta": 1163.8241, "W": 814.3622, "Re": 836.3310, "Os": 297.8338,
    "Ir": 566.0660, "Pt": 391.1448, "Au": 342.3526, "Hg": 362.0755, "Tl": 792.2378,
    "Pb": 738.8156, "Bi": 617.5296, "Po": 562.6011, "At": 483.6536, "Rn": 412.8275,
    "Fr": 7314.7398, "Ra": 5305.4399, "Ac": 3799.6565, "Th": 2847.2704, "Pa": 2908.9206,
    "U": 2721.5209, "Np": 3032.9760, "Pu": 2815.2366}
# Note: Several elements have two different C6 values listed, for different numbers of unpaired electrons. At present, the value with fewest unpaired electrons is used
# This affects Fe (4 unpaired, 491.3349; 0 unpaired 109.5041), Ru (4 unpaired, 598.1988; 0 unpaired, 239.0071), Os (4 unpaired, 678.5278; 0 unpaired, 297.8338)
# Note also that values are for the element alone, except Ce with data available only for CeH3

# Define a dictionary for the optimised values of the parameter a1 used in C6-only dispersion for different density functional approximations
# Values in atomic units from DOI: 10.1021/acs.jctc.5b00400
CSO_a1 = {
    "BLYP": 1.28, "BP86": 1.01, "PBE": 0.24, "TPSS": 0.72,
    "B3LYP": 0.86, "PBE0": 0.20, "PW6B95": -0.15, "B2PLYP": 0.24}


#############################################################################################################
# Do *not* define constants or conversion factors below here
#############################################################################################################

# Test if the argument is (can be converted to)
# an integer number
def isInt(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


#############################################################################################################
# Potential Fuctions are to be defined below here
#############################################################################################################

def potMorse(r, r0, D, b):
    """
    Morse oscillator potential
    """

    u = D * (1 - math.exp(-b * (r - r0))) ** 2

    return u


def potHarmonic(a, a0, k):
    """"
    Harmonic potential (for stretches and bends)
    """

    u = 0.5 * k * (a - a0) ** 2

    return u


def potSimpleCosine(theta, theta0, k):
    """"
    Extremely simplified cosine potential for torsions
    """

    u = k * (1 + numpy.cos(math.radians(180) + theta - theta0))

    return u


def bond_exp(a, b):
    """
    Exponent a used in the Generalised Lennard-Jones potential for bond stretches
    """
    # Now taking as input the symbols a and b for the two atoms in a bond 
    delta_EN = SymbolToEN[a] - SymbolToEN[b]
    bond_exp = (k_a[a] * k_a[b]) + (k_EN * (delta_EN ** 2))

    return bond_exp


def exp_1_3(a, b):
    """
    Exponent a used in the Generalised Lennard-Jones potential for 1,3-stretches
    """
    # Now taking atomic symbols as input 
    exp_1_3 = k_a13 + (k_b13 * k_a[a] * k_a[b])
    # Special case: if both atoms are in a ring, then k_13r is added to exp_1_3. Detection of rings to be implemented later

    return exp_1_3


def potGLJ(r, r0, k_str, a):
    """"
    Generalised Lennard-Jones potential (for bond and 1,3 stretches)
    """
    # k_str to be set equal to force constant for /that/ bond as read from Hessian as an initial guess, fitting implemented later
    u = k_str * (1 + ((r0 / r) ** a) - 2 * ((r0 / r) ** (a / 2)))

    return u


def DampingFunction(a, b, r):
    """
    Distance dependent damping function for atoms with symbols a and b, and separation r
    """
    r_cov = SymbolToRadius[a] + SymbolToRadius[b]
    # Covalent distance defined as sum of covalent radii for the two atoms - will need to check this is correct

    f_dmp = 1 / (1 + k_damping * ((r / r_cov) ** 4))

    return f_dmp


def potBendNearLinear(a, a0, k_bnd, f_dmp):
    """
    Bending potential for equilibrium angles close to linearity
    """

    u = k_bnd * f_dmp * ((a0 - a) ** 2)

    return u


def potAnyBend(a, a0, k_bnd, f_dmp):
    """
    Double minimum bending potential
    """

    u = k_bnd * f_dmp * ((math.cos(a0) - math.cos(a)) ** 2)

    return u


def ChiralityFunction(theta):
    """
    Function used in torsion potential to give correct mirror symmetry
    """

    f_chiral = 0.5 * (1 - math.erf(theta - math.pi))

    return f_chiral


def potTorsion(theta, theta0, f_dmp, k_tors):
    """
    Torsion potential
    """

    f_chiral = ChiralityFunction(theta)

    u = 0.0
    # Sort out where n comes from in the following sum, check if k_tors ** n or k_tors_n
    # for n in # Range to be determined:
    #    inner_sum = (f_chiral * (1 + math.cos(n * (theta - theta0) + math.pi))) + ((1 - f_chiral) * (1 + math.cos(n * (theta + theta0 - (2 * math.pi)) + math.pi)))
    #    u = u + ((k_tors ** n) * inner_sum)
    # u = u * f_dmp

    return u


def AngleDamping(theta):
    """
    Angle dependent damping function for hydrogen bonding interactions
    """

    f_dmp_a = (0.5 * (math.cos(theta) + 1)) ** 6

    return f_dmp_a


def HBondDamping(a, b, r_ab):
    """
    Distance dependent damping function for hydrogen bonding with donor/acceptor atoms having symbols a and b, separation r_ab
    """
    r_cov = SymbolToRadius[a] + SymbolToRadius[b]
    # Covalent distance defined as sum of covalent radii for the two atoms - will need to check this is correct

    f_dmp_hbnd = 1 / (1 + k_damping * ((r_ab / r_cov) ** 12))

    return f_dmp_hbnd


def AtomicHBondFactor(sym, chg):
    """
    Atomic factor c_hbnd for an atom with symbol sym and atomic charge chg, used in hydrogen bond interaction strength
    """

    c_hbnd = k_hbnd[sym] * (math.exp(-1 * k_q1["hbond"] * chg) / (math.exp(-1 * k_q1["hbond"] * chg) + k_q2["hbond"]))

    return c_hbnd


def HBondStrengthFactor(sym_a, chg_a, r_ah, sym_b, chg_b, r_bh):
    """
    Modifier for the strength of a hydrogen bonding interaction where donor/acceptor atoms a and b are at distances r_ah and r_bh from hydrogen
    """

    c_hbnd_a = AtomicHBondFactor(sym_a, chg_a)
    c_hbnd_b = AtomicHBondFactor(sym_b, chg_b)
    c_hbnd_ab = (c_hbnd_a * (r_ah ** 2) + c_hbnd_b * (r_bh ** 2)) / (r_ah ** 2 + r_bh ** 2)

    return c_hbnd_ab


def potHBond(f_dmp_a, f_dmp_hbnd, c_hbnd_ab, r_ab):
    """
    Function for the hydrogen bonding potential over a given atom triple with donor/acceptor atoms a and b at distance r_ab apart, and calculated damping and strength factors
    """

    u = f_dmp_a * f_dmp_hbnd * (c_hbnd_a / (r_ab ** 3))

    return u


def AtomicXBondFactor(sym, chg):
    """
    Atomic factor c_xbnd_x for a halogen acceptor atom with symbol sym and charge chg, used in halogen bonding potential
    """

    c_xbnd_x = k_X[sym] * (math.exp(-1 * k_q1["xbond"]) / (math.exp(-1 * k_q1["xbond"]) + k_q2["xbond"]))

    return c_xbnd_x


def potXBond(f_dmp_theta, f_dmp_xbnd, c_xbnd_x, r_dx):
    """
    Function for the halogen bonding potential over a given atom triple with donor atom d and halogen x at a distance r_dx from one another, and calculated damping and strength factors
    """

    u = f_dmp_theta * f_dmp_xbnd * (c_xbnd_x / (r_dx ** 2))

    return u


def RadiusFromCn(C6_AB, C8_AB):
    """
    Function for the internuclear distance R0_AB used in London dispersion energy calculations
    """
    R0_AB = C6_AB / C8_AB
    R0_AB = math.sqrt(R0_AB)

    return R0_AB


def VdWCutoffRadius(symA, symB):
    """
    Function to calculate the cutoff radius for a pair of atoms with symbols symA and symB as the average of van der Waals radii
    """
    R_A = SymbolToVdWRadius[symA]
    R_B = SymbolToVdWRadius[symB]
    R0_AB = (R_A + R_B) / 2

    return R0_AB


def potPauliRep(rep_disp_AB, symA, symB, r_AB, C6_AB, C8_AB, typ=1):
    """
    Pairwise formula for the Pauli repulsion between two atoms
    """
    # Calculation of valence electron numbers, and where to implement determination of Cn_AB, still to be worked out
    z_eff_A = SymbolToValenceE[symA] * k_z[symA]
    z_eff_B = SymbolToValenceE[symB] * k_z[symB]
    if typ == 1:
        R_0D3 = VdWCutoffRadius(symA, symB)
    elif typ == 2:
        R_0D3 = RadiusFromCn(C6_AB, C8_AB)  # Note still need to confirm whether R_0D3 and R0_AB are actually equivalent
    u = rep_disp_AB * (z_eff_A * z_eff_B / r_AB) * math.exp(-1 * beta_rep * r_AB / (R_0D3 ** (3 / 2)))

    return u


def BJdamping(a, b, symA, symB, typ=1):
    """
    Function for the Becke-Johnson rational damping for atoms a and b
    """
    # Calculation of cutoff radius is optionally performed using either van der Waals radii or C6 and C8 coefficients
    if typ == 1:
        R_0AB = VdWCutoffRadius(symA, symB)
    elif typ == 2:
        C6_AB = (C6[a] + C6[b]) / 2  # Check this is correct for C6
        C8_AB = 1.0  # 1 used as a placeholder until C8 values can be added as a dictionary
        R_0AB = RadiusFromCn(C6_AB, C8_AB)
    damp = a1 * R_0AB + a2

    return damp


def potLondonDisp(rep_disp_AB, C6_AB, C8_AB, BJdamp_AB, r_AB):
    """
    Function for the London dispersion energy under the D3 scheme employing Becke-Johnson rational damping via BJdamp_AB
    """
    # Calculation of rep_disp_AB, and correct values for the other arguments in this function, to be worked out
    sixterm = C6_AB / (r_AB ** 6 + BJdamp_AB ** 6)
    eightterm = C8_AB / (r_AB ** 8 + BJdamp_AB ** 8)
    u = rep_disp_AB * (sixterm + s8 * eightterm)

    return u


def potElectrostatic(elstat_AB, chg_A, chg_B, r_AB):
    """
    Function for the electrostatic potential between atoms A and B
    """
    # Determination of the screening parameter elstat_AB still to be implemented
    u = elstat_AB * (chg_A * chg_B / r_AB)

    return u


def potCSODisp(C6_AB, R_AB, R0_AB):
    """
    Function for the contribution to London dispersion energy from a single pair of atoms, AB, following the D3(CSO) scheme
    """
    # Note there are different values for a1 depending on functional, default used here is B3LYP
    # Also s6 can have other values, and should if aligning with the B2PLYP functional
    a1 = CSO_a1["B3LYP"]
    s6 = 1
    exparg = (R_AB - (2.5 * R0_AB))
    # print(" \nFor C6-only dispersion calculation, R_AB = " + str(R_AB) + ", R0_AB = " + str(R0_AB))
    if exparg < math.log(sys.float_info.min):
        print("Exponential cannot be computed in C6-only dispersion potential")
        print("Problem value of R_AB -(2.5 * R0_AB) = " + str(exparg))
        print("Arising from: R_AB = " + str(R_AB) + ", R0_AB = " + str(R0_AB))
        ProgramAbort()
    elif exparg > math.log(sys.float_info.max):
        # Internuclear distance R_AB is larger than 709, exponential greater than 1.8e308, and computing it would give an error
        # The contribution from a1/exparg is negligible, so set C6indep equal to the limit value it approaches, that is the constant parameter s6
        C6indep = s6
        # print("Setting a1/(R_AB - (2.5 * R0_AB)) = 0 ")
    else:
        C6indep = s6 + a1 / (1 + math.exp(exparg))
    C6dep = C6_AB / (R_AB ** 6 + (2.5 ** 2) ** 6)
    u = C6indep * C6dep

    return u


#############################################################################################################
# Classes for Force Field Terms defined below
#############################################################################################################

class FFStretch:
    """ A stretching potential"""

    def __init__(self, a, b, r0, typ, arg):
        """ (FFStretch, int, int, number, int, [number]) -> NoneType
    
    A stretch potential between atoms number a and b with equilibrium
    distance r0, of type typ with arguments [arg] where arg[2] and arg[3]
    are the atomic symbols of atoms number a and b respectively
    """

        self.atom1 = a
        self.atom2 = b
        self.r0 = r0
        self.typ = typ
        self.k_str = 0.0

        if typ == 1:
            self.typ = typ
            self.k_str = arg[0]
        elif typ == 2:
            self.D = arg[0]
            self.b = arg[1]
        elif typ == 3:
            self.exp_a = bond_exp(arg[2], arg[3])
            # check whether an Atom can be passed to a function this way
            self.k_str = arg[0]
        elif typ == 4:
            self.exp_a = exp_1_3(arg[2], arg[3])
            self.k_str = arg[0]
        else:
            self.typ = 1
            self.k_str = arg[0]

    def __str__(self):
        """ (FFStretch) -> str
    
    Return a string representation of the stretching potential in this format:
    
    (atom1, atom2, r0, type, arguments)
    
    """

        s = '({0}, {1}, {2}, {3}, '.format(self.atom1, self.atom2, self.r0, self.typ)
        r = ''

        if self.typ == 1:
            r = '{0})'.format(self.k_str)
        elif self.typ == 2:
            r = '{0}, {1})'.format(self.D, self.b)
        elif self.typ == 3 or self.typ == 4:
            r = '{0}, {1}'.format(self.k_str, self.exp_a)

        return s + r

    def __repr__(self):
        """ (FFStretch) -> str
    
    Return a string representation of the stretching potential in this format:
    
    (atom1, atom2, r0, type, arguments)
    
    """

        s = '({0}, {1}, {2}, {3}, '.format(self.atom1, self.atom2, self.r0, self.typ)
        r = ''

        if self.typ == 1:
            r = '{0})'.format(self.k_str)
        elif self.typ == 2:
            r = '{0}, {1})'.format(self.D, self.b)

        return s + r

    def setk(self, k):
        """ (FFStretch) -> NoneType

    Set the force constant k_str equal to k
    """
        self.k_str = k

    def energy(self, r):
        """ Returns the energy of this stretching potential at distance r"""

        energy = 0.0
        if self.typ == 1:
            #      print("Using Harmonic potential for stretch") # REMOVE ONCE FIXED
            print("With r = " + str(r) + ", r0 = " + str(self.r0) + ", k = " + str(self.k))
            energy = potHarmonic(r, self.r0, self.k_str)
        elif self.typ == 2:
            #      print("Using Morse potential for stretch") # REMOVE ONCE FIXED
            #      print("With r = " + str(r) + ", r0 = " + str(self.r0) + ", D = " + str(self.D) + ", b = " + str(self.b)) # REMOVE ONCE FIXED 
            energy = potMorse(r, self.r0, D, b)
        elif self.typ == 3 or self.typ == 4:
            #      print("Using GLJ potential for stretch") # REMOVE ONCE FIXED 
            #      print("With r = " + str(r) + ", r0 = " + str(self.r0) + ", k = " + str(self.k_str) + ", a = " + str(self.exp_a)) # REMOVE ONCE FIXED
            energy = potGLJ(r, self.r0, self.k_str, self.exp_a)

        return energy


class FFBend:
    """ A bending potential"""

    def __init__(self, a, b, c, a0, typ, arg):
        """ (FFStretch, int, int, int, number, int, [number]) -> NoneType
    
    A bending potential between atoms number a, b and c with equilibrium
    angle a0, of type typ with arguments [arg] comprising angle force 
    constant, atomic symbols of atoms a, b, and c, and the distances 
    between atoms a and b and b and c
    """

        self.atom1 = a
        self.atom2 = b
        self.atom3 = c
        self.a0 = a0
        if typ == 1:
            self.typ = typ
            self.k = arg[0]
        elif typ == 2:
            self.typ = 2
            self.k_bnd = arg[0]
            r_12 = arg[4]
            r_23 = arg[5]
            f_dmp_12 = DampingFunction(arg[1], arg[2], r_12)
            f_dmp_23 = DampingFunction(arg[2], arg[3], r_23)
            self.f_dmp = f_dmp_12 * f_dmp_23
            self.typ = 1
            self.k = arg[0]

    def __str__(self):
        """ (FFBend) -> str
    
    Return a string representation of the bending potential in this format:
    
    (atom1, atom2, atom3, a0, type, arguments)
    
    """

        s = '({0}, {1}, {2}, {3}, '.format(self.atom1, self.atom2, self.atom3, self.a0, self.typ)

        if self.typ == 1:
            r = '{0})'.format(self.k)

        return s + r

    def __repr__(self):
        """ (FFBend) -> str
    
    Return a string representation of the bending potential in this format:
    
    (atom1, atom2, atom3, a0, type, arguments)
    
    """

        s = '({0}, {1}, {2}, {3}, '.format(self.atom1, self.atom2, self.atom3, self.a0, self.typ)

        if self.typ == 1:
            r = '{0})'.format(self.k)

        return s + r

    def setk(self, newk):
        """ (FFBend) -> NoneType

    Set the bending force constant k equal to newk
    """
        self.k = newk

    def energy(self, a):
        """ Returns the energy of this bending potential at angle a"""

        energy = 0.0
        if self.typ == 1:
            energy = potHarmonic(a, self.a0, self.k)
        elif self.typ == 2:
            if (math.pi - 0.01) <= self.a0 <= (math.pi + 0.01):
                # Tolerance used here is essentially a placeholder, may need changing in either direction
                energy = potBendNearLinear(a, self.a0, self.k_bnd, self.f_dmp)
            else:
                energy = potAnyBend(a, self.a0, self.k_bnd, self.f_dmp)

        return energy


class FFTorsion:
    """ A torsion potential"""

    def __init__(self, a, b, c, d, theta0, typ, arg):
        """ (FFTorsion, int, int, int, int, number, int, [number]) -> NoneType
    
    A torsion potential between atoms number a, b, c and d with equilibrium
    angle theta0, of type typ with arguments [arg] comprising the dihedral
    force constant, the atomic symbols of atoms a, b, c and d, and the 
    ab, bc and cd bond lengths
    """

        self.atom1 = a
        self.atom2 = b
        self.atom3 = c
        self.atom4 = d
        self.theta0 = theta0
        if typ == 1:
            self.typ = typ
            self.k = arg[0]
        elif typ == 2:
            self.typ = typ
            r_12 = arg[5]
            r_23 = arg[6]
            r_34 = arg[7]
            f_dmp_12 = DampingFunction(arg[1], arg[2], r_12)
            f_dmp_23 = DampingFunction(arg[2], arg[3], r_23)
            f_dmp_34 = DampingFunction(arg[3], arg[4], r_34)
            self.f_dmp = f_dmp_12 * f_dmp_23 * f_dmp_34
            self.k = arg[0]
        else:
            self.typ = 1
            self.k = arg[0]

    def __str__(self):
        """ (FFTorsion) -> str
    
    Return a string representation of the torsion potential in this format:
    
    (atom1, atom2, atom3, atom4, theta0, type, arguments)
    
    """

        s = '({0}, {1}, {2}, {3}, {4}, '.format(self.atom1, self.atom2, self.atom3, self.atom4, self.theta0, self.typ)

        if self.typ == 1:
            r = '{0})'.format(self.k)

        return s + r

    def __repr__(self):
        """ (FFTorsion) -> str
    
    Return a string representation of the torsion potential in this format:
    
    (atom1, atom2, atom3, atom4, theta0, type, arguments)
    
    """

        s = '({0}, {1}, {2}, {3}, {4}, '.format(self.atom1, self.atom2, self.atom3, self.atom4, self.theta0, self.typ)

        if self.typ == 1:
            r = '{0})'.format(self.k)

        return s + r

    def energy(self, theta):
        """ Returns the energy of this torsion potential at angle theta"""

        energy = 0.0
        if self.typ == 1:
            energy = potSimpleCosine(theta, self.theta0, self.k)
        elif typ == 2:
            energy = potTorsion(theta, self.theta0, self.f_dmp, self.k)
        # Will need two cases, one for non-rotatable bonds, the other for rotatable bonds.
        # Probably best to implement via types
        return energy


class FFInversion:
    """ An inversion potential for 3-fold coordinate atoms"""

    def __init__(self, a, b, c, d, phi0, typ, arg):
        """ (FFInversion, int, int, int, int, number, int, [number]) -> NoneType

    An inversion potential between atoms number a, b, c and d (where a is the
    central atom bonded to the other three) with equilibrium out of plane
    angle phi0, of type typ with arguments [arg]
    """

        self.atom1 = a
        self.atom2 = b
        self.atom3 = c
        self.atom4 = d
        self.phi0 = phi0
        if typ == 1:
            self.typ = 1
            self.k_inv = arg[0]
        elif typ == 2:
            self.typ = 2
            self.k_inv = arg[
                0]  # Will need to check there is an appropriate force constant locatable for use here, eventually will require fitting to Hessian
            # Damping constant needs to be checked - set up below on the assumption that a product of distance dependent damping functions for the three bonds will do
            r_12 = arg[5]
            r_13 = arg[6]
            r_14 = arg[7]
            f_dmp_12 = DampingFunction(arg[1], arg[2], r_12)
            f_dmp_13 = DampingFunction(arg[1], arg[3], r_13)
            f_dmp_14 = DampingFunction(arg[1], arg[4], r_14)
            self.f_dmp = f_dmp_12 * f_dmp_13 * f_dmp_14
        else:
            self.typ = 1
            self.k_inv = arg[0]

    def __str__(self):
        """ (FFInversion) -> str

    Return a string representation of the bending potential in this format:

    (atom1, atom2, atom3, atom4, phi0, type, arguments)

    """

        s = '({0}, {1}, {2}, {3}, {4}, '.format(self.atom1, self.atom2, self.atom3, self.atom4, self.phi0, self.typ)

        if self.typ == 1:
            r = '{0})'.format(self.k_inv)

        return s + r

    def __repr__(self):
        """ (FFInversion) -> str

    Return a string representation of the bending potential in this format:

    (atom1, atom2, atom3, atom4, phi0, type, arguments)

    """

        s = '({0}, {1}, {2}, {3}, {4}, '.format(self.atom1, self.atom2, self.atom3, self.atom4, self.phi0, self.typ)

        if self.typ == 1:
            r = '{0})'.format(self.k_inv)

        return s + r

    def setk(self, newk):
        """ (FFInversion) -> NoneType

    Set the inversion force constant for this potential equal to newk
    """
        self.k_inv = newk

    def energy(self, phi):
        """ Returns the energy of this inversion potential at out of plane angle phi"""

        energy = 0.0
        if self.typ == 1:
            energy = potHarmonic(phi, phi0, self.k_inv)  # or use another sutiable simple potential here
        elif self.typ == 2:
            if (math.pi - 0.01) <= self.phi0 <= (math.pi + 0.01):
                # Tolerance used here is essentially a placeholder, may need changing in either direction
                energy = potBendNearLinear(phi, self.phi0, self.k_inv, self.f_dmp)
            else:
                energy = potAnyBend(phi, self.phi0, self.k_inv, self.f_dmp)

        return energy


class FFHBond:
    """ A hydrogen bonding potential for an atom triple """

    def __init__(self, a, b, c, theta, typ, arg):
        """ (FFHBond, int, int, int, number, int, [number]) -> NoneType

    A hydrogen bonding potential between atoms number a, b and c, where b is
    hydrogen and the abc angle is theta, of type typ with arguments [arg]
    comprising the atomic symbols and charges of atoms a, b and c, and distances
    r_ab, r_bc, r_ab
    """

        # Including type may not be necessary, has been done here just for consistency

        self.atomA = a
        self.atomH = b
        self.atomB = c
        self.theta = theta
        if typ == 1:
            self.typ = typ
            self.symA = arg[0]
            self.chgA = arg[1]
            self.symH = arg[2]
            self.chgH = arg[3]
            self.symB = arg[4]
            self.chgB = arg[5]
            self.r_AH = arg[6]
            self.r_BH = arg[7]
            self.r_AB = arg[8]

    def __str__(self):
        """ (FFHbond) -> str

    Return a string representation of the hydrogen bonding potential in this format:

    (atomA, atomH, atomB, theta, type, arguments)

    """

        s = '({0}, {1}, {2}, {3), '.format(self.atomA, self.atomH, self.atomB, self.theta, self.typ)

        if self.typ == 1:
            r = ')'
        # Written now to just close the list, should ultimately expand to allow printing of all arguments in line with other classes

        return s + r

    def __repr__(self):
        """ (FFHbond) -> str

    Return a string representation of the hydrogen bonding potential in this format:

    (atomA, atomH, atomB, theta, type, arguments)

    """

        s = '({0}, {1}, {2}, {3), '.format(self.atomA, self.atomH, self.atomB, self.theta, self.typ)

        if self.typ == 1:
            r = ')'

        return s + r

    def energy(self):
        """ Returns the energy for this hydrogen bonding potential """

        c_hbnd_AB = HBondStrengthFactor(self.symA, self.chgA, self.r_AH, self.symB, self.chgB, self.r_BH)

        f_dmp_theta = AngleDamping(self.theta)

        f_dmp_hbnd = HBondDamping(self.symA, self.symB, self.r_AB)

        energy = potHBond(f_dmp_theta, f_dmp_hbnd, c_hbnd_AB, self.r_AB)

        return energy


# Note - will need to check in adding hbond triples that the indices used above for arguments come out correct 


#############################################################################################################
# Atom class and class methods to be defined below
#############################################################################################################

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

    def __str__(self):
        """ (Atom) -> str
    
    Return a string representation of this Atom in this format:
    
      (SYMBOL, X, Y, Z)
    """

        return '({0}, {1}, {2}, {3})'.format(self.symbol, self.coord[0], self.coord[1], self.coord[2])

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


#############################################################################################################
# Molecule class and class methods to be defined below
#############################################################################################################

class Molecule:
    """A molecule with a name, charge and a list of atoms"""

    def __init__(self, name, charge):
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
        self.stretch = []
        self.str13 = []
        self.bend = []
        self.tors = []
        self.inv = []
        self.hbonds = []
        self.hatoms = []
        self.highENatoms = []
        self.halogens = []
        self.H_QM = numpy.zeros((3, 3))  # Array size arbitrary, just a placeholder for type 

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

    def __str__(self):
        """ (Molecule) -> str
    
    Return a string representation of this Molecule in this format:
    (NAME, CHARGE, MULT, (ATOM1, ATOM2, ...))
    """

        res = ''
        for atom in self.atoms:
            res = res + str(atom) + ', '
        res = res[:-2]
        return '({0}, {1}, {2}, ({3}))'.format(self.name, self.charge, self.mult, res)

    def __repr__(self):
        """ (Molecule) -> str
    
    Return a string representation of this Molecule in this format:
    (NAME, CHARGE, MULT, (ATOM1, ATOM2, ...))
    """

        res = ''
        for atom in self.atoms:
            res = res + str(atom) + ', '
        res = res[:-2]
        return '({0}, {1}, {2}, ({3}))'.format(self.name, self.charge, self.mult, res)

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
        print(str(n))
        if len(cartCoordinates) != 3 * n:
            ProgramWarning()
            print("Cannot update geometry of " + str(self.name) + ", supplied coordinates do not match number of atoms")
        else:
            for i in range(n):
                self.movatom(i, cartCoordinates[3 * i], cartCoordinates[3 * i + 1], cartCoordinates[3 * i + 2])

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
        theta = numpy.arccos(argument)

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
        theta = numpy.arccos(argument)

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
        vnormal_1 = numpy.cross(end_1, bridge)
        vnormal_2 = numpy.cross(bridge, end_2)

        # Construct a set of orthogonal basis vectors to define a frame with vnormal_2 as the x axis
        vcross = numpy.cross(vnormal_2, bridge)
        norm_vn2 = numpy.linalg.norm(vnormal_2)
        norm_b = numpy.linalg.norm(bridge)
        norm_vc = numpy.linalg.norm(vcross)
        basis_vn2 = [vnormal_2[i] / norm_vn2 for i in range(3)]
        basis_b = [bridge[i] / norm_b for i in range(3)]
        basis_cv = [vcross[i] / norm_vc for i in range(3)]

        # Find the signed angle between vnormal_1 and vnormal_2 in the new frame
        vn1_coord_n2 = numpy.dot(vnormal_1, basis_vn2)
        vn1_coord_vc = numpy.dot(vnormal_1, basis_cv)
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
        vnormal_1 = numpy.cross(end_1, bridge)
        vnormal_2 = numpy.cross(bridge, end_2)

        # Construct a set of orthogonal basis vectors to define a frame with vnormal_2 as the x axis
        vcross = numpy.cross(vnormal_2, bridge)
        norm_vn2 = numpy.linalg.norm(vnormal_2)
        norm_b = numpy.linalg.norm(bridge)
        norm_vc = numpy.linalg.norm(vcross)
        basis_vn2 = [vnormal_2[i] / norm_vn2 for i in range(3)]
        basis_b = [bridge[i] / norm_b for i in range(3)]
        basis_cv = [vcross[i] / norm_vc for i in range(3)]

        # Find the signed angle between vnormal_1 and vnormal_2 in the new frame
        vn1_coord_n2 = numpy.dot(vnormal_1, basis_vn2)
        vn1_coord_vc = numpy.dot(vnormal_1, basis_cv)
        psi = math.atan2(vn1_coord_vc, vn1_coord_n2)

        return psi

    def outofplaneangle(self, i):
        """ (Molecule) -> number (in radians)

    Report the out of plane angle described by a set of four atoms in the inv list
    """

        # Calculate the vectors along bonds, and construct a vector plane_norm orthogonal to the plane of end atoms
        threefold = self.threefolds[i]
        atom_c = self.atoms[threefold[0]]
        atom_e1 = self.atoms[threefold[1]]
        atom_e2 = self.atoms[threefold[2]]
        atom_e3 = self.atoms[threefold[3]]
        bond_1 = [atom_e1.coord[i] - atom_c.coord[i] for i in range(3)]
        bond_2 = [atom_e2.coord[i] - atom_c.coord[i] for i in range(3)]
        bond_3 = [atom_e3.coord[i] - atom_c.coord[i] for i in range(3)]
        inplane_12 = [atom_e2.coord[i] - atom_e1.coord[i] for i in range(3)]
        inplane_13 = [atom_e3.coord[i] - atom_e1.coord[i] for i in range(3)]
        plane_norm = numpy.cross(inplane_12, inplane_13)

        # Construct vectors between end atoms and the projection of the central atom on the end-atom plane
        cross_1 = numpy.cross(bond_1, plane_norm)
        cross_2 = numpy.cross(bond_2, plane_norm)
        cross_3 = numpy.cross(bond_2, plane_norm)
        inplane_1 = numpy.cross(plane_norm, cross_1) / numpy.dot(plane_norm, plane_norm)
        inplane_2 = numpy.cross(plane_norm, cross_2) / numpy.dot(plane_norm, plane_norm)
        inplane_3 = numpy.cross(plane_norm, cross_3) / numpy.dot(plane_norm, plane_norm)

        # Calculate the out of plane angle for each of the three bonds
        cos_phi1 = numpy.dot(bond_1, inplane_1) / (numpy.linalg.norm(bond_1) * numpy.linalg.norm(inplane_1))
        cos_phi2 = numpy.dot(bond_2, inplane_2) / (numpy.linalg.norm(bond_2) * numpy.linalg.norm(inplane_2))
        cos_phi3 = numpy.dot(bond_3, inplane_3) / (numpy.linalg.norm(bond_3) * numpy.linalg.norm(inplane_3))

        if (1.0 - (10 ** -15)) <= cos_phi1 and cos_phi1 <= (1.0 + (10 ** -15)):
            phi1 = 0.0
        else:
            phi1 = numpy.arccos(cos_phi1)
        if (1.0 - (10 ** -15)) <= cos_phi2 and cos_phi2 <= (1.0 + (10 ** -15)):
            phi2 = 0.0
        else:
            phi2 = numpy.arccos(cos_phi2)
        if (1.0 - (10 ** -15)) <= cos_phi3 and cos_phi3 <= (1.0 + (10 ** -15)):
            phi3 = 0.0
        else:
            phi3 = numpy.arccos(cos_phi3)

        # Take the numerical average of the three out of plane angles
        # Note - other schemes for obtaining a single out of plane angle could be investigated
        phi = (phi1 + phi2 + phi3) / 3

        return phi

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
            i.mass * ((Ang2Bohr(i.coord[1]) * Ang2Bohr(i.coord[1])) + (Ang2Bohr(i.coord[2]) * Ang2Bohr(i.coord[2]))))
            Ixy = Ixy - i.mass * Ang2Bohr(i.coord[0]) * Ang2Bohr(i.coord[1])
            Ixz = Ixz - i.mass * Ang2Bohr(i.coord[0]) * Ang2Bohr(i.coord[2])
            Iyx = Iyx - i.mass * Ang2Bohr(i.coord[1]) * Ang2Bohr(i.coord[0])
            Iyy = Iyy + (
            i.mass * ((Ang2Bohr(i.coord[0]) * Ang2Bohr(i.coord[0])) + (Ang2Bohr(i.coord[2]) * Ang2Bohr(i.coord[2]))))
            Iyz = Iyz - i.mass * Ang2Bohr(i.coord[1]) * Ang2Bohr(i.coord[2])
            Izx = Izx - i.mass * Ang2Bohr(i.coord[2]) * Ang2Bohr(i.coord[0])
            Izy = Izy - i.mass * Ang2Bohr(i.coord[2]) * Ang2Bohr(i.coord[1])
            Izz = Izz + (
            i.mass * ((Ang2Bohr(i.coord[0]) * Ang2Bohr(i.coord[0])) + (Ang2Bohr(i.coord[1]) * Ang2Bohr(i.coord[1]))))
        inertiaTensor.append([Ixx, Ixy, Ixz])
        inertiaTensor.append([Iyx, Iyy, Iyz])
        inertiaTensor.append([Izx, Izy, Izz])
        inertiaTensor = numpy.matrix(inertiaTensor)

        # Diagonalise inertia tensor
        inertiaMoments, inertialAxes = numpy.linalg.eig(inertiaTensor)

        # Orthogonalise eigenvectors (only sometimes necessary)...
        inertialAxes, r = numpy.linalg.qr(inertialAxes)

        # Sort moments from highest to lowest
        idx = inertiaMoments.argsort()[::-1]
        inertiaMoments = inertiaMoments[idx]
        inertialAxes = inertialAxes[:, idx]

        # Transform molecular coordinates into new frame of principal axes of inertia
        for i in self.atoms:
            vector = [i.coord[0], i.coord[1], i.coord[2]]
            vector = numpy.matrix(vector)
            vector = numpy.matrix.transpose(inertialAxes).dot(numpy.matrix.transpose(vector))
            vector = numpy.array(vector).flatten().tolist()
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

    def setHessian(self, H):
        """ (Molecule) -> NoneType

    Set the Quantum Mechanically calculated Hessian, H_QM, equal to H
    """
        self.H_QM = H

    def setQMenergy(self, E):
        """ (Molecule) -> NoneType

    Set the Quantum Mechanically calculated equilibrium energy, Ee_QM, equal to E
    """
        self.Ee_QM = E

    def addFFStretch(self, a, b, r0, typ, arg):
        """ (Molecule) -> NoneType

    Adds a stretching potential between atoms a and b to the list of stretches
    """
        # Make sure a < b
        if a < b:
            c = a
            d = b
        else:
            c = b
            d = a

        # Note: There's no check if the stretch already exists. Mainly because there's no reason not to have
        # two different functions adding energy to the same "stretch"

        # Append stretch to list if doesn't exist and is plausible
        if a >= 0 and b >= 0 and a <= len(self.atoms) and b <= len(self.atoms) and c != d:
            self.stretch.append(FFStretch(c, d, r0, typ, arg))
        #    print("Adding FFStretch with bond " + str([a, b]) + " and fc = " + str(arg[0])) # REMOVE ONCE FIXED

    def addFFStr13(self, a, b, r0, typ, arg):
        """ (Molecule) -> NoneType

    Adds a stretching potential between atoms a and b to the list of stretches
    """
        # Make sure a < b
        if a < b:
            c = a
            d = b
        else:
            c = b
            d = a

        # Note: There's no check if the stretch already exists. Mainly because there's no reason not to have
        # two different functions adding energy to the same "stretch"

        # Append stretch to list if doesn't exist and is plausible
        if a >= 0 and b >= 0 and a <= len(self.atoms) and b <= len(self.atoms) and c != d:
            self.str13.append(FFStretch(c, d, r0, typ, arg))
        #    print("Adding FFStretch with bond " + str([a, b]) + " and fc = " + str(arg[0])) # REMOVE ONCE FIXED

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

    def addFFBend(self, a, b, c, a0, typ, arg):
        """ (Molecule) -> NoneType

    Adds an angle bend potential between atoms a, b and c to the list of angle bends.
    """

        # Note: There's no check if the bend already exists. Mainly because there's no reason not to have
        # two different functions adding energy to the same "bend"

        # Append bend to list if doesn't exist and is plausible
        # (it's a bit unsatisfactory, but I don't think there are
        # better sanity checks)
        if a >= 0 and b >= 0 and c >= 0 and a <= len(self.atoms) and b <= len(self.atoms) and c <= len(
                self.atoms) and a != b and a != c and b != c:
            self.bend.append(FFBend(a, b, c, a0, typ, arg))

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

    def addFFTorsion(self, a, b, c, d, theta0, typ, arg):
        """ (Molecule) -> NoneType

    Adds a dihedral torsion potential between atoms a, b, c and d to the list of dihedral torsions.
    """

        # Note: There's no check if the bend already exists. Mainly because there's no reason not to have
        # two different functions adding energy to the same "torsion"

        # Append torsion to list if doesn't exist and is plausible
        # (it's a bit unsatisfactory, but I don't think there are
        # better sanity checks)
        if a >= 0 and b >= 0 and c >= 0 and d >= 0 and a <= len(self.atoms) and b <= len(self.atoms) and c <= len(
                self.atoms) and d <= len(self.atoms) and a != b and a != c and a != d and b != c and b != d and c != d:
            self.tors.append(FFTorsion(a, b, c, d, theta0, typ, arg))

    def addThreefold(self, a, b, c, d):
        """ (Molecule) -> NoneType

    Adds the atom a and its bonding partners b, c and d to the list of 3-fold coordinated atom groups
    """

        # Check if these atoms are already in the list inv
        exists = False
        for i in self.threefolds:
            if i == [a, b, c, d]:
                exists = True

        # Append this three-fold coordinated group to list if it doesn't exist and is plausible
        # Check currently based upon those for angles and dihedrals, for lack of a better alternative
        if exists == False and a >= 0 and b >= 0 and c >= 0 and d >= 0 and a <= len(self.atoms) and b <= len(
                self.atoms) and c <= len(self.atoms) and d <= len(
                self.atoms) and a != b and a != c and a != d and b != c and b != d and c != d:
            self.threefolds.append([a, b, c, d])

    def addFFInversion(self, a, b, c, d, phi0, typ, arg):
        """ (Molecule) -> NoneType

    Adds an out-of-plane inversion type potential between atoms a, b, c and d to the list of inversion
    """

        # Note: As for other potentials, no check to see if this inversion exists already since there's no
        # reason not to allow two different functions to contribute energy to the same 'inversion'

        # Append inversion to list if it doesn't exist and is plausible
        if a >= 0 and b >= 0 and c >= 0 and d >= 0 and a <= len(self.atoms) and b <= len(self.atoms) and c <= len(
                self.atoms) and d <= len(self.atoms) and a != b and a != c and a != d and b != c and b != d and c != d:
            self.inv.append(FFInversion(a, b, c, d, phi0, typ, arg))

    def addFFHBond(self, a, b, c, theta, typ, arg):
        """ (Molecule) -> NoneType

    Adds a hydrogen bonding interaction potential between donor/acceptor atoms a and c and hydrogen atom b to the list of hbonds
    """

        # Note: No check to see if this interaction exists already, consistent with the treatment of other potentials

        # Append new hydrogen bonding interaction to the list if it's plausible
        # Where plausibility is checked by distinct atoms only
        # (Checks for distance over which hydrogen bonds could exist and the electronegativity of atoms involved could be added, but will be used to identify interactions to add in any case)
        if a >= 0 and b >= 0 and c >= 0 and a <= len(self.atoms) and b <= len(self.atoms) and c <= len(
                self.atoms) and a != b and a != c and b != c:
            self.hbonds.append(FFHBond(a, b, c, theta, typ, arg))

    def addHAtom(self, a):
        """ (Molecule) -> NoneType

    Adds the atom a to the list of hydrogen atoms present in the molecule
    """

        # Note currently no check for whether this atom is already in the list
        # May be worth implementing such a check. Could also check atomic symbol, though that will be done before calling this method
        if a >= 0 and a <= len(self.atoms):
            self.hatoms.append(a)

    def addhighENatom(self, a):
        """ (Molecule) -> NoneType

    Adds the atom a to the list of atoms sufficiently electronegative to be involved in hydrogen bonding
    """

        # Note currently no check for whether the atom is already in the list
        # As with list of hydrogens, could implement atomic symbol check for extra care
        if a >= 0 and a <= len(self.atoms):
            self.highENatoms.append(a)

    def addXatom(self, a):
        """ (Molecule) -> NoneType

    Adds the atom a to the list of halogens present in the molecule
    """

        # Currently no check to see if this atom already belongs to the list
        # Also no check on atomic symbol at present, but one could be implemented
        if a >= 0 and a <= len(self.atoms):
            self.halogens.append(a)

    def screen_ES(self, a, b):
        """ (Molecule) -> Number
    
    Calculates the number of bonds between atoms a and b, and returns the appropriate value of the topological screening parameter elstat_AB
    """
        # Determine the number of bonds separating atoms a and b
        #    print("Determining number of bonds between atoms: " + str((a, b))) # REMOVE ONCE FIXED
        n_bonds = 4
        if a == b:
            n_bonds = 0
        elif [a, b] in self.bonds or [b, a] in self.bonds:
            n_bonds = 1
        else:
            bonds_a = []
            bonds_b = []
            # Find bonds from atom a to any other atom, and from atom b to any other atom
            for bond in self.bonds:
                if bond[0] == a:
                    bonds_a.append(bond)
                elif bond[1] == a:
                    bonds_a.append([bond[1], bond[0]])
                elif bond[0] == b:
                    bonds_b.append([bond[1], bond[0]])
                elif bond[1] == b:
                    bonds_b.append(bond)
                    # Check whether there is a common third atom bonded to both a and b, stopping if one such atom is found 
                #      print("Bonds from atom " + str(a) +": " + str(bonds_a)) # REMOVE ONCE FIXED
                #      print("Bonds from atom " + str(b) +": " + str(bonds_b)) # REMOVE ONCE FIXED
            for bd_a in bonds_a:
                for bd_b in bonds_b:
                    if bd_a[1] == bd_b[0]:
                        n_bonds = 2
                        break
                if n_bonds == 2:
                    break
            # If not, find all bonds 1 bond away from each of atom a and atom b
            if n_bonds > 2:
                bonds_a1 = []
                bonds_b1 = []
                for bond in self.bonds:
                    for bd_a in bonds_a:
                        if bd_a[1] == bond[0]:
                            bonds_a1.append(bond)
                        elif bd_a[1] == bond[1]:
                            bonds_a1.append([bond[1], bond[0]])
                    for bd_b in bonds_b:
                        if bd_b[0] == bond[0]:
                            bonds_b1.append([bond[1], bond[0]])
                        elif bd_b[0] == bond[1]:
                            bonds_b1.append(bond)
                            # Then check whether there is one of those bonds in common to both a and b, ending the loop if so
                        #        print("Bonds one bond away from atom " + str(a) +": " + str(bonds_a1)) # REMOVE ONCE FIXED
                        #        print("Bonds one bond away from atom " + str(b) +": " + str(bonds_b1)) # REMOVE ONCE FIXED
                for bd_a1 in bonds_a1:
                    for bd_b1 in bonds_b1:
                        if bd_a1 == bd_b1:
                            n_bonds = 3
                            break
                    if n_bonds == 3:
                        break
                        # All cases with more than three bonds between atoms a and b are treated identically, so do not check for 4 bonds or more explicitly  
                        # Determine the value of the screening parameter for the number of covalent bonds separating a and b
                    #    print("Number of bonds between atoms " + str((a, b)) + " = " + str(n_bonds)) # REMOVE ONCE FIXED
        if n_bonds <= 2:
            elstat_AB = 0
        elif n_bonds == 3:
            elstat_AB = E_ES_14
        elif n_bonds >= 4:
            elstat_AB = 1
        #    print("Screening parameter value for atoms " + str((a, b)) + " = " + str(elstat_AB)) # REMOVE ONCE FIXED
        return elstat_AB

    def screen_RepDisp(self, a, b):
        """ (Molecule) -> Number
    
    Calculates the number of bonds between atoms a and b, and returns the appropriate value of the topological screening parameter rep_disp_AB 
    """
        # Determine the number of bonds separating atoms a and b
        n_bonds = 5
        if a == b:
            n_bonds = 0
        elif [a, b] in self.bonds or [b, a] in self.bonds:
            n_bonds = 1
        else:
            bonds_a = []
            bonds_b = []
            # Find bonds from atom a to any other atom, and from atom b to any other atom
            for bond in self.bonds:
                if bond[0] == a:
                    bonds_a.append(bond)
                elif bond[1] == a:
                    bonds_a.append([bond[1], bond[0]])
                elif bond[0] == b:
                    bonds_b.append([bond[1], bond[0]])
                elif bond[1] == b:
                    bonds_b.append(bond)
            # Check whether there is a common third atom bonded to both a and b, stopping if one such atom is found 
            for bd_a in bonds_a:
                for bd_b in bonds_b:
                    if bd_a[1] == bd_b[0]:
                        n_bonds = 2
                        break
                if n_bonds == 2:
                    break
            # If not, find all bonds 1 bond away from each of atom a and atom b
            if n_bonds > 2:
                bonds_a1 = []
                bonds_b1 = []
                for bond in self.bonds:
                    for bd_a in bonds_a:
                        if bd_a[1] == bond[0]:
                            bonds_a1.append(bond)
                        elif bd_a[1] == bond[1]:
                            bonds_a1.append([bond[1], bond[0]])
                    for bd_b in bonds_b:
                        if bd_b[0] == bond[0]:
                            bonds_b1.append([bond[1], bond[0]])
                        elif bd_b[0] == bond[1]:
                            bonds_b1.append(bond)
                # Then check whether there is one of those bonds in common to both a and b, ending the loop if so
                for bd_a1 in bonds_a1:
                    for bd_b1 in bonds_b1:
                        if bd_a1 == bd_b1:
                            n_bonds = 3
                            break
                    if n_bonds == 3:
                        break
            # If not, find all bonds 2 bonds away from each of atom a and atom b
            if n_bonds > 3:
                bonds_a2 = []
                bonds_b2 = []
                for bond in self.bonds:
                    for bd_a1 in bonds_a1:
                        if bd_a1[1] == bond[0]:
                            bonds_a2.append(bond)
                        elif bd_a1[1] == bond[1]:
                            bonds_a2.append([bond[1], bond[0]])
                    for bd_b1 in bonds_b1:
                        if bd_b1[0] == bond[0]:
                            bonds_b2.append([bond[1], bond[0]])
                        elif bd_b1[0] == bond[1]:
                            bonds_b2.append(bond)
                # Then check whether one of these bonds is common to both a and b, ending the loop if so
                for bd_a2 in bonds_a2:
                    for bd_b2 in bonds_b2:
                        if bd_a2 == bd_b2:
                            n_bonds = 4
                            break
                    if n_bonds == 4:
                        break
                        # All cases with more than four bonds between atoms a and b are treated identically, so do not check for 1,6- and higher cases explicitly
        # Determine the value of the screening parameter for the number of covalent bonds separating a and b
        if n_bonds <= 2:
            rep_disp_AB = 0
        elif n_bonds == 3 or n_bonds == 4:
            rep_disp_AB = E_disp_rep
        elif n_bonds >= 5:
            rep_disp_AB = 1

        return rep_disp_AB

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

    def FFEnergy(self, cartCoordinates, verbosity=0, dtyp=1):
        """ (Molecule) -> number (Force Field energy)

      Returns a number containing the molecular energy according to the current Force Field definition at structure
      specified by the provided cartesian coordinates.
      The dispersion correction used is specified by dtyp, with 1 for C6-only calculating cutoff radius from van der Waals radii, 2 for full D3 using C6 and C8 coefficients
    """

        energy = 0.0
        if verbosity >= 1:
            print("Initial energy for calculation = " + str(energy))
        energy = energy + self.Ee_QM
        if verbosity >= 1:
            print("With QM energy of equilibrium structure, energy = " + str(energy))

            #  if verbosity >= 1: # REMOVE ONCE FIXED
            #    print("Omitting stretching energy") # REMOVE ONCE FIXED

        for i in self.stretch:
            distance = (cartCoordinates[3 * i.atom1] - cartCoordinates[3 * i.atom2]) ** 2
            distance += (cartCoordinates[3 * i.atom1 + 1] - cartCoordinates[3 * i.atom2 + 1]) ** 2
            distance += (cartCoordinates[3 * i.atom1 + 2] - cartCoordinates[3 * i.atom2 + 2]) ** 2
            distance = math.sqrt(distance)
            #      print("Adding energy for bond " + str([i.atom1, i.atom2]) + ", distance = " + str(distance) + " energy = " + str(i.energy(distance)) + " to total")
            energy = energy + i.energy(distance)
        if verbosity >= 1:
            print("With bond stretches, energy = " + str(energy))

        for i in self.str13:
            distance = (cartCoordinates[3 * i.atom1] - cartCoordinates[3 * i.atom2]) ** 2
            distance += (cartCoordinates[3 * i.atom1 + 1] - cartCoordinates[3 * i.atom2 + 1]) ** 2
            distance += (cartCoordinates[3 * i.atom1 + 2] - cartCoordinates[3 * i.atom2 + 2]) ** 2
            distance = math.sqrt(distance)
            energy = energy + i.energy(distance)
        if verbosity >= 1:
            print("With 1,3-stretches, energy = " + str(energy))

        for i in self.bend:
            d_bond_1 = (cartCoordinates[3 * i.atom1] - cartCoordinates[3 * i.atom2]) ** 2
            d_bond_1 += (cartCoordinates[3 * i.atom1 + 1] - cartCoordinates[3 * i.atom2 + 1]) ** 2
            d_bond_1 += (cartCoordinates[3 * i.atom1 + 2] - cartCoordinates[3 * i.atom2 + 2]) ** 2
            d_bond_1 = math.sqrt(d_bond_1)

            d_bond_2 = (cartCoordinates[3 * i.atom2] - cartCoordinates[3 * i.atom3]) ** 2
            d_bond_2 += (cartCoordinates[3 * i.atom2 + 1] - cartCoordinates[3 * i.atom3 + 1]) ** 2
            d_bond_2 += (cartCoordinates[3 * i.atom2 + 2] - cartCoordinates[3 * i.atom3 + 2]) ** 2
            d_bond_2 = math.sqrt(d_bond_2)

            d_non_bond = (cartCoordinates[3 * i.atom1] - cartCoordinates[3 * i.atom3]) ** 2
            d_non_bond += (cartCoordinates[3 * i.atom1 + 1] - cartCoordinates[3 * i.atom3 + 1]) ** 2
            d_non_bond += (cartCoordinates[3 * i.atom1 + 2] - cartCoordinates[3 * i.atom3 + 2]) ** 2
            d_non_bond = math.sqrt(d_non_bond)

            # Use those distances and the cosine rule to calculate bond angle theta
            numerator = d_bond_1 ** 2 + d_bond_2 ** 2 - d_non_bond ** 2
            denominator = 2 * d_bond_1 * d_bond_2
            argument = numerator / denominator
            theta = numpy.arccos(argument)
            energy = energy + i.energy(theta)
        if verbosity >= 1:
            print("With bends, energy = " + str(energy))

        for i in self.tors:
            # Calculate the vectors lying along bonds, and their cross products
            atom_e1 = [cartCoordinates[3 * i.atom1], cartCoordinates[3 * i.atom1 + 1], cartCoordinates[3 * i.atom1 + 2]]
            atom_b1 = [cartCoordinates[3 * i.atom2], cartCoordinates[3 * i.atom2 + 1], cartCoordinates[3 * i.atom2 + 2]]
            atom_b2 = [cartCoordinates[3 * i.atom3], cartCoordinates[3 * i.atom3 + 1], cartCoordinates[3 * i.atom3 + 2]]
            atom_e2 = [cartCoordinates[3 * i.atom4], cartCoordinates[3 * i.atom4 + 1], cartCoordinates[3 * i.atom4 + 2]]
            end_1 = [atom_e1[i] - atom_b1[i] for i in range(3)]
            bridge = [atom_b1[i] - atom_b2[i] for i in range(3)]
            end_2 = [atom_b2[i] - atom_e2[i] for i in range(3)]
            vnormal_1 = numpy.cross(end_1, bridge)
            vnormal_2 = numpy.cross(bridge, end_2)

            # Construct a set of orthogonal basis vectors to define a frame with vnormal_2 as the x axis
            vcross = numpy.cross(vnormal_2, bridge)
            norm_vn2 = numpy.linalg.norm(vnormal_2)
            norm_b = numpy.linalg.norm(bridge)
            norm_vc = numpy.linalg.norm(vcross)
            basis_vn2 = [vnormal_2[i] / norm_vn2 for i in range(3)]
            basis_b = [bridge[i] / norm_b for i in range(3)]
            basis_cv = [vcross[i] / norm_vc for i in range(3)]

            # Find the signed angle between vnormal_1 and vnormal_2 in the new frame
            vn1_coord_n2 = numpy.dot(vnormal_1, basis_vn2)
            vn1_coord_vc = numpy.dot(vnormal_1, basis_cv)
            psi = math.atan2(vn1_coord_vc, vn1_coord_n2)
            energy = energy + i.energy(psi)
        if verbosity >= 1:
            print("With torsion, energy = " + str(energy))

        for i in self.inv:
            # Calculate the vectors along bonds, and construct a vector plane_norm orthogonal to the plane of end ato
            atom_c = [cartCoordinates[3 * i.atom1], cartCoordinates[3 * i.atom1 + 1], cartCoordinates[3 * i.atom1 + 2]]
            atom_e1 = [cartCoordinates[3 * i.atom2], cartCoordinates[3 * i.atom2 + 1], cartCoordinates[3 * i.atom2 + 2]]
            atom_e2 = [cartCoordinates[3 * i.atom3], cartCoordinates[3 * i.atom3 + 1], cartCoordinates[3 * i.atom3 + 2]]
            atom_e3 = [cartCoordinates[3 * i.atom4], cartCoordinates[3 * i.atom4 + 1], cartCoordinates[3 * i.atom4 + 2]]
            bond_1 = [atom_e1[i] - atom_c[i] for i in range(3)]
            bond_2 = [atom_e2[i] - atom_c[i] for i in range(3)]
            bond_3 = [atom_e3[i] - atom_c[i] for i in range(3)]
            inplane_12 = [atom_e2[i] - atom_e1[i] for i in range(3)]
            inplane_13 = [atom_e3[i] - atom_e1[i] for i in range(3)]
            plane_norm = numpy.cross(inplane_12, inplane_13)

            # Construct vectors between end atoms and the projection of the central atom on the end-atom pla
            cross_1 = numpy.cross(bond_1, plane_norm)
            cross_2 = numpy.cross(bond_2, plane_norm)
            cross_3 = numpy.cross(bond_2, plane_norm)
            inplane_1 = numpy.cross(plane_norm, cross_1) / numpy.dot(plane_norm, plane_norm)
            inplane_2 = numpy.cross(plane_norm, cross_2) / numpy.dot(plane_norm, plane_norm)
            inplane_3 = numpy.cross(plane_norm, cross_3) / numpy.dot(plane_norm, plane_norm)

            # Calculate the out of plane angle for each of the three bonds
            cos_phi1 = numpy.dot(bond_1, inplane_1) / (numpy.linalg.norm(bond_1) * numpy.linalg.norm(inplane_1))
            cos_phi2 = numpy.dot(bond_2, inplane_2) / (numpy.linalg.norm(bond_2) * numpy.linalg.norm(inplane_2))
            cos_phi3 = numpy.dot(bond_3, inplane_3) / (numpy.linalg.norm(bond_3) * numpy.linalg.norm(inplane_3))

            if (1.0 - (10 ** -15)) <= cos_phi1 and cos_phi1 <= (1.0 + (10 ** -15)):
                phi1 = 0.0
            else:
                phi1 = numpy.arccos(cos_phi1)
            if (1.0 - (10 ** -15)) <= cos_phi2 and cos_phi2 <= (1.0 + (10 ** -15)):
                phi2 = 0.0
            else:
                phi2 = numpy.arccos(cos_phi2)
            if (1.0 - (10 ** -15)) <= cos_phi3 and cos_phi3 <= (1.0 + (10 ** -15)):
                phi3 = 0.0
            else:
                phi3 = numpy.arccos(cos_phi3)

            # Take the numerical average of the three out of plane angles
            # Note - other schemes for obtaining a single out of plane angle could be investigated
            phi = (phi1 + phi2 + phi3) / 3

            energy = energy + i.energy(phi)
        if verbosity >= 1:
            print("With inversion, energy = " + str(energy))
            # Don't forget to add non-bonded interactions here

        #    print("Omitting all non-covalent interactions") # REMOVE ONCE FIXED
        #    print("Omitting hydrogen bonding interactions") # REMOVE ONCE FIXED
        #    print("Omitting all non-covalent interactions except hydrogen bonding")
        #    print("Omitting all non-covalent interactions expcept halogen bonding")
        #    print("Omitting all non-covalent interactions except Pauli repulsion")         
        #    print("Omitting all non-covalent interactions except electrostatics")
        #    print("Omitting all non-covalent interactions except dispersion")

        e_hbnd = 0.0
        for i in self.hatoms:
            atH = [cartCoordinates[3 * i], cartCoordinates[3 * i + 1], cartCoordinates[3 * i + 2]]
            for j in self.highENatoms:
                atA = [cartCoordinates[3 * j], cartCoordinates[3 * j + 1], cartCoordinates[3 * j + 2]]
                # Calculate distance between hydrogen i and electronegative atom j
                # Determine whether they are (likely) joined by a covalent bond
                dist_HA = (atH[0] - atA[0]) ** 2
                dist_HA += (atH[1] - atA[1]) ** 2
                dist_HA += (atH[2] - atA[2]) ** 2
                dist_HA = math.sqrt(dist_HA)
                bond_dist_HA = SymbolToRadius[self.atoms[i].symbol] + SymbolToRadius[self.atoms[j].symbol]
                if dist_HA <= bond_dist_HA:
                    for k in self.highENatoms:
                        atB = [cartCoordinates[3 * j], cartCoordinates[3 * j + 1], cartCoordinates[3 * j + 2]]
                        # Calculate distance between hydrogen i and electronegative atom k
                        # Determine whether they are close enough for a hydrogen bonding interaction
                        dist_HB = (atH[0] - atB[0]) ** 2
                        dist_HB += (atH[1] - atB[1]) ** 2
                        dist_HB += (atH[2] - atB[2]) ** 2
                        dist_HB = math.sqrt(dist_HB)
                        Hbond_dist_HB = SymbolToVdWRadius[self.atoms[i].symbol] + SymbolToVdWRadius[
                            self.atoms[k].symbol]
                        # If so, and if A and B are distinct, take the triple AHB to be involved in hydrogen bonding and use to calculate energy
                        if dist_HB <= Hbond_dist_HB and atA != atB:
                            dist_AB = (atA[0] - atB[0]) ** 2
                            dist_AB += (atA[1] - atB[1]) ** 2
                            dist_AB += (atA[2] - atB[2]) ** 2
                            dist_AB = math.sqrt(dist_AB)

                            # Use the calculated distances and the cosine rule to calculate AHB andgle theta
                            numerator = dist_HA ** 2 + dist_HB ** 2 - dist_AB ** 2
                            denominator = 2 * dist_HA * dist_HB
                            argument = numerator / denominator
                            theta = numpy.arccos(argument)

                            # Calculate the relevant constants for a hydrogen bonding potential
                            c_hbnd_AB = HBondStrengthFactor(j.symbol, j.charge, dist_HA, k.symbol, k.charge, dist_HB)
                            f_dmp_theta = AngleDamping(theta)
                            f_dmp_hbnd = HBondDamping(j.symbol, k.symbol, dist_AB)

                            # Calculate the energy of this interaction, and add to the total hydrogen bonding contribution
                            e_hbnd = e_hbnd + potHBond(f_dmp_theta, f_dmp_hbnd, c_hbnd_AB, dist_AB)
        # Calculate total hydrogen bonding contribution from the sum over AHB triples, and add to energy
        e_hbnd = -1 * e_hbnd
        energy = energy + e_hbnd
        if verbosity >= 1:
            print("With hydrogen bonding, energy = " + str(energy))
        #    print("Omitting halogen bonding interactions")

        e_xbnd = 0.0
        for i in self.halogens:
            atX = [cartCoordinates[3 * i], cartCoordinates[3 * i + 1], cartCoordinates[3 * i + 2]]
            # Locate bonding partner(s), Y, for X by searching through the list of bonds in the molecule
            bondsX = []
            for j in range(len(self.bonds)):
                if self.bonds[j][0] == i:
                    bondsX.append([self.bonds[j][1], self.bonds[j][0]])
                elif self.bonds[j][1] == i:
                    bondsX.append(self.bonds[j])
            # Locate donor atoms D that could participate in halogen bonding (check for not being bonded to X implemented later)
            # Sum of van der Waals radii currently employed as a check for this
            donors = []
            for k in range(len(self.atoms)):
                sym = self.atoms[k].symbol
                donorsyms = ["N", "O", "F", "P", "S", "Cl", "As", "Se", "Br", "Sb", "Te", "I", "Bi", "Bi", "Po", "At",
                             "Uup", "Lv", "Uus"]
                if sym in donorsyms:
                    atD = [cartCoordinates[3 * k], cartCoordinates[3 * k + 1], cartCoordinates[3 * k + 2]]
                    dist_XD = (atX[0] - atD[0]) ** 2
                    dist_XD += (atX[1] - atD[1]) ** 2
                    dist_XD += (atX[2] - atD[2]) ** 2
                    dist_XD = math.sqrt(dist_XD)
                    Xbond_dist_XD = SymbolToVdWRadius[self.atoms[k].symbol] + SymbolToVdWRadius[self.atoms[i].symbol]
                    if dist_XD <= Xbond_dist_XD:
                        donors.append([k, dist_XD])
            # For each triple formed by a bond YX and donor D, calculate halogen bonding potential
            for bond in bondsX:
                for donor in donors:
                    if donor[0] != bond[0] and donor[0] != bond[1]:
                        symY = bond[0].symbol
                        coordY = [cartCoordinates[3 * bond[0]], cartCoordinates[3 * bond[0] + 1],
                                  cartCoordinates[3 * bond[0] + 2]]
                        symX = self.atoms[i].symbol
                        coordX = atX
                        symD = donor[0].symbol
                        coordD = [cartCoordinates[3 * donor[0]], cartCoordinates[3 * donor[0] + 1],
                                  cartCoordinates[3 * donor[0] + 2]]
                        r_XD = self.atoms[donor[1]].symbol

                        # Calculate the DXY angle
                        r_XY = (coordY[0] - coordX[0]) ** 2
                        r_XY += (coordY[1] - coordX[1]) ** 2
                        r_XY += (coordY[2] - coordX[2]) ** 2
                        r_XY = math.sqrt(r_XY)

                        r_DY = (coordY[0] - coordD[0]) ** 2
                        r_DY += (coordY[1] - coordD[1]) ** 2
                        r_DY += (coordY[2] - coordD[2]) ** 2
                        r_DY = math.sqrt(r_DY)

                        numerator = r_XD ** 2 + r_XY ** 2 - r_DY ** 2
                        denominator = 2 * r_XD * r_XY
                        argument = numerator / denominator
                        theta = numpy.arccos(argument)

                        # Calculate the relevant constants for a halogen bonding potential
                        f_dmp_theta = AngleDamping(theta)
                        f_dmp_xbnd = HBondDamping(symX, symD, r_XD)
                        c_xbnd = AtomicXBondFactor(symX, self.atoms[i].charge)

                        # Calculate the energy of this halogen bonding interaction and add to the total
                        e_xbnd = e_xbnd + potXBond(f_dmp_theta, f_dmp_xbnd, c_xbnd_x, r_XD)
                        # Calculate total halogen bonding contribution from the sum over DXY triples, and add to energy of the molecule
        e_xbnd = -1 * e_xbnd
        energy = energy + e_xbnd
        if verbosity >= 1:
            print("With halogen bonding, energy = " + str(energy))
        #    print("Omitting all Pauli repulsion interactions")

        e_Pauli = 0.0
        for i in range(len(self.atoms)):
            for j in range(i + 1, len(self.atoms)):
                symA = self.atoms[i].symbol
                symB = self.atoms[j].symbol
                # Calculate the distance between atoms i and j
                coordA = [cartCoordinates[3 * i], cartCoordinates[3 * i + 1], cartCoordinates[3 * i + 2]]
                coordB = [cartCoordinates[3 * j], cartCoordinates[3 * j + 1], cartCoordinates[3 * j + 2]]
                distance = (coordA[0] - coordB[0]) ** 2
                distance += (coordA[1] - coordB[1]) ** 2
                distance += (coordA[2] - coordB[2]) ** 2
                distance = math.sqrt(distance)
                # Calculate the required screening parameter, then the energy for this pair, and add to the total
                # Note that proper calculation will require the D3 cutoff radii R_0D3 which are yet to be worked in
                rep_disp_AB = self.screen_RepDisp(i, j)
                C6_A = C6[symA]
                C6_B = C6[symB]
                C6_AB = (C6_A + C6_B) / 2
                C8_AB = C6_AB  # Correct value to be implemented later, set equal to C6_AB solely to test
                energy_AB = potPauliRep(rep_disp_AB, symA, symB, distance, C6_AB, C8_AB)
                e_Pauli = e_Pauli + energy_AB
        energy = energy + e_Pauli
        if verbosity >= 1:
            print("With Pauli repulsion, energy = " + str(energy))
        #    print("Omitting all electrostatic interactions")

        e_ES = 0.0
        for i in range(len(self.atoms)):
            for j in range(i + 1, len(self.atoms)):
                # Note QM computed atomic charges at equilibrium structure should be used as per QMDFF
                # Currently the charge in class atom is used, which comes from atomic number
                chgA = self.atoms[i].QMcharge
                chgB = self.atoms[j].QMcharge

                # Calculate the distance between atoms i and j
                coordA = [cartCoordinates[3 * i], cartCoordinates[3 * i + 1], cartCoordinates[3 * i + 2]]
                coordB = [cartCoordinates[3 * j], cartCoordinates[3 * j + 1], cartCoordinates[3 * j + 2]]
                distance = (coordA[0] - coordB[0]) ** 2
                distance += (coordA[1] - coordB[1]) ** 2
                distance += (coordA[2] - coordB[2]) ** 2
                distance = math.sqrt(distance)

                # Calculate the required screening parameter, and the energy for this paiwise interaction, then add to the total
                elstat_AB = self.screen_ES(i, j)
                energy_AB = potElectrostatic(elstat_AB, chgA, chgB, distance)
                e_ES = e_ES + energy_AB
            #        print("Adding ES energy for atoms " + str([i, j]) + " with charges " + str([chgA, chgB]) + ", distance " + str(distance) + ", screening parameter " + str(elstat_AB) + ", giving energy = " + str(energy_AB))
        energy = energy + e_ES
        if verbosity >= 1:
            print("With electrostatic interactions, energy = " + str(energy))
        #    print("Omitting all dispersion interactions")

        e_disp = 0.0
        for i in range(len(self.atoms)):
            for j in range(len(self.atoms)):
                symA = self.atoms[i].symbol
                symB = self.atoms[j].symbol
                # Calculate the distance between atoms i and j
                coordA = [cartCoordinates[3 * i], cartCoordinates[3 * i + 1], cartCoordinates[3 * i + 2]]
                coordB = [cartCoordinates[3 * j], cartCoordinates[3 * j + 1], cartCoordinates[3 * j + 2]]
                distance = (coordA[0] - coordB[0]) ** 2
                distance += (coordA[1] - coordB[1]) ** 2
                distance += (coordA[2] - coordB[2]) ** 2
                distance = math.sqrt(distance)
                # Calculate the required paramenters and thence the energy for this pairwise interaction, then add to the total
                # Note that this is incomplete until the D3 cutoff radii R_0D3, as well as the coefficients C6_AB and C8_AB, are incorporated properly
                rep_disp_AB = self.screen_RepDisp(i, j)
                C6_A = C6[symA]
                C6_B = C6[symB]
                C6_AB = (C6_A + C6_B) / 2
                C8_AB = C6_AB  # To be completed - temporarily set equal to C6 for test run only
                BJdamp_AB = BJdamping(i, j, symA, symB,
                                      dtyp)  # Calculation of cutoff radii is incorporated in this function
                if dtyp == 1:
                    R0_AB = VdWCutoffRadius(symA, symB)
                    energy_AB = potCSODisp(C6_AB, distance, R0_AB)
                elif dtyp == 2:
                    energy_AB = potLondonDisp(rep_disp_AB, C6_AB, C8_AB, BJdamp_AB, distance)
                e_disp = e_disp + energy_AB
        energy = energy + e_disp
        if verbosity >= 1:
            print("With London dispersion, energy = " + str(energy))

        # Calculation of polarisation energy (for solute-solvent) to go here in future
        # Left out for version 1 as optional, only important as intermolecular interactions

        if verbosity >= 1:
            print("Total energy: " + str(energy))
        return energy

    def kdepFFEnergy(self, cartCoordinates, ForceConstants, verbosity=0, dtyp=1):
        """ (Molecule) -> number (Force Field energy)

      Returns a number containing the molecular energy according to the current Force Field definition at a fixed structure specified by cartCoordinates
      and using the stretch, bend and inversion force constants specified by ForceConstants
      The contribution from non-covalent interactions, and the torsional force constants, are fixed
      The dispersion potential
      The dispersion correction used is specified by dtyp, with 1 for C6-only calculating cutoff radius from van der Waals radii, 2 for full D3 using C6 and C8 coefficients
    """
        # Note the function fed to the optimiser will need to have only force constants as variables, so must fix cartesian coordinates somehow for the molecule.
        energy = 0.0
        #    print("Initial force constants for k-dependent energy calculation:") # REMOVE ONCE FIXED
        #    print(ForceConstants) # REMOVE ONCE FIXED
        #    print("Input coordinates for k-dependent energy calculation:") # REMOVE ONCE FIXED
        #    print(cartCoordinates) #REMOVE ONCE FIXED
        if verbosity >= 1:
            print("Initial energy for calculation = " + str(energy))
        energy = energy + self.Ee_QM
        if verbosity >= 1:
            print("With QM energy of equilibrium structure, energy = " + str(energy))

        #    if verbosity >= 1: # REMOVE ONCE FIXED
        #      print("Omitting stretching energy") # REMOVE ONCE FIXED

        for j in range(len(self.stretch)):
            i = self.stretch[j]
            k_str0 = i.k_str  # Store the value of k_str originally associated with this stretching potential
            i.setk(
                ForceConstants[j])  # Set k_str equal to the value specified in the force constant list for this stretch
            distance = (cartCoordinates[3 * i.atom1] - cartCoordinates[3 * i.atom2]) ** 2
            distance += (cartCoordinates[3 * i.atom1 + 1] - cartCoordinates[3 * i.atom2 + 1]) ** 2
            distance += (cartCoordinates[3 * i.atom1 + 2] - cartCoordinates[3 * i.atom2 + 2]) ** 2
            distance = math.sqrt(distance)
            #      print("Adding energy for bond " + str([i.atom1, i.atom2]) + ", distance = " + str(distance) + " energy = " + str(i.energy(distance)) + " to total")
            energy = energy + i.energy(distance)
            i.setk(
                k_str0)  # Restore the original value of k_str so this stretching potential is not permanently modified by the calculation
        if verbosity >= 1:
            print("With bond stretches, energy = " + str(energy))

        for j in range(len(self.str13)):
            i = self.str13[j]
            k_str0 = i.k_str  # Store the value of k_str originally associated with this 1,3-stretching potential
            i.setk(ForceConstants[len(
                self.stretch) + j])  # Set k_str equal to the value specified in the force constant list for this 1,3-stretch
            distance = (cartCoordinates[3 * i.atom1] - cartCoordinates[3 * i.atom2]) ** 2
            distance += (cartCoordinates[3 * i.atom1 + 1] - cartCoordinates[3 * i.atom2 + 1]) ** 2
            distance += (cartCoordinates[3 * i.atom1 + 2] - cartCoordinates[3 * i.atom2 + 2]) ** 2
            distance = math.sqrt(distance)
            energy = energy + i.energy(distance)
            i.setk(
                k_str0)  # Restore the original value of k_str so this 1,3-stretching potential is not permanently modified by the energy calculation
        if verbosity >= 1:
            print("With 1,3-stretches, energy = " + str(energy))

        for j in range(len(self.bend)):
            i = self.bend[j]
            k_bnd0 = i.k_bnd  # Store the value of k_bnd originally associated with this bending potential
            i.setk(ForceConstants[len(self.stretch) + len(
                self.str13) + j])  # Se k_bnd equal to the value specified in the force constant list for this bend 
            d_bond_1 = (cartCoordinates[3 * i.atom1] - cartCoordinates[3 * i.atom2]) ** 2
            d_bond_1 += (cartCoordinates[3 * i.atom1 + 1] - cartCoordinates[3 * i.atom2 + 1]) ** 2
            d_bond_1 += (cartCoordinates[3 * i.atom1 + 2] - cartCoordinates[3 * i.atom2 + 2]) ** 2
            d_bond_1 = math.sqrt(d_bond_1)

            d_bond_2 = (cartCoordinates[3 * i.atom2] - cartCoordinates[3 * i.atom3]) ** 2
            d_bond_2 += (cartCoordinates[3 * i.atom2 + 1] - cartCoordinates[3 * i.atom3 + 1]) ** 2
            d_bond_2 += (cartCoordinates[3 * i.atom2 + 2] - cartCoordinates[3 * i.atom3 + 2]) ** 2
            d_bond_2 = math.sqrt(d_bond_2)

            d_non_bond = (cartCoordinates[3 * i.atom1] - cartCoordinates[3 * i.atom3]) ** 2
            d_non_bond += (cartCoordinates[3 * i.atom1 + 1] - cartCoordinates[3 * i.atom3 + 1]) ** 2
            d_non_bond += (cartCoordinates[3 * i.atom1 + 2] - cartCoordinates[3 * i.atom3 + 2]) ** 2
            d_non_bond = math.sqrt(d_non_bond)

            # Use those distances and the cosine rule to calculate bond angle theta
            numerator = d_bond_1 ** 2 + d_bond_2 ** 2 - d_non_bond ** 2
            denominator = 2 * d_bond_1 * d_bond_2
            argument = numerator / denominator
            theta = numpy.arccos(argument)
            energy = energy + i.energy(theta)
            i.setk(
                k_bnd0)  # Restore the original value of the bending force constant so that this potentially is not permanently changed by the energy calculation
        if verbosity >= 1:
            print("With bends, energy = " + str(energy))

        for i in self.tors:
            # Calculate the vectors lying along bonds, and their cross products
            atom_e1 = [cartCoordinates[3 * i.atom1], cartCoordinates[3 * i.atom1 + 1], cartCoordinates[3 * i.atom1 + 2]]
            atom_b1 = [cartCoordinates[3 * i.atom2], cartCoordinates[3 * i.atom2 + 1], cartCoordinates[3 * i.atom2 + 2]]
            atom_b2 = [cartCoordinates[3 * i.atom3], cartCoordinates[3 * i.atom3 + 1], cartCoordinates[3 * i.atom3 + 2]]
            atom_e2 = [cartCoordinates[3 * i.atom4], cartCoordinates[3 * i.atom4 + 1], cartCoordinates[3 * i.atom4 + 2]]
            end_1 = [atom_e1[i] - atom_b1[i] for i in range(3)]
            bridge = [atom_b1[i] - atom_b2[i] for i in range(3)]
            end_2 = [atom_b2[i] - atom_e2[i] for i in range(3)]
            vnormal_1 = numpy.cross(end_1, bridge)
            vnormal_2 = numpy.cross(bridge, end_2)

            # Construct a set of orthogonal basis vectors to define a frame with vnormal_2 as the x axis
            vcross = numpy.cross(vnormal_2, bridge)
            norm_vn2 = numpy.linalg.norm(vnormal_2)
            norm_b = numpy.linalg.norm(bridge)
            norm_vc = numpy.linalg.norm(vcross)
            basis_vn2 = [vnormal_2[i] / norm_vn2 for i in range(3)]
            basis_b = [bridge[i] / norm_b for i in range(3)]
            basis_cv = [vcross[i] / norm_vc for i in range(3)]

            # Find the signed angle between vnormal_1 and vnormal_2 in the new frame
            vn1_coord_n2 = numpy.dot(vnormal_1, basis_vn2)
            vn1_coord_vc = numpy.dot(vnormal_1, basis_cv)
            psi = math.atan2(vn1_coord_vc, vn1_coord_n2)
            energy = energy + i.energy(psi)
        if verbosity >= 1:
            print("With torsion, energy = " + str(energy))

        for j in range(len(self.inv)):
            i = self.inv[j]
            k_inv0 = i.k_inv  # Store the force constant k originally associated with this inversion potential 
            i.setk(ForceConstants[len(self.stretch) + len(self.str13) + len(
                self.bend) + j])  # Set k_inv equal to the value specified for this inversion potential in the force constants list
            # Calculate the vectors along bonds, and construct a vector plane_norm orthogonal to the plane of end ato
            atom_c = [cartCoordinates[3 * i.atom1], cartCoordinates[3 * i.atom1 + 1], cartCoordinates[3 * i.atom1 + 2]]
            atom_e1 = [cartCoordinates[3 * i.atom2], cartCoordinates[3 * i.atom2 + 1], cartCoordinates[3 * i.atom2 + 2]]
            atom_e2 = [cartCoordinates[3 * i.atom3], cartCoordinates[3 * i.atom3 + 1], cartCoordinates[3 * i.atom3 + 2]]
            atom_e3 = [cartCoordinates[3 * i.atom4], cartCoordinates[3 * i.atom4 + 1], cartCoordinates[3 * i.atom4 + 2]]
            bond_1 = [atom_e1[i] - atom_c[i] for i in range(3)]
            bond_2 = [atom_e2[i] - atom_c[i] for i in range(3)]
            bond_3 = [atom_e3[i] - atom_c[i] for i in range(3)]
            inplane_12 = [atom_e2[i] - atom_e1[i] for i in range(3)]
            inplane_13 = [atom_e3[i] - atom_e1[i] for i in range(3)]
            plane_norm = numpy.cross(inplane_12, inplane_13)

            # Construct vectors between end atoms and the projection of the central atom on the end-atom pla
            cross_1 = numpy.cross(bond_1, plane_norm)
            cross_2 = numpy.cross(bond_2, plane_norm)
            cross_3 = numpy.cross(bond_2, plane_norm)
            inplane_1 = numpy.cross(plane_norm, cross_1) / numpy.dot(plane_norm, plane_norm)
            inplane_2 = numpy.cross(plane_norm, cross_2) / numpy.dot(plane_norm, plane_norm)
            inplane_3 = numpy.cross(plane_norm, cross_3) / numpy.dot(plane_norm, plane_norm)

            # Calculate the out of plane angle for each of the three bonds
            cos_phi1 = numpy.dot(bond_1, inplane_1) / (numpy.linalg.norm(bond_1) * numpy.linalg.norm(inplane_1))
            cos_phi2 = numpy.dot(bond_2, inplane_2) / (numpy.linalg.norm(bond_2) * numpy.linalg.norm(inplane_2))
            cos_phi3 = numpy.dot(bond_3, inplane_3) / (numpy.linalg.norm(bond_3) * numpy.linalg.norm(inplane_3))

            if (1.0 - (10 ** -15)) <= cos_phi1 and cos_phi1 <= (1.0 + (10 ** -15)):
                phi1 = 0.0
            else:
                phi1 = numpy.arccos(cos_phi1)
            if (1.0 - (10 ** -15)) <= cos_phi2 and cos_phi2 <= (1.0 + (10 ** -15)):
                phi2 = 0.0
            else:
                phi2 = numpy.arccos(cos_phi2)
            if (1.0 - (10 ** -15)) <= cos_phi3 and cos_phi3 <= (1.0 + (10 ** -15)):
                phi3 = 0.0
            else:
                phi3 = numpy.arccos(cos_phi3)

            # Take the numerical average of the three out of plane angles
            # Note - other schemes for obtaining a single out of plane angle could be investigated
            phi = (phi1 + phi2 + phi3) / 3

            energy = energy + i.energy(phi)
            i.setk(
                k_inv0)  # Restore the original value of k_inv so that this inversion potential is not permanently modified by the energy caculation
        if verbosity >= 1:
            print("With inversion, energy = " + str(energy))

        #    if verbosity >= 1: # REMOVE ONCE FIXED
        #      print("Omitting all non-covalent interactions") # REMOVE ONCE FIXED

        #    if verbosity >=1: # REMOVE ONCE FIXED
        #      print("Total energy:") # REMOVE ONCE FIXED
        #    return energy # REMOVE ONCE FIXED

        #    print("Omitting hydrogen bond interactions") # REMOVE ONCE FIXED
        #    print("Omitting all non-covalent interactions except hydrogen bonding")
        #    print("Omitting all non-covalent interactions expcept halogen bonding")   
        #    print("Omitting all non-covalent interactions except Pauli repulsion")   
        #    print("Omitting all non-covalent interactions except electrostatics")
        #    print("Omitting all non-covalent interactions except dispersion")

        e_hbnd = 0.0
        for i in self.hatoms:
            atH = [cartCoordinates[3 * i], cartCoordinates[3 * i + 1], cartCoordinates[3 * i + 2]]
            for j in self.highENatoms:
                atA = [cartCoordinates[3 * j], cartCoordinates[3 * j + 1], cartCoordinates[3 * j + 2]]
                # Calculate distance between hydrogen i and electronegative atom j
                # Determine whether they are (likely) joined by a covalent bond
                dist_HA = (atH[0] - atA[0]) ** 2
                dist_HA += (atH[1] - atA[1]) ** 2
                dist_HA += (atH[2] - atA[2]) ** 2
                dist_HA = math.sqrt(dist_HA)
                bond_dist_HA = SymbolToRadius[self.atoms[i].symbol] + SymbolToRadius[self.atoms[j].symbol]
                if dist_HA <= bond_dist_HA:
                    for k in self.highENatoms:
                        atB = [cartCoordinates[3 * j], cartCoordinates[3 * j + 1], cartCoordinates[3 * j + 2]]
                        # Calculate distance between hydrogen i and electronegative atom k
                        # Determine whether they are close enough for a hydrogen bonding interaction
                        dist_HB = (atH[0] - atB[0]) ** 2
                        dist_HB += (atH[1] - atB[1]) ** 2
                        dist_HB += (atH[2] - atB[2]) ** 2
                        dist_HB = math.sqrt(dist_HB)
                        Hbond_dist_HB = SymbolToVdWRadius[self.atoms[i].symbol] + SymbolToVdWRadius[
                            self.atoms[k].symbol]
                        # If so, and if A and B are distinct, take the triple AHB to be involved in hydrogen bonding and use to calculate energy
                        if dist_HB <= Hbond_dist_HB and atA != atB:
                            dist_AB = (atA[0] - atB[0]) ** 2
                            dist_AB += (atA[1] - atB[1]) ** 2
                            dist_AB += (atA[2] - atB[2]) ** 2
                            dist_AB = math.sqrt(dist_AB)

                            # Use the calculated distances and the cosine rule to calculate AHB andgle theta
                            numerator = dist_HA ** 2 + dist_HB ** 2 - dist_AB ** 2
                            denominator = 2 * dist_HA * dist_HB
                            argument = numerator / denominator
                            theta = numpy.arccos(argument)

                            # Calculate the relevant constants for a hydrogen bonding potential
                            c_hbnd_AB = HBondStrengthFactor(j.symbol, j.charge, dist_HA, k.symbol, k.charge, dist_HB)
                            f_dmp_theta = AngleDamping(theta)
                            f_dmp_hbnd = HBondDamping(j.symbol, k.symbol, dist_AB)

                            # Calculate the energy of this interaction, and add to the total hydrogen bonding contribution
                            e_hbnd = e_hbnd + potHBond(f_dmp_theta, f_dmp_hbnd, c_hbnd_AB, dist_AB)
        # Calculate total hydrogen bonding contribution from the sum over AHB triples, and add to energy
        e_hbnd = -1 * e_hbnd
        energy = energy + e_hbnd
        if verbosity >= 1:
            print("With hydrogen bonding, energy = " + str(energy))
        #    print("Omitting halogen bonding interactions")

        e_xbnd = 0.0
        for i in self.halogens:
            atX = [cartCoordinates[3 * i], cartCoordinates[3 * i + 1], cartCoordinates[3 * i + 2]]
            # Locate bonding partner(s), Y, for X by searching through the list of bonds in the molecule
            bondsX = []
            for j in range(len(self.bonds)):
                if self.bonds[j][0] == i:
                    bondsX.append([self.bonds[j][1], self.bonds[j][0]])
                elif self.bonds[j][1] == i:
                    bondsX.append(self.bonds[j])
            # Locate donor atoms D that could participate in halogen bonding (check for not being bonded to X implemented later)
            # Sum of van der Waals radii currently employed as a check for this
            donors = []
            for k in range(len(self.atoms)):
                sym = self.atoms[k].symbol
                donorsyms = ["N", "O", "F", "P", "S", "Cl", "As", "Se", "Br", "Sb", "Te", "I", "Bi", "Bi", "Po", "At",
                             "Uup", "Lv", "Uus"]
                if sym in donorsyms:
                    atD = [cartCoordinates[3 * k], cartCoordinates[3 * k + 1], cartCoordinates[3 * k + 2]]
                    dist_XD = (atX[0] - atD[0]) ** 2
                    dist_XD += (atX[1] - atD[1]) ** 2
                    dist_XD += (atX[2] - atD[2]) ** 2
                    dist_XD = math.sqrt(dist_XD)
                    Xbond_dist_XD = SymbolToVdWRadius[self.atoms[k].symbol] + SymbolToVdWRadius[self.atoms[i].symbol]
                    if dist_XD <= Xbond_dist_XD:
                        donors.append([k, dist_XD])
            # For each triple formed by a bond YX and donor D, calculate halogen bonding potential
            for bond in bondsX:
                for donor in donors:
                    if donor[0] != bond[0] and donor[0] != bond[1]:
                        symY = bond[0].symbol
                        coordY = [cartCoordinates[3 * bond[0]], cartCoordinates[3 * bond[0] + 1],
                                  cartCoordinates[3 * bond[0] + 2]]
                        symX = self.atoms[i].symbol
                        coordX = atX
                        symD = donor[0].symbol
                        coordD = [cartCoordinates[3 * donor[0]], cartCoordinates[3 * donor[0] + 1],
                                  cartCoordinates[3 * donor[0] + 2]]
                        r_XD = self.atoms[donor[1]].symbol

                        # Calculate the DXY angle
                        r_XY = (coordY[0] - coordX[0]) ** 2
                        r_XY += (coordY[1] - coordX[1]) ** 2
                        r_XY += (coordY[2] - coordX[2]) ** 2
                        r_XY = math.sqrt(r_XY)

                        r_DY = (coordY[0] - coordD[0]) ** 2
                        r_DY += (coordY[1] - coordD[1]) ** 2
                        r_DY += (coordY[2] - coordD[2]) ** 2
                        r_DY = math.sqrt(r_DY)

                        numerator = r_XD ** 2 + r_XY ** 2 - r_DY ** 2
                        denominator = 2 * r_XD * r_XY
                        argument = numerator / denominator
                        theta = numpy.arccos(argument)

                        # Calculate the relevant constants for a halogen bonding potential
                        f_dmp_theta = AngleDamping(theta)
                        f_dmp_xbnd = HBondDamping(symX, symD, r_XD)
                        c_xbnd = AtomicXBondFactor(symX, self.atoms[i].charge)

                        # Calculate the energy of this halogen bonding interaction and add to the total
                        e_xbnd = e_xbnd + potXBond(f_dmp_theta, f_dmp_xbnd, c_xbnd_x, r_XD)
                        # Calculate total halogen bonding contribution from the sum over DXY triples, and add to energy of the molecule
        e_xbnd = -1 * e_xbnd
        energy = energy + e_xbnd
        if verbosity >= 1:
            print("With halogen bonding, energy = " + str(energy))
        #    print("Omitting all Pauli repulsion interactions")

        e_Pauli = 0.0
        for i in range(len(self.atoms)):
            for j in range(i + 1, len(self.atoms)):
                symA = self.atoms[i].symbol
                symB = self.atoms[j].symbol
                # Calculate the distance between atoms i and j
                coordA = [cartCoordinates[3 * i], cartCoordinates[3 * i + 1], cartCoordinates[3 * i + 2]]
                coordB = [cartCoordinates[3 * j], cartCoordinates[3 * j + 1], cartCoordinates[3 * j + 2]]
                distance = (coordA[0] - coordB[0]) ** 2
                distance += (coordA[1] - coordB[1]) ** 2
                distance += (coordA[2] - coordB[2]) ** 2
                distance = math.sqrt(distance)
                # Calculate the required screening parameter, then the energy for this pair, and add to the total
                # Note that proper calculation will require the D3 cutoff radii R_0D3 which are yet to be worked in
                rep_disp_AB = self.screen_RepDisp(i, j)
                C6_A = C6[symA]
                C6_B = C6[symB]
                C6_AB = (C6_A + C6_B) / 2
                C8_AB = C6_AB  # Correct value to be implemented later, set equal to C6_AB solely to test
                energy_AB = potPauliRep(rep_disp_AB, symA, symB, distance, C6_AB, C8_AB)
                e_Pauli = e_Pauli + energy_AB
        energy = energy + e_Pauli
        if verbosity >= 1:
            print("With Pauli repulsion, energy = " + str(energy))
        #    print("Omitting all electrostatic interactions")

        e_ES = 0.0
        for i in range(len(self.atoms)):
            for j in range(i + 1, len(self.atoms)):
                # Note QM computed atomic charges at equilibrium structure should be used as per QMDFF
                # Currently the charge in class atom is used, which comes from atomic number
                chgA = self.atoms[i].QMcharge
                chgB = self.atoms[j].QMcharge

                # Calculate the distance between atoms i and j
                coordA = [cartCoordinates[3 * i], cartCoordinates[3 * i + 1], cartCoordinates[3 * i + 2]]
                coordB = [cartCoordinates[3 * j], cartCoordinates[3 * j + 1], cartCoordinates[3 * j + 2]]
                distance = (coordA[0] - coordB[0]) ** 2
                distance += (coordA[1] - coordB[1]) ** 2
                distance += (coordA[2] - coordB[2]) ** 2
                distance = math.sqrt(distance)

                # Calculate the required screening parameter, and the energy for this paiwise interaction, then add to the total
                elstat_AB = self.screen_ES(i, j)
                energy_AB = potElectrostatic(elstat_AB, chgA, chgB, distance)
                e_ES = e_ES + energy_AB
            #        print("Adding ES energy for atoms " + str([i, j]) + " with charges " + str([chgA, chgB]) + ", distance " + str(distance) + ", screening parameter " + str(elstat_AB) + ", giving energy = " + str(energy_AB))
        energy = energy + e_ES
        if verbosity >= 1:
            print("With electrostatic interactions, energy = " + str(energy))
        #    print("Omitting all dispersion interactions")

        e_disp = 0.0
        for i in range(len(self.atoms)):
            for j in range(len(self.atoms)):
                symA = self.atoms[i].symbol
                symB = self.atoms[j].symbol
                # Calculate the distance between atoms i and j
                coordA = [cartCoordinates[3 * i], cartCoordinates[3 * i + 1], cartCoordinates[3 * i + 2]]
                coordB = [cartCoordinates[3 * j], cartCoordinates[3 * j + 1], cartCoordinates[3 * j + 2]]
                distance = (coordA[0] - coordB[0]) ** 2
                distance += (coordA[1] - coordB[1]) ** 2
                distance += (coordA[2] - coordB[2]) ** 2
                distance = math.sqrt(distance)
                # Calculate the required paramenters and thence the energy for this pairwise interaction, then add to the total
                # Note that this is incomplete until the D3 cutoff radii R_0D3, as well as the coefficients C6_AB and C8_AB, are incorporated properly
                rep_disp_AB = self.screen_RepDisp(i, j)
                C6_A = C6[symA]
                C6_B = C6[symB]
                C6_AB = (C6_A + C6_B) / 2
                C8_AB = C6_AB  # To be completed - temporarily set equal to C6 for test run only
                BJdamp_AB = BJdamping(i, j, symA, symB, dtyp)  # Calculation of cutoff radii incorporated here 
                if dtyp == 1:
                    R0_AB = VdWCutoffRadius(symA, symB)
                    energy_AB = potCSODisp(C6_AB, distance, R0_AB)
                elif dtyp == 2:
                    energy_AB = potLondonDisp(rep_disp_AB, C6_AB, C8_AB, BJdamp_AB, distance)
                e_disp = e_disp + energy_AB
        energy = energy + e_disp
        if verbosity >= 1:
            print("With London dispersion, energy = " + str(energy))

        # Calculation of polarisation energy (for solute-solvent) to go here in future
        # Left out for version 1 as optional, only important as intermolecular interactions

        if verbosity >= 1:
            print("Total energy:")
        return (energy)

    def kdepHessian(self, ForceConstants):
        """ (Molecule) -> 3N x 3N matrix

    Returns the Hessian matrix for the molecule as calculated numerically from the Force field energy with the specified list of force constants
    """
        # For greater flexibility in usage, a list of Cartesian Coordinates could also be given as an argument if dependence upon force constants alone were not desired
        # Take the original Cartesian coordinates of the molecule as initial geometry
        coords = self.cartesianCoordinates()
        #    print("Initial coordinates for FF Hessian calculation:") # REMOVE ONCE FIXED
        #    print(coords) # REMOVE ONCE FIXED
        epsilon = 1 * (10 ** -5)
        #    print("For FF Hessian calculation, epsilon = " + str(epsilon)) # REMOVE ONCE FIXED
        # Use the finite difference approximation to calculate first derivatives at the initial geometry
        #    print("Calculating approximate first derivatives:") # REMOVE ONCE FIXED
        deriv1 = scipy.optimize.approx_fprime(coords, self.kdepFFEnergy, epsilon,
                                              ForceConstants)  # Note verbosity option not passed as an argument, so cannot be used from kdepFFEnergy at present except by modifying default values
        # Check whether the syntax for additional arguments in approx_fprime is correct here or whether they should be in a list
        # Also whether a approx_fprime works as well with a class method as with an independently defined function
        # If not, may need to write out in full
        # Set up a zero matrix to become the Force Field Hessian
        n = len(coords)
        H_FF = numpy.zeros((n, n))
        # For each coordinate, displace by epsilon and calculate the second derivatives using another finite difference approximation
        for i in range(n):
            x0 = coords[i]
            coords[i] = x0 + epsilon
            deriv2 = scipy.optimize.approx_fprime(coords, self.kdepFFEnergy, epsilon,
                                                  ForceConstants)  # Same comments on verbosity and checks apply as above
            # Place the calculated second derivatives for coordinate i into the ith column of the Hessian matrix
            H_FF[:, i] = (deriv2 - deriv1) / epsilon
            coords[i] = x0
        return H_FF

    def HessianDiffSquared(self, ForceConstants):
        """
    Objective function to be minimised in the Hessian fit
    Gives squared deviation between QM Hessian H_QM and Force Field Hessian H_FF
    """
        # Take the QM calculated Hessian stored as an attribute of the molecule
        H_QM = self.H_QM
        #    print("QM Hessian used for Hessian difference:") # REMOVE ONCE FIXED
        #    print(H_QM) # REMOVE ONCE FIXED
        # Calculate the Force Field Hessian for the given force constants
        #    print("Calculated FF Hessian:")
        H_FF = self.kdepHessian(ForceConstants)  # Temporary print
        #    print(H_FF)

        sqdev = 0.0
        # Given H_QM and H_FF as arrays of equal size and shape, iterate over the individual entries of each
        for i in range(int(numpy.sqrt(H_QM.size))):
            for j in range(int(numpy.sqrt(H_QM.size))):
                # Take the difference between entries and square it
                diff = H_QM[i, j] - H_FF[i, j]
                diff = diff ** 2
                # Add this to the total squared deviation
                sqdev = sqdev + diff

        return sqdev


#############################################################################################################
# Most important function so far: Read Quantum Chemistry output file and construct WellFaRe Molecule from it
#############################################################################################################

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
            print("\nReading of Mulliken charges finished. \nAdding charges to atoms in WellFARe molecule: ",
                  molecule.name)
        for i in charges:
            readBuffer = i.split()
            n = int(readBuffer[0]) - 1
            molecule.atoms[n].setq(float(readBuffer[2]))
            if verbosity >= 2:
                print(molecule.atoms[n].__repr__())
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
                print(molecule.atoms[n].__repr__())
        f.close()

    # EQUILIBRIUM ENERGY READING SECTION
    QM_energies = []
    # Read through Gaussian file, read energies from SCF in A.U.
    if program == "g09":
        f = open(filename, 'r')
        for line in f:
            if line.find("SCF Done:") != -1:
                if verbosity >= 3:
                    print("\nEnergy from SCF cycle found")
                    print(str(line))
                QM_energies.append(line)
        # Take the last SCF energy from the file and assign that value as the QM equilibrium energy of the molecule
        i = QM_energies[len(QM_energies) - 1]
        readBuffer = i.split()
        QMenergy = float(readBuffer[4])
        molecule.setQMenergy(QMenergy)
        if verbosity >= 1:
            print("\nReading of QM equilibrium energy complete")
            if verbosity >= 2:
                print("Ee_QM = " + str(QMenergy))
        f.close()
    # Read through ORCA file and assign the final single point energy reported as the QM equilibrium energy of the molecule
    if program == "orca":
        f = open(filename, 'r')
        for line in f:
            if line.find("FINAL SINGLE POINT ENERGY") != -1:
                if verbosity >= 3:
                    print("\nSingle point energy found")
                    print(str(line))
                QM_energies.append(line)
        i = QM_energies[len(QM_energies) - 1]
        readBuffer = i.split()
        #    print(readBuffer) # REMOVE AFTER TESTING
        QMenergy = float(readBuffer[4])
        molecule.setQMenergy(QMenergy)
        if verbosity >= 1:
            print("\nReading of QM equilibrium energy complete")
            if verbosity >= 2:
                print("Ee_QM = " + str(QMenergy))
        f.close()

    # BOND ORDER READING SECTION
    bo = []
    bo = numpy.zeros((molecule.numatoms(), molecule.numatoms()))
    if program == "g09":
        f = open(filename, 'r')
        for line in f:
            if line.find("Atomic Valencies and Mayer Atomic Bond Orders:") != -1:
                if verbosity >= 2:
                    print("\nAtomic Valencies and Mayer Atomic Bond Orders found, reading data")
                bo = numpy.zeros((molecule.numatoms(), molecule.numatoms()))
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
            numpy.set_printoptions(suppress=True)
            numpy.set_printoptions(formatter={'float': '{: 0.3f}'.format})
            print(bo)
    if program == "orca":
        f = open(filename, 'r')
        for line in f:
            if line.find("Mayer bond orders larger than 0.1") != -1:
                if verbosity >= 2:
                    print("\nMayer bond orders larger than 0.1 found, reading data")
                bo = numpy.zeros((molecule.numatoms(), molecule.numatoms()))
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
            numpy.set_printoptions(suppress=True)
            numpy.set_printoptions(formatter={'float': '{: 0.3f}'.format})
            print(bo)

    # FORCE CONSTANT READING SECTION
    H = []
    H = numpy.zeros((3 * molecule.numatoms(), 3 * molecule.numatoms()))
    if program == "g09":
        f = open(filename, 'r')
        for line in f:
            if line.find("Force constants in Cartesian coordinates") != -1:
                if verbosity >= 2:
                    print("\nForce constants in Cartesian coordinates, reading data")
                H = numpy.zeros((3 * molecule.numatoms(), 3 * molecule.numatoms()))
                while True:
                    readBuffer = f.__next__()
                    # Check if the whole line is integers only (Header line)
                    if isInt("".join(readBuffer.split())) == True:
                        # And use this information to label the columns
                        columns = readBuffer.split()
                    # Once we find the FormGI statement, we're done reading
                    elif readBuffer.find("FormGI is forming") != -1 or readBuffer.find(
                            "Cartesian forces in FCRed") != -1:
                        break
                    else:
                        row = readBuffer.split()
                        for i in range(0, len(row) - 1):
                            H[int(row[0]) - 1][int(columns[i]) - 1] = row[i + 1].replace('D', 'E')
                            H[int(columns[i]) - 1][int(row[0]) - 1] = row[i + 1].replace('D', 'E')
        molecule.setHessian(H)  # Store H as the QM calculated Hessian for this molecule  
        if verbosity >= 3:
            print("\nForce constants in Cartesian coordinates (Input orientation):")
            # numpy.set_printoptions(suppress=True)
            # numpy.set_printoptions(formatter={'float': '{: 0.3f}'.format})
            print(H)
        f.close()

    # Test if we actually have Mayer Bond orders
    if numpy.count_nonzero(bo) != 0:
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
                    SymbolToRadius[molecule.atoms[i].symbol] + SymbolToRadius[molecule.atoms[j].symbol]) * distfactor:
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

    # Same for threefolds: Use angles to determine where they are
    if verbosity >= 2:
        print("\nAdding threefolds to WellFARe molecule: ", molecule.name)
    for i in range(0, len(molecule.angles)):
        for j in range(i + 1, len(molecule.angles)):
            for k in range(j + 1, len(molecule.angles)):
                if molecule.angles[i][1] == molecule.angles[j][1] == molecule.angles[k][1]:
                    if molecule.angles[i][0] == molecule.angles[j][0] and molecule.angles[j][2] == molecule.angles[k][
                        2] and molecule.angles[i][2] == molecule.angles[k][0]:
                        molecule.addThreefold(molecule.angles[i][1], molecule.angles[i][0], molecule.angles[j][2],
                                              molecule.angles[i][2])
                        if verbosity >= 2:
                            print(
                                " {:<3} ({:3d}), {:<3} ({:3d}), {:<3} ({:3d}) and {:<3} ({:3d}) ({: 7.2f} deg)".format(
                                    molecule.atoms[molecule.angles[i][1]].symbol, molecule.angles[i][1],
                                    molecule.atoms[molecule.angles[i][0]].symbol, molecule.angles[i][0],
                                    molecule.atoms[molecule.angles[j][2]].symbol, molecule.angles[j][2],
                                    molecule.atoms[molecule.angles[i][2]].symbol, molecule.angles[i][2],
                                    math.degrees(molecule.outofplaneangle(len(molecule.threefolds) - 1))))
                    if molecule.angles[i][0] == molecule.angles[j][0] and molecule.angles[j][2] == molecule.angles[k][
                        0] and molecule.angles[i][2] == molecule.angles[k][2]:
                        molecule.addThreefold(molecule.angles[i][1], molecule.angles[i][0], molecule.angles[j][2],
                                              molecule.angles[i][2])
                        if verbosity >= 2:
                            print(
                                " {:<3} ({:3d}), {:<3} ({:3d}), {:<3} ({:3d}) and {:<3} ({:3d}) ({: 7.2f} deg)".format(
                                    molecule.atoms[molecule.angles[i][1]].symbol, molecule.angles[i][1],
                                    molecule.atoms[molecule.angles[i][0]].symbol, molecule.angles[i][0],
                                    molecule.atoms[molecule.angles[j][2]].symbol, molecule.angles[j][2],
                                    molecule.atoms[molecule.angles[i][2]].symbol, molecule.angles[i][2],
                                    math.degrees(molecule.outofplaneangle(len(molecule.threefolds) - 1))))
                    if molecule.angles[i][0] == molecule.angles[j][2] and molecule.angles[j][0] == molecule.angles[k][
                        0] and molecule.angles[i][2] == molecule.angles[k][2]:
                        molecule.addThreefold(molecule.angles[i][1], molecule.angles[i][0], molecule.angles[j][0],
                                              molecule.angles[i][2])
                        if verbosity >= 2:
                            print(
                                " {:<3} ({:3d}), {:<3} ({:3d}), {:<3} ({:3d}) and {:<3} ({:3d}) ({: 7.2f} deg)".format(
                                    molecule.atoms[molecule.angles[i][1]].symbol, molecule.angles[i][1],
                                    molecule.atoms[molecule.angles[i][0]].symbol, molecule.angles[i][0],
                                    molecule.atoms[molecule.angles[j][0]].symbol, molecule.angles[j][0],
                                    molecule.atoms[molecule.angles[i][2]].symbol, molecule.angles[i][2],
                                    math.degrees(molecule.outofplaneangle(len(molecule.threefolds) - 1))))
                    if molecule.angles[i][0] == molecule.angles[j][2] and molecule.angles[j][0] == molecule.angles[k][
                        2] and molecule.angles[i][2] == molecule.angles[k][0]:
                        molecule.addThreefold(molecule.angles[i][1], molecule.angles[i][0], molecule.angles[j][0],
                                              molecule.angles[i][2])
                        if verbosity >= 2:
                            print(
                                " {:<3} ({:3d}), {:<3} ({:3d}), {:<3} ({:3d}) and {:<3} ({:3d}) ({: 7.2f} deg)".format(
                                    molecule.atoms[molecule.angles[i][1]].symbol, molecule.angles[i][1],
                                    molecule.atoms[molecule.angles[i][0]].symbol, molecule.angles[i][0],
                                    molecule.atoms[molecule.angles[j][0]].symbol, molecule.angles[j][0],
                                    molecule.atoms[molecule.angles[i][2]].symbol, molecule.angles[i][2],
                                    math.degrees(molecule.outofplaneangle(len(molecule.threefolds) - 1))))
                    if molecule.angles[i][2] == molecule.angles[j][0] and molecule.angles[j][2] == molecule.angles[k][
                        0] and molecule.angles[i][0] == molecule.angles[k][2]:
                        molecule.addThreefold(molecule.angles[i][1], molecule.angles[i][2], molecule.angles[j][2],
                                              molecule.angles[i][0])
                        if verbosity >= 2:
                            print(
                                " {:<3} ({:3d}), {:<3} ({:3d}), {:<3} ({:3d}) and {:<3} ({:3d}) ({: 7.2f} deg)".format(
                                    molecule.atoms[molecule.angles[i][1]].symbol, molecule.angles[i][1],
                                    molecule.atoms[molecule.angles[i][2]].symbol, molecule.angles[i][2],
                                    molecule.atoms[molecule.angles[j][2]].symbol, molecule.angles[j][2],
                                    molecule.atoms[molecule.angles[i][0]].symbol, molecule.angles[i][0],
                                    math.degrees(molecule.outofplaneangle(len(molecule.threefolds) - 1))))
                    if molecule.angles[i][2] == molecule.angles[j][0] and molecule.angles[j][2] == molecule.angles[k][
                        2] and molecule.angles[i][0] == molecule.angles[k][0]:
                        molecule.addThreefold(molecule.angles[i][1], molecule.angles[i][2], molecule.angles[j][2],
                                              molecule.angles[i][0])
                        if verbosity >= 2:
                            print(
                                " {:<3} ({:3d}), {:<3} ({:3d}), {:<3} ({:3d}) and {:<3} ({:3d}) ({: 7.2f} deg)".format(
                                    molecule.atoms[molecule.angles[i][1]].symbol, molecule.angles[i][1],
                                    molecule.atoms[molecule.angles[i][2]].symbol, molecule.angles[i][2],
                                    molecule.atoms[molecule.angles[j][2]].symbol, molecule.angles[j][2],
                                    molecule.atoms[molecule.angles[i][0]].symbol, molecule.angles[i][0],
                                    math.degrees(molecule.outofplaneangle(len(molecule.threefolds) - 1))))
                    if molecule.angles[i][2] == molecule.angles[j][2] and molecule.angles[j][0] == molecule.angles[k][
                        0] and molecule.angles[i][0] == molecule.angles[k][2]:
                        molecule.addThreefold(molecule.angles[i][1], molecule.angles[i][2], molecule.angles[j][0],
                                              molecule.angles[i][0])
                        if verbosity >= 2:
                            print(
                                " {:<3} ({:3d}), {:<3} ({:3d}), {:<3} ({:3d}) and {:<3} ({:3d}) ({: 7.2f} deg)".format(
                                    molecule.atoms[molecule.angles[i][1]].symbol, molecule.angles[i][1],
                                    molecule.atoms[molecule.angles[i][2]].symbol, molecule.angles[i][2],
                                    molecule.atoms[molecule.angles[j][0]].symbol, molecule.angles[j][0],
                                    molecule.atoms[molecule.angles[i][0]].symbol, molecule.angles[i][0],
                                    math.degrees(molecule.outofplaneangle(len(molecule.threefolds) - 1))))
                    if molecule.angles[i][2] == molecule.angles[j][2] and molecule.angles[j][0] == molecule.angles[k][
                        2] and molecule.angles[i][0] == molecule.angles[k][0]:
                        molecule.addThreefold(molecule.angles[i][1], molecule.angles[i][2], molecule.angles[j][0],
                                              molecule.angles[i][0])
                        if verbosity >= 2:
                            print(
                                " {:<3} ({:3d}), {:<3} ({:3d}), {:<3} ({:3d}) and {:<3} ({:3d}) ({: 7.2f} deg)".format(
                                    molecule.atoms[molecule.angles[i][1]].symbol, molecule.angles[i][1],
                                    molecule.atoms[molecule.angles[i][2]].symbol, molecule.angles[i][2],
                                    molecule.atoms[molecule.angles[j][0]].symbol, molecule.angles[j][0],
                                    molecule.atoms[molecule.angles[i][0]].symbol, molecule.angles[i][0],
                                    math.degrees(molecule.outofplaneangle(len(molecule.threefolds) - 1))))

    # Now that we know bonds, angles, dihedrals and threefolds we determine the corresponding force constants
    # Bonds first:
    if verbosity >= 2:
        print("\nAdding Force Field bond stretching terms to WellFARe molecule: ", molecule.name)
    for i in range(0, len(molecule.bonds)):
        # print(molecule.atoms[molecule.bonds[i][0]].coord[1])
        a = numpy.array([molecule.atoms[molecule.bonds[i][0]].coord[0], molecule.atoms[molecule.bonds[i][0]].coord[1],
                         molecule.atoms[molecule.bonds[i][0]].coord[2]])
        b = numpy.array([molecule.atoms[molecule.bonds[i][1]].coord[0], molecule.atoms[molecule.bonds[i][1]].coord[1],
                         molecule.atoms[molecule.bonds[i][1]].coord[2]])
        c1 = (a - b)
        c2 = (b - a)
        c = numpy.zeros(molecule.numatoms() * 3)
        c[3 * molecule.bonds[i][0]] = c1[0]
        c[3 * molecule.bonds[i][0] + 1] = c1[1]
        c[3 * molecule.bonds[i][0] + 2] = c1[2]
        c[3 * molecule.bonds[i][1]] = c2[0]
        c[3 * molecule.bonds[i][1] + 1] = c2[1]
        c[3 * molecule.bonds[i][1] + 2] = c2[2]
        c = c / numpy.linalg.norm(c)
        fc = numpy.dot(numpy.dot(c, H), numpy.transpose(c))
        if fc < 0.002:
            ProgramWarning()
            print(" This force constant is smaller than 0.002")
        if verbosity >= 2:
            print(" {:<3} ({:3d}) and {:<3} ({:3d}) (Force constant: {: .3f})".format(
                molecule.atoms[molecule.bonds[i][0]].symbol, molecule.bonds[i][0],
                molecule.atoms[molecule.bonds[i][1]].symbol, molecule.bonds[i][1], fc))
        molecule.addFFStretch(molecule.bonds[i][0], molecule.bonds[i][1],
                              molecule.atmatmdist(molecule.bonds[i][0], molecule.bonds[i][1]), 3,
                              [fc, "b", molecule.atoms[molecule.bonds[i][0]].symbol,
                               molecule.atoms[molecule.bonds[i][1]].symbol])
    # Note "b" as an argument in the previouw line is a placeholder so that indices are consistent in the FFstretch class
    # it  would need replacing with the appropriate value to make using the  Morse potential an option

    # Then 1,3-stretches:
    if verbosity >= 2:
        print("\nAdding Force Field 1,3-bond stretching terms to WellFARe molecule: ", molecule.name)
    for i in range(0, len(molecule.angles)):
        a = numpy.array([molecule.atoms[molecule.angles[i][0]].coord[0], molecule.atoms[molecule.angles[i][0]].coord[1],
                         molecule.atoms[molecule.angles[i][0]].coord[2]])
        b = numpy.array([molecule.atoms[molecule.angles[i][2]].coord[0], molecule.atoms[molecule.angles[i][2]].coord[1],
                         molecule.atoms[molecule.angles[i][2]].coord[2]])
        c1 = (a - b)
        c2 = (b - a)
        c = numpy.zeros(molecule.numatoms() * 3)
        c[3 * molecule.angles[i][0]] = c1[0]
        c[3 * molecule.angles[i][0] + 1] = c1[1]
        c[3 * molecule.angles[i][0] + 2] = c1[2]
        c[3 * molecule.angles[i][1]] = c2[0]
        c[3 * molecule.angles[i][1] + 1] = c2[1]
        c[3 * molecule.angles[i][1] + 2] = c2[2]
        c = c / numpy.linalg.norm(c)
        fc = numpy.dot(numpy.dot(c, H), numpy.transpose(c))
        if fc < 0.002:
            ProgramWarning()
            print(" This force constant is smaller than 0.002")
        if verbosity >= 2:
            print(" {:<3} ({:3d}) and {:<3} ({:3d}) (Force constant: {: .3f})".format(
                molecule.atoms[molecule.angles[i][0]].symbol, molecule.angles[i][0],
                molecule.atoms[molecule.angles[i][0]].symbol, molecule.angles[i][0], fc))
        molecule.addFFStr13(molecule.angles[i][0], molecule.angles[i][2],
                            molecule.atmatmdist(molecule.angles[i][0], molecule.angles[i][2]), 4,
                            [fc, "b", molecule.atoms[molecule.angles[i][0]].symbol,
                             molecule.atoms[molecule.angles[i][2]].symbol])

    # Then angle bends:
    if verbosity >= 2:
        print("\nAdding Force Field angle bending terms to WellFARe molecule: ", molecule.name)
    for i in range(0, len(molecule.angles)):
        a = numpy.array([molecule.atoms[molecule.angles[i][0]].coord[0], molecule.atoms[molecule.angles[i][0]].coord[1],
                         molecule.atoms[molecule.angles[i][0]].coord[2]])
        b = numpy.array([molecule.atoms[molecule.angles[i][1]].coord[0], molecule.atoms[molecule.angles[i][1]].coord[1],
                         molecule.atoms[molecule.angles[i][1]].coord[2]])
        c = numpy.array([molecule.atoms[molecule.angles[i][2]].coord[0], molecule.atoms[molecule.angles[i][2]].coord[1],
                         molecule.atoms[molecule.angles[i][2]].coord[2]])
        aprime = a - b
        bprime = c - b
        p = numpy.cross(aprime, bprime)
        adprime = numpy.cross(p, aprime)
        bdprime = numpy.cross(bprime, p)
        c = numpy.zeros(molecule.numatoms() * 3)
        c[3 * molecule.angles[i][0]] = adprime[0]
        c[3 * molecule.angles[i][0] + 1] = adprime[1]
        c[3 * molecule.angles[i][0] + 2] = adprime[2]
        c[3 * molecule.angles[i][2]] = bdprime[0]
        c[3 * molecule.angles[i][2] + 1] = bdprime[1]
        c[3 * molecule.angles[i][2] + 2] = bdprime[2]
        # Temporary fix to avoid divide-by-zero errors follows, may be replaced by better check in future
        if c.all() == numpy.zeros(molecule.numatoms() * 3).all():
            print("Zero vector returned while extracting angle bend force constants, skipping normalisation")
        else:
            c = c / numpy.linalg.norm(c)
        fc = numpy.dot(numpy.dot(c, H), numpy.transpose(c))
        if fc < 0.002:
            ProgramWarning()
            print(" This force constant is smaller than 0.002")
        if verbosity >= 2:
            print(" {:<3} ({:3d}), {:<3} ({:3d}), {:<3} ({:3d}) and {:<3} ({:3d}) (Force constant: {: .3f})".format(
                molecule.atoms[molecule.angles[i][0]].symbol, molecule.angles[i][0],
                molecule.atoms[molecule.angles[i][1]].symbol, molecule.angles[i][1],
                molecule.atoms[molecule.angles[i][1]].symbol, molecule.angles[i][1],
                molecule.atoms[molecule.angles[i][2]].symbol, molecule.angles[i][2], fc))
        molecule.addFFBend(molecule.angles[i][0], molecule.angles[i][1], molecule.angles[i][2], molecule.bondangle(i),
                           2, [fc, molecule.atoms[molecule.angles[i][0]].symbol,
                               molecule.atoms[molecule.angles[i][1]].symbol,
                               molecule.atoms[molecule.angles[i][2]].symbol,
                               molecule.atmatmdist(molecule.angles[i][0], molecule.angles[i][1]),
                               molecule.atmatmdist(molecule.angles[i][1], molecule.angles[i][2])])
    # currently initiating bends with extra information in arguments list to avoid calling molecule or atom class methods inside FFBend.
    #  These quantities might ultimately be better included explicitly. 

    # Then dihedral torsions:
    if verbosity >= 2:
        print("\nAdding Force Field torsion terms to WellFARe molecule: ", molecule.name)
    for i in range(0, len(molecule.dihedrals)):
        a = numpy.array(
            [molecule.atoms[molecule.dihedrals[i][0]].coord[0], molecule.atoms[molecule.dihedrals[i][0]].coord[1],
             molecule.atoms[molecule.dihedrals[i][0]].coord[2]])
        b = numpy.array(
            [molecule.atoms[molecule.dihedrals[i][1]].coord[0], molecule.atoms[molecule.dihedrals[i][1]].coord[1],
             molecule.atoms[molecule.dihedrals[i][1]].coord[2]])
        c = numpy.array(
            [molecule.atoms[molecule.dihedrals[i][2]].coord[0], molecule.atoms[molecule.dihedrals[i][2]].coord[1],
             molecule.atoms[molecule.dihedrals[i][2]].coord[2]])
        d = numpy.array(
            [molecule.atoms[molecule.dihedrals[i][3]].coord[0], molecule.atoms[molecule.dihedrals[i][3]].coord[1],
             molecule.atoms[molecule.dihedrals[i][3]].coord[2]])
        aprime = a - b
        dprime = d - c
        c1prime = c - b
        c2prime = b - c
        p1 = numpy.cross(aprime, c1prime)
        p2 = numpy.cross(dprime, c2prime)
        c = numpy.zeros(molecule.numatoms() * 3)
        c[3 * molecule.dihedrals[i][0]] = p1[0]
        c[3 * molecule.dihedrals[i][0] + 1] = p1[1]
        c[3 * molecule.dihedrals[i][0] + 2] = p1[2]
        c[3 * molecule.dihedrals[i][2]] = p2[0]
        c[3 * molecule.dihedrals[i][2] + 1] = p2[1]
        c[3 * molecule.dihedrals[i][2] + 2] = p2[2]
        if c.all() == numpy.zeros(molecule.numatoms() * 3).all():
            print(
                "Zero vector returned while extracting force constants, skipping normalisation")  # Avoids fc=nan error for cases where all three atoms lie in the plane of two coordinate axes
        else:
            c = c / numpy.linalg.norm(c)
        # Note that the above is just an initial fix for cases where there would otherwise be a divide-by-zero error
        # These arise where several atoms have 0 in one coordinate which propogates through cross-products and are a 
        # side-effect of orienting molecule along principal axes from the centre of mass as bonds lie in a coordinate plane
        # Better fix may be possible/necessary - could translate molecule, for instance
        fc = numpy.dot(numpy.dot(c, H), numpy.transpose(c))
        if fc < 0.002:
            ProgramWarning()
            print(" This force constant is smaller than 0.002")
        if verbosity >= 2:
            print(" {:<3} ({:3d}), {:<3} ({:3d}), {:<3} ({:3d}) and {:<3} ({:3d}) (Force constant: {: .3f})".format(
                molecule.atoms[molecule.dihedrals[i][0]].symbol, molecule.dihedrals[i][0],
                molecule.atoms[molecule.dihedrals[i][1]].symbol, molecule.dihedrals[i][1],
                molecule.atoms[molecule.dihedrals[i][2]].symbol, molecule.dihedrals[i][2],
                molecule.atoms[molecule.dihedrals[i][3]].symbol, molecule.dihedrals[i][3], fc))
        molecule.addFFTorsion(molecule.dihedrals[i][0], molecule.dihedrals[i][1], molecule.dihedrals[i][2],
                              molecule.dihedrals[i][3], molecule.dihedralangle(i), 1,
                              [fc, molecule.atoms[molecule.dihedrals[i][0]].symbol,
                               molecule.atoms[molecule.dihedrals[i][1]].symbol,
                               molecule.atoms[molecule.dihedrals[i][2]].symbol,
                               molecule.atoms[molecule.dihedrals[i][3]].symbol,
                               molecule.atmatmdist(molecule.dihedrals[i][0], molecule.dihedrals[i][1]),
                               molecule.atmatmdist(molecule.dihedrals[i][1], molecule.dihedrals[i][2]),
                               molecule.atmatmdist(molecule.dihedrals[i][2], molecule.dihedrals[i][3])])
    # As for bends, arg list now includes atom symbols and bond lengths, which could be separated out later

    # Threefold inversions last
    if verbosity >= 2:
        print("\nAdding Force Field inversion terms ro WellFARe molecule: ", molecule.name)
    # (Extracting force constants to be implemented later)
    for i in range(0, len(molecule.threefolds)):
        a = numpy.array(
            [molecule.atoms[molecule.threefolds[i][0]].coord[0], molecule.atoms[molecule.threefolds[i][0]].coord[1],
             molecule.atoms[molecule.threefolds[i][0]].coord[2]])
        b = numpy.array(
            [molecule.atoms[molecule.threefolds[i][1]].coord[0], molecule.atoms[molecule.threefolds[i][1]].coord[1],
             molecule.atoms[molecule.threefolds[i][1]].coord[2]])
        c = numpy.array(
            [molecule.atoms[molecule.threefolds[i][2]].coord[0], molecule.atoms[molecule.threefolds[i][2]].coord[1],
             molecule.atoms[molecule.threefolds[i][2]].coord[2]])
        d = numpy.array(
            [molecule.atoms[molecule.threefolds[i][3]].coord[0], molecule.atoms[molecule.threefolds[i][3]].coord[1],
             molecule.atoms[molecule.threefolds[i][3]].coord[2]])
        ba = a - b
        cb = b - c
        db = b - d
        dc = c - d
        bprime = numpy.cross(-cb, -db)
        cprime = numpy.cross(-dc, cb)
        dprime = numpy.cross(db, dc)
        aprime = numpy.cross(bprime, numpy.cross(-ba, bprime)) / numpy.dot(bprime, bprime)
        c = numpy.zeros(molecule.numatoms() * 3)
        c[3 * molecule.threefolds[i][0]] = aprime[0]
        c[3 * molecule.threefolds[i][0] + 1] = aprime[1]
        c[3 * molecule.threefolds[i][0] + 2] = aprime[2]
        c[3 * molecule.threefolds[i][1]] = bprime[0]
        c[3 * molecule.threefolds[i][1] + 1] = bprime[1]
        c[3 * molecule.threefolds[i][1] + 2] = bprime[2]
        c[3 * molecule.threefolds[i][2]] = cprime[0]
        c[3 * molecule.threefolds[i][2] + 1] = cprime[1]
        c[3 * molecule.threefolds[i][2] + 2] = cprime[2]
        c[3 * molecule.threefolds[i][3]] = dprime[0]
        c[3 * molecule.threefolds[i][3] + 1] = dprime[1]
        c[3 * molecule.threefolds[i][3] + 2] = dprime[2]
        # Temporary fix to avoid divide-by-zero errors follows, may be replaced by better check in future
        if c.all() == numpy.zeros(molecule.numatoms() * 3).all():
            print("Zero vector returned while extracting inversion force constants, skipping normalisation")
        else:
            c = c / numpy.linalg.norm(c)
        fc = numpy.dot(numpy.dot(c, H), numpy.transpose(c))
        if fc < 0.002:
            ProgramWarning()
            print(" This force constant is smaller than 0.002")
        if verbosity >= 2:
            print(" {:<3} ({:3d}), {:<3} ({:3d}), {:<3} ({:3d}) and {:<3} ({:3d}) (Force constant: {: .3f})".format(
                molecule.atoms[molecule.threefolds[i][0]].symbol, molecule.threefolds[i][0],
                molecule.atoms[molecule.threefolds[i][1]].symbol, molecule.threefolds[i][1],
                molecule.atoms[molecule.threefolds[i][2]].symbol, molecule.threefolds[i][2],
                molecule.atoms[molecule.threefolds[i][3]].symbol, molecule.threefolds[i][3], fc))
        molecule.addFFInversion(molecule.threefolds[i][0], molecule.threefolds[i][1], molecule.threefolds[i][2],
                                molecule.threefolds[i][3], molecule.outofplaneangle(i), 2,
                                [fc, molecule.atoms[molecule.threefolds[i][0]].symbol,
                                 molecule.atoms[molecule.threefolds[i][1]].symbol,
                                 molecule.atoms[molecule.threefolds[i][2]].symbol,
                                 molecule.atoms[molecule.threefolds[i][3]].symbol,
                                 molecule.atmatmdist(molecule.threefolds[i][0], molecule.threefolds[i][1]),
                                 molecule.atmatmdist(molecule.threefolds[i][0], molecule.threefolds[i][2]),
                                 molecule.atmatmdist(molecule.threefolds[i][0], molecule.threefolds[i][3])])

    # Moving to noncovalent interactions
    # Locate hydrogen atoms and add them to the list hatoms
    if verbosity >= 2:
        print("\nListing hydrogen atoms in WellFARe molecule: ", molecule.name)
    for i in range(0, len(molecule.atoms)):
        if molecule.atoms[i].symbol == "H":
            molecule.addHAtom(i)
            if verbosity >= 2:
                print(molecule.atoms[i].symbol, molecule.atoms[i].coord)

    # Locate high electronegatvity atoms and add them to the list highENatoms
    if verbosity >= 2:
        print("Listing highly electronegative atoms in WellFARe molecule: ", molecule.name)
    for i in range(0, len(molecule.atoms)):
        sym = molecule.atoms[i].symbol
        if sym == "N" or sym == "O" or sym == "F" or sym == "S" or sym == "Cl":
            molecule.addhighENatom(i)
            if verbosity >= 2:
                print(sym, molecule.atoms[i].coord)

    # Locate halogen atoms and add them to the list halogens
    if verbosity >= 2:
        print("Listing non-F halogen atoms in WellFARe molecule: ", molecule.name)
    for i in range(0, len(molecule.atoms)):
        sym = molecule.atoms[i].symbol
        if sym == "Cl" or sym == "Br" or sym == "I" or sym == "At":
            molecule.addXatom(i)
            if verbosity >= 2:
                print(sym, molecule.atoms[i].coord)


    # Locate hydrogen bonding triples AHB and create instances of FFHBond
    if verbosity >= 2:
        print("\nAdding hydrogen bonds to WellFARe molecule: ", molecule.name)
    for i in range(0, len(molecule.bonds)):
        sym1 = molecule.atoms[molecule.bonds[i][0]].symbol
        sym2 = molecule.atoms[molecule.bonds[i][1]].symbol
        if sym1 == "H":
            atH = molecule.bonds[i][0]
            if sym2 == "N" or sym2 == "O" or sym2 == "F" or sym2 == "S" or sym2 == "Cl":
                atA = molecule.bonds[i][1]
                for j in range(0, len(molecule.atoms)):
                    sym3 = molecule.atoms[j].symbol
                    if sym3 == "N" or sym3 == "O" or sym3 == "F" or sym3 == "S" or sym3 == "Cl":
                        r = molecule.atmatmdist(atH, j)
                        r_check = SymbolToVdWRadius[molecule.atoms[atH].symbol] + SymbolToVdWRadius[
                            molecule.atoms[j].symbol]  # Using sum of van der Waals radii
                        if r <= r_check and j != atA:
                            theta = molecule.anybondangle(atA, atH, j)
                            molecule.addFFHBond(atA, atH, j, theta, 1,
                                                [sym2, molecule.atoms[atA].charge, sym1, molecule.atoms[atH].charge,
                                                 sym3, molecule.atoms[j].charge, molecule.atmatmdist(atA, atH),
                                                 molecule.atmatmdist(j, atH), molecule.atmatmdist(atA, j)])
                            if verbosity >= 2:
                                print(
                                    " ({:<3}, {:<3}, {:<3} {:3.2f} deg), {:<3}, [{:<3}, {:3.2f}, {:<3}, {:3.2f}, {:<3}, {:3.2f}, {:3.2f}, {:3.2f}, {:3.2f}]".format(
                                        atA, atH, j, theta, 1, sym2, molecule.atoms[atA].charge, sym1,
                                        molecule.atoms[atH].charge, sym3, molecule.atoms[j].charge,
                                        molecule.atmatmdist(atA, atH), molecule.atmatmdist(j, atH),
                                        molecule.atmatmdist(atA, j)))
        elif sym2 == "H":
            atH = molecule.bonds[i][1]
            if sym1 == "N" or sym1 == "O" or sym1 == "F" or sym1 == "S" or sym1 == "Cl":
                atA = molecule.bonds[i][0]
                for j in range(0, len(molecule.atoms)):
                    sym3 = molecule.atoms[j].symbol
                    if sym3 == "N" or sym3 == "O" or sym3 == "F" or sym3 == "S" or sym3 == "Cl":
                        r = molecule.atmatmdist(atH, j)
                        r_check = SymbolToVdWRadius[sym1] + SymbolToVdWRadius[
                            sym2]  # Sum of van der Waals radii again used as check
                        if r <= r_check and j != atA:
                            theta = molecule.anybondangle(atA, atH, j)
                            molecule.addFFHBond(atA, atH, j, theta, 1,
                                                [sym1, molecule.atoms[atA].charge, sym2, molecule.atoms[atH].charge,
                                                 sym3, molecule.atoms[j].charge, molecule.atmatmdist(atA, atH),
                                                 molecule.atmatmdist(j, atH), molecule.atmatmdist(atA, j)])
                            if verbosity >= 2:
                                print(
                                    " ({:<3}, {:<3}, {:<3} {:3.2f} deg), {:<3}, [{:<3}, {:3.2f}, {:<3}, {:3.2f}, {:<3}, {:3.2f}, {:3.2f}, {:3.2f}, {:3.2f}]".format(
                                        atA, atH, j, theta, 1, sym1, molecule.atoms[atA].charge, sym2,
                                        molecule.atoms[atH].charge, sym3, molecule.atoms[j].charge,
                                        molecule.atmatmdist(atA, atH), molecule.atmatmdist(j, atH),
                                        molecule.atmatmdist(atA, j)))


# End of routine

################################################################################
# Last function needed to set up force field: carry out Hessian fit to determine
# force constants for stretches, bends and inversion potentials
################################################################################

def fitForceConstants(molecule, verbosity=0):
    if verbosity >= 1:
        print("\nFitting force constants for WellFARe molecule: ", molecule.name)
    # Construct a list of the force constants initially assigned to the molecule
    ForceConstants = []
    for i in range(len(molecule.stretch)):
        ForceConstants.append(molecule.stretch[i].k_str)
    for i in range(len(molecule.str13)):
        ForceConstants.append(molecule.str13[i].k_str)
    for i in range(len(molecule.bend)):
        ForceConstants.append(molecule.bend[i].k)
    for i in range(len(molecule.inv)):
        ForceConstants.append(molecule.inv[i].k_inv)
    InitialFC = ForceConstants
    if verbosity >= 1:
        print("\nForce constants to be optimised:")
        print(ForceConstants)

    # Carry out Hessian fitting procedure to determine the appropriate values of those force constants
    timestamp("Running Optimisation ")  # REMOVE ONCE FIXED
    xopt = scipy.optimize.fmin_bfgs(molecule.HessianDiffSquared, ForceConstants,
                                    gtol=0.01)  # Other optimisers might be more suitable, and extra parameters can be specified if needed
    # Tolerance has been increased from the default 1e-05 in order to speed up the optimisation
    if verbosity >= 1:
        if verbosity >= 2:
            print("\nInitial Force constants:")
            print(InitialFC)
            print("\QM Hessian:")
            print(molecule.H_QM)
            print("\nInitial Force Field Hessian:")
            print(molecule.kdepHessian(InitialFC))
        timestamp("\nFitted Force constants: ")
        print(xopt)
        print("\nForce Field Hessian with those force constants:")
        print(molecule.kdepHessian(xopt))

    # Assign the fitted force constants to the corresponding force field potentials
    # Note: could add verbosity option here to print a representation of each after changing force constant, but probably not required
    for i in range(len(molecule.stretch)):
        molecule.stretch[i].setk(xopt[i])
    for i in range(len(molecule.str13)):
        molecule.str13[i].setk(xopt[len(molecule.stretch) + i])
    for i in range(len(molecule.bend)):
        molecule.bend[i].setk(xopt[len(molecule.stretch) + len(molecule.str13) + i])
    for i in range(len(molecule.inv)):
        molecule.inv[i].setk(xopt[len(molecule.stretch) + len(molecule.str13) + len(molecule.bend) + i])


######################################################################################
# Functions for additional calculations, optional to the program, defined below here #
######################################################################################

def dissociateBond(molecule, atom1, atom2, epsilon, cutoff, verbosity=1):
    """
  Function to progressively increase separation between atoms number a1 and a2 in increments of size epsilon until it exceeds the specified cutoff, and calculate potential energy at each new geometry
  Output can be used to generate a dissociation curve
  """

    if verbosity >= 1:
        print("\nSetting up for bond dissociation")
        print("epsilon = " + str(epsilon))
        print("cutoff = " + str(cutoff))
    # Identify the atoms to be moved apart
    a1 = molecule.atoms[atom1]
    a2 = molecule.atoms[atom2]
    # Check that the two are actually bonded in the molecule, give error/warning if not
    if [atom1, atom2] in molecule.bonds:
        bond = [atom1, atom2]
        if verbosity >= 2:
            print("Bond to dissociate: " + str(bond))
    elif [atom2, atom1] in molecule.bonds:
        bond = [atom2, atom1]
        if verbosity >= 2:
            print("Bond to dissociate: " + str(bond))
    else:
        ProgramWarning()
        print("\nAtoms for dissociation are not bonded in input structure")
    # Determine how many other bonding partners each atom has
    bonded_a1 = []
    bonded_a2 = []
    for i in range(len(molecule.atoms)):
        if atom1 != i and atom2 != i:
            if [atom1, i] in molecule.bonds or [i, atom1] in molecule.bonds:
                bonded_a1.append(i)
            elif [atom2, i] in molecule.bonds or [i, atom2] in molecule.bonds:
                bonded_a2.append(i)
    nbonds_a1 = len(bonded_a1)
    nbonds_a2 = len(bonded_a2)
    # Set the atom with fewest bonding partners to move, the other to remain stationary
    atM = a1
    atomM = atom1
    atS = a2
    atomS = atom2
    if nbonds_a1 <= nbonds_a2:
        if verbosity >= 1:
            print("Fixing atom " + str(atom2))
            print("Moving atom " + str(atom1))
    elif nbonds_a2 < nbonds_a2:
        atM = a2
        atomM = atom2
        atS = a1
        atomS = atom1
        if verbosity >= 1:
            print("Fixing atom " + str(atom1))
            print("Moving atom " + str(atom2))
    # Calculate the vector along which movement should occur
    dvector = [atM.coord[i] - atS.coord[i] for i in range(3)]
    norm_dv = numpy.linalg.norm(dvector)
    unitdv = [dvector[i] / norm_dv for i in range(3)]
    movevector = [unitdv[i] * epsilon for i in range(3)]

    # Calculate initial energy, at half the equilibrium bond length
    DissocEnergies = []
    r0 = molecule.atmatmdist(atom1, atom2)
    ri = r0 / 2
    cartCoords = molecule.cartesianCoordinates()
    initialmv = [unitdv[i] * (r0 / 2) for i in range(3)]
    for i in range(3):
        cartCoords[(atomM * 3) + i] = cartCoords[(atomM * 3) + i] - initialmv[i]
    if verbosity >= 1:
        print("\nCalculating initial energy, at separation " + str(ri) + ":")
        ei = molecule.FFEnergy(cartCoords, verbosity=1)
    else:
        ei = molecule.FFEnergy(cartCoords, verbosity=0)
    DissocEnergies.append([ri, Bohr2Ang(ri), ei, au2kcal_mol(ei)])
    # Iteratively increase separation and re-calculate energy for the distorted geometry until separation exceeds the given cutoff
    r = ri
    while r <= cutoff:
        for i in range(3):
            cartCoords[(atomM * 3) + i] = cartCoords[(atomM * 3) + i] + movevector[i]
        r = r + epsilon
        if verbosity >= 1:
            print("\nCalculating energy at separation " + str(r) + ":")
            e = molecule.FFEnergy(cartCoords, verbosity=1)
        else:
            e = molecule.FFEnergy(cartCoords, verbosity=0)
        DissocEnergies.append([r, Bohr2Ang(r), e, au2kcal_mol(e)])
    else:
        if verbosity >= 1:
            print("\nCutoff reached, dissociation calculation complete")
    # NOTE Currently moving only one atom - full version should move whole bonded fragment
    # Calculate total increase in separation, total energy change
    nsteps = len(DissocEnergies)
    delta_r = DissocEnergies[nsteps - 1][0] - ri
    delta_e = DissocEnergies[nsteps - 1][1] - ei

    # Output results in a nice format
    print("\nNumber of points evaluated: " + str(nsteps))
    print(
        "\nCalculated distance and dissociation energy for " + str(molecule.name) + " atoms " + str(atomS) + ", " + str(
            atomM) + ":")
    print("---------------------------------------------------------------------")
    print('{:<11}   {:<11}   {:<11}      {:<11}'.format("r (AU)", "r (Angstrom)", "E (Hartree)", "E (kcal/mol)"))
    for i in range(len(DissocEnergies)):
        data = DissocEnergies[i]
        print('{:.8f}    {:.8f}     {:.7f}     {:.5f}'.format(data[0], data[1], data[2], data[3]))
    print("\nTotal increase in separation: " + str(delta_r))
    print("\nTotal energy change: " + str(delta_e))


# End of routine

def ObjFuncSEAM(X, reactant, product):
    """
  Objective function employed in the SEAM method to find a minimum on the intersction of reactant and product potential energy surfaces
  Here reactant and product both belong to the class Molecule, and X is a list of cartesian coordinates with the Lagrange multiplier L appended
  """
    cartCoordinates = numpy.zeros(len(X) - 1)
    for i in range(len(X) - 1):
        cartCoordinates[i] = X[i]
    lm = X[len(X) - 1]

    E_r = reactant.FFEnergy(cartCoordinates, verbosity=0, dtyp=1)  # Extra arguments included as a reminder of options
    E_p = product.FFEnergy(cartCoordinates, verbosity=0, dtyp=1)

    L = (E_r + E_p) - lm * (E_r - E_p)

    return L


def dObjFuncSEAM(X, reactant, product):
    """
  Function to calculate the (first) derivate of the objective function for the SEAM method by finite differences
  """
    dL = numpy.zeros(len(X))
    h = 1e-3  # Step size for use in numerical evaluation of gradient, can be tailored

    for i in range(len(X)):
        dX = numpy.zeros(len(X))
        dX[i] = h
        dL[i] = (ObjFuncSEAM((X + dX), reactant, product) - ObjFuncSEAM((X - dX), reactant, product)) / (2 * h)

    return dL


def TSbySEAM(reactant, product, verbosity=1):
    """
  Function to carry out the search for a transition state by the SEAM method
  """

    # Set up the list of variables X to be optimised, taking a 50:50 interpolation of coordinates for reactant and product as the initial guess
    CoordsR = reactant.cartesianCoordinates()
    print(CoordsR)  # Testing only
    CoordsP = product.cartesianCoordinates()
    guessL = 0  # Initial guess for the Lagrange multiplier can be modified as desired
    if len(CoordsR) == len(CoordsP):
        guessCoords = numpy.zeros(len(CoordsR))
        for i in range(len(CoordsR)):
            guessCoords[i] = (CoordsR[i] + CoordsP[i]) / 2
    else:
        print("Molecules provided to SEAM method have incompatible coordinates (length mismatch)")
        ProgramError()  # Consider arranging such that the function does not execute further in this case
    X = numpy.zeros(len(guessCoords) + 1)
    for i in range(len(guessCoords)):
        X[i] = guessCoords[i]
    X[len(guessCoords)] = guessL
    if verbosity >= 1:
        print("Starting SEAM with initial coordinates + multiplier:")
        print(X)

    # Starting from X, optimise to a stationary point of the objective function
    X_opt = scipy.optimize.fsolve(dObjFuncSEAM, X, args=(reactant, product))

    # Check that the stationary point found is a minimum?

    if verbosity >= 1:
        print("Candidate transition state + multiplier located by SEAM:")
        print(X_opt)

    # Calculate the energy at the transition state (should be equal using either force field)
    TS = numpy.zeros(len(X_opt) - 1)
    for i in range((len(X_opt) - 1)):
        TS[i] = X_opt[i]
    E_TS = reactant.FFEnergy(TS)

    # Print out the SEAM transition state in a format readable by Gaussian 
    print("\nSEAM Transition state in Gaussian format:")
    print(reactant.gaussStringatX(TS))
    print("\nEnergy at Transition state: " + str(E_TS))

    return TS


################################################################################
#                                                                              #
# This is the part of the program where the command line arguments are defined #
#                                                                              #
################################################################################

parser = argparse.ArgumentParser(
    description="WellFAReFF: Wellington Fast Assessment of Reactions - Force Field",
    epilog="recognised filetypes: g09, orca")
parser.add_argument("-r", "--reactant", metavar='file', help="input file with qc data of the reactant",
                    default="g09-dielsalder-r.log")
parser.add_argument("-p", "--product", metavar='file', help="input file with qc data of the product",
                    default="g09-dielsalder-p.log")
parser.add_argument("-v", "--verbosity", help="increase output verbosity", type=int, choices=[0, 1, 2, 3], default=1)

args = parser.parse_args()

###############################################################################
#                                                                             #
# The main part of the program starts here                                    #
#                                                                             #
###############################################################################

# Print GPL v3 statement and program header
ProgramHeader()

# print("Number of Atoms: ", reactant_mol.numatoms(), "Multiplicity: ", reactant_mol.mult)

# print(molecule)
# print("Molecular mass = ", reactant_mol.mass())
# reactant_mol.orient()

# print(reactant_mol.gaussString())

# print("Bonds:")
# for i in reactant_mol.bonds:
#   print(i)
#
# print("")
# print("Angles:")
# for i in reactant_mol.angles:
#   print(i)
#
# print("")
# print("Angles in degrees:")
# for i in range(len(reactant_mol.angles)):
#   print(math.degrees(reactant_mol.bondangle(i)))
#
# print("")
# print("Dihedrals:")
# for i in reactant_mol.dihedrals:
#   print(i)
#
# print("")
# print("Dihedral angles in degrees:")
# for i in range(len(reactant_mol.dihedrals)):
#   print(math.degrees(reactant_mol.dihedralangle(i)))

# print("")
# print("Bond Stretches:")
# for i in reactant_mol.stretch:
#   print(i)
#
# print("")
# print("1-3 Bond Stretches:")
# for i in reactant_mol.str13:
#   print(i)
#
# print("")
# print("Angle Bends:")
# for i in reactant_mol.bend:
#   print(i)
#
# print("")
# print("Dihedral Torsions:")
# for i in reactant_mol.tors:
#   print(i)

reactant_mol = Molecule("Reactant", 0)
extractCoordinates(args.reactant, reactant_mol, verbosity=args.verbosity)
fitForceConstants(reactant_mol, verbosity=args.verbosity)

print("\nCartesian Coordinates of Reactant (as one list):")
print(reactant_mol.cartesianCoordinates())

print("\nForce Field Energy of Reactant:")
print(reactant_mol.FFEnergy(reactant_mol.cartesianCoordinates(), verbosity=args.verbosity))

print("\nGeometry Optimizer (Reactant):")
initialcoords2optimiseR = reactant_mol.cartesianCoordinates()
xopt = scipy.optimize.fmin_bfgs(reactant_mol.FFEnergy, initialcoords2optimiseR, gtol=0.00005)
print("\Optimized Geometry coordinates (Reactant):")
print(xopt)

reactant_mol.setGeometry(xopt)
print("\nOptimized Geometry in Gaussian format (Reactant):")
print(reactant_mol.gaussString())


# product_mol = Molecule("Product",0)
# extractCoordinates(infile2, product_mol, verbosity = 2)
# fitForceConstants(product_mol, verbosity = 2)

# print("\nCartesian Coordinates of Product (as one list):")
# print(product_mol.cartesianCoordinates())

# print("\nForce Field Energy of Product:")
# print(product_mol.FFEnergy(product_mol.cartesianCoordinates(), verbosity = 1))

# print("\nGeometry Optimizer (Product):")
# initialcoords2optimiseP = product_mol.cartesianCoordinates()
# xopt = scipy.optimize.fmin_bfgs(product_mol.FFEnergy, initialcoords2optimiseP, gtol=0.00005)
# print("\Optimized Geometry coordinates (Product):")
# print(xopt)

# product_mol.setGeometry(xopt)
# print("\nOptimized Geometry in Gaussian format (Product):")
# print(product_mol.gaussString())


# print("\nDistort Geometry by interpolation and print energy again:")
# coordinates2optimiseR = reactant_mol.cartesianCoordinates()
# coordinates2optimiseP = product_mol.cartesianCoordinates()

# coordinates2optimiseR = (numpy.array(coordinates2optimiseR)+(numpy.array(coordinates2optimiseP))/2.0)

# print(reactant_mol.FFEnergy(coordinates2optimiseR, verbosity = 1))

# print("\nGeometry Optimizer:")
# xopt = scipy.optimize.fmin_bfgs(reactant_mol.FFEnergy, coordinates2optimiseR, gtol=0.00005)
# print("\nOptimized Geometry coordinates:")
# print(xopt)

# print("\nBond Dissociation:")
# dissociateBond(reactant_mol, 0, 1, 10**-3, 14)

# TSbySEAM(reactant_mol, product_mol, verbosity = 1)

# Test the screening procedure used to determine appropriate values of elstat_AB
# print("\nTesting electrostatic topological screening parameter procedure\n")
# for i in range(len(reactant_mol.atoms)):
#  for j in range(len(reactant_mol.atoms)):
#    reactant_mol.screen_ES(i, j)

ProgramFooter()
