#!/usr/bin/python

import sys
import getopt
import math
import numpy

#--------------------------------------------------------------
# Define some handy functions that will come in useful
#--------------------------------------------------------------

def iofiles(argv):
   inputfile = ''
   outputfile = ''
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
   except getopt.GetoptError:
      print 'test.py -i <inputfile> -o <outputfile>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'test.py -i <inputfile> -o <outputfile>'
         sys.exit()
      elif opt in ("-i", "--ifile"):
         inputfile = arg
      elif opt in ("-o", "--ofile"):
         outputfile = arg
   if inputfile == '':
      print "Input file not specified on command line. Will use default."
      inputfile="g09-ethane.log"
   print "Input file is : ", inputfile
   if outputfile == "":
      print "Output file not specified on command line. Will use default."
      outputfile="WellFAReFF.log"
   print "Output file is: ", outputfile
   return (inputfile, outputfile)

# Define dictionary to convert atomic symbols to atomic numbers
SymbolToNumber = {
"H"  :1,
"He" :2,
"Li" :3,
"Be" :4,
"B"  :5,
"C"  :6,
"N"  :7,
"O"  :8,
"F"  :9,
"Ne" :10,
"Na" :11,
"Mg" :12,
"Al" :13,
"Si" :14,
"P"  :15,
"S"  :16,
"Cl" :17,
"Ar" :18,
"K"  :19,
"Ca" :20,
"Sc" :21,
"Ti" :22,
"V"  :23,
"Cr" :24,
"Mn" :25,
"Fe" :26,
"Co" :27,
"Ni" :28,
"Cu" :29,
"Zn" :30,
"Ga" :31,
"Ge" :32,
"As" :33,
"Se" :34,
"Br" :35,
"Kr" :36,
"Rb" :37,
"Sr" :38,
"Y"  :39,
"Zr" :40,
"Nb" :41,
"Mo" :42,
"Tc" :43,
"Ru" :44,
"Rh" :45,
"Pd" :46,
"Ag" :47,
"Cd" :48,
"In" :49,
"Sn" :50,
"Sb" :51,
"Te" :52,
"I"  :53,
"Xe" :54,
"Cs" :55,
"Ba" :56,
"La" :57,
"Ce" :58,
"Pr" :59,
"Nd" :60,
"Pm" :61,
"Sm" :62,
"Eu" :63,
"Gd" :64,
"Tb" :65,
"Dy" :66,
"Ho" :67,
"Er" :68,
"Tm" :69,
"Yb" :70,
"Lu" :71,
"Hf" :72,
"Ta" :73,
"W"  :74,
"Re" :75,
"Os" :76,
"Ir" :77,
"Pt" :78,
"Au" :79,
"Hg" :80,
"Tl" :81,
"Pb" :82,
"Bi" :83,
"Po" :84,
"At" :85,
"Rn" :86,
"Fr" :87,
"Ra" :88,
"Ac" :89,
"Th" :90,
"Pa" :91,
"U"  :92,
"Np" :93,
"Pu" :94,
"Am" :95,
"Cm" :96,
"Bk" :97,
"Cf" :98,
"Es" :99,
"Fm" :100,
"Md" :101,
"No" :102,
"Lr" :103,
"Rf" :104,
"Db" :105,
"Sg" :106,
"Bh" :107,
"Hs" :108,
"Mt" :109,
"Ds" :110,
"Rg" :111,
"Cn" :112,
"Uut":113,
"Fl" :114,
"Uup":115,
"Lv" :116,
"Uus":117,
"Uuo":118}

# Invert the above: atomic numbers to atomic symbols
NumberToSymbol = {v: k for k, v in SymbolToNumber.items()}

# Define dictionary to convert atomic symbols to masses
SymbolToMass = {
"H" : 1.00794,
"He": 4.002602,
"Li": 6.941,
"Be": 9.012182,
"B": 10.811,
"C": 12.0107,
"N": 14.0067,
"O": 15.9994,
"F": 18.9984032,
"Ne": 20.1797,
"Na": 22.98976928,
"Mg": 24.3050,
"Al": 26.9815386,
"Si": 28.0855,
"P": 30.973762,
"S": 32.065,
"Cl": 35.453,
"Ar": 39.948,
"K": 39.0983,
"Ca": 40.078,
"Sc": 44.955912,
"Ti": 47.867,
"V": 50.9415,
"Cr": 51.9961,
"Mn": 54.938045,
"Fe": 55.845,
"Co": 58.933195,
"Ni": 58.6934,
"Cu": 63.546,
"Zn": 65.38,
"Ga": 69.723,
"Ge": 72.64,
"As": 74.92160,
"Se": 78.96,
"Br": 79.904,
"Kr": 83.798,
"Rb": 85.4678,
"Sr": 87.62,
"Y": 88.90585,
"Zr": 91.224,
"Nb": 92.90638,
"Mo": 95.96,
"Tc": 98.0,
"Ru": 101.07,
"Rh": 102.90550,
"Pd": 106.42,
"Ag": 107.8682,
"Cd": 112.411,
"In": 114.818,
"Sn": 118.710,
"Sb": 121.760,
"Te": 127.60,
"I": 126.90447,
"Xe": 131.293,
"Cs": 132.9054519,
"Ba": 137.327,
"La": 138.90547,
"Ce": 140.116,
"Pr": 140.90765,
"Nd": 144.242,
"Pm": 145.0,
"Sm": 150.36,
"Eu": 151.964,
"Gd": 157.25,
"Tb": 158.92535,
"Dy": 162.500,
"Ho": 164.93032,
"Er": 167.259,
"Tm": 168.93421,
"Yb": 173.054,
"Lu": 174.9668,
"Hf": 178.49,
"Ta": 180.94788,
"W": 183.84,
"Re": 186.207,
"Os": 190.23,
"Ir": 192.217,
"Pt": 195.084,
"Au": 196.966569,
"Hg": 200.59,
"Tl": 204.3833,
"Pb": 207.2,
"Bi": 208.98040,
"Po": 209.0,
"At": 210.0,
"Rn": 222.0,
"Fr": 223.0,
"Ra": 226.0,
"Ac": 227.0,
"Th": 232.03806,
"Pa": 231.03588,
"U": 238.02891,
"Np": 237.0,
"Pu": 244.0,
"Am": 243.0,
"Cm": 247.0,
"Bk": 247.0,
"Cf": 251.0,
"Es": 252.0,
"Fm": 257.0,
"Md": 258.0,
"No": 259.0,
"Lr": 262.0,
"Rf": 267.0,
"Db": 268.0,
"Sg": 271.0,
"Bh": 272.0,
"Hs": 270.0,
"Mt": 276.0,
"Ds": 281.0,
"Rg": 280.0,
"Cn": 285.0,
"Uut": 284.0,
"Uuq": 289.0,
"Uup": 288.0,
"Uuh": 293.0,
"Uuo": 294.0}

# Define dictionary to convert atomic symbols to covalent radii (in Angstrom)
SymbolToRadius = {
"H"  : 0.37,
"He" : 0.32,
"Li" : 1.34,
"Be" : 0.90,
"B"  : 0.82,
"C"  : 0.77,
"N"  : 0.75,
"O"  : 0.73,
"F"  : 0.71,
"Ne" : 0.69,
"Na" : 1.54,
"Mg" : 1.30,
"Al" : 1.18,
"Si" : 1.11,
"P"  : 1.06,
"S"  : 1.02,
"Cl" : 0.99,
"Ar" : 0.97,
"K"  : 1.96,
"Ca" : 1.74,
"Sc" : 1.44,
"Ti" : 1.36,
"V"  : 1.25,
"Cr" : 1.27,
"Mn" : 1.39,
"Fe" : 1.25,
"Co" : 1.26,
"Ni" : 1.21,
"Cu" : 1.38,
"Zn" : 1.31,
"Ga" : 1.26,
"Ge" : 1.22,
"As" : 1.19,
"Se" : 1.16,
"Br" : 1.14,
"Kr" : 1.10,
"Rb" : 2.11,
"Sr" : 1.92,
"Y"  : 1.62,
"Zr" : 1.48,
"Nb" : 1.37,
"Mo" : 1.45,
"Tc" : 1.56,
"Ru" : 1.26,
"Rh" : 1.35,
"Pd" : 1.31,
"Ag" : 1.53,
"Cd" : 1.48,
"In" : 1.44,
"Sn" : 1.41,
"Sb" : 1.38,
"Te" : 1.35,
"I"  : 1.33,
"Xe" : 1.30,
"Cs" : 2.25,
"Ba" : 1.98,
"La" : 1.69,
"Ce" : 1.70,
"Pr" : 1.70,
"Nd" : 1.70,
"Pm" : 1.70,
"Sm" : 1.70,
"Eu" : 1.70,
"Gd" : 1.70,
"Tb" : 1.70,
"Dy" : 1.70,
"Ho" : 1.70,
"Er" : 1.70,
"Tm" : 1.70,
"Yb" : 1.70,
"Lu" : 1.60,
"Hf" : 1.50,
"Ta" : 1.38,
"W"  : 1.46,
"Re" : 1.59,
"Os" : 1.28,
"Ir" : 1.37,
"Pt" : 1.28,
"Au" : 1.44,
"Hg" : 1.49,
"Tl" : 1.48,
"Pb" : 1.47,
"Bi" : 1.46,
"Po" : 1.50,
"At" : 1.50,
"Rn" : 1.45,
"Fr" : 1.50,
"Ra" : 1.50,
"Ac" : 1.50,
"Th" : 1.50,
"Pa" : 1.50,
"U"  : 1.50,
"Np" : 1.50,
"Pu" : 1.50,
"Am" : 1.50,
"Cm" : 1.50,
"Bk" : 1.50,
"Cf" : 1.50,
"Es" : 1.50,
"Fm" : 1.50,
"Md" : 1.50,
"No" : 1.50,
"Lr" : 1.50,
"Rf" : 1.50,
"Db" : 1.50,
"Sg" : 1.50,
"Bh" : 1.50,
"Hs" : 1.50,
"Mt" : 1.50,
"Ds" : 1.50,
"Rg" : 1.50,
"Cn" : 1.50,
"Uut" : 1.50,
"Uuq" : 1.50,
"Uup" : 1.50,
"Uuh" : 1.50,
"Uus" : 1.50,
"Uuo" : 1.50}

def Ang2Bohr(ang):
    return ang*1.889725989

def Bohr2Ang(bohr):
    return bohr/1.889725989

def extractCoordinates(filename):
  f = open(filename,'r')
  program = "N/A"
  for line in f:
  	if line.find("Entering Gaussian System, Link 0=g09") != -1:
  	  program = "g09"
  	  break
  	elif line.find("* O   R   C   A *") != -1:
  	  program = "orca"
  	  break
  f.close()
  geom = []
  geom2 = []
  if program == "g09":
    f = open(filename,'r')
    for line in f:
      if line.find("Standard orientation:") != -1:
        for i in range(0,4):
          readBuffer = f.next()
        while True:
          readBuffer = f.next()
          if readBuffer.find("-----------") == -1:
            geom.append(readBuffer)
          else:
            break
        break
    for i in geom:
      readBuffer=i.split()
      geom2.append([NumberToSymbol[int(readBuffer[1])],float(readBuffer[3]),float(readBuffer[4]),float(readBuffer[5])])
    f.close()
  elif program == "orca":
    f = open(filename,'r')
    for line in f:
      if line.find("CARTESIAN COORDINATES (ANGSTROEM)") != -1:
        readBuffer = f.next()
        while True:
          readBuffer = f.next()
          if readBuffer and readBuffer.strip():
            geom.append(readBuffer)
          else:
            break
        break
    for i in geom:
      readBuffer=i.split()
      geom2.append([readBuffer[0],float(readBuffer[1]),float(readBuffer[2]),float(readBuffer[3])])
    f.close()
  return geom2

def calcDistance(atom1, atom2):
    distance=(atom1[1]-atom2[1])*(atom1[1]-atom2[1])
    distance=distance+((atom1[2]-atom2[2])*(atom1[2]-atom2[2]))
    distance=distance+((atom1[3]-atom2[3])*(atom1[3]-atom2[3]))
    return math.sqrt(distance)

def calcAngle(atom1,atom2,atom3):
    v1 = [atom1[1]-atom2[1],atom1[2]-atom2[2],atom1[3]-atom2[3]]
    v2 = [atom3[1]-atom2[1],atom3[2]-atom2[2],atom3[3]-atom2[3]]
    v1mag = math.sqrt((v1[0]*v1[0])+(v1[1]*v1[1])+(v1[2]*v1[2]))
    v1norm = [v1[0]/v1mag, v1[1]/v1mag, v1[2]/v1mag]
    v2mag = math.sqrt((v2[0]*v2[0])+(v2[1]*v2[1])+(v2[2]*v2[2]))
    v2norm = [v2[0]/v2mag, v2[1]/v2mag, v2[2]/v2mag]
    res = (v1norm[0]*v2norm[0])+(v1norm[1]*v2norm[1])+(v1norm[2]*v2norm[2])
    return round(math.degrees(math.acos(res)),3)

def calcDihedral(atom1,atom2,atom3,atom4):
    v1 = [atom1[1]-atom2[1],atom1[2]-atom2[2],atom1[3]-atom2[3]]
    v1 = v1/numpy.linalg.norm(v1)
    v2 = [atom3[1]-atom2[1],atom3[2]-atom2[2],atom3[3]-atom2[3]]
    v2 = v2/numpy.linalg.norm(v2)
    v3 = [atom4[1]-atom3[1],atom4[2]-atom3[2],atom4[3]-atom3[3]]
    v3 = v3/numpy.linalg.norm(v3)
    n1 = numpy.cross(v1,v2)
    n1 = n1/numpy.linalg.norm(n1)
    n2 = numpy.cross(v2,v3)
    n2 = n2/numpy.linalg.norm(n2)
    res1 = numpy.dot(n1,n2)
    res2 = numpy.cross(n1,v2)
    res3 = numpy.dot(res2,n2)
    return round(math.degrees(math.atan2(res3,res1)),3)

# Routine to find bonds based on several alternatives
# method: 0 = within 110% of covalent radius, 1 = Wiberg Bond Index > 0.5 (TBI)
def findBonds(geom,method):
  bonds = []
  if method == 0:
    for i in range(0,len(geom)):
       for j in range(i+1,len(geom)):
           if calcDistance(geom[i],geom[j])<=(SymbolToRadius[geom[i][0]]+SymbolToRadius[geom[j][0]])*1.1:
              bonds.append([i,j])
  elif method == 1:
    print "Not implemented yet!"
  return bonds

# Routine to find angles based on the bonds already identified
def findAngles(bonds):
  angles = []
  for i in range(0,len(bonds)):
      for j in range(i+1,len(bonds)):
          if bonds[i][0]==bonds[j][0]:
              angles.append([bonds[i][1],bonds[i][0],bonds[j][1]])
          if bonds[i][0]==bonds[j][1]:
              angles.append([bonds[i][1],bonds[i][0],bonds[j][0]])
          if bonds[i][1]==bonds[j][0]:
              angles.append([bonds[i][0],bonds[i][1],bonds[j][1]])
          if bonds[i][1]==bonds[j][1]:
              angles.append([bonds[i][0],bonds[i][1],bonds[j][0]])
  return angles

# Routine to find dihedrals based on the angles already identified
def findDihedrals(angles):
  dihedrals = []
  for i in range(0,len(angles)):
      for j in range(i+1,len(angles)):
          if angles[i][1]==angles[j][0] and angles[i][2]==angles[j][1]:
              dihedrals.append([angles[i][0],angles[i][1],angles[i][2],angles[j][2]])
          if angles[i][1]==angles[j][2] and angles[i][2]==angles[j][1]:
              dihedrals.append([angles[i][0],angles[i][1],angles[i][2],angles[j][0]])
          if angles[i][1]==angles[j][0] and angles[i][0]==angles[j][1]:
              dihedrals.append([angles[i][2],angles[j][0],angles[j][1],angles[j][2]])
          if angles[i][1]==angles[j][2] and angles[i][0]==angles[j][1]:
              dihedrals.append([angles[i][2],angles[j][2],angles[j][1],angles[j][0]])
  return dihedrals

#--------------------------------------------------------------
# The main part of the program starts here
#--------------------------------------------------------------

# Print GPL v3 statement
print "WellFAReFF  Copyright (C) 2015 Matthias Lein"
print "This program comes with ABSOLUTELY NO WARRANTY"
print "This is free software, and you are welcome to redistribute it"
print "under certain conditions."
print

# Determine the name of the file to be read (outfile currently unused)
infile, outfile = iofiles(sys.argv[1:])

# Extract coordinates from file, then find bonds (based on distance or bond order)
# then determine angles and dihedrals based on bonds and angles respectively.
geometry = extractCoordinates(infile)
bonds = findBonds(geometry,0)
angles = findAngles(bonds)
dihedrals = findDihedrals(angles)

# Print coordinates in a nice, readable way
print
for i in range(0,len(geometry)):
  print "{:<3} {: .8f} {: .8f} {: .8f}".format(geometry[i][0], geometry[i][1], geometry[i][2], geometry[i][3])

# Print bonds in a nice, readable way
print
for i in range(0,len(bonds)):
    print geometry[bonds[i][0]][0],"-",geometry[bonds[i][1]][0]," bond (%.4f)" % calcDistance(geometry[bonds[i][0]],geometry[bonds[i][1]])," between atoms ",bonds[i][0]," and ", bonds[i][1]

# Print angles in a nice, readable way
print
for i in range(0,len(angles)):
    print geometry[angles[i][0]][0],"-",geometry[angles[i][1]][0],"-",geometry[angles[i][2]][0]," angle between atoms ",angles[i][0],angles[i][1],angles[i][2], " ({:7.3f})".format(calcAngle(geometry[angles[i][0]],geometry[angles[i][1]],geometry[angles[i][2]]))

# Print dihedrals in a nice, readable way
print
for i in range(0,len(dihedrals)):
    print geometry[dihedrals[i][0]][0],"-",geometry[dihedrals[i][1]][0],"-",geometry[dihedrals[i][2]][0],"-",geometry[dihedrals[i][3]][0]," dihedral between atoms ",dihedrals[i][0],dihedrals[i][1],dihedrals[i][2],dihedrals[i][3]," ({:8.3f})".format(calcDihedral(geometry[dihedrals[i][0]],geometry[dihedrals[i][1]],geometry[dihedrals[i][2]],geometry[dihedrals[i][3]]))
