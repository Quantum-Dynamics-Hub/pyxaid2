#***********************************************************
# * Copyright (C) 2013 Alexey V. Akimov
# * This file is distributed under the terms of the
# * GNU General Public License as published by the
# * Free Software Foundation; either version 3 of the
# * License, or (at your option) any later version.
# * http://www.gnu.org/copyleft/gpl.txt
#***********************************************************/

import os
import sys
import math

def cryst2cart(a1,a2,a3,r):
# auxiliary function
# convert vertor <r> in crystal (alat) coordinates with the cell defined by
# vectors a1, a2, a3, to the Cartesian coordinate xyz
    x = [0.0,0.0,0.0]
    for i in [0,1,2]:
        x[i] = a1[i]*r[0] + a2[i]*r[1] + a3[i]*r[2]

    return x


def convert(out_filename,T,dt,pdb_prefix):
# out_filename - name of the file which contains the MD trajectory
# this file is the QE MD trajectory output
# The function will convert this output file into pdb file (containing trajectory)
# No more than T steps from the out_filename file will be used
# dt - the difference of indexes of the frames which are written consequetively
# such that if you dt = 5 it will write frames 0,5,10,15,etc. with 0 - being the input configuration

    dt = dt - 1
    verbose = 0

    # Parameters
    nat = 0   # Number of atoms per unit cell
    is_nat = 0# flag defining if nat variable has been set
    is_a1 = 0 # flag defining if the cell variables are set
    is_a2 = 0 # flag defining if the cell variables are set
    is_a3 = 0 # flag defining if the cell variables are set
    is_alat=0 # flag defining scaling constant for the cell
    start = 0 # flag to start reading the coordinates
    t = 0     # index of the input file
    at_line = 0 # line with the coordinates of at_line-th atom is being written

    a1 = [0.0, 0.0, 0.0]
    a2 = [0.0, 0.0, 0.0]
    a3 = [0.0, 0.0, 0.0]

    # Read the file
    if verbose==1:
        print "Reading file", out_filename
    f = open(out_filename,"r")

    fr = open("tmp","w")
    fr.close()

    t = 0  # time (configuration)

    nframe = dt

    for a in f:
        if (start==1 or start==2) and t<T and nframe==dt:
            if at_line<nat:
                tmp  = a.split()
                symb = tmp[0]
                r    = [0.0, 0.0, 0.0]
                if(start==1):
                    r    = [float(tmp[1]),float(tmp[2]),float(tmp[3])]
                    symb = tmp[0]
                elif(start==2):
                    r    = [float(tmp[6]),float(tmp[7]),float(tmp[8])]
                    symb = tmp[1]

                scl  = alat * 0.52918  # alat in a.u., coordinates will be in Angstrom
                name = ""
                elt = ""
                if len(symb)==1:
                    name = "   "+symb
                    elt = " "+symb
                elif len(symb)==2:
                    name = "  "+symb
                    elt = symb
                else:
                    name = "    "
             
                fr.write("ATOM  %5.d %s MOL M%4.d    %8.3f%8.3f%8.3f%6.2f%6.2f      XYZ1%s  \n"  % (at_line+1,name,1,scl*r[0],scl*r[1],scl*r[2],1.0,0.0,elt))
            else:
                fr.close()
                t = t + 1
                start = 0
                nframe = 0
            at_line = at_line + 1

        elif start==0:
            if is_nat==0:            
                if a.find("number of atoms/cell")!=-1:
                    tmp = a.split()
                    nat = int(float(tmp[4]))
                    #print "nat = %5d" % nat
                    is_nat = 1
            if is_a1==0:
                if a.find("a(1) =")!=-1:
                    tmp = a.split()
                    a1 = [float(tmp[3]),float(tmp[4]),float(tmp[5])]
                    #print "a1 = ", a1
                    is_a1 = 1
            if is_a2==0:
                if a.find("a(2) =")!=-1:
                    tmp = a.split()
                    a2 = [float(tmp[3]),float(tmp[4]),float(tmp[5])]
                    #print "a2 = ", a2
                    is_a2 = 1
            if is_a3==0:
                if a.find("a(3) =")!=-1:
                    tmp = a.split()
                    a3 = [float(tmp[3]),float(tmp[4]),float(tmp[5])]
                    #print "a3 = ", a3
                    is_a3 = 1
            if is_alat==0:
                if a.find("lattice parameter (alat)  =")!=-1:
                    tmp = a.split()
                    alat = float(tmp[4])
                    #print "alat = ", alat
                    is_alat = 1

            else:
                if a.find("ATOMIC_POSITIONS")!=-1:
                    # Set flag to write positions
                    #print "Atomic_positions record is found! nframe=%d t=%d start=%d" % (nframe,t,start)
                    if t<T:
                        if nframe==dt:
                            fr = open("%s%d.pdb" % (pdb_prefix,t), "w")
                            A = math.sqrt(a1[0]**2 + a1[1]**2 + a1[2]**2)
                            B = math.sqrt(a2[0]**2 + a2[1]**2 + a2[2]**2)
                            C = math.sqrt(a3[0]**2 + a3[1]**2 + a3[2]**2)
                            AB = a1[0]*a2[0] + a1[1]*a2[1] + a1[2]*a2[2]
                            AC = a1[0]*a3[0] + a1[1]*a3[1] + a1[2]*a3[2]
                            BC = a2[0]*a3[0] + a2[1]*a3[1] + a2[2]*a3[2]
                            gam = math.degrees(math.acos(AB/(A*B)))
                            bet = math.degrees(math.acos(AC/(A*C)))
                            alp = math.degrees(math.acos(BC/(B*C)))

                            A = A*alat*0.52918
                            B = B*alat*0.52918
                            C = C*alat*0.52918
                            fr.write("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n" % (A,B,C,alp,bet,gam))

                            start = 1
                            at_line = 0
                        else:
                            start = 0
                            nframe = nframe + 1
                            t = t + 1
                    else:
                        break

                elif a.find("site n.     atom                  positions (alat units)")!=-1:
                    # Set flag to write positions
                    #print "Atomic_positions record is found! nframe=%d t=%d start=%d" % (nframe,t,start)
                    if t<T:
                        if nframe==dt:
                            fr = open("%s%d.pdb" % (pdb_prefix,t), "w")
                            A = math.sqrt(a1[0]**2 + a1[1]**2 + a1[2]**2)
                            B = math.sqrt(a2[0]**2 + a2[1]**2 + a2[2]**2)
                            C = math.sqrt(a3[0]**2 + a3[1]**2 + a3[2]**2)
                            AB = a1[0]*a2[0] + a1[1]*a2[1] + a1[2]*a2[2]
                            AC = a1[0]*a3[0] + a1[1]*a3[1] + a1[2]*a3[2]
                            BC = a2[0]*a3[0] + a2[1]*a3[1] + a2[2]*a3[2]
                            gam = math.degrees(math.acos(AB/(A*B)))
                            bet = math.degrees(math.acos(AC/(A*C)))
                            alp = math.degrees(math.acos(BC/(B*C)))

                            A = A*alat*0.52918
                            B = B*alat*0.52918
                            C = C*alat*0.52918
                            fr.write("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n" % (A,B,C,alp,bet,gam))

                            start = 2
                            at_line = 0
                        else:
                            start = 0
                            nframe = nframe + 1
                            t = t + 1
                    else:
                        break
           
    f.close()
   
# Example of usage
#> from PYXAID import*
#> out2pdb.convert("x.md.out",250,25,"snaps/snap_")
# This will create MD snapshots at times 0 (input configuration), 25 (25-th nuclear configuration), 50, etc.
# the files will be collcted in the folder /snaps and named snap_0.pdb, snap_25.pdb, snap_50.pdb, etc.

