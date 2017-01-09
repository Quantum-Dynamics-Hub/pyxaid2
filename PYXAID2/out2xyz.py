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

def cryst2cart(a1,a2,a3,r):
# auxiliary function
# convert vertor <r> in crystal (alat) coordinates with the cell defined by
# vectors a1, a2, a3, to the Cartesian coordinate xyz
    x = [0.0,0.0,0.0]
    for i in [0,1,2]:
        x[i] = a1[i]*r[0] + a2[i]*r[1] + a3[i]*r[2]

    return x


def convert(out_filename,T,dt,xyz_filename):
# out_filename - name of the file which contains the MD trajectory
# this file is the QE MD trajectory output
# The function will convert this output file into xyz file (containing trajectory)
# No more than T steps from the out_filename file will be used
# dt - number of steps between output frames, so dt = 5 will output frames 0, 5, 10, 15, etc.
# xyz_filename - is the prefix of the file to which the result is written

#    dt = dt - 1
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

    fr = open(xyz_filename,"w")
    fr.close()
    fr = open(xyz_filename,"w")
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
                fr.write("%s   %5.3f  %5.3f  %5.3f \n"  %(symb,scl*r[0],scl*r[1],scl*r[2]))
            else:
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
#                   print "Atomic_positions record is found! nframe=%d t=%d start=%d" % (nframe,t,start)
                    if t<T:
                        if nframe==dt:
#                            print "Start of output"
                            fr.write("%d\n" % nat)
                            fr.write("molecule\n")

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
#                            print "Start of output"
                            fr.write("%d\n" % nat)
                            fr.write("molecule\n")

                            start = 2
                            at_line = 0
                        else:
                            start = 0
                            nframe = nframe + 1
                            t = t + 1
                    else:
                        break


    fr.close()
    f.close()
   

# Example of usage
#> from PYXAID import*
#> out2xyz.convert("x.md.out",250,25,"snaps/traj.xyz")
# This will create the MD trajectory file in .xyz format with the snapshots takes at times 0
# (input configuration), 25 (25-th nuclear configuration), 50, etc.
# the snapshots will written in the file trahj.xyz in the folder /snaps 


