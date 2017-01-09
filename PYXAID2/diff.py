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


def diff(out_filename,t0,T,at,res_filename,opt):
# out_filename - name of the file which contains the MD trajectory
# at - is an array (list) of atom lists for which we want to compute diffusion constants (starting from 1, not from 0)
# t0 initial frames are skipped
# T - is the maximal time frame which will be counted (if if exists in the MD trajectory)
# res_filename - is where the results are written
# opt - option: = 0 - the diffusion coefficient is calculated
#               = 1 - the "reaction-coordinate" is calculates (mean-square displacement per time step)

    print "Running diff function in module diff"
    print "Arguments are:"
    print " out_filename=",out_filename
    print " t0 = ",t0
    print " T = ", T
    print " at = ", at
    print " res_filename = ", res_filename
    print " opt = ", opt

    sz = len(at)  # number of groups of atoms
    Rprev = []
    Rt = []
    R0 = []
    D  = []
    for i in range(0,sz):
        rprev = []
        rt = []
        r0 = []
        d  = []
        for j in range(0,len(at[i])):
            rprev.append([0.0,0.0,0.0])
            rt.append([0.0,0.0,0.0])
            r0.append([0.0,0.0,0.0])
            d.append([0.0,0.0,0.0])
        Rprev.append(rprev)
        Rt.append(rt)
        R0.append(r0)
        D.append(d)
         

    #print D
    Dtot = range(0,sz)  # total diffusion coefficient for each group
    i = 0
    while i<sz:
        Dtot[i] = 0.0
        i = i + 1

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
    f = open(out_filename,"r")

    t = 0  # time (configuration)

    convert = 0.1  # 1 A^2/fs = convert cm^2/s
    convert = convert * 100000  # in units 10^-5

    indx = range(0,sz)  # index of the currently processing position in the internal indexing of small arrays R0,Rt,D
                        # it is an analog of the at_line in global indexing, for i-th group of atoms
    i = 0
    while i<sz:
        indx[i] = 0
        i = i + 1

    res = open(res_filename,"w")
    res.close()
    res = open(res_filename,"w")

    _D = []  # Diffusion coefficients to return

    for a in f:
        if t>=T:
            break
        if start==1 and t<T:
            if at_line<nat:
                scl  = alat * 0.52918  # alat in a.u., coordinates will be in Angstrom
                tmp  = a.split()
                symb = tmp[0]
                r    = [scl*float(tmp[1]),scl*float(tmp[2]),scl*float(tmp[3])]
                if t<=t0:
                    i = 0
                    while i<sz:  # run through all groups of atoms
                        if (at_line+1) in at[i]:  # atoms in i-th group
                            Rprev[i][indx[i]] = r[:]
                            Rt[i][indx[i]] = r[:]
                            R0[i][indx[i]] = r[:]
                            indx[i] = indx[i] + 1
                        i = i + 1
                else:
                    i = 0
                    while i<sz: # run through all groups of atoms
                        if (at_line+1) in at[i]:  # Found specific atoms
                            Rt[i][indx[i]] = r[:]

                            for j in [0,1,2]:  # all components x, y, z
                                if opt==0:
                                    D[i][indx[i]][j] = 0.5*convert*((Rt[i][indx[i]][j] - R0[i][indx[i]][j])**2 / float(t-t0) ) 
                                elif opt==1:                        
                                    D[i][indx[i]][j] = 0.5*convert*((Rt[i][indx[i]][j] - Rprev[i][indx[i]][j])**2 / float(1.0) )
                                Dtot[i] = Dtot[i] + D[i][indx[i]][j]

                            indx[i] = indx[i] + 1
                        i = i + 1
            else:
                # Print D:
                if t-t0>=0.0:
                    line = str(t-t0)+"  "
                    i = 0
                    _d = range(0,sz)
                    while i<sz:
                        _d[i] = Dtot[i]/(3.0*float(len(at[i])))
                        line = line + str( _d[i] )+"  "
                        i = i + 1
                    line = line + "\n"
                    res.write(line)
#                    res.write("%15.10f    %15.10f\n" % (t-t0, Dtot/(3.0*float(sz))) )
                    _D.append(_d)
                i = 0
                while i<sz:
                    Dtot[i] = 0.0
                    # Current becomes previous:
                    Rprev[i] = list(Rt[i])
                    i = i + 1

                t = t + 1
                start = 0
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
                    start = 1
                    at_line = 0
                    i = 0
                    while i<sz:
                        indx[i] = 0
                        i = i + 1
                    if t<T:
                        pass
                    else:
                        break

    res.close()
    f.close()

    return list(_D)


def diff_ave(out_filename,T0,Tmax,dT,at,res_filename,opt):
# This is the same function as diff, but it also averages the comuted value
# over different starting points. Basically this is to split a single 
# trajectory into pieces and average over initial conditions

    D = []
    t0 = T0
    T = t0 + dT
    count = 0
    sz0 = 0
    ngrp = 0

    while T<=Tmax:
        _D = diff(out_filename,t0,T,at,"tmp.diff",opt)
        if count==0:
            sz0 = len(_D)  # length of this set
            ngrp = len(_D[0]) # number of groups
            D = list(_D)
            count = 1
        else:
            sz = len(_D)
            if sz==sz0:
                i = 0
                while i<sz0:
                    j = 0
                    while j<ngrp:
                        D[i][j] = D[i][j] + _D[i][j]
                        j = j + 1
                    i = i + 1
                count = count + 1

        t0 = t0 + dT
        T = T + dT

    res = open(res_filename,"w")
    i = 0
    while i<sz0:
        line = str(i)+"  "
        j = 0
        while j<ngrp:
            D[i][j] = D[i][j] / float(count)
            line = line + str(D[i][j])+"  "
            j = j + 1
        line = line +"\n"
        res.write(line)
#        res.write("%d  %15.10f\n" % (i,D[i]))
        i = i + 1

    res.close()




# Example of usage
#at = []
#for i in range(0,12):
#    at.append(26+3*i)
#    at.append(27+3*i)
# At = [at]
#diff("x.md.out",100,1000,At,"diff.dat",1)
#diff_ave("x.md.out",100,50000,1000,At,"diff0_ave.dat",0)
