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

def xdatcar2inp(out_filename,templ_filename,wd,prefix,t0,tmax,dt):
# out_filename - name of the file which contains the MD trajectory
# templ_filename - name of the template file for input generation, should not contain atomic positions!
# prefix - is the prefix of the files generated at output
# wd - working directory - will be created 
# t0 and t1 - define the starting and final frames
# dt - defines the spacing between frames which are written
# this is defined as a difference between written configuration indexes:
# so if dt = 5, the following frames will be written: 0,5,10,15, etc...

    verbose = 0
    # Read the template file
    f_templ = open(templ_filename,"r")
    T = f_templ.readlines()
    f_templ.close()
    T_sz = len(T)

    if os.path.isdir(wd):
        pass
    else:
        # Create working directory and generate the files
        os.system("mkdir %s" % wd)

    # Parameters
    nat = 0   # Number of atoms per unit cell
    is_nat = 0# flag defining if nat variable has been set
    start = 0 # flag to start printing coordinates
    t = 0     # index of the input file
    at_line = 0 # line with the coordinates of at_line-th atom is being written

    # Read the file
    if verbose==1:
        print "Reading file", out_filename

    out_file = open(out_filename,"r")
    f = out_file.readlines()
    out_file.close()


    f_t = open("%s/tmp" % wd, "w")
    f_t.close()


    # The format of XDATCAR is straight, so we can extract info sequentially
    scl = float(f[1].split()[0])   # scaling factor

    # Unit cell
    tmp = f[2].split()
    tv1 = [scl*float(tmp[0]),scl*float(tmp[1]),scl*float(tmp[2])]   

    tmp = f[3].split()
    tv2 = [scl*float(tmp[0]),scl*float(tmp[1]),scl*float(tmp[2])]   

    tmp = f[4].split()
    tv3 = [scl*float(tmp[0]),scl*float(tmp[1]),scl*float(tmp[2])]   


    types = f[5].split()   # atomic types
    ntyp = len(types)   # number of atomic types
    print types
      

    tmp = f[6].split()
    nat = []  # number of atoms of different type
    totnat = 0  # total number of atoms (all types)
    for elt in tmp:
        nelt = int(float(elt))
        nat.append(nelt)
        totnat = totnat + nelt

    print nat


    x = []
    for a in range(0,totnat):
        x.append([0.0, 0.0, 0.0])



    start = 0
    cnt = 0 # count of how many atoms have been read for given frame
    t = 0

    for a in f[8:]:
#        print a
        #======= Try to read coordinates first ======
        if start==0:

            tmp = a.split()
            if len(tmp)==3:
                x[cnt] = [float(tmp[0]),float(tmp[1]),float(tmp[2])]
                cnt = cnt + 1

            if cnt==totnat:
                start = 1
                cnt = 0        
                  

        #======= Printing files ==========
        elif start==1:
            if t>=t0 and t<=tmax:
                if math.fmod(t-t0,dt)==0:
                    f_t = open("%s/%s.%d.in" % (wd,prefix,t), "w")
                    # Write the template header
                    for b in T:
                        f_t.write(b)
        

                    # Write the header for positions        
                    f_t.write("ATOMIC_POSITIONS (alat)\n")

                    # Write positions
                    i = 0
                    for nt in range(0,ntyp):
                        for na in range(0,nat[nt]):

                            X = x[i][0]*tv1[0] + x[i][1]*tv2[0] + x[i][2]*tv3[0]
                            Y = x[i][0]*tv1[1] + x[i][1]*tv2[1] + x[i][2]*tv3[1]
                            Z = x[i][0]*tv1[2] + x[i][1]*tv2[2] + x[i][2]*tv3[2]
                            S = types[nt]
                            i = i + 1

                            f_t.write("%3s  %8.5f  %8.5f  %8.5f\n" % (S,X,Y,Z))
                    
                    f_t.write("\n")

                
            elif t>tmax:
                break;
        
            t = t + 1            
            start = 0  # printing is over, accumulate data


   

# Example of usage
#xdatcar2inp("xdatacar","x0.scf.in", "wd", "x0.scf",0,5,1)

