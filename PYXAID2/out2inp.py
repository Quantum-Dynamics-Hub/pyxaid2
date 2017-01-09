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

def out2inp(out_filename,templ_filename,wd,prefix,t0,tmax,dt):
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
    start = 0 # flag to start reading the coordinates
    t = 0     # index of the input file
    at_line = 0 # line with the coordinates of at_line-th atom is being written

    # Read the file
    if verbose==1:
        print "Reading file", out_filename
    f = open(out_filename,"r")

    f_t = open("%s/tmp" % wd, "w")
    f_t.close()

    for a in f:
        #print t
#        print "is_nat=%d, nat=%d, at_line=%d, start=%d, t=%d" % (is_nat,nat,at_line,start,t)
        if start==1:
            if at_line<=nat:
                f_t.write(a)
            else:
                # A blank line just in case
                f_t.write("\n")
                f_t.close()
                t = t + 1
                start = 0
            at_line = at_line + 1

        if start==2:
            if at_line<nat:
                a_tmp = a.split()
                f_t.write(a_tmp[1]+"   "+a_tmp[6]+"   "+a_tmp[7]+"  "+a_tmp[8]+"\n")
            else:
                # A blank line just in case
                f_t.write("\n")
                f_t.close()
                t = t + 1
                start = 0
            at_line = at_line + 1


        elif start==0:
            if is_nat==0:            
                if a.find("number of atoms/cell")!=-1:
                    tmp = a.split()
                    nat = int(float(tmp[4]))
                    is_nat = 1
            else:
                if a.find("ATOMIC_POSITIONS")!=-1:
                    if t>=t0 and t<=tmax:
                        if math.fmod(t-t0,dt)==0:
                            f_t = open("%s/%s.%d.in" % (wd,prefix,t), "w")
                            # Write the template header
                            i = 0
                            while i<T_sz:
                                f_t.write(T[i])
                                i = i + 1
                            # Write the header for positions
                            f_t.write(a)
                            # Set flag to write positions
                            start = 1
                            at_line = 0
                        else:
                            t = t + 1
                    elif t>tmax:
                        break;
                    else:
                        t = t + 1

                elif a.find("site n.     atom                  positions (alat units)")!=-1:
                    # Set flag to write positions
                    #print "Atomic_positions record is found! nframe=%d t=%d start=%d" % (nframe,t,start)
                    if t>=t0 and t<=tmax:
                        if math.fmod(t-t0,dt)==0:
                            f_t = open("%s/%s.%d.in" % (wd,prefix,t), "w")
                            # Write the template header
                            i = 0
                            while i<T_sz:
                                f_t.write(T[i])
                                i = i + 1
                            # Write the header for positions
                            f_t.write("ATOMIC_POSITIONS (alat)\n")
                            # Set flag to write positions
                            start = 2
                            at_line = 0

                        else:
                            t = t + 1
                    elif t>tmax:
                        break;
                    else:
                        t = t + 1


    f.close()
   


