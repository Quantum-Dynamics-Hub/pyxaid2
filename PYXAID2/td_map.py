#***********************************************************
# * Copyright (C) 2013 Alexey V. Akimov
# * This file is distributed under the terms of the
# * GNU General Public License as published by the
# * Free Software Foundation; either version 3 of the
# * License, or (at your option) any later version.
# * http://www.gnu.org/copyleft/gpl.txt
#***********************************************************/

def make_map(enfile,popfile,outfile,do_sort):
# enfile  - the na-md output file, containing the energies of the states
# popfile -    - // -            , containing the populations (either se or sh) of the states
# outfile - the name of the output file, where the result (map) is written
# do_sort - flag, defining if the energies must be sorted (= 1) or not (otherwise)
#           should be set to 1 in order to produce reasonable map when result is plotted

    # Reading input files
    f = open(enfile,"r")
    A = f.readlines()
    f.close()

    f1 = open(popfile,"r")
    B = f1.readlines()
    f1.close()


    # Determining dimensions
    T = len(A) # number of time steps
    x = len(A[0].split())  # number of columns
    Nst = (x-4)/2  # number of states

    print "T = ", T
    print "Number of states = ", Nst

    # Processing input files
    E = []  # energies
    P = []  # populations
    t = 0
    while t<T:  # loop over time
        e = []
        p = []

        st = 0
        tmpe = A[t].split()
        tmpp = B[t].split()
        while st<Nst:  # loop over states
            val_e = float(tmpe[5+2*st])
            val_p = float(tmpp[3+2*st])
            e.append(val_e)
            p.append(val_p)
            st = st + 1
        E.append(e)
        P.append(p)
        
        t = t + 1

    # Placeholder for some math
    # Sort energy for each time step (simple,but inefficient sorting scheme)
    if(do_sort==1):
        t = 0
        while t<T:
            i = 0
            while i<(Nst-1):               
                emin = E[t][i]
                indx = i

                j = i+1
                while j<Nst:
                    if E[t][j]<emin:
                        emin = E[t][j]
                        indx = j
                    j = j + 1
                # Swap elements i and indx
                # don't forget to swap populations simulataneousely with the swaps of the energy scale
                se = E[t][i]
                E[t][i] = E[t][indx]
                E[t][indx] = se

                sp = P[t][i]
                P[t][i] = P[t][indx]
                P[t][indx] = sp

                i = i + 1
            t = t + 1
        


    # Print the 2D map
    f2 = open(outfile,"w")
    t = 0
    while t<T:
        st = 0
        while st<Nst:
            f2.write("%f  %f  %10.8f \n" % (float(t),E[t][st],P[t][st]))
            st = st + 1
        f2.write("\n")
        t = t + 1
    f2.close()


# Example of usage
#make_map("me_energies0","se_pop_ex78","se.dat",1)
#make_map("me_energies0","sh_pop_ex78","sh.dat",1)


