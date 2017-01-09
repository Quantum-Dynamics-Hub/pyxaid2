#***********************************************************
# * Copyright (C) 2013 Alexey V. Akimov
# * This file is distributed under the terms of the
# * GNU General Public License as published by the
# * Free Software Foundation; either version 3 of the
# * License, or (at your option) any later version.
# * http://www.gnu.org/copyleft/gpl.txt
#***********************************************************/

# This module defines a set of auxiliary functions on arrays
# Arrays are represented as nested lists

import sys
import os
import os.path

def add_arrays(A, B):
# Here we assume that A and B have the same dimensions
# A, B = T x X - dimensions
    T = len(A)
    X = len(A[0])

    C = []
    t = 0
    while t<T:
        c = []
        j = 0
        while j<X:
            c.append( A[t][j] + B[t][j] )
            j = j + 1
        C.append(c)
        t = t + 1
    return C

def mult_arrays(A, B):
# Here we assume that A and B have the same dimensions
# A, B = T x X - dimensions
    T = len(A)
    X = len(A[0])

    C = []
    t = 0
    while t<T:
        c = []
        j = 0
        while j<X:
            c.append( A[t][j] * B[t][j] )
            j = j + 1
        C.append(c)
        t = t + 1
    return C

def sum_mult_arrays(A, B):
# Here we assume that A and B have the same dimensions
# A, B = T x X - dimensions
    T = len(A)
    X = len(A[0])

    C = []
    t = 0
    while t<T:
        c = [0.0]
        j = 0
        while j<X:
            c[0]  = c[0] + ( A[t][j] * B[t][j] )
            j = j + 1
        C.append(c)
        t = t + 1
    return C


def contract_array(A,Templ):
# A - input, uncontracted, of size T x X
# Templ - template for contractions, 2D list of possible contractions
# each element of T (list by itself) is a list of column indices 
# of A to be summed together in one column


    C = [] # Contracted array
    T = len(A)        # Time
    X = len(A[0])     # Number of uncontracted states
    X_contr = len(Templ)  # Number of contracted states

    t = 0
    while t<T:
        c = []
        j = 0
        while j<X_contr:  # Each possible contraction
            cc = 0.0
            i = 0
            while i<len(Templ[j]):      # Each element in given contraction j
                k = Templ[j][i]         # Index of uncontracted column (microstate to include in macrostate)
                cc = cc + A[t][k]
                i = i + 1
            c.append(cc)
            j = j + 1
        C.append(c)
        t = t + 1
    return C


def get_file(prefix,indx,T,X,shift,N):
# This function reads the file <prefix><indx> as a matrix of dimension
# T by X (this is what it returns), shift - gives the column offset to
# start with
    filename = prefix+str(indx)
    print filename
    Arr = []
    if os.path.exists(filename)==1:

        print "adding file "+filename
        f = open(filename,"r")
        B = f.readlines()
        f.close()

        t = 0
        while t<T:
            tmp = B[t].split()
            j = 0
            arr = []
            while j<X:
                arr.append( float(tmp[shift+N*j]) ) # reads each N-th column
                j = j + 1
            Arr.append(arr)
            t = t + 1
    return Arr


def write_array(prefix,in_ex,T,X,Arr,denom):
# Write Arr divided by denom
    print "writing "+prefix+" file..."
    f = open(prefix+str(in_ex),"w")
    t = 0
    line = ""
    while t<T:
        line = line + "time "+str(t)+" "
        j = 0
        while j<X:
            line = line + "P("+str(j)+")= "+str(Arr[t][j]/denom)+" "
            j = j + 1
        line = line +"\n"
        t = t + 1
    f.write(line)
    f.close()


