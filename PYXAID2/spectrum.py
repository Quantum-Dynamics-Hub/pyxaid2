#***********************************************************
# * Copyright (C) 2013 Alexey V. Akimov
# * This file is distributed under the terms of the
# * GNU General Public License as published by the
# * Free Software Foundation; either version 3 of the
# * License, or (at your option) any later version.
# * http://www.gnu.org/copyleft/gpl.txt
#***********************************************************/

from numpy import *
import math


def get_data(filename,line_start,line_end,col):
# filename - is the name of the file condaining the data
# line_start - the number of the first line to include in data set (1 - for the first line)
# line_end - the number of the last line to include in data set
# col - is the number (1 -for the first column) of the column that contains C values

    # Read the file
    f = open(filename,"r")
    A = f.readlines()
    f.close()

    # convert to indices
    line_start = line_start - 1
    line_end = line_end - 1
    col = col - 1

    # Calculate the average of the input
    X = []
    i = line_start
    while i<=line_end:
        tmp = A[i].split()
        x = float(tmp[col])
        X.append(x)
        i = i + 1

    return list(X)



def do_fft(x,dt,filename):
# Do the FFT of the data x with the time step dt
# The result is written to the <filename>

    f = open(filename,"w")
    signal = abs(fft.fft(x))
    freq   = fft.fftfreq(len(x),dt)

    A = []
    B = []
    C = []

    sz = len(freq)/2
    i = 0
    while i<sz:
        a = 2.0*math.pi*freq[i]
        b = 2.0*math.pi*freq[i]*(100000.0/3.0)
        c = signal[i]
        A.append(a)
        B.append(b)
        C.append(c)
        if i==0:
            c = 0
        f.write(str(a)+" 1/fs "+str(b)+" 1/cm "+str(c)+"\n")
        # transform Nyquist frequency to angular frequency (factor 2)
        # and then transform angular freq. to normal freq. (factor math.pi)
        # and then show it in units of 1/cm
        i = i + 1

    f.close()

    # Return data written to the file
    return A, B, C


# Example of usage
#from pyxais import *
# Read data from the file
#x = spectrum.get_data("autocorr.dat",1)
# Compute spectrum - fft of the signal
#spectrum.do_fft(x,1.0,"autocorr_fft.dat")

