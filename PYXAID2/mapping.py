#*********************************************************************************
#* Copyright (C) 2017 Alexey V. Akimov
#*
#* This file is distributed under the terms of the GNU General Public License
#* as published by the Free Software Foundation, either version 2 of
#* the License, or (at your option) any later version.
#* See the file LICENSE in the root directory of this distribution
#* or <http://www.gnu.org/licenses/>.
#*
#*********************************************************************************/

import os
import sys
import math

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *


def elementary_overlap(psi1, psi2):
    """
    This function creates a matrix of elementary overlaps - the overlaps
    of the 1-electron 2-component spin-orbitals. The spin-orbitals can have an
    arbitrary meaning: either spin-adiabatic (then both alpha and beta spatial
    parts are non-zero) or spin-diabatic (then one of the components is zero)
     
    psi1, psi2 - are the lists of 2 elements, each element is CMATRIX(npw, ndim) with the
    plane-wave coefficients for the spatial part for all spin-orbitals
    ndim may be different for psi1 and psi2

    psi1[0], psi2[0] - alpha
    psi1[1], psi2[1] - beta

    """
    n, m = psi1[0].num_of_cols, psi2[0].num_of_cols

    ovlp = CMATRIX(n,m)   # <psi1 | psi2 >

    ovlp = psi1[0].H() * psi2[0] +  psi1[1].H() * psi2[1]   # Eq. 17

    return ovlp    
    
    


def adi2indx(inp):
# This mapping is almost 1-to-1, just change the indexing convention
# inp - start indexing from 1
# out - start indexing from 0

    sz = len(inp)
    out = [0] * sz

    for i in xrange(sz):
        out[i] = inp[i] - 1

    return out


def dia2indx(inp):
# This function maps a list of integers defining the spin-diabatic configuration
# on the list of indices of the corresponding spin-orbitals
# inp - indices start from 1, positive - for alpha, negative for beta
# out - indices start from 0, even - for alpha, odd -for beta
# Example: [1, -1, 2, -2, 3, -3, etc.] ->   [0, 1, 2, 3, 4, 5, etc.]
    sz = len(inp)
    out = [0] * sz

    for i in xrange(sz):
        # alpha
        if inp[i] > 0: 
            out[i] = 2 * (inp[i] - 1 )
        else:
            out[i] = 2 * (abs(inp[i]) - 1) + 1

    return out



def ovlp_arb(SD1, SD2, S):
# Overlap of two generic SDs: <SD1|SD2>
# SD1, SD2 - are the lists of indices of the spin-orbitals (in the corresponding
# active spaces) that are included in the two configurations
# S - is the matrix in the spaces of 1-el spin-orbitals (either spin-diabatic or 
# spin-adiabatic or both) 
# See Eq. 16 in the write up

    sz1,sz2 = len(SD1), len(SD2)

    s = CMATRIX(sz1,sz2)
    pop_submatrix(S, s, SD1, SD2)

    return det(s)    # Eq. 16



def ovlp_aa(SD1, SD2, S):
# Overlap of two spin-adiabatic SDs
    sd1 = adi2indx(SD1)
    sd2 = adi2indx(SD2)

    return ovlp_arb(sd1, sd2, S)

def ovlp_dd(SD1, SD2, S):
# Overlap of two spin-adiabatic SDs
    sd1 = dia2indx(SD1)
    sd2 = dia2indx(SD2)

    return ovlp_arb(sd1, sd2, S)

def ovlp_da(SD1, SD2, S):
# Overlap of two spin-adiabatic SDs
    sd1 = dia2indx(SD1)
    sd2 = adi2indx(SD2)

    return ovlp_arb(sd1, sd2, S)



# Example
#print dia2indx([1,-1,2,-2,3,-3])
