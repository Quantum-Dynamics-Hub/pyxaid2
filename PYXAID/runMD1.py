#***********************************************************
# * Copyright (C) 2013-2015 Alexey V. Akimov
# * This file is distributed under the terms of the
# * GNU General Public License as published by the
# * Free Software Foundation; either version 3 of the
# * License, or (at your option) any later version.
# * http://www.gnu.org/copyleft/gpl.txt
#***********************************************************/

import os
import sys
from pyxaid_core import *

# Fisrt, we add the location of the library to test to the PYTHON path
if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *







def get_value(params,key,default,typ):
# Function to extract parameter from the dictionary
    # Try to get value from the dictionary
    str_val = "None"
    if key in params:
        if params[key]!=None:
            str_val = params[key]

    # If nothing found - use default value
    if str_val!="None":
        pass  # All is ok
    else: 
        str_val = default
        print "Warning: Parameter with key = %s does not exist in dictionary" % key
        print "Using the default value of %s" % default

    # Convert string to desired data type
    if typ=="s":
        return str_val
    elif typ=="f":
        return float(str_val)
    elif typ=="i":
        return int(float(str_val))


def runMD(params):
#------------ Read the parameters -----------------
# Parameters meaning
# pp_type - pseudopotential type: US - ultra-soft, NC - norm-conserving, PAW - projector-augmented waves
# wd - working directory, where all output (working) files will be written
# rd - results directory, where all final results (energy, NAC, H', etc.) will be written by default it will be set to wd
# This MD uses corrected NAC method

    print "Starting runMD"

    # Now try to get parameters from the input
    BATCH_SYSTEM = get_value(params,"BATCH_SYSTEM","srun","s")  # either "srun" (for SLURM) or "mpirun" (for PBS)
    NP = get_value(params,"NP","1","i")
    EXE = get_value(params,"EXE","","s")
    EXE_EXPORT = get_value(params,"EXE_EXPORT","","s")
    EXE_CONVERT = get_value(params,"EXE_CONVERT","","s")  # this is the path to iotk executable
    start_indx = get_value(params,"start_indx","0","i")
    stop_indx = get_value(params,"stop_indx","1","i")
    dt = get_value(params,"dt","1.0","f") # time step in fs - rescale NAC if actual dt is different
    dt = 41.34145 * dt # convert to a.u., so the NACs are in a.u.
    pp_type = get_value(params,"pp_type","NC","s")
    wd = get_value(params,"wd","wd","s")
    rd = get_value(params,"rd",wd,"s")
    minband = get_value(params,"minband",1,"i")
    maxband = get_value(params,"maxband",2,"i")
#    nocc = get_value(params,"nocc",1,"i")
#    nac_method = get_value(params,"nac_method",0,"i")  # choose what method for NAC calculations to use: 0 -standard, 1-corrected
    prefix0 = get_value(params,"prefix0","x0.scf","s")
#    prefix1 = get_value(params,"prefix1","x1.scf","s")

#    wfc_preprocess = get_value(params,"wfc_preprocess","restore","s") # variants are: normalize, complete, restore
#    do_complete = get_value(params,"do_complete",1,"i") # argument for option "restore"
    
    compute_Hprime = get_value(params,"compute_Hprime",0,"i") # transition dipole moments


    # Sanity/Convention check
    if(minband<=0): 
        print "Error: minband should be >0, current value of minband = ",minband
        sys.exit(0)
    if(minband>maxband):
        print "Error: minband must be smaller or equal to maxband. Current values: minband = ",minband," maxband = ",maxband
        sys.exit(0)
#    if(nocc>maxband):
#        print "Error: nocc must be smaller or equal to maxband. Current values: nocc = ",nocc," maxband = ",maxband
#        sys.exit(0)
#    if(nocc<minband):
#        print "Error: nocc must be larger or equal to minband. Current values: nocc = ",nocc," minband = ",minband
#        sys.exit(0)

    # Convert minband, maxband and nocc from external (QE-consistent, starting from 1) convetion
    # to the internal one (starting from 0)    
#    minband = minband - 1
#    maxband = maxband - 1
#    nocc = nocc - 1

    # Use this for nspin = 1 or 2
    act_sp1 = range(minband, maxband+1)     # min = 1, max = 2 => range(1,3) = [1,2]

    # Use this for nspin = 4
    act_sp2 = range(2*minband-1, 2*(maxband+1)-1 ) # min =1, max = 2 => range(1,5) = [1,2,3,4]


    # Initialize variables
    curr_index = start_indx - 1
    t = start_indx

    print "In runMD: current working directory for python: ",os.getcwd()
    print "In runMD: current working directory for sh:",os.system("echo pwd")

    os.system("mkdir %s" % wd)  # Create the working directory where all output files will be written
                                # results directory should already exist

    while t<=stop_indx:
        print "t= ", t

        dirname = ""
        if t==start_indx:
           print "Starting first point in this batch"
           dirname = "curr0"

        if t>start_indx:
           print "Continuing with other points in this batch"
           dirname = "next0"


        # A common block
        # Run calculations
        os.system( "%s -n %s %s < %s.%d.in > %s.%d.out" % (BATCH_SYSTEM,NP,EXE,prefix0,t,prefix0,t) )
        os.system( "%s -n %s %s < x0.exp.in > x0.exp.out" % (BATCH_SYSTEM,NP,EXE_EXPORT) )

        # Create temporary directory
        os.system("mkdir %s/%s" % (wd, dirname) )

        # Copy some results to that directory
        os.system( "mv %s.%d.out %s/%s" % (prefix0,t, wd, dirname) )
        os.system( "mv *.wfc* %s/%s" % (wd, dirname) )
        os.system( "mv x0.export %s/%s" % (wd, dirname) ) # "x0" - corresponds to x0 as a prefix in input files



        # Now general part - from current and next wavefunctions calculate NACs:
        if curr_index>=start_indx:
            print "Generate NAC from WFCs at two adjacent points"

            # Convert binary files to xml - this is needed in newest version of QE
            # becuase pw_export will produce binary files in any case (bug?)
            os.system("%s convert %s/curr0/x0.export/wfc.1 %s/curr0/x0.export/wfc.1.xml" % (EXE_CONVERT,wd,wd))
            os.system("%s convert %s/next0/x0.export/wfc.1 %s/next0/x0.export/wfc.1.xml" % (EXE_CONVERT,wd,wd))

            ngw, nbnd, nspin, gamma_only = QE_methods.read_qe_wfc_info("%s/curr0/x0.export/wfc.1" % wd , "Kpoint.1")
           
            act_space = act_sp1
            if nspin==4: 
                act_space = act_sp2

            wfc_curr = QE_methods.read_qe_wfc("%s/curr0/x0.export/wfc.1" % wd , "Kpoint.1", act_space)
            wfc_next = QE_methods.read_qe_wfc("%s/next0/x0.export/wfc.1" % wd , "Kpoint.1", act_space)

            e_curr = QE_methods.read_qe_index("%s/curr0/x0.export/index.xml" % wd, act_space )
            e_next = QE_methods.read_qe_index("%s/next0/x0.export/index.xml" % wd, act_space )


            os.system("-rf %s/curr0/x0.export/wfc.1.xml" % wd)
            os.system("-rf %s/next0/x0.export/wfc.1.xml" % wd)


            #-----------------------------------------------------------------
            # Finally compute Hamiltonian and the overlap matrix
            # In this case - any reasonable value for nocc leads to the same results
            # Keep in mind, the curr_wfc0 - is already only a subset of the whole wfc that has been 
            # computed, so all bands - from 0 to maxband - minband will be active!

            S, H = None, None

            if nspin==1 or nspin==2:

                ovlp  = wfc_curr.H() * wfc_next
                H = 0.5*(e_curr + e_next) - (0.5j/dt)*(ovlp - ovlp.H())
                S = 0.5 *(wfc_curr.H() * wfc_curr + wfc_next.H() * wfc_next) # for debug


            elif nspin==4:

                S = wfc_curr.H() * wfc_curr
                St  = wfc_curr.H() * wfc_next  # overlap of wfc at different times

                sx = S.num_of_cols
                ovlp = CMATRIX(sx/2, sx/2)
        ##### Pauli matrices ###
        #
        #        | 0  1 |         | 0  -i |         | 1   0 |
        #  sig1 =|      |  sig2 = |       |  sig3 = |       | 
        #        | 1  0 |         | i   0 |         | 0  -1 |
        #
        ######
 
                sig1 = CMATRIX(sx/2, sx/2)
                sig2 = CMATRIX(sx/2, sx/2)
                sig3 = CMATRIX(sx/2, sx/2)

                nac = CMATRIX(sx/2, sx/2)
                ec = CMATRIX(sx/2, sx/2)
                en = CMATRIX(sx/2, sx/2)

                for n in xrange(sx/2):
                    for k in xrange(sx/2):
                        ovlp.set(n,k, S.get(2*n,2*k) + S.get(2*n+1,2*k+1) )
                        sig1.set(n,k, S.get(2*n,2*k+1) + S.get(2*n+1,2*k) )
                        sig2.set(n,k, (-S.get(2*n,2*k+1) + S.get(2*n+1,2*k))*(1.0j+0.0) )
                        sig3.set(n,k, S.get(2*n,2*k) - S.get(2*n+1,2*k+1) )

                        nac.set(n,k, St.get(2*n,2*k).real + St.get(2*n+1,2*k+1).real, 0.0 )
                    ec.set(n,n, 0.5*(e_curr.get(2*n, 2*n)+e_curr.get(2*n+1, 2*n+1)) )
                    en.set(n,n, 0.5*(e_next.get(2*n, 2*n)+e_next.get(2*n+1, 2*n+1)) )

                H = 0.5*(ec + en) - (0.5j/dt)*(nac - nac.H())
                S = ovlp

                sig1.real().show_matrix("%s/0_sig1_%d_re" % (rd, curr_index) )
                sig1.imag().show_matrix("%s/0_sig1_%d_im" % (rd, curr_index) )
                sig2.real().show_matrix("%s/0_sig2_%d_re" % (rd, curr_index) )
                sig2.imag().show_matrix("%s/0_sig2_%d_im" % (rd, curr_index) )
                sig3.real().show_matrix("%s/0_sig3_%d_re" % (rd, curr_index) )
                sig3.imag().show_matrix("%s/0_sig3_%d_im" % (rd, curr_index) )



 
            H.real().show_matrix("%s/0_Ham_%d_re" % (rd, curr_index) )
            H.imag().show_matrix("%s/0_Ham_%d_im" % (rd, curr_index) )    
            S.real().show_matrix("%s/0_S_%d_re" % (rd, curr_index) )
            S.imag().show_matrix("%s/0_S_%d_im" % (rd, curr_index) )

            
#            ovlp.real().show_matrix("%s/0_ovlp_%d_re" % (rd, curr_index) )
#            ovlp.imag().show_matrix("%s/0_ovlp_%d_im" % (rd, curr_index) )



  
            if compute_Hprime==1:
#                os.system("%s convert %s/curr0/x0.export/grid.1 %s/curr0/x0.export/grid.1.xml" % (EXE_CONVERT,wd,wd))
#                curr_wfc0.QE_read_acsii_grid("%s/curr0/x0.export/grid.1.xml" % wd)
#                curr_wfc0.compute_Hprime(0,maxband-minband,"%s/0_Hprime_%d" % (rd, curr_index) )
                os.system("-rf %s/curr0/x0.export/grid.1.xml" % wd)


            #-----------------------------------------------------------------

            # Remove current run, make next run to be the current one
            os.system("rm -rf %s/curr0" % wd )
            os.system("mv %s/next0 %s/curr0" % (wd, wd) )

            print "old files deleted, new have become old"


# ACHTUNG!!! Restoring wfc makes some complications, so we might need to destroy wfc objects
# after each round of operations and create new objects from the beginning - thia may be safer!

        curr_index = curr_index + 1
               
        print "End of step t=", t
        t = t + 1

#================= End of runMD function =============================
