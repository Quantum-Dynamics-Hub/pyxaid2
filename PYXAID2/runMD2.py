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



def compute_overlap(info, act_space):
# \param[in] info (Python dictionary) - contains all the essential information about the calculations
# \param[in] act_space (Python list of integers) - defines the orbitals belonging to the active space
#

    nstates = len(act_space)

    if info["nspin"]==4:  # non-collinear magnetism
        nstates = nstates / 2

    S = CMATRIX(info["nk"]*nstates, info["nk"]*nstates)

    for ik1 in xrange(info["nk"]):
        for ik2 in range(ik1, info["nk"]):

            s = pw_overlap(info["k"][ik1], info["k"][ik2], coeff[ik1], coeff[ik2], grid[ik1], grid[ik2])
      
        if ik2!=ik1:
            push_submatrix(S, s, range(ik1*2, (ik1+1)*2), range(ik2*2, (ik2+1)*2))
            push_submatrix(S, s.H(), range(ik2*2, (ik2+1)*2), range(ik1*2, (ik1+1)*2))

            print ik1, ik2;  s.show_matrix()
            print ik2, ik1;  s.H().show_matrix()

        else:
            push_submatrix(S, s, range(ik1*2, (ik1+1)*2), range(ik2*2, (ik2+1)*2))

            print ik1, ik2;  s.show_matrix()
        
    print "The overall overlap matrix: "; S.show_matrix()

    return S


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
    nac_method = get_value(params,"nac_method",0,"i")  # choose what method for NAC calculations to use: 0 -standard, 1-corrected
    prefix0 = get_value(params,"prefix0","x0.scf","s")
    prefix1 = get_value(params,"prefix1","x1.scf","s")

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
           dirname1 = "curr1"

        if t>start_indx:
           print "Continuing with other points in this batch"
           dirname = "next0"
           dirname1 = "next1"


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


        if nac_method==2:  # In addition perform the spin-polarized calculation
            os.system( "%s -n %s %s < %s.%d.in > %s.%d.out" % (BATCH_SYSTEM,NP,EXE,prefix1,t,prefix1,t) )
            os.system( "%s -n %s %s < x1.exp.in > x1.exp.out" % (BATCH_SYSTEM,NP,EXE_EXPORT) )

            os.system("mkdir %s/%s" % (wd,dirname1) )

            os.system( "mv %s.%d.out %s/%s" % (prefix1,t, wd, dirname1) )
            os.system( "mv *.wfc* %s/%s" % (wd, dirname1) )
            os.system( "mv x1.export %s/%s" % (wd, dirname1) ) # "x1" - corresponds to x1 as a prefix in input files




        # Now general part - from current and next wavefunctions calculate NACs:
        if curr_index>=start_indx:
            print "Generate NAC from WFCs at two adjacent points"

            info, all_e_dum = QE_methods.read_qe_index("%s/curr0/x0.export/index.xml" % wd, [], 0)



            print "The total # of k-points is: ", info["nk"]

            # read the coefficients and energies for the mluti k-points cases, even if soc requires gamma only

            #====== Current electronic structure ===========
            # we do know the first job is soc case
            dum, e_curr = QE_methods.read_qe_index("%s/curr0/x0.export/index.xml" % wd, act_sp2, 0)
            coeff_curr = []
            grid_curr = []

            for ik in xrange(info["nk"]):
                print ik, info["k"][ik]
                coeff_curr.append( QE_methods.read_qe_wfc("%s/curr0/x0.export/wfc.%i" % (wd, ik+1), act_sp2, 0))
                grid_curr.append( QE_methods.read_qe_wfc_grid("%s/curr0/x0.export/grid.%i" % (wd, ik+1) , 0) )


            #====== Next electronic structure ===========
            dum, e_next = QE_methods.read_qe_index("%s/next0/x0.export/index.xml" % wd, act_sp2, 0)
            coeff_next = []
            grid_next = []

            for ik in xrange(info["nk"]):
                print ik, info["k"][ik]
                coeff_next.append( QE_methods.read_qe_wfc("%s/next0/x0.export/wfc.%i" % (wd, ik+1), act_sp2, 0))
                grid_next.append( QE_methods.read_qe_wfc_grid("%s/next0/x0.export/grid.%i" % (wd, ik+1) , 0) )

            if nac_method==2:
                # read coefficients and energies for spin-polarized case
                # we know wfc.1 and wfc.2 are up and down coefficients

                #====== Current electron electructure ========
                dum, e_curr_spin = QE_methods.read_qe_index("%s/curr1/x1.export/index.xml" % wd, act_sp1, 0)
                e_curr_alpha = e_curr_spin[0]
                e_curr_beta = e_curr_spin[1]

                coeff_curr_alpha = []
                coeff_curr_beta = []
                coeff_curr_alpha.append(QE_methods.read_qe_wfc("%s/curr1/x1.export/wfc.1" % wd,act_sp1,0)) 
                coeff_curr_beta.append(QE_methods.read_qe_wfc("%s/curr1/x1.export/wfc.2" % wd,act_sp1,0)) 

                #====== Next electronic structure ===========
                dum, e_next_spin = QE_methods.read_qe_index("%s/next1/x1.export/index.xml" % wd, act_sp1, 0)
                e_next_alpha = e_next_spin[0]
                e_next_beta = e_next_spin[1]

                coeff_next_alpha = []
                coeff_next_beta = []
                coeff_next_alpha.append( QE_methods.read_qe_wfc("%s/next1/x1.export/wfc.1" % wd, act_sp1, 0))
                coeff_next_beta.append( QE_methods.read_qe_wfc("%s/next1/x1.export/wfc.2" % wd, act_sp1, 0))

 

            S, H, S_aa, S_bb, S_ab, H_aa, H_ab, H_bb = None, None, None, None, None, None,None, None

            if info["nspin"]==4:

                if info["nk"]==1: # Only one k-point

                    S = coeff_curr[0].H() * coeff_curr[0]
                    St = coeff_curr[0].H() * coeff_next[0]  # overlap of wfc at different times

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
                        ec.set(n,n, 0.5*(e_curr[0].get(2*n, 2*n)+e_curr[0].get(2*n+1, 2*n+1)) )
                        en.set(n,n, 0.5*(e_next[0].get(2*n, 2*n)+e_next[0].get(2*n+1, 2*n+1)) )

                    H = 0.5*(ec + en) - (0.5j/dt)*(nac - nac.H())
                    S = ovlp

                    sig1.real().show_matrix("%s/0_sig1_%d_re" % (rd, curr_index) )
                    sig1.imag().show_matrix("%s/0_sig1_%d_im" % (rd, curr_index) )
                    sig2.real().show_matrix("%s/0_sig2_%d_re" % (rd, curr_index) )
                    sig2.imag().show_matrix("%s/0_sig2_%d_im" % (rd, curr_index) )
                    sig3.real().show_matrix("%s/0_sig3_%d_re" % (rd, curr_index) )
                    sig3.imag().show_matrix("%s/0_sig3_%d_im" % (rd, curr_index) )

                else:
                    print "Multiple k-points scheme with SOC is not yet implemented"
                    sys.exit(0)


            H.real().show_matrix("%s/0_Ham_%d_re" % (rd, curr_index) )
            H.imag().show_matrix("%s/0_Ham_%d_im" % (rd, curr_index) )    
            S.real().show_matrix("%s/0_S_%d_re" % (rd, curr_index) )
            S.imag().show_matrix("%s/0_S_%d_im" % (rd, curr_index) )

            

            info_spin = QE_methods.read_qe_wfc_info("%s/curr1/x1.export/wfc.1" % wd,0)
            if info_spin["nspin"]==2:

                print 'begin the projection of adiabatic/diabatic basis'
                
                # check the spin-polarized case has the same number of plane waves as soc
                ngw = info_spin["ngw"]
                norbs = len(act_sp1)

                coeff_curr_alpha_patch_zeros = CMATRIX(ngw,2*norbs)
                coeff_curr_beta_patch_zeros = CMATRIX(ngw,2*norbs)
                coeff_next_alpha_patch_zeros = CMATRIX(ngw,2*norbs)
                coeff_next_beta_patch_zeros = CMATRIX(ngw,2*norbs)

                # adapt the Phi_i (NO-SOC) to Phi_i (SOC) 
                # the soc coefficient are in  pw i_alp i_beta j_alp j_beta etc.
                # so for every alp, set beta to zero; for every beta, set alp to zero
                for i in xrange(norbs):
                    for pw in xrange(ngw):
                        coeff_curr_alpha_patch_zeros.set(pw,2*i,coeff_curr_alpha[0].get(pw,i))
                        coeff_curr_alpha_patch_zeros.set(pw,2*i+1,0)  #patch beta with zero

                        coeff_curr_beta_patch_zeros.set(pw,2*i,0)   #patch alpha with zero
                        coeff_curr_beta_patch_zeros.set(pw,2*i,coeff_curr_beta[0].get(pw,i))

                        coeff_next_alpha_patch_zeros.set(pw,2*i,coeff_next_alpha[0].get(pw,i))
                        coeff_next_alpha_patch_zeros.set(pw,2*i+1,0)   #patch beta with zero

                        coeff_next_beta_patch_zeros.set(pw,2*i,0) #patch alpha with zero
                        coeff_next_beta_patch_zeros.set(pw,2*i,coeff_next_beta[0].get(pw,i))

                # project the spin-adiabatic basis to spin-diabatic basis
                alpha_proj_curr = coeff_curr_beta_patch_zeros.H() * coeff_curr[0]
                beta_proj_curr = coeff_curr[0].H()*coeff_curr_beta_patch_zeros

                alpha_proj_next = coeff_next_beta_patch_zeros.H() * coeff_next[0]
                beta_proj_next = coeff_next[0].H()*coeff_next_beta_patch_zeros

                S_aa = 0.5*(alpha_proj_curr.H()*alpha_proj_curr + alpha_proj_next.H()*alpha_proj_next)  #for debug
                S_ab = 0.5*(alpha_proj_curr.H()*beta_proj_curr + alpha_proj_next.H()*beta_proj_next)
                S_bb = 0.5*(beta_proj_curr.H()*beta_proj_curr + beta_proj_next.H()*beta_proj_next)

                St_aa = alpha_proj_curr.H()*alpha_proj_next  # overlap of wfc at different times
                St_ab = alpha_proj_curr.H()*beta_proj_next
                St_bb = beta_proj_curr.H()*beta_proj_next

                H_aa = 0.5*(e_curr_alpha + e_next_alpha) - (0.5j/dt)*(St_aa - St_aa.H())
                H_ab = 0.5*(e_curr_alpha + e_next_beta) - (0.5j/dt)*(St_ab - St_ab.H())
                H_bb = 0.5*(e_curr_beta + e_next_beta) - (0.5j/dt)*(St_bb - St_bb.H())

                H_aa.real().show_matrix("%s/0_Ham_aa_%d_re" % (rd, curr_index) )
                H_aa.imag().show_matrix("%s/0_Ham_aa_%d_im" % (rd, curr_index) ) 
                H_ab.real().show_matrix("%s/0_Ham_ab_%d_re" % (rd, curr_index) )
                H_ab.imag().show_matrix("%s/0_Ham_ab_%d_im" % (rd, curr_index) ) 
                H_bb.real().show_matrix("%s/0_Ham_bb_%d_re" % (rd, curr_index) )
                H_bb.imag().show_matrix("%s/0_Ham_bb_%d_im" % (rd, curr_index) )    
                S_aa.real().show_matrix("%s/0_S_aa_%d_re" % (rd, curr_index) )
                S_aa.imag().show_matrix("%s/0_S_aa_%d_im" % (rd, curr_index) )
                S_ab.real().show_matrix("%s/0_S_ab_%d_re" % (rd, curr_index) )
                S_ab.imag().show_matrix("%s/0_S_ab_%d_im" % (rd, curr_index) )
                S_bb.real().show_matrix("%s/0_S_bb_%d_re" % (rd, curr_index) )
                S_bb.imag().show_matrix("%s/0_S_bb_%d_im" % (rd, curr_index) )
#            ovlp.real().show_matrix("%s/0_ovlp_%d_re" % (rd, curr_index) )
#            ovlp.imag().show_matrix("%s/0_ovlp_%d_im" % (rd, curr_index) )


  
#            if compute_Hprime==1:
#                os.system("%s convert %s/curr0/x0.export/grid.1 %s/curr0/x0.export/grid.1.xml" % (EXE_CONVERT,wd,wd))
#                curr_wfc0.QE_read_acsii_grid("%s/curr0/x0.export/grid.1.xml" % wd)
#                curr_wfc0.compute_Hprime(0,maxband-minband,"%s/0_Hprime_%d" % (rd, curr_index) )
#                os.system("-rf %s/curr0/x0.export/grid.1.xml" % wd)


            #-----------------------------------------------------------------

            # Remove current run, make next run to be the current one
            os.system("rm -rf %s/curr0" % wd )
            os.system("mv %s/next0 %s/curr0" % (wd, wd) )

            if nac_method==2:
                os.system("rm -rf %s/curr1" % wd )
                os.system("mv %s/next1 %s/curr1" % (wd, wd) )

            print "old files deleted, new have become old"


# ACHTUNG!!! Restoring wfc makes some complications, so we might need to destroy wfc objects
# after each round of operations and create new objects from the beginning - thia may be safer!

        curr_index = curr_index + 1
               
        print "End of step t=", t
        t = t + 1

#================= End of runMD function =============================
