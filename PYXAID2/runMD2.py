#***********************************************************
# * Copyright (C) 2017 Alexey V. Akimov
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


    # nac_method = 0 : non spin polarized/SOC, multiple k-point
    # nac_method = 1 : SOC
    # nac_method = 2 : spin polarized
    # nac_method = 3 : perform 1 and 2 at the same time 
    if nac_method == 0:
        print "you are doing a non spin polarized calculation for NAC"
    elif nac_method == 1:
        print "you are doing a SOC calculation for NAC"
    elif nac_method == 2:
        print "you are doing a spin polarized calculation for NAC"
    elif nac_method == 3:
        print "you are doing adiabatic/diabatic projrction with SOC for NAC, \
              it will perform the SOC and spin polarized calculation \
              at the same time"
    else:
        print "Error: nac_method must be one of the values in [0,1,2,3]"
        sys.exit(0)

    # Use this for nspin = 1 or 2
    act_sp1 = range(minband, maxband+1)     # min = 1, max = 2 => range(1,3) = [1,2]

    # Use this for nspin = 4
    act_sp2 = range(2*minband-1, 2*(maxband+1)-1 ) # min =1, max = 2 => range(1,5) = [1,2,3,4]

    # Initialize variables
    curr_index = start_indx - 1
    t = start_indx

    print "In runMD: current working directory for python: ",os.getcwd()
    print "In runMD: current working directory for sh:",os.system("echo $(pwd)")

    os.system("mkdir %s" % wd)  # Create the working directory where all output files will be written
                                # results directory should already exist

    while t<=stop_indx:
        print "t= ", t

        dirname = ""
        if t==start_indx:
           print "Starting first point in this batch"
           dirname0 = "curr0"
           dirname1 = "curr1"

        if t>start_indx:
           print "Continuing with other points in this batch"
           dirname0 = "next0"
           dirname1 = "next1"


        # A common block
        # Run calculations
        # A regular calculation anyway, no mattter whether it includes the SOC effect
        if nac_method == 0 or nac_method == 1 or nac_method == 3:
            os.system( "%s -n %s %s < %s.%d.in > %s.%d.out" % (BATCH_SYSTEM,NP,EXE,prefix0,t,prefix0,t) )
            os.system( "%s -n %s %s < x0.exp.in > x0.exp.out" % (BATCH_SYSTEM,NP,EXE_EXPORT) )

            dirname = dirname0

            # Create temporary directory
            os.system("mkdir %s/%s" % (wd, dirname) )

            # Copy some results to that directory
            os.system( "mv %s.%d.out %s/%s" % (prefix0,t, wd, dirname) )
            os.system( "mv *.wfc* %s/%s" % (wd, dirname) )
            os.system( "mv x0.export %s/%s" % (wd, dirname) ) # "x0" - corresponds to x0 as a prefix in input files


        if nac_method == 2 or nac_method == 3:  # In addition perform the spin-polarized calculation
            os.system( "%s -n %s %s < %s.%d.in > %s.%d.out" % (BATCH_SYSTEM,NP,EXE,prefix1,t,prefix1,t) )
            os.system( "%s -n %s %s < x1.exp.in > x1.exp.out" % (BATCH_SYSTEM,NP,EXE_EXPORT) )

            dirname = dirname1

            os.system("mkdir %s/%s" % (wd,dirname) )

            os.system( "mv %s.%d.out %s/%s" % (prefix1,t, wd, dirname) )
            os.system( "mv *.wfc* %s/%s" % (wd, dirname) )
            os.system( "mv x1.export %s/%s" % (wd, dirname) ) # "x1" - corresponds to x1 as a prefix in input files



        # Now general part - from current and next wavefunctions calculate NACs:
        # First see wther the calculation is what we wanted
        if curr_index>=start_indx:
            print "Generate NAC from WFCs at two adjacent points"

            # some checks
            # for non soc/spin-polarized cases
            if nac_method == 0:
                info, all_e_dum = QE_methods.read_qe_index("%s/curr0/x0.export/index.xml" % wd, [], 0)

                if info["nspin"] != 1:
                    print "Error,you are not running the non spin polarized calculation \
                           check your setting with nspin"

                    sys.exit(0)

                print "The total # of k-points (non spin polarized calculation) is: ", info["nk"]

            # for soc cases
            if nac_method == 1 or nac_method == 3:
                info, all_e_dum = QE_methods.read_qe_index("%s/curr0/x0.export/index.xml" % wd, [], 0)

                if info["nspin"] != 4:
                    print "Error,you are not running SOC calculation (generating the spin-adiabatic basis) \
                           check you setting with nspin"

                    sys.exit(0)

                print "The total # of k-points (soc) is: ", info["nk"]

            # for the spin polarized case 

            if nac_method == 2 or nac_method==3:
                info1, all_e_dum = QE_methods.read_qe_index("%s/curr1/x1.export/index.xml" % wd, [], 0)

                if info1["nspin"] != 2:
                    print "Error, you are not running spin polarized calc (generating the diabatic basis),\
                           check you settings with nspin"

                    sys.exit(0)

                print "The total # of k-points (spin-polarized) including up and down components is: ", info1["nk"]

            # read the coefficients and energies for the mluti k-points cases, even if some cases require gamma only
            
            
            if nac_method == 0 or nac_method == 1 or nac_method == 3:
                # read the coefficients anyway
                # no matter it includes the soc

                if info["nspin"] == 1:
                    act_space = act_sp1
                elif info["nspin"] == 4:
                    act_space = act_sp2

                #====== Current electronic structure ===========
                dum, e_curr = QE_methods.read_qe_index("%s/curr0/x0.export/index.xml" % wd, act_space, 0)
                coeff_curr = []
                grid_curr = []

                for ik in xrange(info["nk"]):
                    print ik, info["k"][ik]
                    coeff_curr.append( QE_methods.read_qe_wfc("%s/curr0/x0.export/wfc.%i" % (wd, ik+1), act_space, 0))
                    grid_curr.append( QE_methods.read_qe_wfc_grid("%s/curr0/x0.export/grid.%i" % (wd, ik+1) , 0) )

                #====== Next electronic structure ===========
                dum, e_next = QE_methods.read_qe_index("%s/next0/x0.export/index.xml" % wd, act_space, 0)
                coeff_next = []
                grid_next = []

                for ik in xrange(info["nk"]):
                    print ik, info["k"][ik]
                    coeff_next.append( QE_methods.read_qe_wfc("%s/next0/x0.export/wfc.%i" % (wd, ik+1), act_space, 0))
                    grid_next.append( QE_methods.read_qe_wfc_grid("%s/next0/x0.export/grid.%i" % (wd, ik+1) , 0) )


            if nac_method == 2 or nac_method == 3:
                # read coefficients and energies for spin-polarized case
                # we know wfc.nk+1 and wfc.nk+2 are up and down coefficients, for nk in range(0,info["nk"])

                act_space = act_sp1

                #====== Current electron electructure ========
                # we do use the info1 dict here, since the info1 contains the key for k-vector
                dum, e_curr1 = QE_methods.read_qe_index("%s/curr1/x1.export/index.xml" % wd, act_space, 0)
                coeff_curr1 = []
                grid_curr1 = []

                for ik in xrange(info1["nk"]):
                    print ik, info1["k"][ik]
                    coeff_curr1.append(QE_methods.read_qe_wfc("%s/curr1/x1.export/wfc.%i" % (wd,ik+1),act_space,0)) 
                    grid_curr1.append( QE_methods.read_qe_wfc_grid("%s/curr1/x1.export/grid.%i" % (wd, ik+1) , 0) )
                    

                #====== Next electronic structure ===========
                dum, e_next1 = QE_methods.read_qe_index("%s/next1/x1.export/index.xml" % wd, act_space, 0)
                coeff_next1 = []
                grid_next1 = []


                for ik in xrange(info1["nk"]):
                    print ik, info1["k"][ik]
                    coeff_next1.append(QE_methods.read_qe_wfc("%s/next1/x1.export/wfc.%i" % (wd,ik+1),act_space,0))
                    grid_next1.append( QE_methods.read_qe_wfc_grid("%s/next1/x1.export/grid.%i" % (wd, ik+1) , 0) )
                
                #sys.exit(0)
 


            ######################################### NAC calculation #######################################
            # Finally compute Hamiltonian and the overlap matrix
            S, H, S_dia, H_dia, S_aa, S_bb, S_ab, H_aa, H_ab, H_bb = None, None, None, None, None, None,None, None, None, None

            # non spin-polarized case
            if nac_method == 0 or nac_method == 1 or nac_method == 3:

                if info["nspin"]==1:  # non soc case
                    if info["nk"]==1: # Only one k-point
                        ovlp  = coeff_curr[0].H() * coeff_next[0]
                        H = 0.5*(e_curr[0] + e_next[0]) - (0.5j/dt)*(ovlp - ovlp.H())
                        S = 0.5 *(coeff_curr[0].H() * coeff_curr[0] + coeff_next[0].H() * coeff_next[0]) # for debug

                    else: 
                        print "you are deal with many kpoint"

                        as_sz = len(act_sp1)
                        H = CMATRIX(info["nk"]*as_sz, info["nk"]*as_sz )
                        S = CMATRIX(info["nk"]*as_sz, info["nk"]*as_sz )
                        for ik1 in xrange(info["nk"]):
                            for ik2 in range(ik1, info["nk"]):
                                ovlp_cc = pw_overlap(info["k"][ik1], info["k"][ik2], coeff_curr[ik1], coeff_curr[ik2], grid_curr[ik1], grid_curr[ik2])
                                ovlp_nn = pw_overlap(info["k"][ik1], info["k"][ik2], coeff_next[ik1], coeff_next[ik2], grid_next[ik1], grid_next[ik2])
                                ovlp_cn = pw_overlap(info["k"][ik1], info["k"][ik2], coeff_curr[ik1], coeff_next[ik2], grid_curr[ik1], grid_next[ik2])
                    
                                h_cc = CMATRIX(as_sz, as_sz)
                                h_nn = CMATRIX(as_sz, as_sz)

                                for i1 in xrange(as_sz):
                                    for j1 in xrange(as_sz):
                                        h_cc.set(i1, j1, 0.5*(e_curr[ik1].get(i1,i1) + e_curr[ik2].get(j1,j1))*ovlp_cc.get(i1,j1)) 
                                        h_nn.set(i1, j1, 0.5*(e_next[ik1].get(i1,i1) + e_next[ik2].get(j1,j1))*ovlp_nn.get(i1,j1))

                                h = 0.5*(h_cc + h_nn)  -(0.5j/dt)*(ovlp_cn - ovlp_cn.H()) 
                                s = 0.5*(ovlp_cc + ovlp_nn)
      
                                if ik2!=ik1:
                                    push_submatrix(S, s, range(ik1*as_sz, (ik1+1)*as_sz), range(ik2*as_sz, (ik2+1)*as_sz))
                                    push_submatrix(S, s.H(), range(ik2*as_sz, (ik2+1)*as_sz), range(ik1*as_sz, (ik1+1)*as_sz))

                                    push_submatrix(H, h, range(ik1*as_sz, (ik1+1)*as_sz), range(ik2*as_sz, (ik2+1)*as_sz))
                                    push_submatrix(H, h.H(), range(ik2*as_sz, (ik2+1)*as_sz), range(ik1*as_sz, (ik1+1)*as_sz))


                                else:
                                    push_submatrix(S, s, range(ik1*as_sz, (ik1+1)*as_sz), range(ik2*as_sz, (ik2+1)*as_sz))
                                    push_submatrix(H, h, range(ik1*as_sz, (ik1+1)*as_sz), range(ik2*as_sz, (ik2+1)*as_sz))
                        


                elif info["nspin"]==4:
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

            
            # spin-polarized case
            if nac_method == 2 or nac_method == 3:
                if info1["nk"] == 2:  #gamma point only

                    # assume we use only the alpha coefficient, similar as PYXAID1, ham() function
                    # H_dia =  Eii - i*hbar*(<i(t)|j(t+dt)> - <i(t+dt)|j(t)>)

                    ovlp  = coeff_curr1[0].H() * coeff_next1[0]   
                    H_dia = 0.5*(e_curr1[0] + e_next1[0]) - (0.5j/dt)*(ovlp - ovlp.H())
                    S_dia = 0.5 *(coeff_curr1[0].H() * coeff_curr1[0] + coeff_next1[0].H() * coeff_next1[0])

                    print 'begin the projection of adiabatic/diabatic basis'

                    # check whether the adiabatic and diabatic basis have the same number of plane waves
                    # the reason why I used the read_qe_wfc_info is because I will need the ngw 
                    # for constructing the matrix store the diabatic wavefunction. 
                    # But the read_qe_index does not read it, so in order to avoid the changes in the Libra code, 
                    # I use the read_qe_wfc_info.

                    info_wfc = QE_methods.read_qe_wfc_info("%s/curr0/x0.export/wfc.1" % wd,0)
                    info_wfc1 = QE_methods.read_qe_wfc_info("%s/curr1/x1.export/wfc.1" % wd,0)

                    if info_wfc["ngw"] != info_wfc1["ngw"]:

                        print "Error: the number of plane waves between diabatic and adiabatic does not equal"
                        sys.exit(0)
                
                    # project the spin-adiabatic basis onto spin-diabatic basis
                    # |i> => <a|i>    |i> => <b|i>
                    # i is the spin adiabatic basis which mixes the alp(a) and bet(b) spin from SOC calc,NPW*M matrix; 
                    # alp (a) or bet (b) is the spin diabatic basis from spin-polarized calc, NPW*N matrix
                    # NPW-the number of plane wave, N-the number of diabatic state, M-the number of adiabatic state
                    # Then we will get a N*M matrix for <a|i> and <b|i>
                    # coeff_curr1[0] - the current spin diabatic state for alp, NPW*N matrix
                    # coeff_curr1[1] - the current spin diabatic state for bet, NPW*N matrix
                    # coeff_next1[0] - the next spin diabatic state for alp, NPW*N matrix
                    # coeff_next1[1] - the next spin diabatic state for bet, NPW*N matrix
                    # coeff_curr[0] - the current spin adiabatic state which mix alp and bet, NPW*M matrix
                    # e_curr[0] - the current adiabatic state energy, a M*M matrix
                    # e_curr[1] - the next adiabatic state energy, a M*M matrix

                    alpha_proj_curr = coeff_curr1[0].H() * coeff_curr[0]
                    beta_proj_curr = coeff_curr1[1].H()*coeff_curr[0]

                    alpha_proj_next = coeff_next1[0].H() * coeff_next[0]
                    beta_proj_next = coeff_next1[1].H()*coeff_next[0]


                    # electronic Ham in the form Hab = <a|i> Ei <i|b>, 
                    # sum over all adiabatic state, M
                    # Ei is the spin adiabatic state energy
                    # finally we get a N*N matrix for the electronic Ham

                    N = len(act_sp1)
                    M = len(act_sp2)

                    h_aa_cc = CMATRIX(N, N)
                    h_ab_cc = CMATRIX(N, N)
                    h_bb_cc = CMATRIX(N, N)

                    h_aa_nn = CMATRIX(N, N)
                    h_ab_nn = CMATRIX(N, N)
                    h_bb_nn = CMATRIX(N, N)

                    for n1 in range(0,N):
                        for n2 in range(0,N):

                            tmp_aa_cc,tmp_ab_cc,tmp_bb_cc = 0.0,0.0,0.0
                            tmp_aa_nn,tmp_ab_nn,tmp_bb_nn = 0.0,0.0,0.0
                            for m in range(0,M):
                                tmp_aa_cc += alpha_proj_curr.get(n1,m)*alpha_proj_curr.H().get(m,n2)*e_curr[0].get(m, m)
                                tmp_ab_cc += alpha_proj_curr.get(n1,m)*beta_proj_curr.H().get(m,n2)*e_curr[0].get(m, m)
                                tmp_bb_cc += beta_proj_curr.get(n1,m)*beta_proj_curr.H().get(m,n2)*e_curr[0].get(m, m)

                                tmp_aa_nn += alpha_proj_next.get(n1,m)*alpha_proj_next.H().get(m,n2)*e_next[0].get(m, m)
                                tmp_ab_nn += alpha_proj_next.get(n1,m)*beta_proj_next.H().get(m,n2)*e_next[0].get(m, m)
                                tmp_bb_nn += beta_proj_next.get(n1,m)*beta_proj_next.H().get(m,n2)*e_next[0].get(m, m)

                            h_aa_cc.set(n1,n2,tmp_aa_cc)
                            h_ab_cc.set(n1,n2,tmp_ab_cc)
                            h_bb_cc.set(n1,n2,tmp_bb_cc)

                            h_aa_nn.set(n1,n2,tmp_aa_nn)
                            h_ab_nn.set(n1,n2,tmp_ab_nn)
                            h_bb_nn.set(n1,n2,tmp_bb_nn)                         

                    # overlap of spin-diabatic states at different times
                    # <a_curr|a_next>, <a_curr|b_next>, <b_curr|b_next>
                    #     St_aa              St_ab          St_bb
                    St_aa = coeff_curr1[0].H() * coeff_next1[0]   
                    St_ab = coeff_curr1[0].H() * coeff_next1[1]
                    St_bb = coeff_curr1[1].H() * coeff_next1[1]


                    # invoke the standard form of H_vib_ab
                    # H_vib_ab = H_ab - i*hbar*d_ab
                    H_aa = 0.5*(h_aa_cc + h_aa_nn) - (0.5j/dt)*(St_aa - St_aa.H())
                    H_ab = 0.5*(h_ab_cc + h_ab_nn) - (0.5j/dt)*(St_ab - St_ab.H())
                    H_bb = 0.5*(h_bb_cc + h_bb_nn) - (0.5j/dt)*(St_bb - St_bb.H())
 
                    # H_aa = h_aa_cc - (0.5j/dt)*(St_aa - St_aa.H())
                    # H_ab = h_ab_cc - (0.5j/dt)*(St_ab - St_ab.H())
                    # H_bb = h_bb_cc - (0.5j/dt)*(St_bb - St_bb.H())

                    H_aa.real().show_matrix("%s/0_Ham_aa_%d_re" % (rd, curr_index) )
                    H_aa.imag().show_matrix("%s/0_Ham_aa_%d_im" % (rd, curr_index) ) 
                    H_ab.real().show_matrix("%s/0_Ham_ab_%d_re" % (rd, curr_index) )
                    H_ab.imag().show_matrix("%s/0_Ham_ab_%d_im" % (rd, curr_index) ) 
                    H_bb.real().show_matrix("%s/0_Ham_bb_%d_re" % (rd, curr_index) )
                    H_bb.imag().show_matrix("%s/0_Ham_bb_%d_im" % (rd, curr_index) )    


                    H_dia.real().show_matrix("%s/0_Ham_dia_%d_re" % (rd, curr_index) )
                    H_dia.imag().show_matrix("%s/0_Ham_dia_%d_im" % (rd, curr_index) )    
                    S_dia.real().show_matrix("%s/0_S_dia_%d_re" % (rd, curr_index) )
                    S_dia.imag().show_matrix("%s/0_S_dia_%d_im" % (rd, curr_index) )

                else:

                    print "multiple k-point for spin-polarized case is not yet implemented"
                    sys.exit(0)


            H.real().show_matrix("%s/0_Ham_%d_re" % (rd, curr_index) )
            H.imag().show_matrix("%s/0_Ham_%d_im" % (rd, curr_index) )    
            S.real().show_matrix("%s/0_S_%d_re" % (rd, curr_index) )
            S.imag().show_matrix("%s/0_S_%d_im" % (rd, curr_index) )


            #-----------------------------------------------------------------

            # Remove current run, make next run to be the current one
            if nac_method == 0 or nac_method == 1 or nac_method==3:
                os.system("rm -rf %s/curr0" % wd )
                os.system("mv %s/next0 %s/curr0" % (wd, wd) )

            if nac_method==2 or nac_method==3:
                os.system("rm -rf %s/curr1" % wd )
                os.system("mv %s/next1 %s/curr1" % (wd, wd) )

            print "old files deleted, new have become old"


# ACHTUNG!!! Restoring wfc makes some complications, so we might need to destroy wfc objects
# after each round of operations and create new objects from the beginning - thia may be safer!

        curr_index = curr_index + 1
               
        print "End of step t=", t
        t = t + 1

#================= End of runMD function =============================
