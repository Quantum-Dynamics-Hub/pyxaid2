#***********************************************************
# * Copyright (C) 2017 Wei Li and Alexey V. Akimov
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





def merge_orbitals(Ca, Cb):
    """
    This function puts two matrices together into a single matrix
    """

    npw_a = Ca.num_of_rows
    N_mo_a = Ca.num_of_cols
 
    npw_b = Cb.num_of_rows
    N_mo_b = Cb.num_of_cols

    if npw_a != npw_b:
        print "Error: The number of rows of the two matrices should be equal"
        sys.exit(0)

    C = CMATRIX(npw_a, N_mo_a + N_mo_b)

    push_submatrix(C, Ca, range(0,npw_a), range(0,N_mo_a));
    push_submatrix(C, Cb, range(0,npw_a), range(N_mo_a, N_mo_a + N_mo_b));

    return C


def split_orbitals_energies(C, E):
    """ In SOC, non-collinear case, the orbitals are 2-component spinors:
             | psi_i^alp |          | E_i_alp          |
     psi_i = |           |,  so E = |                  |
             | psi_i^bet |          |          E_i_bet |
  
    So, the wfc we read from the QE calculations, in this case is composed of such
    pairs, going in this order, so:

    psi_i^alp, psi_i^bet, psi_{i+1}^alp, psi_{i+1}^bet, ...

    Thus, there are 2*N_mo_adi columns, so we need to extract spin-components 
    into 2 matrices

    """

    N_pw = C.num_of_rows
    N_mo_adi = C.num_of_cols/2   # C.num_of_cols has to be even

    C_alp = CMATRIX(N_pw, N_mo_adi)
    C_bet = CMATRIX(N_pw, N_mo_adi)
    E_alp = CMATRIX(N_mo_adi, N_mo_adi)
    E_bet = CMATRIX(N_mo_adi, N_mo_adi)


    stenc1, stenc2 = [], []

    for i in xrange(N_mo_adi):
        stenc1.append(2*i)
        stenc2.append(2*i+1)

    pop_submatrix(C, C_alp, range(0, N_pw), stenc1 )
    pop_submatrix(C, C_bet, range(0, N_pw), stenc2 )

    pop_submatrix(E, E_alp, stenc1, stenc1)
    pop_submatrix(E, E_bet, stenc2, stenc2)

    return C_alp, C_bet, E_alp, E_bet

    


def orthogonalize_orbitals(C):
    """
    C = N_pw x N_mo   (just alpha or beta orbitals)
    flag == 2:   C = N_pw x 2*N_mo (both alpha and beta orbitals)

    This function takes an input of orbitals (C), which may not
    be rigorously orthogonal, finds a suitable transformation (U)
    and converts them into rigorously orthogonal orbitals (C_tilda)

    C_tilda = C * U, so if you want

    C_tilda^+  * C_tilda = I, you'll find that

    U = S^{-1/2}, where S = C^+ * C

    """

    S = C.H() * C  # overlap matrix

    S_half = CMATRIX(S.num_of_rows, S.num_of_cols)
    S_i_half = CMATRIX(S.num_of_rows, S.num_of_cols)

    sqrt_matrix(S, S_half, S_i_half)

    C_tilda = C * S_i_half

    return C_tilda


def orthogonalize_orbitals2(Ca,Cb):
    """
    Ca and Cb = N_pw x N_mo   - represent the spin-components
    of the adiabatic states

    This function takes an input of orbitals (C), which may not
    be rigorously orthogonal, finds a suitable transformation (U)
    and converts them into rigorously orthogonal orbitals (C_tilda)

    For each channel:

    C_tilda = C * U, so if you want

    C_tilda^+  * C_tilda = I, you'll find that

    U = S^{-1/2}, where S = Ca^+ * Ca + Cb^+ * Cb

    """

    S = Ca.H() * Ca + Cb.H() * Cb  # overlap matrix

    S_half = CMATRIX(S.num_of_rows, S.num_of_cols)
    S_i_half = CMATRIX(S.num_of_rows, S.num_of_cols)

    sqrt_matrix(S, S_half, S_i_half)

    Ca_tilda = Ca * S_i_half
    Cb_tilda = Cb * S_i_half

    return Ca_tilda, Cb_tilda



def compute_H_el(C_dia_a, C_dia_b, C_adi_a, C_adi_b, E_adi_a, E_adi_b):
    """
    C_dia_ - spin-diabatic wavefunctions. CMATRIX(N_pw, N_mo_dia)
    C_adi_ - spin-adiabatic wavefunctions (with SOC). CMATRIX(N_pw, N_mo_adi)
    E_adi_ - spin-adiabatic energies (with SOC). CMATRIX(N_mo_adi, N_mo_adi)

    _a - alph,  _b - beta

    Here, N_pw, N_mo_dia, N_mo_adi are the numbers of planewaves (basis),
    diabatic and adiabatic states, respectively.

    The diabatic orbitals are obtained from spin-polarized calculations, but without 
    SOC
    The spin-adiabatic orbitals are obtained from non-collinear calculations (with 
    SOC) that also mix alpha and beta components. 

    To compute the matrix elements of electronic Hamiltonian in spin-diabatic 
    basis, We use the following formula:

    ^
    H  =  sum  {  |i> E_i <i| }
           i

    here |i> is spin-orbital  |i> = |i_a> * |a> + |i_b> * |b>


    <x|H|y> = sum { <x|i> * E_i * <i|y>  }
               i

    where i runs over all spin-adiabatic states. 
    
    Because |i> mixes the components, the couplings will be non-zero

    Also, E_adi_a and E_adi_b are the same, so below we'll be using only one

    """

    E_adi = E_adi_a

# Matrix of projection of spin-adiabatic states onto spin-diabatic
# (N_mo_dia x N_mo_adi) = (N_mo_dia x N_pw) x (N_pw x N_mo_adi)
    Sai = C_dia_a.H() * C_adi_a  #  <a|i>
    Sbi = C_dia_b.H() * C_adi_b  #  <b|i>

# Now compute the Hamiltonian
# (N_mo_dia x N_mo_dia) = (N_mo_dia x N_mo_adi) x (N_mo_adi x N_mo_adi) x (N_mo_adi x N_mo_dia)
#  
    H_aa = Sai * E_adi * Sai.H()
    H_bb = Sbi * E_adi * Sbi.H()
    H_ab = Sai * E_adi * Sbi.H()


    return H_aa, H_ab, H_bb



def compute_Hvib(coeff_curr1, coeff_next1, coeff_curr0, coeff_next0, E_adi_curr1, E_adi_next1,e_curr0, e_next0, params):
    """ 
    The first set of coefficients corresponds to the diabatic orbitals:
    coeff_curr and coeff_next - [0] - alpha, [1] - beta
    The second - to adiabatic, energies are also adiabatic [0] - includes both
    spin channels, but the number of orbitals is doubled
    """

    do_orth = params["do_orth"]   # whether to correct orbitals to make them orthogonal
    rd = params["root_directory"] # of where the files will be printed out
    curr_index = params["curr_index"] # a counter for the files
    print_overlaps = params["print_overlaps"] # whether we need overlaps (debug)
    dt = params["dt"]


    #========== Orthogonalize orbitals =========================
    # spin-adiabatic orbitals
    c_adi_curr_a, c_adi_curr_b, e_adi_curr_a, e_adi_curr_b = None, None, None, None
    c_adi_next_a, c_adi_next_b, e_adi_next_a, e_adi_next_b = None, None, None, None

    if do_orth:
        a, b, e_adi_curr_a, e_adi_curr_b = split_orbitals_energies(coeff_curr1[0], E_adi_curr1[0])
        c_adi_curr_a, c_adi_curr_b = orthogonalize_orbitals2(a, b)

        a, b, e_adi_next_a, e_adi_next_b = split_orbitals_energies(coeff_next1[0], E_adi_next1[0])
        c_adi_next_a, c_adi_next_b = orthogonalize_orbitals2(a, b)
    else:
        # Here, we copy matrices by references, but that is okay
        # since we don't modify the final matrix
        c_adi_curr_a, c_adi_curr_b, e_adi_curr_a, e_adi_curr_b = split_orbitals_energies(coeff_curr1[0], E_adi_curr1[0])
        c_adi_next_a, c_adi_next_b, e_adi_next_a, e_adi_next_b = split_orbitals_energies(coeff_next1[0], E_adi_next1[0])
        #c_adi_curr_a, c_adi_curr_b, e_adi_curr_a, e_adi_curr_b = split_orbitals(coeff_curr[0])
        #c_adi_next_a, c_adi_next_b, e_adi_next_a, e_adi_next_b = split_orbitals(coeff_next[0])



    # spin-diabatic orbitals
    c_dia_curr_a, c_dia_curr_b = None, None
    if do_orth:
        c_dia_curr_a = orthogonalize_orbitals(coeff_curr0[0])
        c_dia_curr_b = orthogonalize_orbitals(coeff_curr0[1]) 
    else:
        c_dia_curr_a = coeff_curr0[0]
        c_dia_curr_b = coeff_curr0[1]

    c_dia_next_a, c_dia_next_b = None, None
    if do_orth:
        c_dia_next_a = orthogonalize_orbitals(coeff_next0[0])
        c_dia_next_b = orthogonalize_orbitals(coeff_next0[1])
    else:
        c_dia_next_a = coeff_next0[0]
        c_dia_next_b = coeff_next0[1]



    #========== Compute the electronic Ham ==============
    H_aa_curr, H_ab_curr, H_bb_curr = compute_H_el(c_dia_curr_a, c_dia_curr_b, c_adi_curr_a, c_adi_curr_b, e_adi_curr_a, e_adi_curr_b)
    H_aa_next, H_ab_next, H_bb_next = compute_H_el(c_dia_next_a, c_dia_next_b, c_adi_next_a, c_adi_next_b, e_adi_next_a, e_adi_next_b)


    # overlap of spin-diabatic states at different times
    # <a_curr|a_next>, <a_curr|b_next>, <b_curr|b_next>
    #     St_aa              St_ab          St_bb
    # The <alp|bet> = 0, because these are spin-diabatic states
    St_aa = c_dia_curr_a.H() * c_dia_next_a
    St_bb = c_dia_curr_b.H() * c_dia_next_b

    # Hvib = H_el - i*hbar*d_ab
    Hvib_aa = 0.5*(H_aa_curr + H_aa_next) - (0.5j/dt)*(St_aa - St_aa.H())
    Hvib_ab = 0.5*(H_ab_curr + H_ab_next)
    Hvib_bb = 0.5*(H_bb_curr + H_bb_next) - (0.5j/dt)*(St_bb - St_bb.H())

    sz=Hvib_aa.num_of_cols
    for i in range(0,sz):
        for j in range(0,sz):
            if i==j:
               Hvib_aa.set(i,i,0.5*(e_curr0[0].get(i,i) + e_next0[0].get(i,i)))
               Hvib_bb.set(i,i,0.5*(e_curr0[1].get(i,i) + e_next0[1].get(i,i)))
            else:
                pass

    # Print out the results
    Hvib_aa.real().show_matrix("%s/0_Hvib_aa_%d_re" % (rd, curr_index) )
    Hvib_aa.imag().show_matrix("%s/0_Hvib_aa_%d_im" % (rd, curr_index) )
    Hvib_ab.real().show_matrix("%s/0_Hvib_ab_%d_re" % (rd, curr_index) )
    Hvib_ab.imag().show_matrix("%s/0_Hvib_ab_%d_im" % (rd, curr_index) )
    Hvib_bb.real().show_matrix("%s/0_Hvib_bb_%d_re" % (rd, curr_index) )
    Hvib_bb.imag().show_matrix("%s/0_Hvib_bb_%d_im" % (rd, curr_index) )

    # For the debug (and quality control) - print out the overlaps
    if print_overlaps:
        ovlp = c_adi_curr_a.H() * c_adi_curr_a + c_adi_curr_b.H() * c_adi_curr_b
        ovlp_aa = c_dia_curr_a.H() * c_dia_curr_a
        ovlp_bb = c_dia_curr_b.H() * c_dia_curr_b
 
        ovlp.real().show_matrix("%s/S_adi_%d_re" % (rd, curr_index) )
        ovlp.imag().show_matrix("%s/S_adi_%d_im" % (rd, curr_index) )

        ovlp_aa.real().show_matrix("%s/S_dia_alp_%d_re" % (rd, curr_index) )
        ovlp_aa.imag().show_matrix("%s/S_dia_alp_%d_im" % (rd, curr_index) )

        ovlp_bb.real().show_matrix("%s/S_dia_bet_%d_re" % (rd, curr_index) )
        ovlp_bb.imag().show_matrix("%s/S_dia_bet_%d_im" % (rd, curr_index) )



def read_all(wd, order, ind, act_space, info):
    """
    This function reads index, wfc and grid files from a given directory
    The number of wfc and grid files may be larger than 1 - this is the
    case of spin-polarized or multiple k-points calculations

    Parameters:
    wd (string) - "working directory", the path to cuur* and next* directories
    order (string) - "curr" or "next"
    ind (integer) - 0, 1 This index is also used to enumerate
                  the prefixes of the export directories (types of calculations)
    act_space (list of ints) - defines the indices of the orbitals we are 
                 interested in, so to minimize the computational burden
    info (dictionary) - contains some info about the calculation (e.g. the # 
                of the k-points to read)

    Return values:
    The function returs lists containing: energies, plane wave coefficients,
    and the grid point vectors
    """

    # Verbosity level in 3 functions for reading.
    # Set to 1 if you debug, 0 in the production runs
    verb0 = 0   # for index
    verb1 = 0   # for wfc
    verb2 = 0   # for grid

    file0 = "%s/%s%i/x%i.export/index.xml" % (wd, order, ind, ind)
    print "Reading index from file ", file0
    dum, e = QE_methods.read_qe_index(file0, act_space, verb0)

    coeff = []
    grid = []

    for ik in xrange(info["nk"]):
        print "Handling the k-point %i with coordinates: %8.5f %8.5f %8.5f " \
         % (ik, info["k"][ik].x, info["k"][ik].y, info["k"][ik].z)

        file1 = "%s/%s%i/x%i.export/wfc.%i" % (wd, order, ind, ind, ik+1)
        print "Reading the wfc from file ",file1
        coeff.append( QE_methods.read_qe_wfc(file1, act_space, verb1))

        file2 = "%s/%s%i/x%i.export/grid.%i" % (wd, order, ind, ind, ik+1)
        print "Reading the grid from file ", file2
        grid.append( QE_methods.read_qe_wfc_grid(file2 , verb2) )

    return e, coeff, grid



def runMD(params):
#------------ Read the parameters -----------------
# Parameters meaning
# pp_type - pseudopotential type: US - ultra-soft, NC - norm-conserving, PAW - projector-augmented waves
# wd - working directory, where all output (working) files will be written
# rd - results directory, where all final results (energy, NAC, H', etc.) will be written by default it will be set to wd
# This MD uses corrected NAC method

    tim = Timer()

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
    minband_soc = get_value(params,"minband_soc",1,"i")
    maxband_soc = get_value(params,"maxband_soc",2,"i")
    nac_method = get_value(params,"nac_method",0,"i")  # choose what method for NAC calculations to use: 0 -standard, 1-corrected
    prefix0 = get_value(params,"prefix0","x0.scf","s")
    prefix1 = get_value(params,"prefix1","x1.scf","s")

    
    compute_Hprime = get_value(params,"compute_Hprime",0,"i") # transition dipole moments


    # Sanity/Convention check
    if(minband<=0): 
        print "Error: minband should be >0, current value of minband = ",minband
        sys.exit(0)
    if(minband>maxband):
        print "Error: minband must be smaller or equal to maxband. Current values: minband = ",minband," maxband = ",maxband
        sys.exit(0)


    # nac_method = 0 : non spin polarized/SOC, multiple k-point
    # nac_method = 1 : spin polarized
    # nac_method = 2 : SOC
    # nac_method = 3 : perform 1 and 2 at the same time 
    if nac_method == 0:
        print "you are doing a non spin-polarized calculation for NAC"
    elif nac_method == 1:
        print "you are doing a spin-polarized calculation for NAC"
    elif nac_method == 2:
        print "you are doing a SOC calculation for NAC"
    elif nac_method == 3:
        print "you are doing adiabatic/diabatic projection with SOC for NAC, \
              it will perform the SOC and spin polarized calculation \
              at the same time"
    else:
        print "Error: nac_method must be one of the values in [0,1,2,3]"
        sys.exit(0)

    # Use this for nspin = 1 or 2
    act_sp1 = range(minband, maxband+1)     # min = 1, max = 2 => range(1,3) = [1,2]

    # Use this for nspin = 4
    act_sp2 = range(2*minband_soc-1, 2*(maxband_soc+1)-1 ) # min =1, max = 2 => range(1,5) = [1,2,3,4]

    # Initialize variables
    curr_index = start_indx - 1
    t = start_indx

    print "In runMD: current working directory for python: ",os.getcwd()
    print "In runMD: current working directory for sh:",os.system("echo $(pwd)")

    os.system("mkdir %s" % wd)  # Create the working directory where all output files will be written
                                # results directory should already exist

    while t<=stop_indx:
        print ">>>>>>>>>>>>>>>>>>>>  t= ", t, " <<<<<<<<<<<<<<<<<<<<<"

        dirname = ""
        if t==start_indx:
           print "Starting first point in this batch"
           dirname0 = "curr0"
           dirname1 = "curr1"

        if t>start_indx:
           print "Continuing with other points in this batch"
           dirname0 = "next0"
           dirname1 = "next1"


        tim.start()
        # A common block
        # Run calculations
        # A regular calculation anyway
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


        if nac_method == 2 or nac_method == 3:  # In addition perform the soc calculation
            os.system( "%s -n %s %s < %s.%d.in > %s.%d.out" % (BATCH_SYSTEM,NP,EXE,prefix1,t,prefix1,t) )
            os.system( "%s -n %s %s < x1.exp.in > x1.exp.out" % (BATCH_SYSTEM,NP,EXE_EXPORT) )

            dirname = dirname1

            os.system("mkdir %s/%s" % (wd,dirname) )

            os.system( "mv %s.%d.out %s/%s" % (prefix1,t, wd, dirname) )
            os.system( "mv *.wfc* %s/%s" % (wd, dirname) )
            os.system( "mv x1.export %s/%s" % (wd, dirname) ) # "x1" - corresponds to x1 as a prefix in input files

        print "Time to run first calculations = ", tim.stop(); 

        # Now general part - from current and next wavefunctions calculate NACs:
        # First see wther the calculation is what we wanted
        if curr_index>=start_indx:
            tim.start()
            print "Generate NAC from WFCs at two adjacent points"

            # some checks
            # for non soc/spin-polarized cases
            if nac_method == 0:
                info0, all_e_dum0 = QE_methods.read_qe_index("%s/curr0/x0.export/index.xml" % wd, [], 0)

                if info0["nspin"] != 1:
                    print "Error,you are not running the non spin polarized calculation \
                           check your setting with nspin"

                    sys.exit(0)

                print "The total # of k-points (non spin polarized calculation) is: ", info0["nk"]

            # for the spin polarized case 
            if nac_method == 1 or nac_method == 3:
                info0, all_e_dum0 = QE_methods.read_qe_index("%s/curr0/x0.export/index.xml" % wd, [], 0)

                if info0["nspin"] != 2:
                    print "Error, you are not running spin polarized calc (generating the diabatic basis),\
                           check you settings with nspin"

                    sys.exit(0)

                print "The total # of k-points (spin-polarized) including up and down components is: ", info0["nk"]


            # for soc cases 
            if nac_method == 2 or nac_method==3:
                info1, all_e_dum1 = QE_methods.read_qe_index("%s/curr1/x1.export/index.xml" % wd, [], 0)

                if info1["nspin"] != 4:
                    print "Error,you are not running SOC calculation (generating the spin-adiabatic basis) \
                           check you setting with nspin"

                    sys.exit(0)

                print "The total # of k-points (soc) is: ", info1["nk"]



            # read the coefficients and energies for the mluti k-points cases, even if some cases require gamma only
                        
            if nac_method == 0 or nac_method == 1 or nac_method == 3:
                # read the coefficients anyway

                act_space = act_sp1
                #====== Current electronic structure ========
                e_curr0, coeff_curr0, grid_curr0 = read_all(wd, "curr", 0, act_space, info0)

                #====== Next electronic structure ===========
                e_next0, coeff_next0, grid_next0 = read_all(wd, "next", 0, act_space, info0)



            if nac_method == 2 or nac_method == 3:
                # includes the soc
                act_space = act_sp2

                #====== Current electron electructure =======
                e_curr1, coeff_curr1, grid_curr1 = read_all(wd, "curr", 1, act_space, info1)

                #====== Next electronic structure ===========
                e_next1, coeff_next1, grid_next1 = read_all(wd, "next", 1, act_space, info1)

                
 
            print "Time to read index, wfc, and wfc grids = ", tim.stop();

            ######################################### NAC calculation #######################################
            # Finally compute Hamiltonian and the overlap matrix
            S, H, S_dia, H_dia, H_soc, S_soc  = None, None, None, None, None, None

            # non spin-polarized case
            if nac_method == 0 or nac_method == 1 or nac_method == 3:
             

                if info0["nspin"]==1:  # non SOC case
                    if info0["nk"]==1: # Only one k-point

                        orthogonalize=1
                        if orthogonalize==1:
                            print "Do internal orbital orthogonalization"
                            coeff_curr0[0] = orthogonalize_orbitals(coeff_curr0[0])
                            coeff_next0[0] = orthogonalize_orbitals(coeff_next0[0])


                        ovlp_cn  = coeff_curr0[0].H() * coeff_next0[0]
                        H = 0.5*(e_curr0[0] + e_next0[0]) - (0.5j/dt)*(ovlp_cn - ovlp_cn.H())
                        S = 0.5 *(coeff_curr0[0].H() * coeff_curr0[0] + coeff_next0[0].H() * coeff_next0[0]) # for debug

                    else: 
                        print "you are dealing with multiple kpoints"

                        as_sz = len(act_sp1)
                        H = CMATRIX(info0["nk"]*as_sz, info0["nk"]*as_sz )
                        S = CMATRIX(info0["nk"]*as_sz, info0["nk"]*as_sz )

                        # optional orthogonalization - to mitigate the round off errors
                        orthogonalize=1
                        if orthogonalize==1:
                            print "Do internal orbital orthogonalization"
                            for ik1 in xrange(info0["nk"]):
                                ovlp_cc = pw_overlap(info0["k"][ik1], info0["k"][ik1], coeff_curr0[ik1], coeff_curr0[ik1], grid_curr0[ik1], grid_curr0[ik1])
                                ovlp_nn = pw_overlap(info0["k"][ik1], info0["k"][ik1], coeff_next0[ik1], coeff_next0[ik1], grid_next0[ik1], grid_next0[ik1])

                                ovlp_cc_half = CMATRIX(ovlp_cc.num_of_rows, ovlp_cc.num_of_cols)
                                ovlp_cc_i_half = CMATRIX(ovlp_cc.num_of_rows, ovlp_cc.num_of_cols)
                                sqrt_matrix(ovlp_cc, ovlp_cc_half, ovlp_cc_i_half)
                                coeff_curr0[ik1] = coeff_curr0[ik1] * ovlp_cc_i_half

                                ovlp_nn_half = CMATRIX(ovlp_nn.num_of_rows, ovlp_nn.num_of_cols)
                                ovlp_nn_i_half = CMATRIX(ovlp_nn.num_of_rows, ovlp_nn.num_of_cols)
                                sqrt_matrix(ovlp_nn, ovlp_nn_half, ovlp_nn_i_half)
                                coeff_next0[ik1] = coeff_next0[ik1] * ovlp_nn_i_half




                        """
                        The convention for the matrices for multiple k-points is:

                              |  x_11  x_12 ... |
                        X =   |  x_21  x_22 ... |
                              }  ...      ...   |

                        here, each x_ij block is a  as_sz x as_sz matrix describing
                        the interactions of the as_sz orbitals of a k-point i and 
                        as_sz orbitals of a k-point j

                        So X is has a k-point first block-structure

                        """                     

                        for ik1 in xrange(info0["nk"]):
                            for ik2 in range(ik1, info0["nk"]):
                                tim.start()
                                ovlp_cc = pw_overlap(info0["k"][ik1], info0["k"][ik2], coeff_curr0[ik1], coeff_curr0[ik2], grid_curr0[ik1], grid_curr0[ik2])
                                ovlp_nn = pw_overlap(info0["k"][ik1], info0["k"][ik2], coeff_next0[ik1], coeff_next0[ik2], grid_next0[ik1], grid_next0[ik2])
                                #ovlp_nc = pw_overlap(info["k"][ik1], info["k"][ik2], coeff_next[ik1], coeff_curr[ik2], grid_next[ik1], grid_curr[ik2])
                                ovlp_cn = pw_overlap(info0["k"][ik1], info0["k"][ik2], coeff_curr0[ik1], coeff_next0[ik2], grid_curr0[ik1], grid_next0[ik2])
                    
                                print "Time to compute 3 overlaps for the pair of k-points ", ik1, " ", ik2," is ", tim.stop()

                                
                                h_cc = CMATRIX(as_sz, as_sz)
                                h_nn = CMATRIX(as_sz, as_sz)

                                tim.start()
                                for i1 in xrange(as_sz):
                                    for j1 in xrange(as_sz):
                                        h_cc.set(i1, j1, 0.5*(e_curr0[ik1].get(i1,i1) + e_curr0[ik2].get(j1,j1))*ovlp_cc.get(i1,j1)) 
                                        h_nn.set(i1, j1, 0.5*(e_next0[ik1].get(i1,i1) + e_next0[ik2].get(j1,j1))*ovlp_nn.get(i1,j1))

                                h = 0.5*(h_cc + h_nn)  - (0.5j/dt)*(ovlp_cn - ovlp_cn.H()) 
                                s = 0.5*(ovlp_cc + ovlp_nn)
      
                                if ik2!=ik1:
                                    push_submatrix(S, s, range(ik1*as_sz, (ik1+1)*as_sz), range(ik2*as_sz, (ik2+1)*as_sz))
                                    push_submatrix(S, s.H(), range(ik2*as_sz, (ik2+1)*as_sz), range(ik1*as_sz, (ik1+1)*as_sz))

                                    push_submatrix(H, h, range(ik1*as_sz, (ik1+1)*as_sz), range(ik2*as_sz, (ik2+1)*as_sz))
                                    push_submatrix(H, h.H(), range(ik2*as_sz, (ik2+1)*as_sz), range(ik1*as_sz, (ik1+1)*as_sz))


                                else:
                                    push_submatrix(S, s, range(ik1*as_sz, (ik1+1)*as_sz), range(ik2*as_sz, (ik2+1)*as_sz))
                                    push_submatrix(H, h, range(ik1*as_sz, (ik1+1)*as_sz), range(ik2*as_sz, (ik2+1)*as_sz))
                                print "Time to push matrices is ", tim.stop()


                elif info0["nspin"]==2:

                    if info0["nk"] == 2:  #single k-point! 

                        # assume we use only the alpha coefficient, similar as PYXAID1, ham() function
                        # H_dia =  Eii - i*hbar*(<i(t)|j(t+dt)> - <i(t+dt)|j(t)>)

                        orthogonalize = 1
                        if orthogonalize==1:
                            print "Do internal orbital orthogonalization"
                            coeff_curr0[0] = orthogonalize_orbitals(coeff_curr0[0])
                            coeff_next0[0] = orthogonalize_orbitals(coeff_next0[0])

                        ovlp_cn  = coeff_curr0[0].H() * coeff_next0[0]   
                        H_dia = 0.5*(e_curr0[0] + e_next0[0]) - (0.5j/dt)*(ovlp_cn - ovlp_cn.H())
                        S_dia = 0.5 *(coeff_curr0[0].H() * coeff_curr0[0] + coeff_next0[0].H() * coeff_next0[0])
                    
                    else:

                        print "multiple k-point for spin-polarized case is not yet implemented"
                        sys.exit(0)



            if nac_method == 2 or nac_method == 3:

                if info1["nk"]==1: # Only one k-point

                    S = coeff_curr1[0].H() * coeff_curr1[0]
                    St = coeff_curr1[0].H() * coeff_next1[0]  # overlap of wfc at different times


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
                        ec.set(n,n, 0.5*(e_curr1[0].get(2*n, 2*n)+e_curr1[0].get(2*n+1, 2*n+1)) )
                        en.set(n,n, 0.5*(e_next1[0].get(2*n, 2*n)+e_next1[0].get(2*n+1, 2*n+1)) )

                    H_soc = 0.5*(ec + en) - (0.5j/dt)*(nac - nac.H())
                    S_soc = ovlp

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
            if  nac_method == 3:
                # check whether the adiabatic and diabatic basis have the same number of plane waves
                # the reason why I used the read_qe_wfc_info is because I will need the ngw 
                # to check the consistency 
                # But the read_qe_index does not read it, so in order to avoid the changes in the Libra code, 
                # I use the read_qe_wfc_info.
                info_wfc0 = QE_methods.read_qe_wfc_info("%s/curr0/x0.export/wfc.1" % wd,0)
                info_wfc1 = QE_methods.read_qe_wfc_info("%s/curr1/x1.export/wfc.1" % wd,0)

                if info_wfc0["ngw"] != info_wfc1["ngw"]:
                    print "Error: the number of plane waves between diabatic and adiabatic does not equal"
                    sys.exit(0)


                prms = {"do_orth": 1, "root_directory" : rd, "curr_index" : curr_index, "print_overlaps" : 1, "dt": dt}
                compute_Hvib(coeff_curr1, coeff_next1, coeff_curr0, coeff_next0, e_curr1, e_next1, e_curr0, e_next0, prms)



            if nac_method == 0:
                H.real().show_matrix("%s/0_Ham_%d_re" % (rd, curr_index) )
                H.imag().show_matrix("%s/0_Ham_%d_im" % (rd, curr_index) )

            if nac_method == 1 or nac_method == 3:
                H_dia.real().show_matrix("%s/0_H_dia_%d_re" % (rd, curr_index) )
                H_dia.imag().show_matrix("%s/0_H_dia_%d_im" % (rd, curr_index) )

            if nac_method == 2 or nac_method == 3:
                H_soc.real().show_matrix("%s/0_H_soc_%d_re" % (rd, curr_index) )
                H_soc.imag().show_matrix("%s/0_H_soc_%d_im" % (rd, curr_index) )



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
