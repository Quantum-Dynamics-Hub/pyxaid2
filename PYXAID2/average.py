#***********************************************************
# * Copyright (C) 2013 Alexey V. Akimov
# * This file is distributed under the terms of the
# * GNU General Public License as published by the
# * Free Software Foundation; either version 3 of the
# * License, or (at your option) any later version.
# * http://www.gnu.org/copyleft/gpl.txt
#***********************************************************/


###################################################################################
#
#  Description:
#  This file performs averaging of results of the cppnamd calculations
#  over different initial conditions (initial times of the MD trajectory).
#  It may also perform contraction of different initial excitations (microstates)
#  into macrostates and compute corresponding populations and energies
#
#  Requirements:
#  It should be run in the directory which contains:
#  1) output files:   out0, out1, etc.
#  2) multi-electron population files: me_pop0, me_pop1, etc.
#  3) energies of the microstates: me_energies0, me_energies1, etc.
#  Thus basically in the same directory in which cppnamd was executed
#
#  Output:
#  For opt==1 or opt==12
#  1) SH-based populations of all microstates for different initial excitations (microstates): 
#     sh_pop_ex0, sh_pop_ex1, etc.
#  2) SE-based populations of all microstates for different initial excitations (microstates):
#     se_pop_ex0, se_pop_ex1, etc.
#  3) SH-based energies of the system for different initial excitations (microstates):
#     sh_en_ex0, sh_en_ex1, etc.
#  4) SE-based energies of the system for different initial excitations (microstates):
#     se_en_ex0, se_en_ex1, etc.

#  For opt==2 or opt==12
#  5) SH-based populations of all macrostates for different initial excitations (macrostates):
#     sh_pop_macro_ex0, sh_pop_macro_ex1, etc.
#  6) SE-based populations of all macrostates for different initial excitations (macrostates):
#     se_pop_macro_ex0, se_pop_macro_ex1, etc.
#
####################################################################################

from aux import *


def average(namdtime,num_states,iconds,opt,MS,inp_dir,res_dir):
#>>>>>>>>>>>>>>>>>> Input description <<<<<<<<<<<<<<<<<<<<<<<<<<<
# namdtime - is a length of na-md trajectories (in # of steps)
# num_states - number of microstates
# iconds - list of initial conditions - iconds[i][0] - time, iconds[i][1] - state
# opt - is a option for execution
#opt = 1  # run only first-part
#opt = 2  # run only seconf part
#opt = 12 # run both first and the second parts
# MS - array of macrostates, e.g.
# MS[0] = [0] - means macrostate 0 contains only the microstate 0
# MS[1] = [1,2,3] - means macrostate 1 contains microstates 1,2,and 3

#>>>>>>>>>>>>>>>>>>> Pre-processing <<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Determine how many and what initial excitations are considered
    dist_ex = []      # Distinct excitations (excite to the same state, may be at different times)
    all_ex = []       # All initial conditions
    niconds = len(iconds)  # Number of initial conditions
    i = 0
    while i<niconds:
        e = iconds[i][1]
        all_ex.append(e)
        if (e in dist_ex):
            pass 
        else:
            dist_ex.append(e)
        i = i + 1

    print "All initial conditions", all_ex
    print "Distinct initial excitations", dist_ex
    print "Number of initial conditions = ", len(all_ex)
    print "Number distinct initial excitations = ", len(dist_ex)


#>>>>>>>>>>>>>>>>>>>>> MAIN PROGRAM <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#>>>>>>>>> Fisrt part microstates populations and energies <<<<<<

# It turns out that for this first part we need to use MS containing element-wise states:
# so:
    ms = []
    for i in dist_ex:
        ms.append([i])

   # and now use ms instead of MS, but only for the first part

    if opt==1 or opt==12:

        for in_ex in dist_ex:
            print "initial excitation = ", in_ex

            # Initialize populations for given initial excitation
            print "initialization..."
            num_runs = 0.0

            P = [] #sh populations
            C = [] #SE populations
            EP = [] # sh-based energies
            EC = [] # se-based energies

            # Warning: The following initialization of 2-D lists is correct
            # We need create separate lists, so the loops should not be decoupled
            t = 0
            while t<namdtime:
                p,c = [],[]
                j = 0
                while j<len(ms):
                    p.append(0.0)
                    c.append(0.0)
                    j = j + 1

                P.append(p)
                C.append(c)
                EP.append([0.0])
                EC.append([0.0])
                t = t + 1

            # Look only for those files which correspond to given initial excitation
            print "calculation..."
            i = 0
            while i<len(all_ex):
                if all_ex[i]==in_ex:
                    num_runs = num_runs + 1.0
                    # Then file out<i> -corresponds to given initial excitation

                    #===================== me_energies file =========================
                    eE = get_file(inp_dir+"/me_energies",i,namdtime,num_states,5,2)

                    #======================= out file ===============================
                    pP = get_file(inp_dir+"/out",i,namdtime,num_states,3,2)
                    pe = sum_mult_arrays(eE,pP)
                    #P = add_arrays(P,pP)
                    EP = add_arrays(EP,pe)
                    #print len(pP),len(pP[0]),len(MS),len(MS[0])
                    cpP = contract_array(pP,ms) # Dimension:  T x num_macro_states
                    P = add_arrays(P,cpP)


                    #===================== me_pop file ==============================
                    cC = get_file(inp_dir+"/me_pop",i,namdtime,num_states,3,2)
                    ce = sum_mult_arrays(eE,cC)
                    #C = add_arrays(C,cC)
                    EC = add_arrays(EC,ce)
                    ccC = contract_array(cC,ms)  # Dimension: T x num_macro_states
                    C = add_arrays(C,ccC)

                    #================================================================

                i = i + 1
            print "Number of runs for this excitations = ", num_runs

#            write_array(res_dir+"/sh_pop_ex",in_ex,namdtime,num_states,P,num_runs)
#            write_array(res_dir+"/se_pop_ex",in_ex,namdtime,num_states,C,num_runs)
            write_array(res_dir+"/sh_pop_ex",in_ex,namdtime,len(ms),P,num_runs)
            write_array(res_dir+"/se_pop_ex",in_ex,namdtime,len(ms),C,num_runs)

            write_array(res_dir+"/sh_en_ex",in_ex,namdtime,1,EP,num_runs)
            write_array(res_dir+"/se_en_ex",in_ex,namdtime,1,EC,num_runs)



#>>>>>>>> Second part Macrostates populations and energies <<<<<<
    if opt==2 or opt==12:
        num_macro_states = len(MS)

        i = 0
        while i<num_macro_states:
            print "initial macro-excitation = ", i

            # Initialize populations for given initial excitation
            print "initialization..."

            P = [] #sh populations of macrostates
            C = [] #SE populations of macrostates
            En = [] # energies of macrostates

            # Warning: The following initialization of 2-D lists is correct
            # We need create separate lists, so the loops should not be decoupled
            t = 0
            while t<namdtime:
                p,c,en = [],[],[]
                j = 0
                while j<num_macro_states:
                    p.append(0.0)
                    c.append(0.0)
                    en.append(0.0)
                    j = j + 1

                P.append(p)
                C.append(c)
                En.append(en)
                t = t + 1

            # Look only for those files which correspond to given initial excitation
            print "calculation..."
            j = 0
            num_runs = len(MS[i])  # Number of microstates corresponding to given macrostate
            print "Number of microstates for this macrostate = ", num_runs

            while j<num_runs:
                # Then file sh_pop<MS[i][j]> -corresponds to given initial excitation
                #======================= energy file ===============================
                pEn = get_file(inp_dir+"/me_energies",MS[i][j],namdtime,num_states,5,2) # Dimension: T x num_states
                cpEn = contract_array(pEn,MS) # Dimension:  T x num_macro_states
                En = add_arrays(En,cpEn)

                #======================= out file ===============================
                pP = get_file(res_dir+"/sh_pop_ex",MS[i][j],namdtime,num_states,3,2) # Dimension: T x num_states
                cpP = contract_array(pP,MS) # Dimension:  T x num_macro_states
                P = add_arrays(P,cpP)
  
                #===================== me_pop file ==============================
                cC = get_file(res_dir+"/se_pop_ex",MS[i][j],namdtime,num_states,3,2) # Dimension: T x num_states
                ccC = contract_array(cC,MS)  # Dimension: T x num_macro_states
                C = add_arrays(C,ccC)

                #================================================================
                j = j + 1

            write_array(res_dir+"/energy_macro_ex",i,namdtime,num_macro_states,En,num_runs)
            write_array(res_dir+"/sh_pop_macro_ex",i,namdtime,num_macro_states,P,num_runs)
            write_array(res_dir+"/se_pop_macro_ex",i,namdtime,num_macro_states,C,num_runs)

            i = i + 1


# End of function
