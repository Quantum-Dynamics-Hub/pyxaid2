#***********************************************************
# * Copyright (C) 2013 Alexey V. Akimov
# * This file is distributed under the terms of the
# * GNU General Public License as published by the
# * Free Software Foundation; either version 3 of the
# * License, or (at your option) any later version.
# * http://www.gnu.org/copyleft/gpl.txt
#***********************************************************/

def ground_state(Nmin,Nmax):
    st = []
    for i in range(Nmin,Nmax+1):
        st.append(i)
        st.append(-i)
    return ["GS",st]


def ground_state1(occ_orb):
# from the list
    st = []
    for i in occ_orb:
        st.append(i)
        st.append(-i)
    return ["GS",st]


def single_excitations(Nmin,Nmax,HOMO,nspin):
# form all single excitations from [Nmin,HOMO] to [LUMO,Nmax]
# where LUMO = HOMO + 1
# nspin - defines if we want all spin orientations(2) or only one(1)

    st = []  # states
    LUMO = HOMO + 1

    cnt = 0
    for j in range(LUMO,Nmax+1):
        for i in range(Nmin,HOMO+1):
            ex1  = [] # i->j
            ex2 = []
            for v in range(Nmin,HOMO+1):
                if v==i: 
                    ex1.append(j)
                    ex1.append(-i)
                    ex2.append(-j)
                    ex2.append(i)
                else:
                    ex1.append(v)
                    ex1.append(-v)
                    ex2.append(v)
                    ex2.append(-v)

            if nspin>=1:
                st.append(["SE"+str(cnt),ex1])
                cnt = cnt + 1
            if nspin>=2:
                st.append(["SE"+str(cnt),ex2])
                cnt = cnt + 1

    return st

def single_excitations1(occ_orb,virt_orb,nspin):
# form all single excitations from occ_orb to virt_orb
# nspin - defines if we want all spin orientations(2) or only one(1)

    st = []  # states
    #LUMO = HOMO + 1

    cnt = 0
    for j in virt_orb:
        for i in occ_orb:
            ex1  = [] # i->j
            ex2 = []
            for v in occ_orb:
                if v==i:
                    ex1.append(j)
                    ex1.append(-i)
                    ex2.append(-j)
                    ex2.append(i)
                else:
                    ex1.append(v)
                    ex1.append(-v)
                    ex2.append(v)
                    ex2.append(-v)

            if nspin>=1:
                st.append(["SE"+str(cnt),ex1])
                cnt = cnt + 1
            if nspin>=2:
                st.append(["SE"+str(cnt),ex2])
                cnt = cnt + 1

    return st


def double_excitations(Nmin,Nmax,HOMO,nspin):
# form all single excitations from [Nmin,HOMO] to [LUMO,Nmax]
# where LUMO = HOMO + 1
# nspin - defines if we want all spin orientations(2) or only one(1)

    st = []  # states
    LUMO = HOMO + 1

    cnt = 0
    for j1 in range(LUMO,Nmax+1):
        for j2 in range(LUMO,Nmax+1):
            for i1 in range(Nmin,HOMO+1):
                for i2 in range(Nmin,HOMO+1):
                    ex1 = [] # (i1,i2)->(j1,j2)
                    ex2 = []
                    ex3 = []
                    ex4 = [] 
                    for v in range(Nmin,HOMO+1):
                        if v==i1:
                            ex1.append(j1)
                            ex1.append(-i1)
                            ex2.append(j1)
                            ex2.append(-i1)

                            ex3.append(-j1)
                            ex3.append(i1)
                            ex4.append(-j1)
                            ex4.append(i1)

                        if v==i2:
                            ex1.append(j2)
                            ex1.append(-i2)
                            ex2.append(j2)
                            ex2.append(-i2)

                            ex3.append(-j2)
                            ex3.append(i2)
                            ex4.append(-j2)
                            ex4.append(i2)

                        else:
                            ex1.append(v)
                            ex1.append(-v)
                            ex2.append(v)
                            ex2.append(-v)
                            ex3.append(v)
                            ex3.append(-v)
                            ex4.append(v)
                            ex4.append(-v)

                    if(nspin>=1):
                    # For now include all configurations
                        st.append(["DE"+str(cnt),ex1])
                        cnt = cnt + 1
                        st.append(["DE"+str(cnt),ex2])
                        cnt = cnt + 1
                        st.append(["DE"+str(cnt),ex3])
                        cnt = cnt + 1
                        st.append(["DE"+str(cnt),ex4])
                        cnt = cnt + 1
                        
    return st


def spin(ex):
    nup = 0
    ndn = 0
    for a in ex:
        if a>0:
            nup = nup + 1
        else:
            ndn = ndn + 1

    return nup-ndn


def num_in(num,A):
# number of occurences of num in list A
    sz = len(A)
    res = 0
    for a in A:
        if a==num:
            res = res + 1
    return res

def is_equal(ex,a):
# ex and a - are lists
# they are equal if each element of ex is contained in a the same number of times 

    res = 1
    sz = len(ex)
    i = 0
    while i<sz:
        if num_in(ex[i],ex)!=num_in(ex[i],a):
            res = 0
        i = i + 1

    return res



def is_included(ex,all):
    res = 0
    for a in all:
        if is_equal(ex,a)==1:
            res = 1

    return res

def double_excitations1(o_orb,v_orb,nspin):
# form all single excitations from occ_orb to virt_orb
# nspin - defines if we want all spin orientations(2) or only one(1)

    st = []  # states
    #LUMO = HOMO + 1

    occ_orb = []
    virt_orb = []

    # Convert to more excplicit ranges
    for j1 in o_orb:
        occ_orb.append(j1)
        occ_orb.append(-j1)
    for j1 in v_orb:
        virt_orb.append(j1)
        virt_orb.append(-j1)


    nvirt = len(virt_orb)
    nocc  = len(occ_orb)
    cnt = 0

    all = []
    # Outer loop goes over all virtual orbitals, for double excitation
    # this is not a problem because there can be up to two electrons
    # with the same spin on the same spatial orbital
    for j1 in v_orb:
        for j2 in v_orb:
            # Portion from the occ orbitals which remain after double excitations    
            i1 = 0
            while i1<nocc:
                i2 = 0
                while i2<nocc:
                    if i2!=i1:
                        # i1, i2 are now proper indexes of the occ orbitals being removed
                        # So now add all occupied orbitals except for those with indexes i1,i2               
                        ex1 = [j1,j2]
                        ex2 = [j1,-j2]
                        ex3 = [-j1,j2]
                        ex4 = [-j1,-j2]
  
                        i = 0
                        while i<nocc:
                            if i!=i1 and i!=i2:
                                ex1.append(occ_orb[i])
                                ex2.append(occ_orb[i])
                                ex3.append(occ_orb[i])
                                ex4.append(occ_orb[i])
                            i = i + 1

                        
                        if spin(ex1)==0:
                            if is_included(ex1,all)==0:
                                st.append(["DE"+str(cnt),ex1])
                                all.append(ex1)
                                cnt = cnt + 1
                        if spin(ex2)==0:
                            if is_included(ex2,all)==0:
                                st.append(["DE"+str(cnt),ex2])
                                all.append(ex2)
                                cnt = cnt + 1
                        if spin(ex3)==0 and j1!=j2:
                            # in case j1==j2 the configuration ex3 is identical to ex2
                            if is_included(ex3,all)==0:
                                st.append(["DE"+str(cnt),ex3])
                                all.append(ex3)
                                cnt = cnt + 1
                        if spin(ex4)==0:
                            if is_included(ex4,all)==0:
                                st.append(["DE"+str(cnt),ex4])
                                all.append(ex4)
                                cnt = cnt + 1




                    i2 = i2 + 1
                i1 = i1 + 1


    return st

