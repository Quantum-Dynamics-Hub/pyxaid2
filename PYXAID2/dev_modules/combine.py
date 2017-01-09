#
# This module converts the output files (energies and "real" files - couplings) computed
# by older scripts to the new format, supported by PYXAID

def split_energy(en_file,res_dir,imin):
# Split a single file containing energies to a number of
# smaller files, each of which contains energies for a given 
# step. Resulting files get indexes from imin to imin + # of lines in en_file

    f = open(en_file,"r")
    A = f.readlines()
    f.close()

    sz = len(A)
    i = 0
    while i<sz:
        f = open("%senergy%i" % (res_dir,imin+i),"w")
        f.write(A[i])
        f.close()
        i = i + 1

    return sz


def combine_old_VASP(en_file,real_file,Ham_re,Ham_im):
# This is for compatibility with group VASP-based codes

    f = open(en_file,"r")
    A = f.readlines()
    f.close()

    xa = A[0].split()
    E = []
    for x in xa:
        E.append(float(x))

    # Number of states
    sz = len(E)

    # Real files
    f = open(real_file,"r")
    B = f.readlines()
    f.close()

    H_re = []
    for i in range(0,sz):
        xb = B[i].split()
        h = []
        for j in range(0,sz):
            h.append(float(xb[j]))
        H_re.append(h)


    # In the instructios below we assume NACs are d_ij = <i|d/dt|j>, so
    # to get Ham we multiply by -hbar, hbar is chosen to give output
    # in Ry units, for simplicity
    # Energies are alreay in Ry, so do not change them

    # Print Ham_re
    # Note here we are doing conversion of im to re and vice versa
    f = open(Ham_re,"w")
    for i in range(0,sz):
        line = ""
        for j in range(0,sz):
            if i==j:
                line = line + "%15.10f " % (E[i] / 13.60569253) # older VASP-based codes report energies in eV
            else:
                line = line + "%15.10f " % ((-1.0/13.60569253) * 0.0 )
        f.write(line+"\n")
    f.close()

    # Print Ham_im
    # Note here we are doing conversion of im to re and vice versa
    f = open(Ham_im,"w")
    for i in range(0,sz):
        line = ""
        for j in range(0,sz):
            if i==j:
                line = line + "%15.10f " % 0.0
            else:
                line = line + "%15.10f " % ((-1.0/13.60569253) * H_re[i][j])
        f.write(line+"\n")
    f.close()





def combine_old_QE(en_file,nac_re_file,nac_im_file,Ham_re,Ham_im):
# This is for compatibility with my older QE-based PYXAID versions

    f = open(en_file,"r")
    A = f.readlines()
    f.close()

    xa = A[0].split()
    E = []
    for x in xa:
        E.append(float(x))

    # Number of states
    sz = len(E)


    # NAC - real part
    f = open(nac_re_file,"r")
    B = f.readlines()
    f.close()
    
    H_re = []
    for i in range(0,sz):
        xb = B[i].split()   
        h = []
        for j in range(0,sz):
            h.append(float(xb[j]))
        H_re.append(h)

    # NAC - imaginary part
    f = open(nac_im_file,"r")
    C = f.readlines()
    f.close()

    H_im = []
    for i in range(0,sz):
        xc = C[i].split()
        h = []
        for j in range(0,sz):
            h.append(float(xc[j]))
        H_im.append(h)


    # In the instructios below we assume NACs are d_ij = <i|d/dt|j>, so
    # to get Ham we multiply by -hbar, hbar is chosen to give output
    # in Ry units, for simplicity
    # Energies are alreay in Ry, so do not change them

    # Print Ham_re
    # Note here we are doing conversion of im to re and vice versa
    f = open(Ham_re,"w")
    for i in range(0,sz):
        line = ""
        for j in range(0,sz):
            if i==j:
                line = line + "%15.10f " % E[i]
            else:
                line = line + "%15.10f " % ((-1.0/13.60569253) * H_im[i][j] )
        f.write(line+"\n")
    f.close()

    # Print Ham_im
    # Note here we are doing conversion of im to re and vice versa
    f = open(Ham_im,"w")
    for i in range(0,sz):
        line = ""
        for j in range(0,sz):
            if i==j:
                line = line + "%15.10f " % 0.0
            else:
                line = line + "%15.10f " % ((-1.0/13.60569253) * H_re[i][j])  
        f.write(line+"\n")
    f.close()


def combine_all_QE(imin,imax,rt):

    i = imin
    while i < imax:
        print i
        combine_old_QE("%senergy%i" % (rt,i), "%snac%i_re" % (rt,i), "%snac%i_im" % (rt,i), "%s0_Ham_%i_re" % (rt,i), "%s0_Ham_%i_im" % (rt,i))
        i = i + 1


def combine_all_VASP(en_file, res_dir,imin):
# res_dir - is the place where resulting 1-line energy files will be written
    sz = split_energy(en_file,res_dir,imin)

    i = imin
    while i < imin + sz:
        print i
        combine_old_VASP("%senergy%i" % (res_dir,i), "%sreal%i" % (res_dir,i), "%s0_Ham_%i_re" % (res_dir,i), "%s0_Ham_%i_im" % (res_dir,i))
        i = i + 1



#combine_all_QE(0,1578,"res/")
combine_all_VASP("energy","REAL_1000_1200/",1001)

