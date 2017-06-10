from PYXAID import *
# Example of usage
opt = 1
scl1 = 13.60569253 # Ry to eV
scl2 = 1000.0 # some scaling to get plottable data
tmin = 0
tmax = 50


# Energy, couplings and H' in space of orbital indexes
excitation_spectrum.ham_map("res/0_Ham_",   tmin,tmax,"_re" ,opt,scl1,"spectr/ave_Ham_re.dat")
excitation_spectrum.ham_map("res/0_Ham_",   tmin,tmax,"_im" ,opt,scl1,"spectr/ave_Ham_im.dat")

# Same, but in space of orbital energies
excitation_spectrum.ham_map1("res/0_Ham_","_re","res/0_Ham_",   tmin,tmax,"_re" ,opt,scl1,scl1,"spectr/1ave_Ham_re.dat")
excitation_spectrum.ham_map1("res/0_Ham_","_re","res/0_Ham_",   tmin,tmax,"_im" ,opt,scl1,scl1,"spectr/1ave_Ham_im.dat")

