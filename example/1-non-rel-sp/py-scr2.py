from PYXAID2 import *
import os

user = 1 #weili or eric
#user = 2 #alexey

nsteps_per_job = 150
tot_nsteps = 300

# Step 1 - split MD trajectory on the time steps
# Provide files listed below: "GaAs-md.out" and "x0.scf.in", "x1.scf.in"
# IMPORTANT: 
# 1) use only ABSOLUTE path for PP in x.scf.in file
# 2) provided input file is just a template, so do not include coordinates
rt="/home/eric/src/Libra/for_develop/pyxaid2/example/"
out2inp.out2inp(rt+"GaAs-md.out","x0.scf.in","wd","x0.scf",0,tot_nsteps,1)  # non-relativistic setup
#out2inp.out2inp(rt+"GaAs-md.out","x1.scf.in","wd","x1.scf",0,tot_nsteps,1)  # relativistic setup


# Step 2 - distribute all time steps into groups(jobs) 
# several time steps per group - this is to accelerate calculations
# creates a "customized" submit file for each job and submit it - run
# a swarm in independent calculations (trajectory pieces)
#(HTC paradigm)
# Provide the files below: 
# submit_templ.pbs - template for submit files - manually edit the variables
# x.exp.in - file for export of the wavefunction

if user==1: 
   os.system("cp submit.pbs wd")
elif user==2:
   os.system("cp submit_templ.slm wd")

os.system("cp x0.exp.in wd")
#os.system("cp x1.exp.in wd")
os.chdir("wd")

if user==1:
   #distribute.distribute(0,tot_nsteps,nsteps_per_job,"submit.pbs",["x0.exp.in","x1.exp.in"],["x0.scf","x1.scf"],1)
   distribute.distribute(0,tot_nsteps,nsteps_per_job,"submit.pbs",["x0.exp.in"],["x0.scf"],1)
   #distribute.distribute(0,tot_nsteps,nsteps_per_job,"submit.pbs",["x1.exp.in"],["x1.scf"],1)
elif user==2:
   #distribute.distribute(0,tot_nsteps,nsteps_per_job,"submit_templ.slm",["x0.exp.in","x1.exp.in"],["x0.scf","x1.scf"],2)
   #distribute.distribute(0,tot_nsteps,nsteps_per_job,"submit_templ.slm",["x0.exp.in"],["x0.scf"],2)
   #distribute.distribute(0,tot_nsteps,nsteps_per_job,"submit_templ.slm",["x1.exp.in"],["x1.scf"],2)



