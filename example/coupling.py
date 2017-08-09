import math
import string

nsteps = 200

energy_file=open("Hab_re","w")
minum_band = 1
maximum_band = 2

for n in range(0,nsteps):
   Re = "res/Hvib_ab_"+str(n)+"_re"  
   FR = open(Re, "r")
   FRr = FR.readlines()
   FR.close()

   energy_file.write("%d   " % n )
   for i in range(minum_band,maximum_band+1):
      for j in range(minum_band,maximum_band+1):
       temp=[]
       temp = FRr[i - 1].split()
       val1 = 13.606*float(temp[j -1])*1000
       energy_file.write("   %f    " % val1)
   energy_file.write(" \n ")





