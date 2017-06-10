import os
import sys
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *

act_space = [91,92,93,94]

#wd = os.getcwd()+"/wd_test"
wd = "wd/job0/wd_test"
coeff_1 = QE_methods.read_qe_wfc("%s/curr0/x0.export/wfc.1" % wd , "Kpoint.1", act_space)

nbnd = coeff_1.num_of_cols
print "The number of bands = ", nbnd
 
ovlp_1  = coeff_1.H() * coeff_1


print "<mix|mix> = "; ovlp_1.show_matrix()
print "\n"

for n in xrange(nbnd):
    print n, ovlp_1.get(n,n)
print "\n"

for n in xrange(nbnd/2):
    print n, ovlp_1.get(2*n,2*n) + ovlp_1.get(2*n+1,2*n+1)
print "\n"

ovlp = CMATRIX(nbnd/2, nbnd/2)
for n in xrange(nbnd/2):
    for k in xrange(nbnd/2):
        ovlp.set(n,k, ovlp_1.get(2*n,2*k) + ovlp_1.get(2*n+1,2*k+1) )

print "spinor overlap matrix = "; ovlp.show_matrix()





