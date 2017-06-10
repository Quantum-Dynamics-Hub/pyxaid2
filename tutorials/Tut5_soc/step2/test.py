import os
import sys
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *

act_space = [6,7,8,9]

#wd = os.getcwd()+"/wd_test"
wd = "wd/job0/wd_test"
wfc_curr = QE_methods.read_qe_wfc("%s/curr0/x0.export/wfc.1" % wd , "Kpoint.1", act_space)
wfc_next = QE_methods.read_qe_wfc("%s/next0/x0.export/wfc.1" % wd , "Kpoint.1", act_space)

#e_curr = QE_methods.read_qe_index("%s/curr0/x0.export/index.xml" % wd , act_space )

e_curr = QE_methods.read_qe_index("%s/curr0/x0.export/index.xml" % wd , act_space )
e_next = QE_methods.read_qe_index("%s/next0/x0.export/index.xml" % wd, act_space )


S = wfc_curr.H() * wfc_curr
print "S"
S.show_matrix()

sx = S.num_of_cols
ovlp = CMATRIX(sx/2, sx/2)
for n in xrange(sx/2):
    for k in xrange(sx/2):
        ovlp.set(n,k, S.get(2*n,2*k) + S.get(2*n+1,2*k+1) )

print "Ovlp"
ovlp.show_matrix()

