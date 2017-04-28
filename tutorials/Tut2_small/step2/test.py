import os
import sys
import math

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
from libra_py import *

act_space = [6,7]

#wd = os.getcwd()+"/wd_test"
wd = "wd_test"
#wfc_curr = QE_methods.read_qe_wfc("%s/curr0/x0.export/wfc.1" % wd , "Kpoint.1", act_space)
#wfc_next = QE_methods.read_qe_wfc("%s/next0/x0.export/wfc.1" % wd , "Kpoint.1", act_space)

#e_curr = QE_methods.read_qe_index("%s/curr0/x0.export/index.xml" % wd , act_space )

e_curr = QE_methods.read_qe_index("%s/curr0/x0.export/index.xml" % wd , act_space )
e_next = QE_methods.read_qe_index("%s/next0/x0.export/index.xml" % wd, act_space )

