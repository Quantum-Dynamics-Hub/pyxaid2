#!/bin/bash
#PBS -l nodes=1:ppn=5
#PBS -N qespresso_NAC
#PBS -j oe
#PBS -q batch
#PBS -l walltime=30:00:00

cd $PBS_O_WORKDIR
echo $PBS_O_WORKDIR

# the below settings works for my laptop, please customize your own environmental variables!
###### some conventional settings ###########
export PYTHONPATH=/home/eric/src/Libra/libra-code/_build/src:$PYTHONPATH
export LD_LIBRARY_PATH=/home/eric/src/Libra/libra-code/_build/src:$LD_LIBRARY_PATH
export PYTHONPATH=/home/eric/src/Libra/pyxaid2:$PYTHONPATH
export LD_LIBRARY_PATH=/home/eric/src/Libra/pyxaid2:/opt/boost/1.55.0/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/opt/intel/composer_xe_2013.1.117/mkl/lib/intel64:$LD_LIBRARY_PATH
##### path for the QE module #####
exe_qespresso=/home/eric/src/qe/espresso-5.1.2/bin/pw.x
exe_export=/home/eric/src/qe/espresso-5.1.2/bin/pw_export.x
exe_convert=/home/eric/src/qe/espresso-5.1.2/bin/iotk
##### MPI #####
MPIRUN=/home/eric/src/ompi/1.8.5-gcc/bin/mpirun
NP=$(wc -l $PBS_NODEFILE | awk '{print $1}')

# These will be assigned automatically, leave them as they are
param1=
param2=

res=/home/eric/src/Libra/for_develop/pyxaid2/example/1-non-rel-sp/res

# This is invocation of the scripts which will further handle NA-MD calclculations
# on the NAC calculation step
python -c "from PYXAID2 import *
params = { }
params[\"BATCH_SYSTEM\"]=\"$MPIRUN\"
params[\"NP\"]=$NP
params[\"EXE\"]=\"$exe_qespresso\"
params[\"EXE_EXPORT\"]=\"$exe_export\"
params[\"EXE_CONVERT\"] =\"$exe_convert\"
params[\"start_indx\"]=$param1
params[\"stop_indx\"]=$param2
params[\"wd\"]=\"wd_test\"
params[\"rd\"]=\"$res\"
params[\"minband\"]=36
params[\"maxband\"]=37
params[\"minband_soc\"]=36
params[\"maxband_soc\"]=37
params[\"nac_method\"]=1
params[\"wfc_preprocess\"]=\"complete\"
params[\"do_complete\"]=1
params[\"prefix0\"]=\"x0.scf\"
params[\"prefix1\"]=\"x1.scf\"
params[\"compute_Hprime\"]=0
params[\"pptype\"]=\"US\"
print params
runMD2.runMD(params)
"




