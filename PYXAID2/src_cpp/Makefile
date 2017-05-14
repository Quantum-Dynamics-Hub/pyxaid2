#FLAGS= -std=c++11 -O2 -fPIC # -fno-for-scope -O2 -fPIC
FLAGS= -fPIC -O2 -fpermissive -w -g
CPP=c++

# BOOST
# UB CCR
p1=/util/academic/boost/v1.57.0/include/boost
pl1=/util/academic/boost/v1.57.0/lib

# PYTHON
# UB CCR
p2=/util/academic/python/2.7.6/include/python2.7
pl2=/util/academic/python/2.7.6/lib

# LIBRA
p3=/user/alexeyak/Soft/libra-code/src
pl3=/user/alexeyak/Soft/libra-code/_build/src
pl3a=/user/alexeyak/Soft/libra-code/_build/src/util


I=-I ${p1} -I ${p2} -I ${p3}
L=-L ${pl1} -L ${pl2} -L ${pl3} -L ${pl3a}

all: pyxaid_core.so

#random.o: random.cpp random.h
#	${CPP} ${FLAGS} ${I} -c random.cpp

#matrix.o: matrix.cpp matrix.h
#	${CPP} ${FLAGS} ${I} -c matrix.cpp

state.o: state.cpp state.h
	${CPP} ${FLAGS} ${I} -c state.cpp

io.o: io.cpp io.h
	${CPP} ${FLAGS} ${I} -c io.cpp

InputStructure.o: InputStructure.cpp InputStructure.h
	${CPP} ${FLAGS} ${I} -c InputStructure.cpp

ElectronicStructure.o: ElectronicStructure.cpp ElectronicStructure.h
	${CPP} ${FLAGS} ${I} -c ElectronicStructure.cpp

namd.o: namd.cpp namd.h
	${CPP} ${FLAGS} ${I} -c namd.cpp

namd_export.o: namd_export.cpp namd_export.h
	${CPP} ${FLAGS} ${I} -c namd_export.cpp

#aux.o: aux.cpp aux.h
#	${CPP} ${FLAGS} ${I} -c aux.cpp

wfc_basic_methods.o: wfc_basic_methods.cpp wfc.h
	${CPP} ${FLAGS} ${I} -c wfc_basic_methods.cpp

wfc_QE_methods.o: wfc_QE_methods.cpp wfc.h 
	${CPP} ${FLAGS} ${I} -c wfc_QE_methods.cpp

wfc_functions.o: wfc_functions.cpp wfc.h
	${CPP} ${FLAGS} ${I} -c wfc_functions.cpp

wfc_export.o: wfc_basic_methods.o wfc_QE_methods.o wfc_functions.o wfc_export.cpp
	${CPP} ${FLAGS} ${I} -c wfc_export.cpp

pyxaid_core.o: pyxaid_core.cpp 
	${CPP} ${FLAGS} ${I} -c pyxaid_core.cpp



pyxaid_core.so: pyxaid_core.o wfc_export.o wfc_functions.o wfc_QE_methods.o wfc_basic_methods.o \
        state.o ElectronicStructure.o namd.o namd_export.o InputStructure.o io.o
	${CPP} ${FLAGS} ${I} -shared -o pyxaid_core.so pyxaid_core.o wfc_export.o wfc_functions.o \
        wfc_QE_methods.o wfc_basic_methods.o state.o ElectronicStructure.o namd.o \
        namd_export.o InputStructure.o io.o ${L} -lboost_python -llibra_core -lutil
	cp pyxaid_core.so ../.

clean:
	rm *.o
	rm pyxaid_core.so
	rm ../pyxaid_core.so

