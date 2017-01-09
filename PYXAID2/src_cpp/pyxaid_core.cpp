/***********************************************************
 * Copyright (C) 2013 PYXAID group
 * This file is distributed under the terms of the
 * GNU General Public License as published by the
 * Free Software Foundation; either version 3 of the
 * License, or (at your option) any later version.
 * http://www.gnu.org/copyleft/gpl.txt
***********************************************************/

#include <boost/python.hpp>
#include <iostream>
#include <iomanip>
#include "wfc_export.h"
#include "namd_export.h"
#include "matrix.h"
using namespace std;
using namespace boost::python;


class info{

  public:

  // Constructor 
  info(){}

  // Class members
  void version(){
    cout<<"================================================================================\n";
    cout<<"PYXAID: PYthon eXtension for Ab Inition Dynamics version 1.0\n";
    cout<<"/***********************************************************\n";
    cout<<" * Copyright (C) 2013 PYXAID group\n";
    cout<<" * This program is free software distributed under the terms of the\n";
    cout<<" * GNU General Public License as published by the\n";
    cout<<" * Free Software Foundation; either version 3 of the\n";
    cout<<" * License, or (at your option) any later version.\n";
    cout<<" * http://www.gnu.org/copyleft/gpl.txt\n";
    cout<<"***********************************************************/\n";
  }

  void developers(){
    cout<<"===== Name ========================== Contact info ==============\n";
    cout<<"Alexey V. Akimov              e-mail: alexvakimov@gmail.com      \n";
    cout<<"Oleg V. Prezhdo               e-mail: oprezhdo@chem.rochester.edu\n";
    cout<<"...\n";
  }

  void documentation(){
    cout<<"Coming soon...\n";
  }
 
};



BOOST_PYTHON_MODULE(pyxaid_core)
{

    class_<info>("info",init<>())
        .def("version",&info::version)
        .def("developers",&info::developers)
        .def("documentation", &info::documentation)
    ;

    export_wfc();
    export_namd();
}

