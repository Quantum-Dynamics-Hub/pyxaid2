/***********************************************************
 * Copyright (C) 2013 Alexey V. Akimov
 * This file is distributed under the terms of the
 * GNU General Public License as published by the
 * Free Software Foundation; either version 3 of the
 * License, or (at your option) any later version.
 * http://www.gnu.org/copyleft/gpl.txt
***********************************************************/

#ifndef state_h
#define state_h

#include <vector>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <boost/python.hpp>
using namespace boost::python;
using namespace std;

int ext2int(int,vector<int>&);
int delta(vector<int>& A,vector<int>& B,int& a,int& b);

  
class me_state{
// Multi-electron state
// Basically it is a Slater product

public:
  std::string name;          // label of the determinant
  // Data
  vector<int> active_space;  // All occupied and virtual orbitals involved in dynamics
  vector<int> actual_state;  // This is an actual(current) state

  double Eshift; // this correction goes to Exc and is read directly from input - this is
                 // to simplify parameter development
  double Exc;    // correlation correction, to better describe the energy
                 // of the state with 2 electrons on the same orbital

  // NAC scaling constants for different pairs of the states
  // the sizes of the below vectors should be equal and elements be in correspondense
  vector<int> nac_scl_indx;
  vector<double> nac_scl;


  // Constructors
  me_state(){}
  me_state(vector<int>& as_,vector<int>& cs_){ active_space = as_; actual_state = cs_; Exc = 0.0; }

  // Basically the constructor
  void set_me_state(vector<int>& as_,vector<int>& cs_){ active_space = as_; actual_state = cs_; Exc = 0.0;}
  
  // Destructor
  ~me_state() { ; ;}

  // Functions
  int calculate_Exc(vector<int>&, vector<int>&, vector<double>&,vector<int>&, vector<double>&);
  void show_state();

};


void input_iconds(boost::python::dict params,int me_numstates,vector<vector<int> >& icond);
void input_states(boost::python::dict params,vector<me_state>& states);



#endif // state_h
