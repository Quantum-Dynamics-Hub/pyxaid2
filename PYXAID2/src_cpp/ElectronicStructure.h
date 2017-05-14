/***********************************************************
 * Copyright (C) 2013 Alexey V. Akimov
 * This file is distributed under the terms of the
 * GNU General Public License as published by the
 * Free Software Foundation; either version 3 of the
 * License, or (at your option) any later version.
 * http://www.gnu.org/copyleft/gpl.txt
***********************************************************/

#ifndef ElectronicStructure_H
#define ElectronicStructure_H

//#include "matrix.h"
#include "state.h"
#include "units_pyxaid.h"
#include "InputStructure.h"
#include <boost/python.hpp>
using namespace boost::python;
using namespace std;

#include "liblibra_core.h"
using namespace liblibra::librandom;
using namespace liblibra::liblinalg;


class ElectronicStructure{

  // For DISH
  void update_decoherence_times(CMATRIX& rates);
  void project_out(int i);
  void decohere(int i); 

  // For integrators
  void rot1(double phi,int i,int j);
  void rot2(double phi,int i,int j);
  void rot(complex<double> Hij,double dt,int i,int j);
  void phase(complex<double> Hii,double dt,int i);

public:

  //=========== Members ===============
  int num_states;                 // number of adiabatic states
  int curr_state;                 // current adiabatic state

  // Wavefunction
  CMATRIX* Ccurr;
  CMATRIX* Cprev;
  CMATRIX* Cnext;
  CMATRIX* A;                      // density matrix - populations and coherences

  // Hamiltonian: meaning of "current" and "next" depends on the algorithm to solve TD-SE
  CMATRIX* Hcurr; // current Hamiltonian                             Hij = H[i*num_states+j]
  CMATRIX* Hprev; // Hamiltonian on previous nuclear time step                 ---
  CMATRIX* Hnext; // Hamiltonian on next nuclear time step                     ---
  CMATRIX* dHdt;  // slope of Hamiltonian - for interpolation scheme           ---

  CMATRIX* Hprimex;
  CMATRIX* Hprimey;
  CMATRIX* Hprimez;

  vector<double> g; // num_states x num_states matrix, reshaped in 1D array

  // DISH variables:
  vector<double> tau_m; // times since last decoherence even for all PES (actually rates, that is inverse times)
  vector<double> t_m;   // time counters for each PES
  


  //========== Methods ================
  // Constructors
  ElectronicStructure(int n){ 
    num_states = n;

    complex<double> tmp(0.0,0.0);
    vector<double> tmp1(n,0.0);
    Ccurr = new CMATRIX(n,1); *Ccurr = tmp;
    Cprev = new CMATRIX(n,1); *Cprev = tmp;
    Cnext = new CMATRIX(n,1); *Cnext = tmp;

    g = std::vector<double>(n*n,0.0);  // g[i*n+j] ~=g[i][j] - probability of i->j transition

    A = new CMATRIX(n,n); *A = tmp;

    Hcurr = new CMATRIX(n,n); *Hcurr = tmp;
    Hprev = new CMATRIX(n,n); *Hprev = tmp;
    Hnext = new CMATRIX(n,n); *Hnext = tmp;
    dHdt  = new CMATRIX(n,n); *dHdt  = tmp;

    Hprimex = new CMATRIX(n,n); *Hprimex = tmp;
    Hprimey = new CMATRIX(n,n); *Hprimey = tmp;
    Hprimez = new CMATRIX(n,n); *Hprimez = tmp;


    tau_m = std::vector<double>(n,0.0);
    t_m = std::vector<double>(n,0.0);
  }

  ElectronicStructure(const ElectronicStructure& es){ // Copy constructor
    int n = es.num_states;
    num_states = n;
    curr_state = es.curr_state;

    complex<double> tmp(0.0,0.0);
    vector<double> tmp1(n,0.0);
    Ccurr = new CMATRIX(n,1);
    Cprev = new CMATRIX(n,1);
    Cnext = new CMATRIX(n,1);

    g = std::vector<double>(n*n,0.0);  // g[i*n+j] ~=g[i][j] - probability of i->j transition
    
    A = new CMATRIX(n,n);

    Hcurr = new CMATRIX(n,n);
    Hprev = new CMATRIX(n,n);
    Hnext = new CMATRIX(n,n);
    dHdt  = new CMATRIX(n,n); 

    Hprimex = new CMATRIX(n,n);
    Hprimey = new CMATRIX(n,n);
    Hprimez = new CMATRIX(n,n);


    tau_m = es.tau_m;
    t_m = es.t_m;

    *Ccurr = *es.Ccurr; *Cprev = *es.Cprev; *Cnext = *es.Cnext;
    g = es.g;  *A = *es.A;
    *Hcurr = *es.Hcurr;  *Hprev = *es.Hprev;  *Hnext = *es.Hnext;  *dHdt  = *es.dHdt;
    *Hprimex = *es.Hprimex; *Hprimey = *es.Hprimey; *Hprimez = *es.Hprimez;
  }
  // Destructor

  ~ElectronicStructure(){
    if(g.size()>0) {g.clear();}
    if(Ccurr!=NULL) { delete Ccurr;}// Ccurr = NULL;}
    if(Cprev!=NULL) {delete Cprev;}// Cprev = NULL;}
    if(Cnext!=NULL) {delete Cnext;}// Cnext = NULL;}
    if(A!=NULL) { delete A;} // A = NULL;}
    if(Hcurr!=NULL){ delete Hcurr;} // Hcurr = NULL;}
    if(Hprev!=NULL){ delete Hprev;}// Hprev = NULL;}
    if(Hnext!=NULL){ delete Hnext;}// Hnext = NULL;}
    if(dHdt!=NULL){ delete dHdt;} // dHdt = NULL;} 
    if(Hprimex!=NULL){ delete Hprimex; }
    if(Hprimey!=NULL){ delete Hprimey; }
    if(Hprimez!=NULL){ delete Hprimez; }
    if(tau_m.size()>0){ tau_m.clear(); }
    if(t_m.size()>0){ t_m.clear(); }
  }


  ElectronicStructure operator=(ElectronicStructure es){
    num_states = es.num_states;
    curr_state = es.curr_state;
   *Ccurr = *es.Ccurr; *Cprev = *es.Cprev; *Cnext = *es.Cnext;
    g = es.g;  *A = *es.A;
    *Hcurr = *es.Hcurr;  *Hprev = *es.Hprev; *Hnext = *es.Hnext;
    *Hprimex = *es.Hprimex; *Hprimey = *es.Hprimey; *Hprimez = *es.Hprimez; 
    *dHdt  = *es.dHdt;
    tau_m = es.tau_m;  t_m = es.t_m;
    return *this;
  }

  ElectronicStructure operator<<(ElectronicStructure es){
// This is basically the same as operator=, but keeps parameters the same
// copied only wavefunction and state

    num_states = es.num_states;
    curr_state = es.curr_state;
    *Cprev = *es.Cprev;
    *Ccurr = *es.Ccurr;
    *Cnext = *es.Cnext;  

// Also need to update the DISH timescales
   tau_m = es.tau_m;
   t_m = es.t_m;

    return *this;
  }

 

  // Other methods
  void set_state(int indx){
    for(int i=0;i<num_states;i++){
      if(i==indx){ Ccurr->M[i] = complex<double>(1.0,0.0); } else{  Ccurr->M[i] = 0.0; }
    }
    curr_state = indx;
  }
  double energy();                // calculate the total energy
  double norm(); // calculate the norm of the wavefunction

  void update_populations();      // update matrix A from Ccurr
  void update_hop_prob(double dt,int is_boltz_flag,double Temp, CMATRIX& Ef);


  void update_hop_prob_fssh(double dt,int is_boltz_flag,double Temp, CMATRIX& Ef,double Ex, CMATRIX&);
  void update_hop_prob_mssh(double dt,int is_boltz_flag,double Temp, CMATRIX& Ef,double Ex, CMATRIX&);
  void update_hop_prob_gfsh(double dt,int is_boltz_flag,double Temp, CMATRIX& Ef,double Ex, CMATRIX&);
  void init_hop_prob1(); 

  void check_decoherence(double dt,int boltz_flag,double Temp, CMATRIX& rates, Random& rnd); // practically DISH correction

  void propagate_coefficients(double dt, CMATRIX& Ef);  // Trotter factorization
  void propagate_coefficients(double dt, CMATRIX& Ef, CMATRIX&);  // Trotter factorization with purostat
  void propagate_coefficients1(double dt,int opt, CMATRIX& Ef); // Finite difference
  void propagate_coefficients2(double dt, CMATRIX& Ef); // "Exact"

};


#endif // ElectronicStructure_H

