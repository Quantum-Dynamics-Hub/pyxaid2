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

#include "matrix.h"
#include "state.h"
#include "units.h"
#include "InputStructure.h"
#include <boost/python.hpp>
using namespace boost::python;
using namespace std;



class ElectronicStructure{

  // For DISH
  void update_decoherence_times(matrix& rates);
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
  matrix* Ccurr;
  matrix* Cprev;
  matrix* Cnext;
  matrix* A;                      // density matrix - populations and coherences

  // Hamiltonian: meaning of "current" and "next" depends on the algorithm to solve TD-SE
  matrix* Hcurr; // current Hamiltonian                             Hij = H[i*num_states+j]
  matrix* Hprev; // Hamiltonian on previous nuclear time step                 ---
  matrix* Hnext; // Hamiltonian on next nuclear time step                     ---
  matrix* dHdt;  // slope of Hamiltonian - for interpolation scheme           ---

  matrix* Hprimex;
  matrix* Hprimey;
  matrix* Hprimez;

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
    Ccurr = new matrix(n,1); *Ccurr = tmp;
    Cprev = new matrix(n,1); *Cprev = tmp;
    Cnext = new matrix(n,1); *Cnext = tmp;

    g = std::vector<double>(n*n,0.0);  // g[i*n+j] ~=g[i][j] - probability of i->j transition

    A = new matrix(n,n); *A = tmp;

    Hcurr = new matrix(n,n); *Hcurr = tmp;
    Hprev = new matrix(n,n); *Hprev = tmp;
    Hnext = new matrix(n,n); *Hnext = tmp;
    dHdt  = new matrix(n,n); *dHdt  = tmp;

    Hprimex = new matrix(n,n); *Hprimex = tmp;
    Hprimey = new matrix(n,n); *Hprimey = tmp;
    Hprimez = new matrix(n,n); *Hprimez = tmp;


    tau_m = std::vector<double>(n,0.0);
    t_m = std::vector<double>(n,0.0);
  }

  ElectronicStructure(const ElectronicStructure& es){ // Copy constructor
    int n = es.num_states;
    num_states = n;
    curr_state = es.curr_state;

    complex<double> tmp(0.0,0.0);
    vector<double> tmp1(n,0.0);
    Ccurr = new matrix(n,1);
    Cprev = new matrix(n,1);
    Cnext = new matrix(n,1);

    g = std::vector<double>(n*n,0.0);  // g[i*n+j] ~=g[i][j] - probability of i->j transition
    
    A = new matrix(n,n);

    Hcurr = new matrix(n,n);
    Hprev = new matrix(n,n);
    Hnext = new matrix(n,n);
    dHdt  = new matrix(n,n); 

    Hprimex = new matrix(n,n);
    Hprimey = new matrix(n,n);
    Hprimez = new matrix(n,n);


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
  void update_hop_prob(double dt,int is_boltz_flag,double Temp,matrix& Ef);


  void update_hop_prob_fssh(double dt,int is_boltz_flag,double Temp,matrix& Ef,double Ex, matrix&);
  void update_hop_prob_mssh(double dt,int is_boltz_flag,double Temp,matrix& Ef,double Ex, matrix&);
  void update_hop_prob_gfsh(double dt,int is_boltz_flag,double Temp,matrix& Ef,double Ex, matrix&);
  void init_hop_prob1(); 

  void check_decoherence(double dt,int boltz_flag,double Temp,matrix& rates); // practically DISH correction

  void propagate_coefficients(double dt,matrix& Ef);  // Trotter factorization
  void propagate_coefficients(double dt,matrix& Ef,matrix&);  // Trotter factorization with purostat
  void propagate_coefficients1(double dt,int opt,matrix& Ef); // Finite difference
  void propagate_coefficients2(double dt,matrix& Ef); // "Exact"

};


#endif // ElectronicStructure_H

