/***********************************************************
 * Copyright (C) 2017 Wei Li and Alexey V. Akimov
 * Copyright (C) 2013-2016 Alexey V. Akimov
 * This file is distributed under the terms of the
 * GNU General Public License as published by the
 * Free Software Foundation; either version 3 of the
 * License, or (at your option) any later version.
 * http://www.gnu.org/copyleft/gpl.txt
***********************************************************/

#include "namd.h"
//#include "aux.h"
#include "io.h"
//#include "random.h"

#include "liblibra_core.h"
using namespace liblibra;
using namespace liblibra::librandom;
using namespace liblibra::libutil;


/*****************************************************************
  Functions implemented in this file:

  void hop(vector<double>& sh_prob,int& hopstate,int numstates)
  void regression(vector<double>& X,vector<double>& Y,int opt,double& a,double& b)
  double decoherence_rates(vector<double>& x,double dt,std::string rt_dir,int regress_mode)
  void Efield(InputStructure& is,double t,matrix& E,double& Eex)
  void propagate_electronic(InputStructure& is,vector<ElectronicStructure>& es,int i, matrix& rates)
  void solve_electronic(InputStructure& is,vector<ElectronicStructure>& es,matrix& rates)
  void run_decoherence_rates(InputStructure& is, vector<ElectronicStructure>& me_es,vector<me_state>& me_states, int icond)
  void run_namd(InputStructure& is, vector<ElectronicStructure>& me_es,vector<me_state>& me_states, int icond) 
  void run_namd1(InputStructure& is, vector<ElectronicStructure>& me_es,vector<me_state>& me_states, int icond)

*****************************************************************/


void hop(vector<double>& sh_prob,int& hopstate,int numstates){
/***********************************************
 sh_prob[i] - is probability to hop from given state  to state i
 hopstate - will contain the state where we actually hopped
************************************************/
  int i;
  double left,right,ksi;

  int in = hopstate; // initial state
  int hstate = -1; // set to an absurd value, so that run fails explicitly if the
                   // surface hopping probabilities are stange
  ksi = rand()/((double)RAND_MAX);

  // But, to avoid the problems, lets renormalize the hopping probabilities
  double nrm = 0.0;
  for(i=0;i<numstates;i++){  nrm += sh_prob[in*numstates+i];  }  

  for(i=0;i<numstates;i++){
    if(i==0){ left = 0.0; right = (sh_prob[in*numstates+i]/nrm); }
    else{ left = right;   right += (sh_prob[in*numstates+i]/nrm); }
    if((left<ksi) && (ksi<=right)){ hstate = i; }
  }
  hopstate = hstate;

  if(hstate==-1){
    std::cout<<"Something is wrong in hop(...) function\nExiting now...\n";
    exit(0);
  }
   

}

void regression(vector<double>& X,vector<double>& Y,int opt,double& a,double& b){
// Linear regression
// opt = 0:   Y =     b*X
// opt = 1:   Y = a + b*X
  if(X.size()!=Y.size()){ cout<<"Error in regression: Size of X array is different from that of Y\n"; exit(0); }

  int sz = X.size();
  double x,y,xy,x2,y2,N;
  // Initialization
  x = 0.0;  y = 0.0;
  x2 = 0.0; y2 = 0.0; xy = 0.0;

//  cout<<"In Regress:\n";  
  for(int i=0;i<sz;i++){
//    cout<<i<<"  "<<X[i]<<"  "<<Y[i]<<endl;
    x += X[i];
    y += Y[i];
    if(opt==1){
      x2 += X[i]*X[i];
      y2 += Y[i]*Y[i];
      xy += X[i]*Y[i];
    }
  }

  if(opt==0){ a = 0.0; b = (y/x); }
  else if (opt==1){ N = sz; b = (N*xy - x*y)/(N*x2 - x*x); a = (y - b*x)/N; }

//  cout<<"In Regress:  y= "<<y<<"  x= "<<x<<" a= "<<a<<"  b= "<<b<<endl;


}

double decoherence_rates(vector<double>& x,double dt,std::string rt_dir,int regress_mode){
/***********************************************
 Computes:
 1) the autocorrelation function of vector x
 Note the size of the autocorr function is 1/2 of 
 the size of vector x
 2) phonon spectrum (FT of the autocorrelation function)
 3) decoherence time
 
 Expected x - fluctuation of the energy difference between two states
***********************************************/
  int len = x.size();
  int sz = (len%2==0)?(len/2):((len-1)/2);
  
  vector<double> C(sz,0.0);  // autocorrelation function
  vector<double> IC(sz,0.0); // first cumulant
  vector<double> IIC(sz,0.0);// second cumulant
  vector<double> D(sz,0.0);  // decoherence function
  vector<double> T,selIIC;  // for regression

  //===== Part 1: Autocorrealtion and decoherence functions ============

  // Normalized autocorrelation functions
  for(int t=0;t<sz;t++){ 
    for(int n=0;n<sz;n++){ C[t] += x[n]*x[n+t];  }
    C[t] /= ((double)sz);
  }//for t

  // Calculate first "cumulants" int_0_t C(t) dt ,for all t
  double sum = 0.0;
  for(t=0;t<sz;t++){ IC[t] = sum;  sum +=  C[t]*(dt/HBAR); }

  // Calculate second "cumulants" int_0_t IC(t) dt ,for all t
  sum = 0.0;
  for(t=0;t<sz;t++){ IIC[t] = sum; sum += IC[t]*(dt/HBAR); }

  // Calculate D(t), see Madrid, et. al
  for(t=0;t<sz;t++){ D[t] = exp(-IIC[t]); }

  // Normalize the autocorrelation function to C[0]
  double nrm = C[0];
  for(t=0;t<sz;t++){ C[t] /= nrm; }

  //===== Part 2: Phonon spectrum (spectral density function) ============
  // Do FT of the normalized autocorrelation function

  // Compute spectral density J
  double dE = 0.0025; // spacing for x (energy) axis for spectral density function = 20 cm^-1
  int Npoints = 400*5; // cover 5 eV range of energies
  vector<double> J(Npoints,0.0);

  for(int w=0;w<Npoints;w++){
    J[w] = 1.0;

    for(int t=1;t<sz;t++){
      double x = (w*dE) * (t*dt);
      J[w] += 2.0*cos(x)*C[t];
    }// for t

    J[w] *= dt;
    J[w] = (J[w]*J[w]/(2.0*M_PI));

  }// for w

  // Output D and its model(based on the fitted parameters)
  ofstream out1((rt_dir+"Spectral_density.txt").c_str(),ios::out);
  for(w=0;w<Npoints;w++){ out1<<"w(eV)= "<<w*dE<<" w(cm^-1)= "<<w*dE*8065.54468111324<<" J= "<<J[w]
                             <<" sqrt(J)= "<<sqrt(J[w])<<endl;
  }
  out1.close();


  //===== Part 3: Decoherence times ============
  // In fact we don't even needed to compute D explicitly
  // Here we take only those pairs of T-D, for which D > eps - some small value
  // => exp(-IIC) > eps => IIC < - ln(eps)
  // If eps = 0.1 => -ln(eps) = 2.3
  //    eps = 0.01 => -ln(eps) = 4.6
  int first = 1;  // this is correction to avoid recurrences!
  for(t=0;t<sz;t++){
    if(first){
      if(IIC[t]<2.3){ 
        T.push_back(t*t*dt*dt); 
        selIIC.push_back(IIC[t]); // sel - selected
      }else{ first = 0; }
    }
    else{ ; ;}
  }
  
  // Do linear regression
  // D(t) = A * exp(-(t/tau)^2), so
  // ln(D(t)) = -IIC(t) = ln(A) -(t/tau)^2, or
  // IIC(t) = -ln(A) + (t/tau)^2
  // linear regression mode is: IIC vs. t^2 with
  // a = -ln(A), b = (1/tau)^2 or sqrt(b) = r_ij - decoherence rate
  double a,b;
  regression(T,selIIC,regress_mode,a,b);
  if(b<0.0){ b = 0.0; }

  // Output D and its model(based on the fitted parameters)
  ofstream out((rt_dir+"Dephasing_function.txt").c_str(),ios::out);
  out<<"Time    D(t)       fitted D(t)     Normalized_autocorrelation_function  Unnormalized_autocorrelation_function   Second cumulant\n";
  for(t=0;t<sz;t++){  out<<t*dt<<"  "<<D[t]<<"  "<<exp(-a) * exp(-b*t*t*dt*dt)<<"  "<<C[t]<<" "<<nrm*C[t]<<"  "<<IIC[t]<<"\n";  }
  out.close();

  return sqrt(b);
}

void Efield(InputStructure& is,double t, CMATRIX& E,double& Eex){
// Field modulation protocol
// is - input parameters
// t - time in fs

  Eex = 0.0;

  if(is.is_field==1){    
    //--------- Direction --------------
    double ix,iy,iz; ix = iy = iz = 0.0;
         if(is.field_dir=="x"){ ix = 1.0; }
    else if(is.field_dir=="y"){ iy = 1.0; }
    else if(is.field_dir=="z"){ iz = 1.0; }
    else if(is.field_dir=="xy"){ ix = iy = 0.5; }
    else if(is.field_dir=="xz"){ ix = iz = 0.5; }
    else if(is.field_dir=="yz"){ iy = iz = 0.5; }
    else if(is.field_dir=="xyz"){ ix = iy = iz = (1.0/3.0); }
    else{ cout<<"Value "<<is.field_dir<<" for the field_dir variable is unknown. Exiting...\n"; exit(0); }

    //--------- Choose modulation protocol ---------
    double Em = 0.0;
    double T = 0.0;  // modulation period
    double Tm = 0.0; // middle of the modulation period
    if(is.field_protocol==1 || is.field_protocol==3){
    /******************************
  Em ^
     |
   1 |        |--------|
     |        |        |
     |        |        |
     |---------------------------->
            Tm-T/2   Tm+T/2       t
    *******************************/
      if(is.is_field_T && is.is_field_Tm){  
        T = is.field_T;
        Tm = is.field_Tm;
        if( ((Tm-0.5*T)<t) && (t<(Tm+0.5*T)) ){  Em = 1.0; }
        else{ Em = 0.0; }
      }else{ // Not defined - then constant for whole time period of simulation
        Em = 1.0;
      }
    }// protocol == 1
    
    else if(is.field_protocol==2){
    /******************************
  Em ^
     |              /\
   1 |             /  \
     |            /    \     
     |           /      \    
     |---------------------------->
              Tm-T/2  Tm+T/2      t
    *******************************/
      if(is.is_field_T && is.is_field_Tm){
        T = is.field_T;
        Tm = is.field_Tm;
        if( ((Tm-0.5*T)<t) && (t<Tm) ){  Em = 2.0*(t-(Tm-0.5*T))/T; }
        else if( (Tm<=t) && (t<=(Tm+0.5*T)) ){ Em = 2.0*((Tm+0.5*T)-t)/T;}
        else{ Em = 0.0; }
      }else{ // Not defined - then constant for whole time period of simulation
        Em = 1.0;
      }

    }// protocol == 2

    //-------- Carrying frequency and amplitude -------
    double omega = 0.0;  // angular frequency [rad/fs]
    double lambda = 1.0; // corresponding wavelength [nm]
    if(is.is_field_freq){
           if(is.field_freq_units=="1/fs"){ omega = 2.0*M_PI*is.field_freq; }     // input is linear frequency
      else if(is.field_freq_units=="rad/fs"){ omega = is.field_freq;  }           // input is angular frequency
      else if(is.field_freq_units=="eV"){ omega = is.field_freq/HBAR; }           // input is energy in eV
      else if(is.field_freq_units=="nm"){ omega = 2.0*M_PI*300.0/is.field_freq; } // input is wavelength in nm
      else{  cout<<"Units of the filed frequency must be specified. Exiting...\n"; exit(0); }
    }// is_field_freq

    lambda = 600.0*M_PI/omega; 

    double Ampl = 0.0; // effective amplitude of the vector potential:
                       // Ampl = e*hbar/m_e * A
    if(is.is_field_fluence){  
      if(is.field_protocol==1){
        Ampl = 0.01038*lambda*sqrt(is.field_fluence/is.field_T);  //result is in eV*Bohr
      }// protocol==1
      else if(is.field_protocol==2){
        Ampl = 0.06857*sqrt(is.field_fluence*lambda/fabs(sin(2.0*omega*is.field_Tm)-cos(2.0*omega*is.field_Tm)));
      }// protocol==2
      else if(is.field_protocol==3){
        Ampl = is.field_fluence;
      }
    }// is_field_fluence


    //============ Now combine all together ==============

    E.M[0] = ix*Em*Ampl*2.0*cos(omega*t);
    E.M[1] = iy*Em*Ampl*2.0*cos(omega*t);
    E.M[2] = iz*Em*Ampl*2.0*cos(omega*t);

    if( ((Tm-0.5*T)<t) && (t<(Tm+0.5*T)) ){  Eex = HBAR*omega; }

  }// is.is_field==1
  else{ E = 0.0; }

}


void propagate_electronic(InputStructure& is,vector<ElectronicStructure>& es,int i, CMATRIX& rates){
  // Propagate coefficients of all adiabatic states
  // Accumulate the surface hopping probabilities for a given nuclear step
  // Yes, we do not make the hops after each electronic step (elec_dt), but
  // The each step contributes to the overall hopping probability (we integrate it
  // over the duration of the nuclear timestep)

  int nel = is.nucl_dt/is.elec_dt; // Number of electronic iterations per 1 nuclear
  int sz = es.size();              // Number of nuclear iterations (ionic steps)
  double tim;                      // time
  double Eex = 0.0;                // bias due to photons
  CMATRIX Ef(3,1);
  

  // Propagate coefficients of all adiabatic states

  // Loop over electronic timesteps withing a single nuclear step
  for(int j=0;j<nel;j++){ 

    // Absolute time
    tim = (i*is.nucl_dt + j*is.elec_dt);

    // Compute field
    Efield(is,tim,Ef,Eex);

    // Propagate coefficients
    if(is.integrator==0){   
      if(is.decoherence==5){   es[i].propagate_coefficients( is.elec_dt,Ef,rates);      } // CPF
      else{                    es[i].propagate_coefficients( is.elec_dt,Ef );       }
    }
    else if(is.integrator==2){
      if(is.decoherence==5){ cout<<"Error: CPF is not available with the integrator=2\nExiting...\n"; exit(0); }
      es[i].propagate_coefficients2(is.elec_dt,Ef);
    }

    // Update time
    es[i].t_m[0] += is.elec_dt; 

    // Update hopping probabilities
    if(is.sh_algo==0){ es[i].update_hop_prob_fssh(is.elec_dt,is.boltz_flag,is.Temp,Ef,Eex,rates);  }
    else if(is.sh_algo==1){  es[i].update_hop_prob_gfsh(is.elec_dt,is.boltz_flag,is.Temp,Ef,Eex,rates);  }
    else if(is.sh_algo==2){  es[i].update_hop_prob_mssh(is.elec_dt,is.boltz_flag,is.Temp,Ef,Eex,rates);  }

  }// for j

}




void run_decoherence_rates(InputStructure& is, vector<ElectronicStructure>& me_es,vector<me_state>& me_states, int icond){
  // The function for computation of the decoherence rates matrix
  cout<<"Entering run_decoherence_rates...\n";

  int sz = me_es.size();              // Number of nuclear iterations (ionic steps)
  int N = me_es[0].num_states;
  CMATRIX rij(N,N);
  ofstream out((is.scratch_dir+"/decoherence_rates_icond"+int2str(icond)+".txt").c_str(),ios::out);

  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      if(i==j){ rij.M[i*N+j]=0.0; }
      else{
        // First lets extract the energy differences of the levels i and j along the trajectory
        vector<double> Eij(sz,0.0);
        double dEij,ave_dEij; ave_dEij = 0.0;
        for(int t=0;t<sz;t++){
          dEij = me_es[t].Hcurr->M[i*N+i].real() - me_es[t].Hcurr->M[j*N+j].real();
          Eij[t] = dEij;
          ave_dEij += dEij;
        }
        ave_dEij /= ((double)sz);
        // Subtract the average value
        for(t=0;t<sz;t++){ Eij[t] -= ave_dEij; }

        // Compute the decoherence rate for pair i,j
        rij.M[i*N+j] = decoherence_rates(Eij,is.nucl_dt,is.scratch_dir+"/icond"+int2str(icond)+"pair"+int2str(i)+"_"+int2str(j),is.regress_mode);
      }
      out<<rij.M[i*N+j].real()<<" ";
    }// for j
    out<<"\n";
  }// for i
  out.close();

}





void run_namd1(InputStructure& is, vector<ElectronicStructure>& me_es,vector<me_state>& me_states, int icond, Random& rnd){
// This version is different from run_namd function in that it does not separate solving TD-SE and computation
// of the surface hopping probabilities. This is because here we inlcude decoherence effects, which effectively
// modify wavefunction (TD-SE solution) along the trajectories stochastically, so it is not possible to separate.
  cout<<"Entering run_namd1 function...\n";

  // Some parameters
  int i,j,n;
  std::string outfile1,outfile2,outfile3,outfile4;
  ofstream out1,out2,out3,out4;
  int nel = is.nucl_dt/is.elec_dt; // Number of electronic iterations per 1 nuclear
  int sz = me_es.size();           // Number of nuclear iterations (ionic steps)
  int nst = me_es[0].num_states;   // Number of electronic states
  int init_state = me_es[0].curr_state;

  // Initialize observables
  int curr_state;  curr_state = me_es[0].curr_state;
  vector<double> tmp(nst,0.0);
  vector<vector<double> > sh_pops(sz,tmp); sh_pops[0][curr_state] = 0.0;
  vector<vector<double> > se_pops(sz,tmp); se_pops[0][curr_state] = 0.0;

  // Decoherence stuff
  vector< vector<double> > r_ij;
  vector< vector<double> > z(nst,std::vector<double>(nst,0.0)); // 2D matrix with all components set to 0.0

  vector<vector<double> > E0(nst,vector<double>(nst,0.0));// average
  vector<vector<double> > d2E_av(nst,vector<double>(nst,0.0)); // average fluctuation of i-j pair
  vector<vector<vector<double> > > d2E(sz,vector<vector<double> >(nst,vector<double>(nst,0.0)));//fluctuation
  CMATRIX rates(nst,nst);

  vector<CMATRIX> sumCij;
  for(int k = 0; k < sz; k++){
      CMATRIX tmp1(nst, nst);
      sumCij.push_back(tmp1);
  }

  if(is.decoherence>0){
    cout<<"Reading decoherence rate matrix for this initial condition...\n";
    std::string filename = is.scratch_dir+"/decoherence_rates_icond"+int2str(icond)+".txt";
    cout<<"Expected filename is: "<<filename<<endl;

    ifstream in;
    in.open(filename.c_str(),ios::in);
    if(in.is_open()){  cout<<"Reading the input from file "<<filename<<endl; }
    else{ cout<<"Error: Can not open file "<<filename<<". Check if this file exists\n"; }
    in.close();

    file2matrix(filename,r_ij);
    rates = CMATRIX(r_ij,z);


    if(is.decoherence==2){ ;; }//is.decoherence == 2 NAC scaling (test method)
    if(is.decoherence==3){ ;; }//is.decoherence == 3 NAC scaling - spectral density variant (test method)
    if(is.decoherence==4){ ;; }//is.decoherence == 4 NAC scaling (same as decoherence==2, but with extrapolation) (test method)
    if(is.decoherence==5){ // Coherence Penalty Functional
     // add nothing special here, but something is done elsewhere
    }    
    if(is.decoherence==6){
    // nothing special:  Reserved for FSSH with wfc collapse - which works fine for ECWR
    }

  }// decoherence > 0

  // Because we will be propagating num_sh_traj independent trajectories at the same time, we will
  // need such an array, initially containing identical entries. Unfortunately this means the memory requirement
  // will increase num_sh_traj_times (hopefully in the future I will optimize such a structure, because the 
  // Hamiltonian(biggest chunk of the memory requirement) is the same everywhere), but for now will do kinda brute force

  /*
  vector<vector<ElectronicStructure> > ME_ES;
  for(int n=0;n<is.num_sh_traj;n++){ 
    vector<ElectronicStructure> ves;
    for(int i=0;i<sz;i++){ ves.push_back(me_es[i]); }
    ME_ES.push_back(ves);
  } // to be sure that we copy objects, not the references
  */
  // Ok. Now, instead of the memory-extensive method we use the following trick:
  // The initial state and coefficients remain unchanged, so we can use the me_es strait away
  // Only we need to make a couple of 2D arrays to store the results.


  //==================== Propagate many-electron orbitals =====================================
  // The outer loop (which calls run_namd1 function) averages over initial conditions

  // Do the hops - averaging over trajectories (stochastic realizations)
  for(n=0;n<is.num_sh_traj;n++){

    me_es[0].set_state(init_state);
    me_es[0].t_m[0] = 0.0; // Time since last hop

    // Loop over time
    for(i=0;i<sz;i++){

      //============ Solve TD-SE and do SH ============
      // Set coefficients and current state index from previous time step to be the current ones
      if(i>0){   me_es[i] << me_es[i-1]; } 

      // Solve TD-SE for i-th time step
      me_es[i].init_hop_prob1();
      propagate_electronic(is,me_es,i,rates);    // update_hop_prob - is called in there 
                                                 // rates are only used if decoherence==5 or decoherence==6
                                                            
      // Calculate the probabilities off all states and hopping probabilities
      me_es[i].update_populations();

      if(is.decoherence==0){  // FSSH
//        curr_state = me_es[i].curr_state;
        hop(me_es[i].g,me_es[i].curr_state,nst);
        curr_state = me_es[i].curr_state;
      }
      else if(is.decoherence==1){  // DISH - currently any value >0
        me_es[i].check_decoherence(is.nucl_dt,is.boltz_flag,is.Temp,rates, rnd);
        curr_state = me_es[i].curr_state;
      }// decoherence == 1

      else if(is.decoherence==2 || is.decoherence==3 || is.decoherence==4){  // NAC scaling
//        curr_state = me_es[i].curr_state;
        hop(me_es[i].g,me_es[i].curr_state,nst);
        curr_state = me_es[i].curr_state;

      }// decoherence == 2
      else if(is.decoherence==5){  // CPF
       // Nothing to do here, because it is MF theory
      }
      else if(is.decoherence==6){  // 
        int st_before = me_es[i].curr_state;

        hop(me_es[i].g,me_es[i].curr_state,nst);

        curr_state = me_es[i].curr_state;
        // Collapse WFC

        if(st_before!=curr_state){ // Hop has happened - collapse wfc

          me_es[i].t_m[0] = 0.0;

          double argg = M_PI*rnd.uniform(-1.0,1.0);
          *me_es[i].Ccurr  = 0.0;
           me_es[i].Ccurr->M[curr_state] = complex<double>( cos(argg), sin(argg) );          
        }

      }// is.decoherence==6

/*  Debug
        cout<<"hop_matrix:\n";
        for(int a=0;a<nst;a++){
          for(int b=0;b<nst;b++){
             cout<<me_es[i].g[a*nst+b]<<"  ";
          }
          cout<<endl;
        }
*/




      // Accumulate SE and SH probabilities for all states
      sh_pops[i][curr_state] += 1.0;
      for(j=0;j<nst;j++){ se_pops[i][j] += me_es[i].A->M[j*nst+j].real(); }
      if(is.td_pop==1){
	      sumCij[i] +=(*(me_es[i].A));
      }
    }// namdtime
  }// for num_sh_traj

  //================ Now output results ======================
  // Output populations as a function of time
  outfile1 = is.scratch_dir+"/me_pop"+int2str(icond);
  out1.open(outfile1.c_str(),ios::out);

  outfile2 = is.scratch_dir+"/out"+int2str(icond);
  out2.open(outfile2.c_str(),ios::out);

  for(i=0;i<sz;i++){
    //---------- SE probabilities ----------
    out1<<"time "<<i<<" "; double tot = 0.0;
    for(j=0;j<nst;j++){
      se_pops[i][j] /= ((double)is.num_sh_traj);
      out1<<"P("<<j<<")= "<<setprecision(10)<<se_pops[i][j]<<"  ";tot += se_pops[i][j];
    } out1<<"Total= "<<tot<<endl;

    //--------- SH probabilities ----------
    out2<<"time "<<i<<" ";
    for(j=0;j<nst;j++){
      sh_pops[i][j] /= ((double)is.num_sh_traj);
      out2<<"P("<<j<<")= "<<setprecision(10)<<sh_pops[i][j]<<" ";
    } out2<<endl;
    
    //------- time-dependent population Cij -----------
    if(is.td_pop==1){
      outfile3 = is.scratch_dir+"/Cij_re_t"+int2str(i)+"_i"+int2str(icond);
      outfile4 = is.scratch_dir+"/Cij_im_t"+int2str(i)+"_i"+int2str(icond);
      out3.open(outfile3.c_str(),ios::out);
      out4.open(outfile4.c_str(),ios::out);
      sumCij[i] = sumCij[i] / is.num_sh_traj;
      for(int i1 = 0; i1 < nst; i1++){
		  for(int j1 = 0; j1 < nst; j1++){
			  out3<<setprecision(10)<<sumCij[i].M[i1*nst+j1].real()<<"   ";
			  out4<<setprecision(10)<<sumCij[i].M[i1*nst+j1].imag()<<"   ";
		  }
		  out3 << endl;
		  out4 << endl;
	  }
     out3.close();
     out4.close();
	}
         
  }

  out1.close();
  out2.close();

  //=======================================================


}

