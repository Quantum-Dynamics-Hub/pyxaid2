/***********************************************************
 * Copyright (C) 2013 Alexey V. Akimov
 * This file is distributed under the terms of the
 * GNU General Public License as published by the
 * Free Software Foundation; either version 3 of the
 * License, or (at your option) any later version.
 * http://www.gnu.org/copyleft/gpl.txt
***********************************************************/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <sstream>
#include <time.h>
#include <stdlib.h>
#include <time.h>
#include "ElectronicStructure.h"
#include "aux.h"
#include "io.h"
#include "namd.h"
#include <boost/python.hpp>
using namespace boost::python;
using namespace std;


int namd(boost::python::dict inp_params){

  time_t t1 = clock();
  complex<double> ihbar(0.0,hbar);

//>>>>>>>>>>>>>>>>>>>>>>>> INITIALIZATION PART <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  // General input parameters from the python dictionary
  InputStructure params(inp_params);

  // Define ME basis
  vector<me_state> me_states;
  input_states(inp_params,me_states); // Convention will be to count from 1, where 1 - is the first orbital from input
                                      // (real and energy files), not necessarily the actual first orbital of the system

  int me_numstates = me_states.size();               // This is a number of true (multi-electron) states to propagate
  int numstates = me_states[0].active_space.size();  // This is a number of 1-electron orbitals to propagate
  int num_elec = me_states[0].actual_state.size();   // This is a number of valent (capable to hop) electrons
  cout<<"Number of 1-electron states(orbitals) = "<<numstates<<endl;
  cout<<"Number of multi-electron states(determinants) = "<<me_numstates<<endl;
  cout<<"Number of electrons in active space = "<<num_elec<<endl;
  
  // Initialize random number generator
  srand(time(0));
  cout<<"RAND_MAX = "<<RAND_MAX<<endl;

  // icond loop - from the input dictionary
  vector< vector<int> > iconds; // iconds[j][0] = init_time[j], iconds[j][1] = init_state[j]
  input_iconds(inp_params,me_numstates,iconds);

  // Set the unit conversion factor 
  double en_scl = 1.0;
  if(params.energy_units=="Ry"){ en_scl = Ry_to_eV; }
  else if(params.energy_units=="Ha"){ en_scl = Ha_to_eV; }

  complex<double> hp_scl(1.0,0.0); // = -ihbar; This is because I initially overlooked i (complex unity) in p operator
                                   // Thus, the field-matter interaction Hamiltonian is a real matrix
                                   // Also: hbar is included in e*hbar/m_e prefactor, coming in Efield function, so for 
                                   // now we don't need to scale values in Hprime_ files

  //-------------------------------------------------------------------
  //------------------- Batch mode preparations -----------------------
  // Find the maximal index of the input files (Hamiltonian), needed for the batch run
  int max_indx = 0;
  for(int icond=0;icond<iconds.size();icond++){
    if(iconds[icond][0]>=max_indx){ max_indx = iconds[icond][0]; }
  }// for icond
  max_indx += params.namdtime;
  cout<<"Maximal Hamiltonian file to read is "<<params.Ham_re_prefix<<(max_indx+1)<<params.Ham_re_suffix<<endl;

  // Read all necessary couplings and transition dipole (if necessary) files - batch mode
  vector< matrix > H_batch;
  vector< matrix > Hprime_x_batch;
  vector< matrix > Hprime_y_batch;
  vector< matrix > Hprime_z_batch;

  if(params.read_couplings=="batch" || params.read_couplings=="batch_all_in_one"){

  // If params.read_couplings=="batch_all_in_one", then the filenames for real and imaginary parts can be the same
  // file - the Hamiltonian will be composed from the diagonal (to become real part) and off-diagonal (to become
  // imaginary part) elements.
  // H_ij = F_ii    if i==j (real, diagonal)
  // H_ij = -i*F_ij   if i!=j (imaginary, off-diagonal, "-" is the convention implying that F_ij = hbar * <i|d/dt|j> )

    for(int j=0;j<max_indx;j++){
      // -------------------- Real part of the Hamiltonian matrix -------------------------------------
      vector< vector<double> > Ham_re, Ham_re_crop;
      std::string Ham_re_file; Ham_re_file = params.Ham_re_prefix + int2string(j) + params.Ham_re_suffix;

      if(params.debug_flag==1){ cout<<"Reading Hamiltonian file(real part) = "<<Ham_re_file<<endl; }
      file2matrix(Ham_re_file,Ham_re);
//      cout<<"Ham_re = :"<<endl; show_2D(Ham_re);
      extract_2D(Ham_re,Ham_re_crop,me_states[0].active_space,-1); // also crop the matrix to the active space
//      cout<<"Ham_re(cropped):"<<endl; show_2D(Ham_re_crop);

      // -------------------- Imaginary part of the Hamiltonian matrix -------------------------------------
      vector< vector<double> > Ham_im, Ham_im_crop;
      std::string Ham_im_file; Ham_im_file = params.Ham_im_prefix + int2string(j) + params.Ham_im_suffix;

      if(params.debug_flag==1){ cout<<"Reading Hamiltonian file(imaginary part) = "<<Ham_im_file<<endl; }
      file2matrix(Ham_im_file,Ham_im);
      extract_2D(Ham_im,Ham_im_crop,me_states[0].active_space,-1); // also crop the matrix to the active space

      //--------------------- Optionally preprocess matrices ----------------------------------------------
      if(params.read_couplings=="batch_all_in_one"){
        for(int a=0;a<Ham_re_crop.size();a++){
          for(int b=0;b<Ham_re_crop.size();b++){
            if(a==b){
              Ham_im_crop[a][a] = 0.0;
            }
            else{ 
              Ham_re_crop[a][b] = 0.0;
              Ham_im_crop[a][b] *= -1.0; 
            }
            
          }// for b
        }// for a
      }// if all_in_one

      // -------------------- Create Hamiltonian matrix ---------------------------------------------------
      matrix Ham(Ham_re_crop,Ham_im_crop);
//      cout<<"Composed Ham = "<<Ham<<endl;

      //--------------------- Scale Ham to units of [Energy] = eV, [time] = fs ----------------------------
      Ham *= en_scl;
      if(params.debug_flag==2){   cout<<"Scaled Ham = "<<Ham<<endl; }

      H_batch.push_back(Ham);


      // -------------------- Real part of the transition dipole matrix -------------------------------------
      if(params.is_field==1){
        vector< vector<double> > Hprime_x, Hprime_x_crop;
        vector< vector<double> > Hprime_y, Hprime_y_crop;
        vector< vector<double> > Hprime_z, Hprime_z_crop;

        std::string Hprime_x_file; Hprime_x_file = params.Hprime_x_prefix + int2string(j) + params.Hprime_x_suffix;
        std::string Hprime_y_file; Hprime_y_file = params.Hprime_y_prefix + int2string(j) + params.Hprime_y_suffix;
        std::string Hprime_z_file; Hprime_z_file = params.Hprime_z_prefix + int2string(j) + params.Hprime_z_suffix;

        if(params.debug_flag==1){ cout<<"Reading transition dipole file in momentum representation (x component) = "
        <<Hprime_x_file<<endl; }
        file2matrix(Hprime_x_file,Hprime_x);
        extract_2D(Hprime_x,Hprime_x_crop,me_states[0].active_space,-1);

        if(params.debug_flag==1){ cout<<"Reading transition dipole file in momentum representation (y component) = "
        <<Hprime_y_file<<endl; }
        file2matrix(Hprime_y_file,Hprime_y);
        extract_2D(Hprime_y,Hprime_y_crop,me_states[0].active_space,-1);

        if(params.debug_flag==1){ cout<<"Reading transition dipole file in momentum representation (z component) = "
        <<Hprime_z_file<<endl; }
        file2matrix(Hprime_z_file,Hprime_z);
        extract_2D(Hprime_z,Hprime_z_crop,me_states[0].active_space,-1);

        vector< vector<double> > tmp(Hprime_x_crop.size(),vector<double>(Hprime_x_crop.size(),0.0));
        //------------ Create matrix --------
        // Compose matrices with complex elements from the matrices with real elements
        // In our case Hprime_ must be purely imaginary, so:
        matrix Hprimex(tmp,Hprime_x_crop);
        matrix Hprimey(tmp,Hprime_y_crop);
        matrix Hprimez(tmp,Hprime_z_crop);

        // Now scale Hprime_ elements by common scaling factor
        Hprimex *= hp_scl; Hprimey *= hp_scl; Hprimez *= hp_scl;
        if(params.debug_flag==2){   
          cout<<"Scaled Hprimex = "<<Hprimex<<endl;
          cout<<"Scaled Hprimey = "<<Hprimey<<endl;
          cout<<"Scaled Hprimez = "<<Hprimez<<endl;
        }

        Hprime_x_batch.push_back(Hprimex);
        Hprime_y_batch.push_back(Hprimey);
        Hprime_z_batch.push_back(Hprimez);


      }// if one wants explicit field effects


    }// for j
    cout<<"end of Hamiltonian files reading\n";
  }// if batch mode and namd

//------------------------------------------------------------------------------

//>>>>>>>>>>>>>>>>>>>>>>>>> MAIN PROGRAM PART <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  cout<<"Starting the program...\n";
  for(icond=0;icond<iconds.size();icond++){  // first_icond may start from 0, not 1

    if(params.debug_flag==2){
      cout<<"Initial condition index = "<<icond<<"     initial_time["<<icond<<"]="<<iconds[icond][0]
                                               <<" initial_me_state["<<icond<<"]="<<iconds[icond][1]<<endl;
    }// debug

    //>>>>> Collect parameters for NA-MD with given initial condition and given trajectory length
    // Electronic structure prototype
    vector<ElectronicStructure> oe_es(params.namdtime,ElectronicStructure(2*numstates));  // one-electron orbitals
    vector<ElectronicStructure> me_es(params.namdtime,ElectronicStructure(me_numstates)); // multi-electron orbitals

 
    for(int j=iconds[icond][0];j<iconds[icond][0]+params.namdtime;j++){ // namdtime also starts from 0
      if(params.debug_flag==1){    cout<<"----------- j = "<<j<<" -------------------"<<endl;}
      int t = (j - iconds[icond][0]);  // Time

      matrix* T;
      matrix *Tx, *Ty, *Tz;

      if(params.read_couplings=="online" || params.read_couplings=="online_all_in_one"){

        // -------------------- Real part of the Hamiltonian matrix -------------------------------------
        vector< vector<double> > Ham_re, Ham_re_crop;
        std::string Ham_re_file; Ham_re_file = params.Ham_re_prefix + int2string(j) + params.Ham_re_suffix;

        if(params.debug_flag==1){ cout<<"Reading Hamiltonian file(real part) = "<<Ham_re_file<<endl; }
        file2matrix(Ham_re_file,Ham_re);
        extract_2D(Ham_re,Ham_re_crop,me_states[0].active_space,-1); // also crop the matrix to the active space

        // -------------------- Imaginary part of the Hamiltonian matrix -------------------------------------
        vector< vector<double> > Ham_im, Ham_im_crop;
        std::string Ham_im_file; Ham_im_file = params.Ham_im_prefix + int2string(j) + params.Ham_im_suffix;

        if(params.debug_flag==1){ cout<<"Reading Hamiltonian file(imaginary part) = "<<Ham_im_file<<endl; }
        file2matrix(Ham_im_file,Ham_im);
        extract_2D(Ham_im,Ham_im_crop,me_states[0].active_space,-1); // also crop the matrix to the active space

        //--------------------- Optionally preprocess matrices ----------------------------------------------
        if(params.read_couplings=="online_all_in_one"){
          for(int a=0;a<Ham_re_crop.size();a++){
            for(int b=0;b<Ham_re_crop.size();b++){
              if(a==b){
                Ham_im_crop[a][a] = 0.0;
              }
              else{
                Ham_re_crop[a][b] = 0.0;
                Ham_im_crop[a][b] *= -1.0;
              }

            }// for b
          }// for a
        }// if all_in_one


        // -------------------- Create Hamiltonian matrix ---------------------------------------------------
        matrix Ham(Ham_re_crop,Ham_im_crop);

        //--------------------- Scale Ham to units of [Energy] = eV, [time] = fs ----------------------------
        Ham *= en_scl;
        if(params.debug_flag==2){   cout<<"Scaled Ham = "<<Ham<<endl; }

        T = new matrix(Ham.n_rows,Ham.n_cols);
        *T = Ham;


        // -------------------- Real part of the transition dipole matrix -------------------------------------
        if(params.is_field==1){
          vector< vector<double> > Hprime_x, Hprime_x_crop;
          vector< vector<double> > Hprime_y, Hprime_y_crop;
          vector< vector<double> > Hprime_z, Hprime_z_crop;

          std::string Hprime_x_file; Hprime_x_file = params.Hprime_x_prefix + int2string(j) + params.Hprime_x_suffix;
          std::string Hprime_y_file; Hprime_y_file = params.Hprime_y_prefix + int2string(j) + params.Hprime_y_suffix;
          std::string Hprime_z_file; Hprime_z_file = params.Hprime_z_prefix + int2string(j) + params.Hprime_z_suffix;

          if(params.debug_flag==1){ cout<<"Reading transition dipole file(x component) = "<<Hprime_x_file<<endl; }
          file2matrix(Hprime_x_file,Hprime_x);
          extract_2D(Hprime_x,Hprime_x_crop,me_states[0].active_space,-1);

          if(params.debug_flag==1){ cout<<"Reading transition dipole file(y component) = "<<Hprime_y_file<<endl; }
          file2matrix(Hprime_y_file,Hprime_y);
          extract_2D(Hprime_y,Hprime_y_crop,me_states[0].active_space,-1);

          if(params.debug_flag==1){ cout<<"Reading transition dipole file(z component) = "<<Hprime_z_file<<endl; }
          file2matrix(Hprime_z_file,Hprime_z);
          extract_2D(Hprime_z,Hprime_z_crop,me_states[0].active_space,-1);

          vector< vector<double> > tmp(Hprime_x_crop.size(),vector<double>(Hprime_x_crop.size(),0.0));
          //------------ Create matrix --------
          matrix Hprimex(Hprime_x_crop,tmp);
          matrix Hprimey(Hprime_y_crop,tmp);
          matrix Hprimez(Hprime_z_crop,tmp);

          Hprimex *= hp_scl; Hprimey *= hp_scl; Hprimez *= hp_scl;
          if(params.debug_flag==2){
            cout<<"Scaled Hprimex = "<<Hprimex<<endl;
            cout<<"Scaled Hprimey = "<<Hprimey<<endl;
            cout<<"Scaled Hprimez = "<<Hprimez<<endl;
          }

          Tx = new matrix(Hprimex.n_rows,Hprimex.n_cols);  *Tx = Hprimex;
          Ty = new matrix(Hprimey.n_rows,Hprimey.n_cols);  *Ty = Hprimey;
          Tz = new matrix(Hprimez.n_rows,Hprimez.n_cols);  *Tz = Hprimez;

        }// if params.is_field==1

      }// "online"
      else if(params.read_couplings=="batch"){ 
        T = new matrix(H_batch[j].n_rows,H_batch[j].n_cols);   *T = H_batch[j];

        if(params.is_field==1){
          Tx = new matrix(Hprime_x_batch[j].n_rows,Hprime_x_batch[j].n_cols); *Tx = Hprime_x_batch[j];
          Ty = new matrix(Hprime_y_batch[j].n_rows,Hprime_y_batch[j].n_cols); *Ty = Hprime_y_batch[j];
          Tz = new matrix(Hprime_z_batch[j].n_rows,Hprime_z_batch[j].n_cols); *Tz = Hprime_z_batch[j];
        }// if params.is_field==1
      }// "batch"

      matrix Hij(T->n_rows,T->n_cols);    Hij = *T; 
      matrix Hij_prime_x(T->n_rows,T->n_cols); Hij_prime_x = 0.0; 
      matrix Hij_prime_y(T->n_rows,T->n_cols); Hij_prime_y = 0.0;
      matrix Hij_prime_z(T->n_rows,T->n_cols); Hij_prime_z = 0.0;
      delete T;

      if(params.is_field==1){   Hij_prime_x = *Tx; Hij_prime_y = *Ty; Hij_prime_z = *Tz;  delete Tx; delete Ty; delete Tz; }


      if(t==0 && params.debug_flag==1){
        
        cout<<"Hij_prime_x  = "<<Hij_prime_x<<endl;
        cout<<"Hij_prime_y  = "<<Hij_prime_y<<endl;
        cout<<"Hij_prime_z  = "<<Hij_prime_z<<endl;

      }

      //Set up properties of the ElectronicStructure objects:    
      //------------------ Common data ----------------------------

        *oe_es[t].Hcurr = 0.0;
        *oe_es[t].Hprimex = 0.0;
        *oe_es[t].Hprimey = 0.0;
        *oe_es[t].Hprimez = 0.0;

        for(int k1=0;k1<numstates;k1++){
          for(int k2=0;k2<numstates;k2++){
            /**************************************
              Setting block matrix:
  
                        i_alp            i_bet
                   __________________________________
            j_alp  | .d[2*k1][2*k2]   .d[2*k1][2*k2+1]
                   |
            j_bet  |.d[2*k1+1][2*k2]  .d[2*k1+1][2*k2+1]

             **************************************/
           //Couplings: dij = <i|d/dt|j> 
           oe_es[t].Hcurr->M[2*k1*(2*numstates)+2*k2] = Hij.M[k1*numstates+k2];
           oe_es[t].Hcurr->M[(2*k1+1)*(2*numstates)+2*k2+1] = Hij.M[k1*numstates+k2];

           //Perturbations
           oe_es[t].Hprimex->M[2*k1*(2*numstates)+2*k2] = Hij_prime_x.M[k1*numstates+k2];
           oe_es[t].Hprimex->M[(2*k1+1)*(2*numstates)+2*k2+1] = Hij_prime_x.M[k1*numstates+k2];

           oe_es[t].Hprimey->M[2*k1*(2*numstates)+2*k2] = Hij_prime_y.M[k1*numstates+k2];
           oe_es[t].Hprimey->M[(2*k1+1)*(2*numstates)+2*k2+1] = Hij_prime_y.M[k1*numstates+k2];

           oe_es[t].Hprimez->M[2*k1*(2*numstates)+2*k2] = Hij_prime_z.M[k1*numstates+k2];
           oe_es[t].Hprimez->M[(2*k1+1)*(2*numstates)+2*k2+1] = Hij_prime_z.M[k1*numstates+k2];


            if(params.alp_bet==0){ // Electrons with a spin, no coupling between alp and bet, default

              oe_es[t].Hcurr->M[(2*k1+1)*(2*numstates)+2*k2] = 
              oe_es[t].Hcurr->M[2*k1*(2*numstates)+2*k2+1] = 0.0;

              oe_es[t].Hprimex->M[(2*k1+1)*(2*numstates)+2*k2] =
              oe_es[t].Hprimex->M[2*k1*(2*numstates)+2*k2+1] = 0.0;

              oe_es[t].Hprimey->M[(2*k1+1)*(2*numstates)+2*k2] =
              oe_es[t].Hprimey->M[2*k1*(2*numstates)+2*k2+1] = 0.0;

              oe_es[t].Hprimez->M[(2*k1+1)*(2*numstates)+2*k2] =
              oe_es[t].Hprimez->M[2*k1*(2*numstates)+2*k2+1] = 0.0;


            }
            else if(params.alp_bet==1){ // Spinless electrons, coupling between alp and bet is !=0, based only
                                        // on spatial part of the wavefunctions
              oe_es[t].Hcurr->M[(2*k1+1)*(2*numstates)+2*k2] = 
              oe_es[t].Hcurr->M[2*k1*(2*numstates)+2*k2+1] = Hij.M[k1*numstates+k2];

              oe_es[t].Hprimex->M[(2*k1+1)*(2*numstates)+2*k2] =
              oe_es[t].Hprimex->M[2*k1*(2*numstates)+2*k2+1] = Hij_prime_x.M[k1*numstates+k2];

              oe_es[t].Hprimey->M[(2*k1+1)*(2*numstates)+2*k2] =
              oe_es[t].Hprimey->M[2*k1*(2*numstates)+2*k2+1] = Hij_prime_y.M[k1*numstates+k2];

              oe_es[t].Hprimez->M[(2*k1+1)*(2*numstates)+2*k2] =
              oe_es[t].Hprimez->M[2*k1*(2*numstates)+2*k2+1] = Hij_prime_z.M[k1*numstates+k2];

            }
          }// for k2
        }// for k1
//      }// if namd

//      cout<<"*(oe_es[t].Hcurr) = "<<*(oe_es[t].Hcurr)<<endl;
    }// namdtime loop - duration of run  - finishes at time init_time[icond]+namdtime
      
    cout<<"One-electron data are set\n";


    //------------------ Not common data ---------------------------
    me_es[0].set_state(iconds[icond][1]); // Coefficients and populations



    //======================== Compute multi-electron Hamiltonians ==========================
    //-------------------- Couplings -----------------------------------------
    // Now we consider all multi-electron states
    for(j=iconds[icond][0];j<iconds[icond][0]+params.namdtime;j++){
      int t = (j - iconds[icond][0]);  // Time
      int I,J;

      *me_es[t].Hcurr = 0.0;
      *me_es[t].Hprimex = 0.0;
      *me_es[t].Hprimey = 0.0;
      *me_es[t].Hprimez = 0.0;

      // Initialize Energies, NACs and Hprime
      for(I=0;I<me_es[t].num_states;I++){
        // This initialization already includes shift of 1-e orbitals and 2-particle corrections
        me_es[t].Hcurr->M[I*me_es[t].num_states+I] = me_states[I].Exc + me_states[I].Eshift;

        for(J=0;J<me_es[t].num_states;J++){
          // off-diagonal elements - NAC
          if(I!=J){ me_es[t].Hcurr->M[I*me_es[t].num_states+J] = 0.0;    }

          // Perturbations
          me_es[t].Hprimex->M[I*me_es[t].num_states+J] = 0.0;
          me_es[t].Hprimey->M[I*me_es[t].num_states+J] = 0.0;
          me_es[t].Hprimez->M[I*me_es[t].num_states+J] = 0.0;

        }// for J
      }// for I



      // Compute many-electron properties from those of the 1-electron
      for(I=0;I<me_es[t].num_states;I++){          // Numerate the me state on which we project.
                                                   // This multi-electron state J is defined by me_states[J]
        for(J=0;J<me_es[t].num_states;J++){

          // Practically there will be only one non-zero contribution, corresponding to
          // the pair of different indexes, for which all other indexes are identical
          // If there are 2 pairs of different indexes (one determinant is doubly-excited
          // with respect to the other) - the corresponding contribution is zero

          int orb_i,orb_j;
          int delt = delta(me_states[I].actual_state,me_states[J].actual_state,orb_i,orb_j);
          if(delt){
            // If we are here - then two configurations differ by not 1 occupied orbital
            if(t==0  &&  params.debug_flag==1){ cout<<"orb_i, orb_j = "<<orb_i<<"  "<<orb_j<<endl;  }
            // Convert external orbital indexes to internal orbital indexes
            orb_i = ext2int(orb_i,me_states[I].active_space); 
            orb_j = ext2int(orb_j,me_states[J].active_space);

            // In the statements below the += operator should be encountered only once over I, J double loop
            // so initialization inside second loop is ok. Also += is effectively = operator.
            // NAC and energy
            me_es[t].Hcurr->M[I*me_es[t].num_states+J] += oe_es[t].Hcurr->M[orb_i*oe_es[t].num_states+orb_j];

            // Perturbations - transition dipole moments
            me_es[t].Hprimex->M[I*me_es[t].num_states+J] += oe_es[t].Hprimex->M[orb_i*oe_es[t].num_states+orb_j];
            me_es[t].Hprimey->M[I*me_es[t].num_states+J] += oe_es[t].Hprimey->M[orb_i*oe_es[t].num_states+orb_j];
            me_es[t].Hprimez->M[I*me_es[t].num_states+J] += oe_es[t].Hprimez->M[orb_i*oe_es[t].num_states+orb_j];


          }// if delt
/*
          if(t==0 &&  params.debug_flag==1){
            // Only for the first time step output info - to check what is the NAC structure of the system
            cout<<"I, J, delt, coupling(not scaled), orb_i, orb_j, Hprimex, Hprimey, Hprimez = "
                <<I<<"  "<<J<<"  "<<delt<<"  "
                <<me_es[t].Hcurr->M[I*me_es[t].num_states+J]<<"  "<<orb_i<<"  "<<orb_j<<"  "
                <<me_es[t].Hprimex->M[I*me_es[t].num_states+J]<<"  "
                <<me_es[t].Hprimey->M[I*me_es[t].num_states+J]<<"  "
                <<me_es[t].Hprimez->M[I*me_es[t].num_states+J]<<"  "

                <<endl;
          }
*/
      }// for J

      // Now scale the coupling!!!
      int sz_scl = me_states[I].nac_scl.size();
      for(int k=0;k<sz_scl;k++){
        J = me_states[I].nac_scl_indx[k];
        me_es[t].Hcurr->M[I*me_es[t].num_states+J] *= me_states[I].nac_scl[k];
      }// for k

      // Compute the energy and the perturbation of the macrostate
      for(int el=0;el<num_elec;el++){
        int orb_i = me_states[I].actual_state[el];        // orbital on which el-th electron sits in current function
        orb_i = ext2int(orb_i,me_states[I].active_space); // internal index of the orbital
        // Energy of I-th basis function (determinant) - contributions of all 1-electron KS orbitals - diagonal terms
        me_es[t].Hcurr->M[I*me_es[t].num_states+I] += oe_es[t].Hcurr->M[orb_i*oe_es[t].num_states+orb_i]; 

        if(params.debug_flag>=1 && t==0){
        cout<<"I= "<<I<<" el= "<<el<<" orb_i= "<<orb_i<<" E_{KS,orb_i}= "
            <<oe_es[t].Hcurr->M[orb_i*oe_es[t].num_states+orb_i]<<" E_{state,I}= "
            << me_es[t].Hcurr->M[I*me_es[t].num_states+I]<<endl;
        }

        me_es[t].Hprimex->M[I*me_es[t].num_states+I] += oe_es[t].Hprimex->M[orb_i*oe_es[t].num_states+orb_i];
        me_es[t].Hprimey->M[I*me_es[t].num_states+I] += oe_es[t].Hprimey->M[orb_i*oe_es[t].num_states+orb_i];
        me_es[t].Hprimez->M[I*me_es[t].num_states+I] += oe_es[t].Hprimez->M[orb_i*oe_es[t].num_states+orb_i];


      }// for el

    }// for I

    for(I=0;I<me_es[t].num_states;I++){          // Numerate the me state on which we project.
                                                 // This multi-electron state J is defined by me_states[J]
      for(J=0;J<me_es[t].num_states;J++){

          if(t==0 &&  params.debug_flag==1){
            // Only for the first time step output info - to check what is the NAC structure of the system
            cout<<"I, J, coupling(scaled), Hprimex, Hprimey, Hprimez = "
                <<I<<"  "<<J<<"  "
                <<me_es[t].Hcurr->M[I*me_es[t].num_states+J]<<"  "
                <<me_es[t].Hprimex->M[I*me_es[t].num_states+J]<<"  "
                <<me_es[t].Hprimey->M[I*me_es[t].num_states+J]<<"  "
                <<me_es[t].Hprimez->M[I*me_es[t].num_states+J]<<"  "
                <<endl;
          }

      }// for J
    }// for I

    //------------------------------------------------------------------------
//      cout<<"me.Hcurr = "<<*(me_es[t].Hcurr)<<endl;

    }// for t = j - ...
    cout<<"Multi-electron couplings and energies are computed\n";



    // Print the energies of multi-electron states
    string outfile = (params.scratch_dir + "/me_energies"+int2string(icond));
    cout<<"The energies of basis  states (with respect to defined ground state) for the this initial condition are written in file "<<outfile<<"\n";
    ofstream out; out.open(outfile.c_str(),ios::out);
    for(j=iconds[icond][0];j<iconds[icond][0]+params.namdtime;j++){
      int t = (j - iconds[icond][0]);  // Time
      out<<"t= "<<j<<"  "<<"E[0]= "<<me_es[t].Hcurr->M[0].real()<<"  ";
      for(int I=0;I<me_es[t].num_states;I++){
        out<<"E["<<I<<"]-E[0]= "<<(me_es[t].Hcurr->M[I*me_es[t].num_states+I]-me_es[t].Hcurr->M[0]).real()<<"  ";
      }// for I
      out<<endl;
    }// for j
    out.close();

    //>>>>> Precompute decoherence rates
    if(params.decoherence>0){
        cout<<"Starting decoherence rates calculation\n";
        run_decoherence_rates(params,me_es,me_states,icond);
    }
    //>>>>> Run NA-MD
//    if(params.runtype=="namd" && params.decoherence==0){
//        cout<<"Starting na-md simulations\n";
//        run_namd(params,me_es,me_states,icond); 
//    }
//    if(params.runtype=="namd" && params.decoherence>0){
        cout<<"Starting na-md simulations with (optional) decoherence\n";
        run_namd1(params,me_es,me_states,icond);
//    }

    oe_es.clear();
    me_es.clear();

  }// icond loop - from which time to start

  
  time_t t2 = clock();
  cout<<"Time in namd is: "<<(t2-t1)/((double)CLOCKS_PER_SEC)<<endl;

  return 0;

}


void export_namd(){
  def("namd",&namd);

}
