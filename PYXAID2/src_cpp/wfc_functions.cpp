/***********************************************************
 * Copyright (C) 2013 Alexey V. Akimov
 * This file is distributed under the terms of the
 * GNU General Public License as published by the
 * Free Software Foundation; either version 3 of the
 * License, or (at your option) any later version.
 * http://www.gnu.org/copyleft/gpl.txt
***********************************************************/

#include "wfc.h"
#include "matrix.h"
#include "units.h"

// Here we define a set of functions which are not the members of wfc, but
// rather take wfc objects as arguments

void overlap(wfc& wfc1,int k1,int minband,int maxband,std::string filename){
// This function computes the overlap of the wavefunctions at given k-point
// and writes the overlap matrix in files filename_re and filename_im
// orbitals are in range [minband,maxband]
  if(minband<0){ cout<<"Error in nac() : minimal band index is 0, given "<<minband<<endl; exit(0); }
  if(maxband>=wfc1.nbands){ cout<<"Error in nac() : maximal band index is "<<wfc1.nbands-1<<", given "<<maxband<<endl; exit(0); }

  ofstream f_re((filename+"_re").c_str(),ios::out);
  ofstream f_im((filename+"_im").c_str(),ios::out);
  int sz = wfc1.nbands;
  
  for(int i=minband;i<=maxband;i++){
    for(int j=minband;j<=maxband;j++){
      complex<double> res = wfc1.kpts[k1].mo[i].conj() * wfc1.kpts[k1].mo[j];   
      f_re<<res.real()<<"  ";
      f_im<<res.imag()<<"  ";
    }// for j
    f_re<<endl;
    f_im<<endl;
  }// for i

  f_re.close();
  f_im.close();

}

void energy(wfc& wfc1,int k,int minband,int maxband,std::string filename){
// This function prints the orbital energies (eigenvalues) in a range of band
// indexes [minmband, maxband] for given wavefunctions
// for given k-point into file filename
  if(minband<0){ cout<<"Error in energy() : minimal band index is 0, given "<<minband<<endl; exit(0); }
  if(maxband>=wfc1.nbands){ cout<<"Error in energy() : maximal band index is "<<wfc1.nbands-1<<", given "<<maxband<<endl; exit(0); }

  ofstream f(filename.c_str(),ios::out);
  for(int i=minband;i<=maxband;i++){ f<<wfc1.kpts[k].mo[i].energy<<"  "; }
  f<<endl;
  f.close();
}


void nac(wfc& wfc1,wfc& wfc2,int k1,int k2,int minband,int maxband,double dt,std::string filename){
// This function computes the non-adiabatic coupling (nac) of the 2 wavefunctions at k-points
// given by indices k1 and k2 (starting from 0)
// and writes the nac matrix in files filename_re and filename_im
// The range of the orbitals for which NACs are computed is: [minband,maxband]
// Here the meaning of wfc1 and wfc2 - is the wavefunction at different time steps
// Note that all wavefunction pereparation steps (normalize,complete,restore) are supposed to be already done

  if(minband<0){ cout<<"Error in nac() : minimal band index is 0, given "<<minband<<endl; exit(0); }
  if(maxband>=wfc1.nbands){ cout<<"Error in nac() : maximal band index is "<<wfc1.nbands-1<<", given "<<maxband<<endl; exit(0); }

  ofstream f_re((filename+"_re").c_str(),ios::out);
  ofstream f_im((filename+"_im").c_str(),ios::out);

  if(wfc1.nbands!=wfc2.nbands){ cout<<"Error in overlap: Wavefunctions have different number of bands\n"; exit(0); }
  for(int i=minband;i<=maxband;i++){
    for(int j=minband;j<=maxband;j++){
      complex<double> res1 = wfc1.kpts[k1].mo[i].conj() * wfc2.kpts[k2].mo[j];
      complex<double> res2 = wfc2.kpts[k2].mo[i].conj() * wfc1.kpts[k1].mo[j];
      complex<double> res = (0.5/dt)*(res1 - res2);
 
      f_re<<res.real()<<"  ";
      f_im<<res.imag()<<"  ";
    }// for j
    f_re<<endl;
    f_im<<endl;
  }// for i

  f_re.close();
  f_im.close();

}

void ham(wfc& wfc1,wfc& wfc2,int k1,int k2,int minband,int maxband,double dt,std::string filename){
// This function computes the Hamiltonian (h) of the 2 wavefunctions at k-points
// given by indices k1 and k2 (starting from 0)
// and writes the H matrix in files filename_re and filename_im - real and imaginary parts 
// The range of the orbitals for which H are computed is: [minband,maxband]
// Here the meaning of wfc1 and wfc2 - is the wavefunction at different time steps
// Note that all wavefunction pereparation steps (normalize,complete,restore) are supposed to be already done
// The non-diagonal elements of the H are the non-adiabatic couplings (with -i*hbar factor), while diagonal elements
// are the orbital energies (pure real)
/******************************************************************
 Below we use formula:
 d_ij(t+dt/2) = (0.5/dt)*(<i(t)|j(t+dt)> - <i(t)|j(t)> + <i(t+dt)|j(t+dt)> - <i(t+dt)|j(t)>) =
 = (0.5/dt)*(<i(t)|j(t+dt)> - <i(t+dt)|j(t)> + <i(t+dt)|j(t+dt)> - <i(t)|j(t)>)
                res1               res2               res3             res4

 In principle if the orbitals are orthonormal the last 2 terms should vanish, but
 practically for plane-wave codes they are not, because the orthogonality is
 given with respect to projector matrix S:  <i|S|j> = delta_ij, but <i|j> !=delta_ij

 If the "restore" (or normalization) operation has been performed on wavefunctions
 the last 2 terms vanish identically.

 Closer look showns that ( res1 - res2 ) is antihermitian, while
 ( res3 - res4 ) is Hermitian, so we do not include the latter in
 calculations (even in case they do not vanish)

 So basically here what is computed:

 H_ij(t+dt/2) = Eii(t+dt/2)*delta_ij - i*hbar*d_ij(t+dt/2)
******************************************************************/

  complex<double> ihbar(0.0,-hbar); // actually minus i*hbar
  if(wfc1.energy_units=="Ry"||wfc1.energy_units=="Rydberg"){ 
    cout<<"Hamiltonian is in Rydberg units\n";
    ihbar /= Ry_to_eV;
  } // all in Ry
  else{  cout<<"Hamiltonian is in eV units\n";  }

  if(minband<0){ cout<<"Error in ham() : minimal band index is 0, given "<<minband<<endl; exit(0); }
  if(maxband>=wfc1.nbands){ cout<<"Error in ham() : maximal band index is "<<wfc1.nbands-1<<", given "<<maxband<<endl; exit(0); }

  ofstream f_re((filename+"_re").c_str(),ios::out);
  ofstream f_im((filename+"_im").c_str(),ios::out);

  if(wfc1.nbands!=wfc2.nbands){ cout<<"Error in overlap: Wavefunctions have different number of bands\n"; exit(0); }
  for(int i=minband;i<=maxband;i++){
    for(int j=minband;j<=maxband;j++){
      complex<double> res1 = wfc1.kpts[k1].mo[i].conj() * wfc2.kpts[k2].mo[j];
      complex<double> res2 = wfc2.kpts[k2].mo[i].conj() * wfc1.kpts[k1].mo[j];
      complex<double> res = ihbar*(0.5/dt)*(res1 - res2);

      if(i==j){ res += 0.5*(wfc1.kpts[k1].mo[i].energy + wfc2.kpts[k2].mo[i].energy); }

      f_re<<res.real()<<"  ";
      f_im<<res.imag()<<"  ";
    }// for j
    f_re<<endl;
    f_im<<endl;
  }// for i

  f_re.close();
  f_im.close();
}



