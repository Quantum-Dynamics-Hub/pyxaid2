/***********************************************************
 * Copyright (C) 2013 Alexey V. Akimov
 * This file is distributed under the terms of the
 * GNU General Public License as published by the
 * Free Software Foundation; either version 3 of the
 * License, or (at your option) any later version.
 * http://www.gnu.org/copyleft/gpl.txt
***********************************************************/

#ifndef wfc_h
#define wfc_h


#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <string>
#include <boost/python.hpp>
#include "matrix.h"
using namespace boost::python;
using namespace std;

// Some useful (potentially) macros
#define SWAP_2(x) ( (((x) & 0xff) << 8) | ((unsigned short)(x) >> 8) )
#define SWAP_4(x) ( ((x) << 24) | (((x) << 8) & 0x00ff0000) | \
         (((x) >> 8) & 0x0000ff00) | ((x) >> 24) )
#define FIX_SHORT(x) (*(unsigned short *)&(x) = SWAP_2(*(unsigned short *)&(x)))
#define FIX_LONG(x) (*(unsigned *)&(x) = SWAP_4(*(unsigned *)&(x)))
#define FIX_FLOAT(x) FIX_LONG(x)


template<class T>
void get_value(T& val,char* memblock,int& pos){
  memcpy((void*)&val, &memblock[pos], sizeof(T)); pos += sizeof(T);
//  cout<<val<<endl;
}


class MO{

public:
  int npw;                        // number of plane waves in MO expansion
  double energy;                  // energy of this orbital (eigenvalue)
  double gweight;
  double fweight;
//  vector< complex<float> > coeff_f; // expansion coefficients
  vector< complex<double> > coeff; // coefficients in double precision

  // Constructor
  MO(){ npw = 0;  energy = 0.0; gweight = 0.0; fweight = 0.0; }
  MO(int _npw){ 
    npw = _npw;
    energy = 0.0; gweight = 0.0; fweight = 0.0;
    coeff = vector<complex<double> >(npw,complex<double>(0.0,0.0));
  }
  // Copy constructor
//  MO(const& MO);

  // Destructor
  ~MO(){ if(coeff.size()>0){ coeff.clear(); } npw = 0; }

  // Operators
  MO operator-();                 // Negation;
  MO operator+(MO ob);
  MO operator-(MO ob);
  void operator+=(MO ob);
  void operator-=(MO ob);
  MO operator/(double num);
  MO operator/(complex<double> num);


  // Methods
  MO conj();
  void normalize();
  void complete();  // add complex conjugate part of the wavefunction

  // Friends
  friend MO operator*(const double& f,  const MO& m1);  // Multiplication of MO and double;
  friend MO operator*(const MO& m1, const double  &f);  // Multiplication of MO and double;
  friend MO operator*(const float& f,  const MO& m1);  // Multiplication of MO and float;
  friend MO operator*(const MO& m1, const float  &f);  // Multiplication of MO and float;
  friend MO operator*(const complex<double>& f,  const MO& m1);  // Multiplication of MO and complex<double>;
  friend MO operator*(const MO& m1, const complex<double>  &f);  // Multiplication of MO and complex<double>;
  friend MO operator*(const complex<float>& f,  const MO& m1);   // Multiplication of MO and complex<float>;
  friend MO operator*(const MO& m1, const complex<float>  &f);   // Multiplication of MO and complex<float>;

  friend complex<double> operator*(const MO& m1,  const MO& m2);  // Multiplication of MO and MO


};


class K_point{

public:
  int kx,ky,kz;  // indexes of this kpoint
  int nbands;    // number of bands (MOs) in this k-point
  int npw;       // number of plane waves in each MO
  vector<MO> mo; // array of MOs

  // Constructor 
  K_point(){  nbands = 0; npw = 0; kx = ky = kz = 0; }
  K_point(int _nbnds){ nbands = _nbnds; npw = 0; kx = ky = kz = 0; mo = vector<MO>(nbands,MO());  }
  K_point(int _nbnds,int _npw){ nbands = _nbnds; npw = _npw; kx = ky = kz = 0; mo = vector<MO>(nbands,MO(npw)); }
  K_point(K_point& k1,int min1,int max1, K_point& k2,int min2,int max2);

  // Copy constructor
  //K_point(const& K_point);
  
  // Destructor
  ~K_point(){ if(mo.size()>0){ mo.clear(); } nbands = 0; npw = 0; }

  // Methods
  void complete();
  void normalize();
  void transform(matrix&);
};


class wfc{

  void aux_line2vec(string line,vector<double>& a);

public:
  // Info
  int nspin;                // type of the spin-polarization used in calculations
  int gamma_only;           // gamma trick:  =1 (true), =0 (false)
  int natoms;               // number of atoms
  double tpiba;             // units of the lattice vectors (reciprocal)
  double alat;              // units of the lattice vectors (real)
  double omega;             // volume of the unit cell
  double efermi;            // fermi energy
  std::string cell_units;   // units of the simulation cell
  std::string energy_units; // units of the energy (eigenvalues)
  vector<double> a1,a2,a3;  // lattice vectors
  vector<double> b1,b2,b3;  // reciprocal lattice vectors

  // Data
  int nkpts;                // number of k-points
  int nbands;               // number of bands - the same for all k-points
  int npw;                  // number of plane waves in each mo - same for all mos
  int is_allocated;         // flag that says which level of memory is allocated 
  vector<K_point> kpts;     // array of k-points
  vector<vector<int> > grid;// internal vector<int> consists of 3 integers

  

  // Constructor
  wfc(){ nkpts = 0; nbands = 0; npw = 0; is_allocated = 0; }
  wfc(int _nkpts){ nkpts = _nkpts; nbands = 0; npw = 0; kpts = vector<K_point>(nkpts,K_point()); is_allocated = 1; }
  wfc(int _nkpts,int _nbnds){ nkpts = _nkpts; nbands = _nbnds; npw = 0; kpts = vector<K_point>(nkpts,K_point(nbands)); is_allocated = 2; }
  wfc(int _nkpts,int _nbnds,int _npw){nkpts = _nkpts; nbands = _nbnds; npw = _npw; kpts = vector<K_point>(nkpts,K_point(nbands,npw)); is_allocated = 3; }
  wfc(wfc& wfc1,int min1,int max1, wfc& wfc2,int min2,int max2);

  // Destructor
  ~wfc(){  if(kpts.size()>0){ kpts.clear(); } if(grid.size()>0){ grid.clear();} nkpts = 0; nbands = 0; npw = 0; is_allocated = 0; }

  // Copy constructor
  //wfc(const wfc&);

  // Methods:
  // QE methods
  void QE_read_binary_wfc(std::string filename,int,int,int);
  void QE_read_acsii_wfc(std::string filename);
  void QE_read_acsii_grid(std::string filename);
  void QE_read_acsii_index(std::string filename);

  // Common methods
  void complete();
  void normalize();
  void transform(int k,matrix&);
  void restore(int k1,int do_complete);
  void set_latt_vectors(boost::python::list);
  void set_reci_vectors(boost::python::list);
  void compute_Hprime(int minband,int maxband,std::string filename);
};


// Functions using the arguments of wfc type
void overlap(wfc& wfc1,int k1,int minband,int maxband,std::string filename);
void energy(wfc& wfc1,int k,int minband,int maxband,std::string filename);
void nac(wfc& wfc1,wfc& wfc2,int k1,int k2,int minband,int maxband,double dt,std::string filename);
void ham(wfc& wfc1,wfc& wfc2,int k1,int k2,int minband,int maxband,double dt,std::string filename);

#endif //wfc_h
