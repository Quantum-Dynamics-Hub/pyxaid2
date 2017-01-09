/***********************************************************
 * Copyright (C) 2013 Alexey V. Akimov
 * This file is distributed under the terms of the
 * GNU General Public License as published by the
 * Free Software Foundation; either version 3 of the
 * License, or (at your option) any later version.
 * http://www.gnu.org/copyleft/gpl.txt
***********************************************************/

#include "wfc.h"

using namespace std;

//======================= MO class methods ========================

MO MO::operator-(){
  MO res; res = *this;
  for(int i=0;i<npw;i++){ res.coeff[i] = -coeff[i]; }
  return res;
}
MO MO::operator+(MO m){
  MO res; res = *this;
  for(int i=0;i<npw;i++){ res.coeff[i] = coeff[i] + m.coeff[i]; }
  return res;
}
MO MO::operator-(MO m){
  MO res; res = *this;
  for(int i=0;i<npw;i++){ res.coeff[i] = coeff[i] - m.coeff[i]; }
  return res;
}
void MO::operator+=(MO m){
  for(int i=0;i<npw;i++){ coeff[i] += m.coeff[i]; }
}
void MO::operator-=(MO m){
  for(int i=0;i<npw;i++){ coeff[i] -= m.coeff[i]; }
}
MO MO::operator/(double num){
  MO res; res = *this;
  for(int i=0;i<npw;i++){ res.coeff[i] /= num;  }
  return res;
}
MO MO::operator/(complex<double> num){
  MO res; res = *this;
  for(int i=0;i<npw;i++){ res.coeff[i] /= num; }
  return res;
}

// Copy constructor
/*
MO::MO(const& MO mo){

  npw = mo.npw;
  energy = mo.energy;
  gweight = mo.gweight;
  fweight = mo.fweight;
  coeff = mo.coeff;
}
*/

MO MO::conj(){
// Return MO which is conjugate to original one
  MO res(npw);
  res.energy = energy; res.gweight = gweight; res.fweight = fweight;
  for(int i=0;i<npw;i++){ res.coeff[i] = std::conj(coeff[i]); }
  return res;
}

void MO::normalize(){
  double norm = 0.0;
  for(int i=0;i<npw;i++){norm += (std::conj(coeff[i]) * coeff[i] ).real();  }
  norm = sqrt(1.0/norm);
  for(i=0;i<npw;i++){ coeff[i] *= norm; }
}

void MO::complete(){
  // Complete the wfc by adding the complex conjugate part
  // First find the index of the G=0
  double tol = 1e-10;
  int indx=0;  
  // The general procedure - find the coefficient for G = (0,0,0)
  // for(int i=0;i<npw;i++){  if( fabs(::imag(coeff[i]))<tol) { indx = i; break; }  }
  // In QE this is simple: indx = 0 corresponds to G=0

  // Now add remaining part
  coeff.resize(2*npw-1);
  double norm = (std::conj(coeff[0]) * coeff[0] ).real();
  for(int i=1;i<npw;i++){
    coeff[npw+i] = std::conj(coeff[i]); 
    norm += 2.0*(std::conj(coeff[i]) * coeff[i] ).real();
  }

  // Update the number of the planewaves in new (completed) wfc
  npw =  2*npw - 1;

  // Finally, normalize the completed wfc
  norm = sqrt(1.0/norm);
  for(i=0;i<npw;i++){ coeff[i] *= norm; }

    
}

template<class T>
MO multiply(T& f,  const MO& m1){
  MO res; res = m1;
  for(int i=0;i<res.npw;i++){ res.coeff[i] *= f; }
  return res;
}

MO operator*(const double& f,  const MO& m1){ return multiply<const double>(f,m1); }
MO operator*(const MO& m1, const double  &f){ return multiply<const double>(f,m1); }
MO operator*(const float& f,  const MO& m1){ return multiply<const float>(f,m1); }
MO operator*(const MO& m1, const float  &f){ return multiply<const float>(f,m1); }
MO operator*(const complex<double>& f,  const MO& m1){ return multiply<const complex<double> >(f,m1); }
MO operator*(const MO& m1, const complex<double>  &f){ return multiply<const complex<double> >(f,m1); }
MO operator*(const complex<float>& f,  const MO& m1){ return multiply<const complex<float> >(f,m1); }
MO operator*(const MO& m1, const complex<float>  &f){ return multiply<const complex<float> >(f,m1); }

complex<double> operator*(const MO& m1,  const MO& m2){
  complex<double> res(0.0,0.0);
  if(m1.npw!=m2.npw){ cout<<"Error: Can not multiply MOs with different basis sizes\n"; exit(0); }
  else{   for(int i=0;i<m1.npw;i++){ res += m1.coeff[i] * m2.coeff[i];}  }
  return res;
}
//==================== K_point methods ================================

K_point::K_point(K_point& k1,int min1,int max1, K_point& k2,int min2,int max2){
// This is a very handy constructor for bouilding multiconfigurational Hamiltonian
// it takes mo in range [min1,max1] from w1
//      and mo in range [min2,max2] from w2

  if(k1.kx!=k2.kx){ cout<<"Error in K_point constructor: Can not merge two k-points with different kx\n"; exit(0); }
  if(k1.ky!=k2.ky){ cout<<"Error in K_point constructor: Can not merge two k-points with different ky\n"; exit(0); }
  if(k1.kz!=k2.kz){ cout<<"Error in K_point constructor: Can not merge two k-points with different kz\n"; exit(0); }

  kx = k1.kx;
  ky = k1.ky;
  kz = k1.kz;

  nbands = max1-min1+1 + max2-min2+1;

  if(k1.npw!=k2.npw){
    cout<<"Error in K_point constructor: Cannot construct a wfc from the wfcs with different basis sizes (npw)\n";
    exit(0);
  }

  npw = k1.npw;

  if(mo.size()>0){ mo.clear(); }

  if(max1>=min2){ cout<<"Warning in K_point constructor: There are several identical bands (MO) in given K_point\n"; }
  for(int i=min1;i<=max1;i++){  mo.push_back(k1.mo[i]); }
  for(    i=min2;i<=max2;i++){  mo.push_back(k2.mo[i]); }

}


// Copy constructor
/*
K_point::K_point(const& K_point k){

  kx = k.kx; ky = k.ky; kz = k.kz;
  nbands = k.nbands;
  npw = k.npw;
  mo = k.mo;
}
*/

void K_point::complete(){
  for(int i=0;i<nbands;i++){  mo[i].complete();  }
  npw = 2*npw - 1;
}

void K_point::normalize(){
  for(int i=0;i<nbands;i++){  mo[i].normalize();  }
}


void K_point::transform(matrix& T){
// T - is nbands x nbands matrix which mixes original MOs to make new LC of MOs
  vector<MO> tmp_mo = mo;
  for(int i=0;i<nbands;i++){
    MO tmp(mo[i].npw); 
    tmp.energy = mo[i].energy;
    tmp.gweight = mo[i].gweight;
    tmp.fweight = mo[i].fweight;

    for(int j=0;j<nbands;j++){
      tmp += T.M[i*nbands+j] * mo[j];
    }
    tmp_mo[i] = tmp;
  }
  mo = tmp_mo;
}

//===================== wfc class methods ============================

wfc::wfc(wfc& w1,int min1,int max1, wfc& w2,int min2,int max2){
// This is a very handy constructor for building of the multiconfigurational Hamiltonian
// it takes mo in range [min1,max1] from w1
//      and mo in range [min2,max2] from w2

// Takes all parameters from the first wfc
  nspin = w1.nspin;
  gamma_only = w1.gamma_only;
  natoms = w1.natoms;
  tpiba = w1.tpiba;
  alat = w1.alat;
  omega = w1.omega;
  efermi = w1.efermi;
  cell_units = w1.cell_units;
  energy_units = w1.energy_units;
  a1 = w1.a1; a2 = w1.a2; a3 = w1.a3;
  b1 = w1.b1; b2 = w1.b2; b3 = w1.b3;

  nkpts = w1.nkpts;

  nbands = max1-min1+1 + max2-min2+1;

  if(w1.npw!=w2.npw){
    cout<<"Error in wfc constructor: Cannot construct a wfc from the wfcs with different basis sizes\n";
    exit(0);
  }
  npw = w1.npw;
  is_allocated = w1.is_allocated;

  // Finally, pick up the mos:
  for(int k=0;k<nkpts;k++){
    K_point kpt(w1.kpts[k],min1,max1,w2.kpts[k],min2,max2); 
    kpts.push_back(kpt);
  }// for k


}

// Copy constructor
/*
wfc::wfc(const wfc& w1){

  nspin = w1.nspin;
  gamma_only = w1.gamma_only;
  natoms = w1.natoms;
  tpiba = w1.tpiba;
  alat = w1.alat;
  omega = w1.omega;
  efermi = w1.efermi;
  cell_units = w1.cell_units;
  energy_units = w1.energy_units;
  a1 = w1.a1; a2 = w1.a2; a3 = w1.a3;
  b1 = w1.b1; b2 = w1.b2; b3 = w1.b3;

  nkpts = w1.nkpts;
  nbands = w1.nbands;
  npw = w1.npw;
  is_allocated = w1.is_allocated;
  kpts = w1.kpts;
  grid = w1.grid;
}
*/

void wfc::complete(){
  for(int i=0;i<nkpts;i++){ kpts[i].complete(); }
  npw = 2*npw - 1;
}

void wfc::normalize(){
  for(int i=0;i<nkpts;i++){ kpts[i].normalize(); }
}

void wfc::transform(int k,matrix& T){
  kpts[k].transform(T);
}

void wfc::restore(int k1,int do_complete){

// This function restores a true wavefunction from the projected one
// Note the projected wavefunction is not orthonormal: <i|j> != delta_ij
// The true wavefunction is orthonormal.
// Some math: phi - projected wfc, psi - true wfc, then:
// phi * phi^+ = S  -iverlap matrix  phi^+ is phi^H - Hermitian conjugate
// we want to find such a transformation T:  psi = T * phi that makes the 
// matrix  psi * psi^+ to be the identity (I) matrix
// psi * psi^+ = T * phi * phi^+ * T^+ = T * S * T^+ = I
// Then T = S^-1/2 if S^-1/2 is Hermitian i.e. S^-1/2 = (S^-1/2)^+ than all is ok
// To find S^-1/2 we do: S = Q * L * Q^-1  - eighenvalue problem for S matrix (Hermitian)
// then S^-1/2 = Q * L^-1/2 * Q^-1 , L is diagonal
// The eigenvectors Q are orthonormal: Q^+ * Q = I so Q^+ = Q^-1
// Thus, the S^-1/2 is Hermitian. Indeed: (S^-1/2)^+ = (Q * L^-1/2 * Q^-1)^+ = 
// =(Q^-1)^+ * (L^-1/2)^+ * Q^+ = Q * L^-1/2 * Q^-1 = S^-1/2

// Parameters: 
// *this - is a function which is to be transformed, it will be overwritten with the true one
// k1 - is a k-point index of the wfc to be transformed
// do_complete - the flag to indicate if we need to complete the input wfc first (=1) or not (=0)


  // Complete 
  if(do_complete==1){ complete(); }

  // Compute overlap matrix
  matrix S(nbands,nbands);

  for(int i=0;i<nbands;i++){
    for(int j=0;j<nbands;j++){
      S.M[i*nbands+j] = kpts[k1].mo[i].conj() * kpts[k1].mo[j];
    }
  }

  //cout<<"S matrix is formed\n";
  //cout<<"S = "<<S<<endl;
  // Find the transformation matrix T
  // which makes S matrix diagonal - eigenvalue problem

  matrix eval(nbands,nbands);
  matrix evec(nbands,nbands);
  matrix eval_sqrt(nbands,nbands);
  //matrix evec_inv(nbands,nbands);

  //cout<<"entering eigen\n";
  S.eigen(1e-12,eval,evec,2); // evect * eval = S * evect
  //cout<<"eval = "<<eval<<endl;
  //cout<<"evec = "<<evec<<endl;
  //evec.inverse(1e-12,evec_inv,1); // We don't need inverse - it is just .H()
  //cout<<"inverse = "<<evec_inv<<endl;
  eval_sqrt = 0.0;
  for(i=0;i<nbands;i++){  eval_sqrt.M[i*nbands+i] = 1.0/ sqrt(eval.M[i*nbands+i]); }

  matrix T(nbands,nbands);
//  T = evec * eval_sqrt * evec_inv; // T = S^-1/2
  T = evec.conj() * eval_sqrt.conj() * evec.T(); // T = (S^-1/2)^*, here evec^-1 = evec.H(), so (evec^-1)^* = evec^T
//  cout<<"T = "<<T<<endl;


  // Finally, restore the real wavefunction (orthonormal)
  transform(k1,T);

}


void wfc::set_latt_vectors(boost::python::list lst){
// Format: lst = a1x, a1y, a1z, a2x, a2y, a2z, a3x, a3y, a3z
 
  int sz = len(lst);
  vector<double> tmp(sz,0.0);
  for(int i=0;i<sz;i++){ tmp[i] = boost::python::extract<double>(lst[i]);  }

  if(sz>=9){
    a1 = vector<double>(3,0.0); a1[0] = tmp[0]; a1[1] = tmp[1]; a1[2] = tmp[2];
    a2 = vector<double>(3,0.0); a2[0] = tmp[3]; a2[1] = tmp[4]; a2[2] = tmp[5];
    a3 = vector<double>(3,0.0); a3[0] = tmp[6]; a3[1] = tmp[7]; a3[2] = tmp[8];

  }else{ cout<<"Error: in set_latt_vectors. Too few items in the list : "<<sz<<" Expected at least 9\n"; exit(0); }
}
 
void wfc::set_reci_vectors(boost::python::list lst){
// Format: lst = b1x, b1y, b1z, b2x, b2y, b2z, b3x, b3y, b3z

  int sz = len(lst);
  vector<double> tmp(sz,0.0);
  for(int i=0;i<sz;i++){ tmp[i] = boost::python::extract<double>(lst[i]);  }

  if(sz>=9){
    b1 = vector<double>(3,0.0); b1[0] = tmp[0]; b1[1] = tmp[1]; b1[2] = tmp[2];
    b2 = vector<double>(3,0.0); b2[0] = tmp[3]; b2[1] = tmp[4]; b2[2] = tmp[5];
    b3 = vector<double>(3,0.0); b3[0] = tmp[6]; b3[1] = tmp[7]; b3[2] = tmp[8];

  }else{ cout<<"Error: in set_latt_vectors. Too few items in the list : "<<sz<<" Expected at least 9\n"; exit(0); }

}



void wfc::compute_Hprime(int minband,int maxband,std::string filename){
// Compute H'_ij, where 0<=minband<=i,j<=maxband<=nbands-1
// Output results in the file: filename

  int g_sz = grid.size(); // grid size
  int sz = maxband - minband + 1;
  int npw = kpts[0].mo[0].npw;
  int is_compl = 0;
  if(minband<0){ cout<<"Error in compute_Hprime: minband<0 : minband = "<<minband<<endl; exit(0); }
  if(maxband>=nbands){ cout<<"Error in compute_Hprime: maxband>=nbands : maxband = "<<maxband<<" nbands = "<<nbands<<endl; exit(0); }
  if(g_sz!=npw){ 
    if(npw==(2*g_sz-1)){
      cout<<"Warning: Using reconstructed (completed) wavefunction\n";
      is_compl = 1;
    }else{
    cout<<"Warning in compute_Hprime: number of plane waves = "
        <<npw<<" is different from the grid size ="<<g_sz<<endl;
    }
    g_sz = min(g_sz,npw);
  }
  if(b1.size()<3){ cout<<"Error in compute_Hprime: b1 vector is not defined\n"; exit(0); }
  if(b2.size()<3){ cout<<"Error in compute_Hprime: b2 vector is not defined\n"; exit(0); }
  if(b3.size()<3){ cout<<"Error in compute_Hprime: b3 vector is not defined\n"; exit(0); }



// Open file
  ofstream outx_re((filename+"x_re").c_str());
  ofstream outx_im((filename+"x_im").c_str());
  ofstream outy_re((filename+"y_re").c_str());
  ofstream outy_im((filename+"y_im").c_str());
  ofstream outz_re((filename+"z_re").c_str());
  ofstream outz_im((filename+"z_im").c_str());



// Compute
  complex<double> scl(1.0,0.0); // scaling factor
  complex<double> scl1(0.0,1.0);
  
  for(int i=minband;i<=maxband;i++){
    for(int j=minband;j<=maxband;j++){

      complex<double> Hx(0.0,0.0),Hy(0.0,0.0),Hz(0.0,0.0);

      for(int g=0;g<g_sz;g++){
        complex<double> tmp,gx,gy,gz;
        tmp = ( (std::conj(kpts[0].mo[i].coeff[g])) * kpts[0].mo[j].coeff[g] );
        gx = scl*(grid[g][0]*b1[0] +  grid[g][1]*b2[0] +  grid[g][2]*b3[0]);
        gy = scl*(grid[g][0]*b1[1] +  grid[g][1]*b2[1] +  grid[g][2]*b3[1]);
        gz = scl*(grid[g][0]*b1[2] +  grid[g][1]*b2[2] +  grid[g][2]*b3[2]);
 
        if(is_compl==0){
          Hx += tmp * gx;  Hy += tmp * gy;  Hz += tmp * gz;
        }
        else if(is_compl==1){
          if(g==0){ Hx += tmp * gx;  Hy += tmp * gy;  Hz += tmp * gz; } // This should give zero!
          // Now the Hprime_ matrices are purely imaginary, for the case of gamma-symmetry.
          else{ Hx += 2.0*scl1*tmp.real() * gx;  Hy += 2.0*scl1*tmp.real() * gy;  Hz += 2.0*scl1*tmp.real() * gz; }
        }
        
      }// for g

      outx_re<<Hx.real()<<"  ";  outx_im<<Hx.imag()<<"  ";
      outy_re<<Hy.real()<<"  ";  outy_im<<Hy.imag()<<"  ";
      outz_re<<Hz.real()<<"  ";  outz_im<<Hz.imag()<<"  ";
    }// for j
    outx_re<<endl;  outx_im<<endl;
    outy_re<<endl;  outy_im<<endl;
    outz_re<<endl;  outz_im<<endl;
  }// for i

  outx_re.close();  outx_im.close();
  outy_re.close();  outy_im.close();
  outz_re.close();  outz_im.close();
}

