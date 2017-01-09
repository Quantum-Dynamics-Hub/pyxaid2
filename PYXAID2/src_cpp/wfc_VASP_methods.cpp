#include "wfc.h"

void wfc::VASP_read_wfc(std::string filename){

  ifstream::pos_type size;
  char * memblock;

  ifstream file;
  file.open(filename.c_str(), ios::in|ios::binary|ios::ate);
  if (file.is_open())
  {
    // First - read all file
    size = file.tellg();
    memblock = new char [size];
    file.seekg (0, ios::beg);
    file.read (memblock, size);
    file.close();

//    cout<<"size of the file = "<<size<<" bytes\n";

    // Read first record
    int pos = 0;  // position in file (in bytes)
    double tmp;
    double rdum,rispin,rtag;

//================== First record ========================
    get_value<double>(rdum,memblock,pos);
    get_value<double>(rispin,memblock,pos);
    get_value<double>(rtag,memblock,pos);

    pos += ((int)rdum - 3*sizeof(double));

//================= Second record ========================
    double inkpt, inband,emax;
    get_value<double>(inkpt,memblock,pos);      nkpts = (int)inkpt;
    get_value<double>(inband,memblock,pos);     nbands = (int)inband;
    get_value<double>(emax,memblock,pos);
    kpts = std::vector<K_point>(nkpts,K_point());

    vector<double> A(9,0.0); for(int i=0;i<9;i++){ get_value<double>(A[i],memblock,pos); }

    pos += ((int)rdum - 12*sizeof(double));

//=============== Other records ==========================
//   1-st k-point
//              1-st  record: npl(1), kpt, eig,gweight,fweight
//     1-st band
//              2-nd  record: coeff for 1-st band (npl entries)
//     2-nd band
//              3-rd  record: coeff for 2-nd band (npl entries)
//      ...
//     nband-th band
//              (nband+1)-th record: coeff for nband-th band (npl entries)
//
//
//  2-nd k-point
//             (nband+2)-th record: npl(2), kpt, eig, gweight, fweight
//     1-st band
//             (nband+3)-th record: coeff for 1-st band (npl(2) entries)
// ... etc.


    double npl;
    for(int ikpt=0;ikpt<nkpts;ikpt++){
      kpts[ikpt].nbands = nbands;
      kpts[ikpt].mo = std::vector<MO>(nbands,MO());

      //============= 1-st sub-record ==========================
      get_value<double>(npl,memblock,pos);
      get_value<double>(tmp,memblock,pos); kpts[ikpt].kx = (int)tmp;
      get_value<double>(tmp,memblock,pos); kpts[ikpt].ky = (int)tmp;
      get_value<double>(tmp,memblock,pos); kpts[ikpt].kz = (int)tmp;

      for(int iband=0;iband<kpts[ikpt].nbands;iband++){
        kpts[ikpt].mo[iband].npw = (int)npl;
        get_value<double>(kpts[ikpt].mo[iband].energy,memblock,pos);
        get_value<double>(kpts[ikpt].mo[iband].gweight,memblock,pos);
        get_value<double>(kpts[ikpt].mo[iband].fweight,memblock,pos);
      }// for bands

      pos += ((int)rdum - (3*kpts[ikpt].nbands+4)*sizeof(double));

      for(iband=0;iband<kpts[ikpt].nbands;iband++){
      //=============== iband-th sub-sub-record ================
        kpts[ikpt].mo[iband].coeff = std::vector< complex<double> >(kpts[ikpt].mo[iband].npw,std::complex<double>(0.0,0.0));
        for(int i=0;i<kpts[ikpt].mo[iband].npw;i++){
          std::complex<float> tmp;
          get_value< std::complex<float> >(tmp,memblock,pos);
          kpts[ikpt].mo[iband].coeff[i] = tmp;
        }
        pos += ((int)rdum - kpts[ikpt].mo[iband].npw*sizeof( std::complex<float>));
      }// for bands

    }// for k-points

    // Free memory
    delete[] memblock;
  }
  else{ cout << "Unable to open file"<<filename<<"\n"; }


}

