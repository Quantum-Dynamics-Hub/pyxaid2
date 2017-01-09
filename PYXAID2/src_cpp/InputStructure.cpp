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
#include "InputStructure.h"
#include "units.h"

using namespace std;

//================= Class InputStructure ==========================

void InputStructure::error(std::string var){
  cout<<"Error: variable "<<var<<" is not defined\n";
}
void InputStructure::warning(std::string var,std::string default_value){
  cout<<"Warning: variable "<<var<<" is not defined. Parameter is set to default value = "<<default_value<<"\n";
}

void InputStructure::init(){
  // Variables are not defined
  is_read_couplings =
//  is_many_electron_algorithm =
  is_namdtime = is_sh_algo = is_num_sh_traj =
  is_boltz_flag = is_debug_flag = is_Temp =
  is_nucl_dt = is_elec_dt = is_integrator =
  is_runtype = 
  is_Ham_re_prefix = is_Ham_re_suffix = 
  is_Ham_im_prefix = is_Ham_im_suffix =
  is_Hprime_x_prefix = is_Hprime_z_prefix = is_Hprime_z_prefix =
  is_Hprime_x_suffix = is_Hprime_z_suffix = is_Hprime_z_suffix = 
//  is_energy_prefix =
//  is_energy_suffix = is_nac_re_prefix = is_nac_re_suffix =
//  is_nac_im_prefix = is_nac_im_suffix = 
  is_energy_in_one_file =
  is_scratch_dir = is_energy_units = is_alp_bet = is_decoherence = 
  is_regress_mode = 
  is_is_field = is_field_dir = is_field_protocol = is_field_Tm = 
  is_field_T = is_field_freq = is_field_freq_units = is_field_fluence = 0;

}

void InputStructure::echo(){

  // Echo of successfully read input variables
  if(is_Ham_re_prefix){ cout<<"Ham_re_prefix = "<<Ham_re_prefix<<endl; }
  if(is_Ham_re_suffix){ cout<<"Ham_re_suffix = "<<Ham_re_suffix<<endl; }
  if(is_Ham_im_prefix){ cout<<"Ham_im_prefix = "<<Ham_im_prefix<<endl; }
  if(is_Ham_im_suffix){ cout<<"Ham_im_suffix = "<<Ham_im_suffix<<endl; }

  if(is_Hprime_x_prefix){ cout<<"Hprime_x_prefix = "<<Hprime_x_prefix<<endl; }
  if(is_Hprime_y_prefix){ cout<<"Hprime_y_prefix = "<<Hprime_y_prefix<<endl; }
  if(is_Hprime_z_prefix){ cout<<"Hprime_z_prefix = "<<Hprime_z_prefix<<endl; }
  if(is_Hprime_x_suffix){ cout<<"Hprime_x_suffix = "<<Hprime_x_suffix<<endl; }
  if(is_Hprime_y_suffix){ cout<<"Hprime_y_suffix = "<<Hprime_y_suffix<<endl; }
  if(is_Hprime_z_suffix){ cout<<"Hprime_z_suffix = "<<Hprime_z_suffix<<endl; }


//  if(is_energy_prefix){ cout<<"energy_prefix = "<<energy_prefix<<endl; }
//  if(is_energy_suffix){ cout<<"energy_suffix = "<<energy_suffix<<endl; }
  if(is_energy_units){ cout<<"energy_units = "<<energy_units<<endl; }
//  if(is_nac_re_prefix){ cout<<"nac_re_prefix = "<<nac_re_prefix<<endl; }
//  if(is_nac_re_suffix){ cout<<"nac_re_suffix = "<<nac_re_suffix<<endl; }
//  if(is_nac_im_prefix){ cout<<"nac_im_prefix = "<<nac_im_prefix<<endl; }
//  if(is_nac_im_suffix){ cout<<"nac_im_suffix = "<<nac_im_suffix<<endl; }
//  if(is_overlap_re_prefix){ cout<<"overlap_re_prefix = "<<overlap_re_prefix<<endl; }
//  if(is_overlap_re_suffix){ cout<<"overlap_re_suffix = "<<overlap_re_suffix<<endl; }
  if(is_energy_in_one_file){ cout<<"energy_in_one_file = "<<energy_in_one_file<<endl; }
  if(is_scratch_dir){ cout<<"scratch_dir = "<<scratch_dir<<endl; }
  if(is_read_couplings){ cout<<"read_couplings = "<<read_couplings<<endl; }
  if(is_read_overlaps){ cout<<"read_overlaps = "<<read_overlaps<<endl; }
//  if(is_many_electron_algorithm){ cout<<"many_electron_algorithm = "<<many_electron_algorithm<<endl; }
  if(is_namdtime){ cout<<"namdtime = "<<namdtime<<endl; }
  if(is_sh_algo){ cout<<"sh_algo = "<<sh_algo<<endl; }
  if(is_num_sh_traj){ cout<<"num_sh_traj = "<<num_sh_traj<<endl; }
  if(is_boltz_flag){ cout<<"boltz_flag = "<<boltz_flag<<endl; }
  if(is_debug_flag){ cout<<"debug_flag = "<<debug_flag<<endl; }
  if(is_Temp){ cout<<"Temp [K] = "<<Temp<<endl; }
  if(is_elec_dt){ cout<<"elec_dt [fs] = "<<elec_dt<<endl; }
  if(is_nucl_dt){ cout<<"nucl_dt [fs] = "<<nucl_dt<<endl; }
  if(is_integrator){ cout<<"integrator = "<<integrator<<endl; }
  if(is_runtype){ cout<<"runtype = "<<runtype<<endl; }
  if(is_alp_bet){ cout<<"alp_bet = "<<alp_bet<<endl; }
  if(is_decoherence){ cout<<"decoherence = "<<decoherence<<endl; }
  if(is_regress_mode){ cout<<"regress_mode = "<<regress_mode<<endl; }
  if(is_is_field){ cout<<"is_field = "<<is_field<<endl; }
  if(is_field_dir){ cout<<"field_dir = "<<field_dir<<endl; }
  if(is_field_protocol){ cout<<"field_protocol = "<<field_protocol<<endl; }
  if(is_field_Tm){ cout<<"field_Tm [fs] = "<<field_Tm<<endl; }
  if(is_field_T){ cout<<"field_T [fs] = "<<field_T<<endl; }
  if(is_field_freq){ cout<<"field_freq = "<<field_freq<<endl; }
  if(is_field_freq_units){ cout<<"field_freq_units = "<<field_freq_units<<endl; }
  if(is_field_fluence){ cout<<"field_fluence [mJ/cm^2] = "<<field_fluence<<endl; }

}

void InputStructure::set_default(){

  // Output warnings and errors for not-read variables
  int err_status = 0;
  int wrn_status = 0;

  if(!is_Ham_re_prefix){ warning("Ham_re_prefix","Ham");    Ham_re_prefix = "Ham";  is_Ham_re_prefix = 1; }
  if(!is_Ham_re_suffix){ warning("Ham_re_suffix","_re");    Ham_re_suffix = "_re";  is_Ham_re_suffix = 1; }
  if(!is_Ham_im_prefix){ warning("Ham_im_prefix","Ham");    Ham_im_prefix = "Ham";  is_Ham_im_prefix = 1; }
  if(!is_Ham_im_suffix){ warning("Ham_im_suffix","_im");    Ham_im_suffix = "_im";  is_Ham_im_suffix = 1; }

  if(!is_Hprime_x_prefix){warning("Hprime_x_prefix","Hprime_x_"); Hprime_x_prefix = "Hprime_x_"; is_Hprime_x_prefix = 1;}
  if(!is_Hprime_y_prefix){warning("Hprime_y_prefix","Hprime_y_"); Hprime_y_prefix = "Hprime_y_"; is_Hprime_y_prefix = 1;}
  if(!is_Hprime_z_prefix){warning("Hprime_z_prefix","Hprime_z_"); Hprime_z_prefix = "Hprime_z_"; is_Hprime_z_prefix = 1;}
  if(!is_Hprime_x_suffix){warning("Hprime_x_suffix","_re"); Hprime_x_suffix = "_re"; is_Hprime_x_suffix = 1;}
  if(!is_Hprime_y_suffix){warning("Hprime_y_suffix","_re"); Hprime_y_suffix = "_re"; is_Hprime_y_suffix = 1;}
  if(!is_Hprime_z_suffix){warning("Hprime_z_suffix","_re"); Hprime_z_suffix = "_re"; is_Hprime_z_suffix = 1;}


//  if(!is_energy_prefix){ warning("energy_prefix","energy"); energy_prefix="energy"; is_energy_prefix = 1; }
//  if(!is_energy_suffix){ warning("energy_suffix","");       energy_suffix = "";     is_energy_suffix = 1; }
  if(!is_energy_units){ warning("energy_units","Ry");       energy_units = "Ry";    is_energy_units = 1; }
//  if(!is_nac_re_prefix){ warning("nac_re_prefix","nac");    nac_re_prefix = "nac";  is_nac_re_prefix = 1; }
//  if(!is_nac_re_suffix){ warning("nac_re_suffix","_re");    nac_re_suffix = "_re";  is_nac_re_suffix = 1; }
//  if(!is_nac_im_prefix){ warning("nac_im_prefix","nac");    nac_im_prefix = "nac";  is_nac_im_prefix = 1; }
//  if(!is_nac_im_suffix){ warning("nac_im_suffix","_im");    nac_im_suffix = "_im";  is_nac_im_suffix = 1; }
//  if(!is_overlap_re_prefix){ warning("overlap_re_prefix","overlap");overlap_re_prefix = "overlap";is_overlap_re_prefix = 1;}
//  if(!is_overlap_re_suffix){ warning("overlap_re_suffix","_re"); overlap_re_suffix = "_re";  is_overlap_re_suffix = 1; }
  if(!is_energy_in_one_file){ warning("energy_in_one_file","false"); energy_in_one_file = "false"; is_energy_in_one_file = 1; }
  if(!is_scratch_dir){ warning("scratch_dir","/"); scratch_dir = "/"; is_scratch_dir = 1; }

  if(!is_read_couplings){ warning("read_couplings","online"); read_couplings = "online"; is_read_couplings = 1; wrn_status++; }
  if(!is_read_overlaps){ warning("read_overlaps","online"); read_overlaps = "online"; is_read_overlaps = 1; wrn_status++; }
  //if(!is_many_electron_algorithm){ warning("many_electron_algorithm","0"); many_electron_algorithm=0; is_many_electron_algorithm=1; wrn_status++; }
  if(!is_namdtime){ warning("namdtime","0"); namdtime = 0; is_namdtime = 1; wrn_status++; }
  if(!is_sh_algo){ warning("sh_algo","0"); sh_algo = 0; is_sh_algo = 1; wrn_status++; }
  if(!is_num_sh_traj){ warning("num_sh_traj","1"); num_sh_traj = 1; is_num_sh_traj = 1; wrn_status++; }
  if(!is_boltz_flag){ warning("boltz_flag","1"); boltz_flag=1; is_boltz_flag = 1; wrn_status++; }
  if(!is_debug_flag){ warning("debug_flag","0"); debug_flag=0; is_debug_flag = 1; wrn_status++; }
  if(!is_Temp){ warning("Temp","300.0"); Temp = 300.0; is_Temp = 1; wrn_status++; }
  if(!is_nucl_dt){ warning("nucl_dt","1.0"); nucl_dt = 1.0; is_nucl_dt = 1; wrn_status++; }
  if(!is_elec_dt){ warning("elec_dt","0.001"); elec_dt = 0.001; is_elec_dt = 1; wrn_status++; }
  if(!is_integrator){ warning("integrator","0"); integrator = 0; is_integrator = 1; wrn_status++; }
  if(!is_runtype){ warning("runtype","namd"); runtype = "namd"; is_runtype = 1; wrn_status++; }
  if(!is_alp_bet){ warning("alp_bet","0"); alp_bet = 0; is_alp_bet = 1; wrn_status++; }
  if(!is_decoherence){ warning("decoherence","0"); decoherence = 0; is_decoherence = 1; wrn_status++; }
  if(!is_regress_mode){ warning("regress_mode","0"); regress_mode = 0; is_regress_mode = 1; wrn_status++; }

  if(!is_is_field){ warning("is_field","0"); is_field = 0; is_is_field = 1; wrn_status++; }
  if(!is_field_dir){ warning("field_dir","xyz"); field_dir = "xyz"; is_field_dir = 1; wrn_status++; }
  if(!is_field_protocol){ warning("field_protocol","2"); field_protocol = 2; is_field_protocol = 1; wrn_status++; }
  if(!is_field_Tm){ warning("field_Tm","5.0"); field_Tm = 5.0; is_field_Tm = 1; wrn_status++; }
  if(!is_field_T){ warning("field_T","5.0"); field_T = 5.0; is_field_T = 1; wrn_status++; }
  if(!is_field_freq){ warning("field_freq","2.0"); field_freq = 2.0; is_field_freq = 1; wrn_status++; }
  if(!is_field_freq_units){ warning("field_freq_units","eV"); field_freq_units = "eV"; is_field_freq_units = 1; wrn_status++; }
  if(!is_field_fluence){ warning("field_fluence","1.0"); field_fluence = 1.0; is_field_fluence = 1; wrn_status++; }



  cout<<"Number of errors = "<<err_status<<endl;
  cout<<"Number of warning = "<<wrn_status<<endl;

  if(err_status>0){ cout<<"Exiting...\n"; exit(0); }
  else{ cout<<"Done with input parameters\n"; }

}

InputStructure::InputStructure(boost::python::dict params){
  init();
  boost::python::list lkeys = params.keys();

  for(int i=0;i<len(lkeys);i++){
    std::string s1;
    s1 = extract<std::string>(lkeys[i]);
//    if(s1=="icondfile"){ icondfile = extract<std::string>(params[s1]); is_icondfile = 1; }
//         if(s1=="energy_prefix"){ energy_prefix = extract<std::string>(params[s1]);  is_energy_prefix = 1; }
//    else if(s1=="energy_suffix"){ energy_suffix = extract<std::string>(params[s1]);  is_energy_suffix = 1; }
         if(s1=="energy_units"){ energy_units = extract<std::string>(params[s1]);  is_energy_units = 1; }

    else if(s1=="Ham_re_prefix"){ Ham_re_prefix = extract<std::string>(params[s1]);  is_Ham_re_prefix = 1; }
    else if(s1=="Ham_re_suffix"){ Ham_re_suffix = extract<std::string>(params[s1]);  is_Ham_re_suffix = 1; }
    else if(s1=="Ham_im_prefix"){ Ham_im_prefix = extract<std::string>(params[s1]);  is_Ham_im_prefix = 1; }
    else if(s1=="Ham_im_suffix"){ Ham_im_suffix = extract<std::string>(params[s1]);  is_Ham_im_suffix = 1; }

    else if(s1=="Hprime_x_prefix"){ Hprime_x_prefix = extract<std::string>(params[s1]); is_Hprime_x_prefix = 1; }
    else if(s1=="Hprime_y_prefix"){ Hprime_y_prefix = extract<std::string>(params[s1]); is_Hprime_y_prefix = 1; }
    else if(s1=="Hprime_z_prefix"){ Hprime_z_prefix = extract<std::string>(params[s1]); is_Hprime_z_prefix = 1; }
    else if(s1=="Hprime_x_suffix"){ Hprime_x_suffix = extract<std::string>(params[s1]); is_Hprime_x_suffix = 1; }
    else if(s1=="Hprime_y_suffix"){ Hprime_y_suffix = extract<std::string>(params[s1]); is_Hprime_y_suffix = 1; }
    else if(s1=="Hprime_z_suffix"){ Hprime_z_suffix = extract<std::string>(params[s1]); is_Hprime_z_suffix = 1; }


//    else if(s1=="nac_re_prefix"){ nac_re_prefix = extract<std::string>(params[s1]);  is_nac_re_prefix = 1; }
//    else if(s1=="nac_re_suffix"){ nac_re_suffix = extract<std::string>(params[s1]);  is_nac_re_suffix = 1; }
//    else if(s1=="nac_im_prefix"){ nac_im_prefix = extract<std::string>(params[s1]);  is_nac_im_prefix = 1; }
//    else if(s1=="nac_im_suffix"){ nac_im_suffix = extract<std::string>(params[s1]);  is_nac_im_suffix = 1; }
//    else if(s1=="overlap_re_prefix"){ overlap_re_prefix = extract<std::string>(params[s1]);  is_overlap_re_prefix = 1; }
//    else if(s1=="overlap_re_suffix"){ overlap_re_suffix = extract<std::string>(params[s1]);  is_overlap_re_suffix = 1; }
    else if(s1=="energy_in_one_file"){ energy_in_one_file = extract<std::string>(params[s1]); is_energy_in_one_file = 1; }
    else if(s1=="scratch_dir"){ scratch_dir = extract<std::string>(params[s1]); is_scratch_dir = 1; }

    else if(s1=="read_couplings") { read_couplings = extract<std::string>(params[s1]); is_read_couplings = 1; }
    else if(s1=="read_overlaps") { read_overlaps = extract<std::string>(params[s1]); is_read_overlaps = 1; }
//    else if(s1=="many_electron_algorithm"){ many_electron_algorithm = extract<int>(params[s1]); is_many_electron_algorithm = 1; }
    else if(s1=="namdtime"){ namdtime = extract<int>(params[s1]); is_namdtime = 1; }
    else if(s1=="sh_algo"){ sh_algo = extract<int>(params[s1]); is_sh_algo = 1; }
    else if(s1=="num_sh_traj"){ num_sh_traj = extract<int>(params[s1]); is_num_sh_traj = 1; }
    else if(s1=="boltz_flag"){ boltz_flag = extract<int>(params[s1]); is_boltz_flag = 1; }
    else if(s1=="debug_flag"){ debug_flag = extract<int>(params[s1]); is_debug_flag = 1; }
    else if(s1=="Temp"){ Temp = extract<double>(params[s1]); is_Temp = 1; }
    else if(s1=="nucl_dt"){ nucl_dt = extract<double>(params[s1]); is_nucl_dt = 1; }
    else if(s1=="elec_dt"){ elec_dt = extract<double>(params[s1]); is_elec_dt = 1; }
    else if(s1=="integrator"){ integrator = extract<int>(params[s1]); is_integrator = 1; }
    else if(s1=="runtype"){ runtype = extract<std::string>(params[s1]); is_runtype = 1; }

    else if(s1=="alp_bet"){ alp_bet = extract<int>(params[s1]); is_alp_bet = 1; }
    else if(s1=="decoherence"){ decoherence = extract<int>(params[s1]); is_decoherence = 1; }
    else if(s1=="regress_mode"){ regress_mode = extract<int>(params[s1]); is_regress_mode = 1; }

    else if(s1=="is_field"){ is_field = extract<int>(params[s1]); is_is_field = 1; }
    else if(s1=="field_dir"){ field_dir = extract<std::string>(params[s1]); is_field_dir = 1; }
    else if(s1=="field_protocol"){ field_protocol = extract<int>(params[s1]); is_field_protocol = 1; }
    else if(s1=="field_Tm"){ field_Tm = extract<double>(params[s1]); is_field_Tm = 1; }
    else if(s1=="field_T"){ field_T = extract<double>(params[s1]); is_field_T = 1; }
    else if(s1=="field_freq"){ field_freq = extract<double>(params[s1]); is_field_freq = 1; }
    else if(s1=="field_freq_units"){ field_freq_units = extract<std::string>(params[s1]); is_field_freq_units = 1; }
    else if(s1=="field_fluence"){ field_fluence = extract<double>(params[s1]); is_field_fluence = 1; }


  }// for i

  echo();
  set_default();
  sanity_check();
}

void InputStructure::sanity_check(){

  // Regression mode
  if(regress_mode==0 || regress_mode==1){ ;; }
  else{
    cout<<"Error: regress_mode = "<<regress_mode<<" is not known\n"; 
    cout<<"Allowed values are:\n";
    cout<<"        0    -  y(t) = a + b*t, with a = 0 and b being fitting parameter\n";
    cout<<"        1    -  y(t) = a + b*t, with both a and b being fitting parameters\n";
    cout<<"Exiting...\n";
    exit(0);
  }

  // File reading-related options
  if(read_couplings=="online" || read_couplings=="batch" || 
     read_couplings=="online_all_in_one" || read_couplings=="batch_all_in_one"){ ;; }
  else{
    cout<<"Error: read_couplings = "<<read_couplings<<" is not known\n";
    cout<<"Allowed values are:\n";
    cout<<"     online   -   read only those files that are needed for each given initial condition"
        <<" at the time the NA-MD with this initial condition is to be computed\n";
    cout<<"     batch    -   read all the files, needed for all initial conditions presented\n";
    cout<<"     online_all_in_one -  same as online, but real and imaginary parts of Hamiltonian are"
        <<" extracted from the same file\n";
    cout<<"     batch_all_in_one -  same as batch, but real and imaginary parts of Hamiltonian are"
        <<" extracted from the same file\n";

    cout<<"Exiting...\n";
    exit(0);
  }

  // Integrator-related options
  if(integrator==0 || integrator==10 || integrator==11 || integrator==2){ ;; }
  else{
    cout<<"Error: integrator = "<<integrator<<" is not known\n";
    cout<<"Allowed values are:\n";
    cout<<"     0   - Trotter splitting scheme (default)\n";
    cout<<"     10  - Finite difference with first order for the dH/dt evaluation\n";
    cout<<"     11  - Finite difference with second order for the dH/dt evaluation\n";
    cout<<"     2   - Exact solution (matrix exponent). May be very slow for big # of states\n";
    cout<<"Exiting...\n";
    exit(0);
  }

  // Time-step related options
  if(nucl_dt<elec_dt){
    cout<<"Error: Nuclear time step nucl_dt = "<<nucl_dt
        <<" can not be smaller that the electronic time step elec_dt = "<<elec_dt<<endl;
    cout<<"Please check these values\n";
    cout<<"Exiting...\n";
    exit(0);
  }

  // Field-related options
  if(is_field){
    if(integrator!=0){
      cout<<"Error: Field is only implemented for integrator = "<<integrator<<endl;
      cout<<"Exiting...\n";
      exit(0);
    }
    if(decoherence>0){
      cout<<"Error: Field is not yet implemented with decoherence\nExiting...\n"; exit(0);
    }
    if(is_field_dir){
      if(field_dir=="x"||field_dir=="y"||field_dir=="z"||
         field_dir=="xy"||field_dir=="xz"||field_dir=="yz"||field_dir=="xyz"
        ){ ;; }
      else{  
        cout<<"Error: field_dir = "<<field_dir<<" is not known\n";
        cout<<"Allowed values are: \"x\", \"y\", \"z\", \"xy\", \"xz\", \"yz\", \"xyz\" \n";
        cout<<"Exiting...\n"; 
        exit(0);
      }
    }
    else{ cout<<"Error: field_dir must be specified\nExiting...\n"; exit(0); }

    if(is_field_protocol){
      if(field_protocol==1|| field_protocol==2 || field_protocol==3){ ;; }
      else{  
        cout<<"Error: field_protocol = "<<field_protocol<<" is not known\n";
        cout<<"Allowed values are:\n";
        cout<<"   1  -  for step function modulation protocol\n";
        cout<<"   2  -  for saw-tooth modulation protocol\n";
        cout<<"   3  -  for constant amplitude (given by fluence value)\n";
        cout<<"Exiting...\n";
        exit(0);
      }
    }
    else{ cout<<"Error: field_protocol must be specified\nExiting...\n"; exit(0); }
  

    // Check 1 - adjust T
    if(is_field_T && is_field_freq && is_field_freq_units){
      double omega = 0.0; // angular frequency
           if(field_freq_units=="1/fs"){ omega = 2.0*M_PI*field_freq; }     // input is linear frequency
      else if(field_freq_units=="rad/fs"){ omega = field_freq;  }           // input is angular frequency
      else if(field_freq_units=="eV"){ omega = field_freq/hbar; }           // input is energy in eV
      else if(field_freq_units=="nm"){ omega = 2.0*M_PI*300.0/field_freq; } // input is wavelength in nm
      else{  cout<<"Units of the field frequency must be specified. Exiting...\n"; exit(0); }

      double oldT,newT;
      oldT = field_T;
      newT = (2.0*M_PI/omega)*floor(oldT*omega*0.5/M_PI);
      cout<<"Duration of the laser pulse is adjusted to the field freqency. Old period = "<<oldT<<" new period = "<<newT<<endl;

      if(is_field_protocol && is_field_Tm){
        // Care in this case, because if tg(omega*Tm) == 1 the amplitude must be infinitely large
        if(field_protocol==2){
          if(tan(omega*field_Tm)==1){ 
            cout<<"Error: field_Tm results in singular vector potential.\n";
            cout<<"Choose another value of field_Tm. Exiting...\n";
            exit(0);
          }
        
        }// protocol==2
      }
    }// Check 1


  }// is_field
}
