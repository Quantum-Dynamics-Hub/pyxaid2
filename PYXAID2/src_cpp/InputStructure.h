/***********************************************************
 * Copyright (C) 2013 Alexey V. Akimov
 * This file is distributed under the terms of the
 * GNU General Public License as published by the
 * Free Software Foundation; either version 3 of the
 * License, or (at your option) any later version.
 * http://www.gnu.org/copyleft/gpl.txt
***********************************************************/

#ifndef InputStructure_H
#define InputStructure_H

#include <string>
#include <boost/python.hpp>
using namespace boost::python;
using namespace std;

class InputStructure{

  void init();
  void echo();
  void set_default();
  void error(std::string);
  void warning(std::string,std::string);
  void sanity_check();

public:
  std::string Ham_re_prefix;    int is_Ham_re_prefix;
  std::string Ham_re_suffix;    int is_Ham_re_suffix;
  std::string Ham_im_prefix;    int is_Ham_im_prefix;
  std::string Ham_im_suffix;    int is_Ham_im_suffix;

  // only real components are non-zero, so below files are for _re
  std::string Hprime_x_prefix;  int is_Hprime_x_prefix;
  std::string Hprime_y_prefix;  int is_Hprime_y_prefix;
  std::string Hprime_z_prefix;  int is_Hprime_z_prefix;
  std::string Hprime_x_suffix;  int is_Hprime_x_suffix;
  std::string Hprime_y_suffix;  int is_Hprime_y_suffix;
  std::string Hprime_z_suffix;  int is_Hprime_z_suffix;

//  std::string energy_prefix;  int is_energy_prefix;
//  std::string energy_suffix;  int is_energy_suffix;
  std::string energy_units;   int is_energy_units;
//  std::string nac_re_prefix;  int is_nac_re_prefix;
//  std::string nac_re_suffix;  int is_nac_re_suffix;
//  std::string nac_im_prefix;  int is_nac_im_prefix;
//  std::string nac_im_suffix;  int is_nac_im_suffix;
//  std::string overlap_re_prefix; int is_overlap_re_prefix;
//  std::string overlap_re_suffix; int is_overlap_re_suffix;
  std::string energy_in_one_file;  int is_energy_in_one_file;
  std::string scratch_dir;    int is_scratch_dir;

  std::string read_couplings; int is_read_couplings;
  std::string read_overlaps;  int is_read_overlaps;
//  int many_electron_algorithm;  int is_many_electron_algorithm;
  int integrator;   int is_integrator;     // choose integration algorithm
  double nucl_dt;   int is_nucl_dt;        // nuclear time step in fs
  double elec_dt;   int is_elec_dt;        // electronic time step in fs
  int namdtime;     int is_namdtime;
  int sh_algo;      int is_sh_algo;        // surface hopping algorithm: 0 = FSSH, 1 = GFSH, 2 = MSSH
  int num_sh_traj;  int is_num_sh_traj;
  int boltz_flag;   int is_boltz_flag;
  double Temp;      int is_Temp;           // Temperature
  int debug_flag;   int is_debug_flag;
  std::string runtype;    int is_runtype;
  int alp_bet;      int is_alp_bet;        // coupling between alpha and beta chanels, 1 - yes, 0 - no
  int decoherence;  int is_decoherence;    // choose the decoherence method to use; 0 - no decoherence
  int regress_mode; int is_regress_mode;   // regression mode used during dephasing times calculations

  // Electromagnetic field
  int is_field;             int is_is_field;       // flag to include explicit field
  std::string field_dir;    int is_field_dir;      // direction of the field
  int field_protocol;       int is_field_protocol; // way the photoexcitation is prepared
  double field_Tm;          int is_field_Tm;       // middle of the excitation period
  double field_T;           int is_field_T;        // excitation period (e.g. duration of the laser pulse)
  double field_freq;        int is_field_freq;     // field frequency 
  std::string field_freq_units; int is_field_freq_units; // units of the excitation frequency
  double field_fluence;     int is_field_fluence;  // fluence of the field in mJ/cm^2


  // Constructor
  InputStructure(boost::python::dict);

};


#endif  // InputStructure_H

