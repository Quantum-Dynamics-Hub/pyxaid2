/***********************************************************
 * Copyright (C) 2013 Alexey V. Akimov
 * This file is distributed under the terms of the
 * GNU General Public License as published by the
 * Free Software Foundation; either version 3 of the
 * License, or (at your option) any later version.
 * http://www.gnu.org/copyleft/gpl.txt
***********************************************************/

#ifndef NAMD_H
#define NAMD_H

#include "aux.h"
#include "InputStructure.h"
#include "ElectronicStructure.h"


void hop(vector<double>& sh_prob,int& state,int num_states);
void solve_electronic(InputStructure& is,vector<ElectronicStructure>& es,matrix&);

void run_decoherence_rates(InputStructure& is, vector<ElectronicStructure>& me_es,vector<me_state>& me_states, int icond);
void run_namd(InputStructure& is, vector<ElectronicStructure>& me_es,vector<me_state>& me_states,int icond);
void run_namd1(InputStructure& is, vector<ElectronicStructure>& me_es,vector<me_state>& me_states, int icond);


#endif // NAMD_H
