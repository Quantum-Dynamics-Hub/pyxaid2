/***********************************************************
 * Copyright (C) 2013 Alexey V. Akimov
 * This file is distributed under the terms of the
 * GNU General Public License as published by the
 * Free Software Foundation; either version 3 of the
 * License, or (at your option) any later version.
 * http://www.gnu.org/copyleft/gpl.txt
***********************************************************/

#ifndef IO_H
#define IO_H

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

int read_file(std::string filename,int verbose,vector<std::string>& A);

void file2matrix(std::string filename,vector< vector<double> >& m);
void file2matrix(std::string filename,vector< vector<double> >& m,double scl);
void file2matrix(std::string filename,vector< vector<int> >& m);

void show_2D(vector< vector<double> >& in);


#endif // IO_H
