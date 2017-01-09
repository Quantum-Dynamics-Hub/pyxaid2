/***********************************************************
 * Copyright (C) 2013 Alexey V. Akimov
 * This file is distributed under the terms of the
 * GNU General Public License as published by the
 * Free Software Foundation; either version 3 of the
 * License, or (at your option) any later version.
 * http://www.gnu.org/copyleft/gpl.txt
***********************************************************/

#include "io.h"
#include <cstdlib>
#include <cstdio>
using namespace std;


int read_file(std::string filename,int verbose,vector<std::string>& A){
/**********************************************************************
  This function reads file <filename> and stores it as a vector of strings
  each string is a line of the file
**********************************************************************/

  // Prepare A
  if(A.size()>0){ A.clear(); }

  // Read the file
  if(verbose==1){ cout<<"Reading file"<<filename<<endl; }
  ifstream f;
  f.open(filename.c_str(),ios::in);
  if(f.is_open()){
    // Estimate the number of lines
    int nlines = 0;
    while(!f.eof()){ std::string line; getline(f,line); nlines++; }
    A.reserve(int(nlines*1.25));
    // Reset file position to beginning
    f.clear();
    f.seekg (0, ios::beg);

    // Actually put the file content in memory
    while(!f.eof()){ std::string line; getline(f,line); A.push_back(line);  }
  }else{ cout<<"Error: Can not open file "<<filename<<endl; }
  f.close();

  return A.size();
}


void file2matrix(std::string filename,vector< vector<double> >& M){
/*****************************************************************
  This function reads the content of the tabular (2D) file into matrix M
*****************************************************************/
  if(M.size()>0){ M.clear(); }
  ifstream in;
  in.open(filename.c_str(),ios::in);
  if(in.is_open()){
    std::string s;
    while(!in.eof()){
      vector<double> line; // line of the file
      getline(in,s);
      stringstream ss(s,stringstream::in|stringstream::out);
      while(ss>>s){  line.push_back(atof(s.c_str()));   }
      if(line.size()>0) { M.push_back(line); }
    }
  }else{ cout<<"Error: Can not open file "<<filename<<". Check if this file exists\n"; exit(0);}
  in.close();
}

void file2matrix(std::string filename,vector< vector<double> >& M,double scl){
/*****************************************************************
  This function reads the contend of the tabular (2D) file into matrix M and scales
  the read data by factor scl
*****************************************************************/
  if(M.size()>0){ M.clear(); }
  ifstream in;
  in.open(filename.c_str(),ios::in);
  if(in.is_open()){
    std::string s;
    while(!in.eof()){
      vector<double> line; // line of the file
      getline(in,s);
      stringstream ss(s,stringstream::in|stringstream::out);
      while(ss>>s){  line.push_back(scl*atof(s.c_str()));   }
      if(line.size()>0) { M.push_back(line); }
    }
  }else{ cout<<"Error: Can not open file "<<filename<<". Check if this file exists\n"; exit(0);}
  in.close();
}


void file2matrix(std::string filename,vector< vector<int> >& M){
/*****************************************************************
  This function reads the content of the tabular (2D) file into matrix M
  Version overloaded for int
*****************************************************************/
  if(M.size()>0){ M.clear(); }
  ifstream in;
  in.open(filename.c_str(),ios::in);
  if(in.is_open()){
    std::string s;
    while(!in.eof()){
      vector<int> line; // line of the file
      getline(in,s);
      stringstream ss(s,stringstream::in|stringstream::out);
      while(ss>>s){  line.push_back(atoi(s.c_str()));   }
      if(line.size()>0) { M.push_back(line);     }
    }
  }else{ cout<<"Error: Can not open file "<<filename<<". Check if this file exists\n"; exit(0);}
  in.close();
}

void show_2D(vector< vector<double> >& in){
/******************************************************************
  This function prints out the matrix in a tabular form
******************************************************************/
  for(int i=0;i<in.size();i++){
    for(int j=0;j<in[i].size();j++){
      cout<<"in["<<i<<"]["<<j<<"]="<<in[i][j]<<" ";
    }
    cout<<endl;
  }
}

