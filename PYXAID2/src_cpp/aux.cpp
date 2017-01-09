/***********************************************************
 * Copyright (C) 2013 Alexey V. Akimov
 * This file is distributed under the terms of the
 * GNU General Public License as published by the
 * Free Software Foundation; either version 3 of the
 * License, or (at your option) any later version.
 * http://www.gnu.org/copyleft/gpl.txt
***********************************************************/


#include "aux.h"

//---------------- Operations on vectors of integers (states) ---------------

void show_vector(vector<int>& A){
  int sz = A.size();
  for(int i=0;i<sz;i++){  cout<<A[i]<<"  "; }
}
int is_in_vector(int a, vector<int>& A){
  int res = 0;
  int sz = A.size();
  for(int i=0;i<sz;i++){ if(a==A[i]){ res=1; break;} }
  return res;
}

int num_in_vector(int a, vector<int>& A, vector<int>& indx){
// Returns how many times a is found in vector A
// indx will contain the positions at which a is found in A
  int res = 0;
  int sz = A.size();
  if(indx.size()>0){ indx.clear(); }
  for(int i=0;i<sz;i++){ if(a==A[i]){ res++; indx.push_back(i); } }
  return res;
}

int find_int(int a,vector<int>& A){
  int res = -1;
  int sz = A.size();
  for(int i=0;i<sz;i++){ if(a==A[i]){ res=i; break;} }
  return res;
}

int is_repeating(vector<int>& A,int& reap){
  int res = 0;
  // Find out if there are repeating elements in vector A
  int sz = A.size();
  for(int i=0;i<sz-1;i++){
    vector<int> tmp = std::vector<int>(A.begin()+i+1,A.end());
    if(is_in_vector(A[i],tmp)){ res = 1; reap = A[i]; break; }
  }
  return res;
}


//-------------------- Operations on strings ----------------------------

void split_line(std::string line, vector<std::string>& arr){

  stringstream ss(line,stringstream::in|stringstream::out);
  std::string s;
  while(ss>>s){ arr.push_back(s); }

}

void split_line2(std::string line,vector<std::string>& arr,char delim){

  std::istringstream f(line);
  std::string s;    
  while(std::getline(f, s, delim)){ arr.push_back(s); }

}

std::string int2str(int inp){
  stringstream ss(stringstream::in | stringstream::out);
  std::string out;
  (ss << inp);  ss >> out;
  return out;
}


int find_section(vector<std::string>& A,std::string marker_beg,std::string marker_end,int min_line,int max_line,int& beg,int& end){ 

  beg = end = -1;
  size_t found;
  int status = 0;

  for(int i=min_line;i<max_line;i++){
    if(beg==-1){
      found = A[i].find(marker_beg);
      if(found!=string::npos){  beg = i;}
    }
    if(end==-1){
      found = A[i].find(marker_end);
      if(found!=string::npos){  end = i;}
    }
    if((beg!=-1) && (end!=-1)){ 
      if(beg<end){status = 1;}
      break;
    }
  }// for i

  return status;
}

std::string extract_s(std::string line, std::string marker){
// line - is the line from which we want to extract some data
// marker - is the keyword before that data, format is:  marker="data"
  size_t pos,pos1,pos2;
  string res;

  pos = line.find(marker);
  if(pos!=string::npos){  // Marker is found
    pos1 = line.find("\"",pos+1);
    pos2 = line.find("\"",pos1+1);
    if(pos1!=string::npos && pos2!=string::npos){
      res = line.substr(pos1+1,pos2-pos1-1); // extract the value
    }
  }
  return res;
}

std::string int2string(int inp){
  stringstream ss(stringstream::in | stringstream::out);
  std::string out;
  (ss << inp);  ss >> out;
  return out;
}


//-------------------- Operations on arrays ----------------------------


void extract_2D(vector< vector<double> >& in, vector< vector<double> >& out, int minx,int maxx, int miny, int maxy ){
/*****************************************************************
  This function extracts 2D sub-array out from 2D array in
  The region is degined by indices: [minx,maxx]x[miny,maxy]
  Axes are: in[x][y], out[x][y]
******************************************************************/
  if(out.size()>0){ out.clear(); }
  for(int x=minx;x<=maxx;x++){
    vector<double> line = std::vector<double>(in[x].begin()+miny,in[x].begin()+maxy+1);
    out.push_back(line);
  }
}

void extract_2D(vector< vector<double> >& in, vector< vector<double> >& out, vector<int>& templ,int shift){
/*****************************************************************
  This function extracts 2D sub-array out from 2D array in
  The region is degined by indices: templ x templ
  temp - is a template for extraction
******************************************************************/
  if(out.size()>0){ out.clear(); }
  int sz = templ.size();
  for(int i=0;i<sz;i++){
    vector<double> line = std::vector<double>(sz,0.0);
    for(int j=0;j<sz;j++){  line[j] = in[templ[i]+shift][templ[j]+shift];   }
    out.push_back(line);
  }
}

void extract_1D(vector<double>& in, vector<double>& out, vector<int>& templ,int shift){
/*****************************************************************
  This function extracts 1D sub-array out from 1D array in
  The region is degined by indices: templ x templ
  temp - is a template for extraction
******************************************************************/
  if(out.size()>0){ out.clear(); }
  int sz = templ.size();
  out = std::vector<double>(sz,0.0);
  for(int i=0;i<sz;i++){  out[i] = in[templ[i]+shift];   }
}
 
