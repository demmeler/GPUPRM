#ifndef UTIL_HPP
#define UTIL_HPP

#include <iostream>
#include <fstream>
#include <string>

#define check(x) if(!(x)){ std::cout<<"check "<<#x<<" failed"<<std::endl;}
#define printvar(x) std::cout<<#x<<"="<<x<<std::endl;
#define printarr(x,n) std::cout<<#x<<"=";for(int i=0;i<n;++i){std::cout<<x[i]<<" ";}std::cout<<std::endl;
#define msg(x) std::cout<<x<<std::endl;


void write_file(std::string path, double* val, int n){
  std::ofstream outfile;
  outfile.open(path.c_str(), std::ios::out | std::ios::binary);
  outfile.write((char*)val, n*sizeof(double));
  outfile.close();
}

void write_file(std::string path, float* val, int n){
  std::ofstream outfile;
  outfile.open(path.c_str(), std::ios::out | std::ios::binary);
  outfile.write((char*)val, n*sizeof(float));
  outfile.close();
}

void write_file(std::string path, int* val, int n){
  std::ofstream outfile;
  outfile.open(path.c_str(), std::ios::out | std::ios::binary);
  outfile.write((char*)val, n*sizeof(int));
  outfile.close();
}

void read_file(std::string path, double* val, int n){
  std::ifstream infile;
  infile.open(path.c_str(), std::ios::in | std::ios::binary);
  infile.read((char*)val, n*sizeof(double));
  infile.close();
}

void read_file(std::string path, float* val, int n){
  std::ifstream infile;
  infile.open(path.c_str(), std::ios::in | std::ios::binary);
  infile.read((char*)val, n*sizeof(float));
  infile.close();
}

void read_file(std::string path, int* val, int n){
  std::ifstream infile;
  infile.open(path.c_str(), std::ios::in | std::ios::binary);
  infile.read((char*)val, n*sizeof(int));
  infile.close();
}






#endif // UTIL_HPP
