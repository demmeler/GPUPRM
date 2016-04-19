#ifndef TICTOC_H
#define TICTOC_H


#include <chrono>
#include <cstdlib>
#include <ctime>

#include "config.h"

#define tick(t1) std::chrono::high_resolution_clock::time_point t1=std::chrono::high_resolution_clock::now();
#define tock(t1) {std::chrono::high_resolution_clock::time_point endtime___=std::chrono::high_resolution_clock::now(); \
                 int time=std::chrono::duration_cast<std::chrono::microseconds>(endtime___-t1).count();                \
                 std::cout << "Since " << #t1 << ": \t" << time << "mcs" << std::endl; }
#define tockval(t1,time__) {std::chrono::high_resolution_clock::time_point endtime___=std::chrono::high_resolution_clock::now(); \
                 time__=std::chrono::duration_cast<std::chrono::microseconds>(endtime___-t1).count();                \
                 std::cout << "Since " << #t1 << ": \t" << time__ << "mcs" << std::endl; }

#define timedelta(t1,t2) {int time=std::chrono::duration_cast<std::chrono::microseconds>(endtime___-t1).count(); \
                          std::cout << "Since " << #t1 << ": \t" << time << "mcs" << std::endl;}




#endif // TICTOC_H
