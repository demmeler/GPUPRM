#ifndef UTIL_HPP
#define UTIL_HPP

#include <iostream>

#define check(x) if(!(x)){ std::cout<<"check "<<#x<<" failed"<<std::endl;}
#define printvar(x) std::cout<<#x<<"="<<x<<std::endl;
#define printarr(x,n) std::cout<<#x<<"=";for(int i=0;i<n;++i){std::cout<<x[i]<<" ";}std::cout<<std::endl;
#define msg(x) std::cout<<x<<std::endl;

#endif // UTIL_HPP
