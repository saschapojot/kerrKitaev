#include <iostream>
#include "solver.hpp"
#include <chrono>

int main() {

//    CONSTS co1;
//    solver sol0=solver(co1);
//    auto start=std::chrono::high_resolution_clock::now();
//    sol0.runCalc();
//    auto stop=std::chrono::high_resolution_clock::now();
//    auto duration=std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
//    std::cout<<duration.count()/1e6<<" seconds"<<std::endl;

//////////////////////////////////////////////////////////////////////////////////////
    int num=17;
     std::vector<CONSTS> consVec;
     float dlmd=1;
     for(int i=-num;i<num;i++){
         CONSTS conTmp;
         conTmp.lmd=(double)dlmd*i;
         consVec.push_back(conTmp);
     }

    for(auto &cParam:consVec){
        solver sol0=solver(cParam);
        auto start=std::chrono::high_resolution_clock::now();
        sol0.runCalc();
        auto stop=std::chrono::high_resolution_clock::now();
        auto duration=std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
        std::cout<<duration.count()/1e6<<" seconds"<<std::endl;
    }




    return 0;
}
