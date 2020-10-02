//
// Created by tusa on 15/9/20.
//

#ifndef KERRKITAEV_CONSTS_HPP
#define KERRKITAEV_CONSTS_HPP

#include <boost/math/constants/constants.hpp>
#include <vector>
#include <complex>
#include <thread>
#include <string>

class CONSTS {

public:
    CONSTS() {

        for (int i = 0; i < this->N; i++) {
            kIndAll.push_back(i);
        }
        /*for (int i = 0; i < this->N / 2; i++) {
            kIndHalf.push_back(i);
        }*/


        int tmp1 = (int) std::thread::hardware_concurrency();
        int tmp2 = 12;
        this->threadNum = tmp1 > tmp2 ? tmp1 : tmp2;
    }

public:


    int N = 500;//fermion number
    double Nd = (double) N;//N in double, to facilitate floating point calculation
    double dk = 2 * M_PI / Nd;

    /*
     * indices of momentum space
     * */
    std::vector<int> kIndAll;//0,1,...,N-1, N in total



    /*
     * Parameters before the quench
     * */
    double mu0 = -6;
    double t0 = 1.0;
    double d0 = -1.0;
    /*
     * Parameters after the quench
     * */
    double mu1 = 0;
    double t1 = t0;
    double d1 = d0;

    double lmd = 0;//nonlinearity strength
    /*
     * consts for time*/

    int R = 160;//small time step number

    int Q = 1000;//large time step number
    double dt = 0.00025;//small time step
    double ds = (double) R * dt;//large time step


    double tol = 1e-16;
    double cutOff = 1.2;

    int threadNum;//= std::thread::hardware_concurrency();

    std::string dir="/home/disk2/Documents/cppCode/kerrKitaev/benchmark/";

    double timeAxisParts=10;

};

#endif //KERRKITAEV_CONSTS_HPP
