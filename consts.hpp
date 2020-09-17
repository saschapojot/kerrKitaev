//
// Created by tusa on 15/9/20.
//

#ifndef KERRKITAEV_CONSTS_HPP
#define KERRKITAEV_CONSTS_HPP

#include <boost/math/constants/constants.hpp>
#include <vector>
#include <complex>
#include <thread>


class CONSTS {

public:
    CONSTS() {

        for (int i = 0; i < this->N + 1; i++) {
            kIndAll.push_back(i);
        }
        /*for (int i = 0; i < this->N / 2; i++) {
            kIndHalf.push_back(i);
        }*/
    }

public:


    int N = 100;//fermion number
    double Nd = (double) N;//N in double, to facilitate floating point calculation
    double dk = 2 * M_PI / Nd;

    /*
     * indices of momentum space
     * */
    std::vector<int> kIndAll;//0,1,...,N, N+1 in total
// std::vector<int> kIndHalf;


    /*
     * Parameters before the quench
     * */
    double mu0 = -1.2;
    double t0 = 1.0;
    double d0 = 1.0;
    /*
     * Parameters after the quench
     * */
    double mu1 = mu0;
    double t1 = t0;
    double d1 = d0;

    double lmd = 1;//nonlinearity strength
    /*
     * consts for time*/

    int R = 20;//small time step number
    int Q = 500;//large time step number
    double dt = 0.01;//small time step
    double ds = (double) R * dt;//large time step


    double tol = 1e-16;
    double cutOff = 0.8;

    int threadNum= std::thread::hardware_concurrency();


};

#endif //KERRKITAEV_CONSTS_HPP
