//
// Created by tusa on 15/9/20.
//

#ifndef KERRKITAEV_CONSTS_HPP
#define KERRKITAEV_CONSTS_HPP

#include <boost/math/constants/constants.hpp>
#include <vector>
#include <complex>



class CONSTS {

public:
    CONSTS(){

        for(int i=0;i<this->N+1;i++){
            kIndAll.push_back(i);
        }
        for(int i=0;i< this->N/2;i++){
            kIndHalf.push_back(i);
        }
    }
public:


    int N = 10;//fermion number
    double Nd = (double) N;//N in double, to facilitate floating point calculation
    double dk = 2 * M_PI / Nd;

    /*
     * indices of momentum space
     * */
    std::vector<int> kIndAll;//0,1,...,N, N+1 in total
    std::vector<int>kIndHalf;


    /*
     * Parameters before the quench
     * */
    double mu0 = -3.0;
    double t0 = 1.0;
    double d0 = 1.0;
    /*
     * Parameters after the quench
     * */
    double mu1 = 0.5;
    double t1 = 1.0;
    double d1 = 1.0;

    double lmd = 0.3;//nonlinearity strength
    /*
     * consts for time*/

    int R=20;//small time step number
    int Q=100;//large time step number
    double dt=0.01;//small time step
    double ds=(double)R*dt;//large time step


    double tol=1e-16;




};

#endif //KERRKITAEV_CONSTS_HPP
