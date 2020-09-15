//
// Created by tusa on 15/9/20.
//

#ifndef KERRKITAEV_CONSTS_HPP
#define KERRKITAEV_CONSTS_HPP

#include <boost/math/constants/constants.hpp>
#include <vector>
const int PI = boost::math::constants::pi<double>();

class CONSTS {

public:
    CONSTS(){

        for(int i=0;i<this->N;i++){
            kIndAll.push_back(i);
        }
        for(int i=0;i< this->N/2;i++){
            kIndHalf.push_back(i);
        }
    }
public:


    int N = 500;//fermion number
    double Nd = (double) N;//N in double, to facilitate floating point calculation
    double dk = 2 * PI / Nd;

    /*
     * indices of momentum space
     * */
    std::vector<int> kIndAll;
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
    double u1 = 0.5;
    double t1 = 1.0;
    double d1 = 1.0;

    double lmd = 0.3;//nonlinearity strength
    /*
     * consts for time*/

    int R=20;//small time step number
    int Q=100;//large time step number
    double dt=0.01;//small time step
    double ds=(double)R*dt;//large time step


};

#endif //KERRKITAEV_CONSTS_HPP
