//
// Created by tusa on 15/9/20.
//

#ifndef KERRKITAEV_SOLVER_HPP
#define KERRKITAEV_SOLVER_HPP

#include <eigen3/Eigen/Dense>
#include <map>
#include "consts.hpp"
#include <complex>
class solver{

public:
    /*
     * Functions to initialize*/
    double b0(const int &k,const CONSTS& con1);
    double c0(const int &k, const CONSTS& con1);
    Eigen::Vector2cd initVec(const int &k, const CONSTS& con1);


    Eigen::Matrix2cd H0(const int& k, const CONSTS& con1);//linear part of the Hamiltonian
    Eigen::Matrix2cd expH0(const Eigen::Matrix2cd & h0val);//exponentiates the H0 matrix in S2
    Eigen::Vector2cd S2(const int& k, const Eigen::Vector2cd& vecStart, const CONSTS & con1);//one step S2
    std::vector<double> phaseD();



public:

    /*
     * Dictionary to hold all solutions for each k, for each small time step mdt
     * */
    std::map<int,std::vector<Eigen::Vector2cd>> solutionAll;
    std::vector<double>tG;//geometric phase at each qds
    std::vector<double> beta;//geometric phase increment

    std::complex<double> integrandD(const int &k);
    std::complex<double> simpsonD(const int &k);


};


#endif //KERRKITAEV_SOLVER_HPP
