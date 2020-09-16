//
// Created by tusa on 15/9/20.
//

#ifndef KERRKITAEV_SOLVER_HPP
#define KERRKITAEV_SOLVER_HPP

#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include <map>
#include "consts.hpp"
#include <complex>
#include <thread>
#include <cmath>
#include <utility>
class solver{
public:
    solver();

public:
    /*
     * Functions to initialize*/
    double b0(const int &k);
    double c0(const int &k);
    Eigen::Vector2cd initVec(const int &k);




    Eigen::Matrix2cd H0(const int& k);//linear part of the Hamiltonian
    Eigen::Matrix2cd expH0(const int &k);//exponentiates the H0 matrix in S2
    Eigen::Vector2cd S2(const int& k, const Eigen::Vector2cd& vecStart);//one step S2
    void calulateVec(const int &k);// calculate the state vectors starting with the kth momentum, write to dict.
    void writeAllVects();//calculate all states with multithreading
   //void toJson();
   // write the state vectors to json
    std::vector<double> phaseD();
    void writeSimpTabOneEntry(const int &k, const int &q);
    void writeSimpTabAllEntries();

    void writeThetaDTabOneEntry(const int&k,const int&q);
    std::complex<double> Jkab(const int&k, const int &a, const int &b);//integrand


public:

    /*
     * Dictionary to hold all solutions for each k, for each small time step mdt
     * */
    CONSTS CON;
    std::map<int,std::vector<Eigen::Vector2cd>> solutionAll;
    std::vector<double>tG;//geometric phase at each qds
    std::vector<double> beta;//geometric phase increment

    std::complex<double> simpsonD(const int &k, const int &a);
   std::vector<std::vector<std::complex<double>>>simpTab;//k=0,1,...,N;a=0,1,...,Q-1
   std::vector<std::vector<std::complex<double>>>thetaDTab;//k=0,1,...,N;q=0,1,...,Q;
};


#endif //KERRKITAEV_SOLVER_HPP
