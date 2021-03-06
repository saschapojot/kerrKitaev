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
#include <fstream>
#include <cmath>
#include <utility>
#include "matplotlib-cpp/matplotlibcpp.h"
#include <memory>
namespace plt=matplotlibcpp;
class solver{
public:
    solver(const CONSTS&);

public:
    /*
     * Functions to initialize*/
    double b0(const int &k);
    double c0(const int &k);
    Eigen::Vector2cd initVec(const int &k);




    Eigen::Matrix2cd H0(const int& k);//linear part of the Hamiltonian
    Eigen::Matrix2cd expH0(const int &k,const double&stept);//exponentiates the H0 matrix in S2
    Eigen::Vector2cd S2(const int& k, const Eigen::Vector2cd& vecStart, const double &stept);//one step S2
    Eigen::Vector2cd S4(const int&k, const Eigen::Vector2cd&vecStart, const double &stept);//one step S4
    void calulateVec(const int &k);// calculate the state vectors starting with the kth momentum.
   //0th
    void writeAllVects();//calculate all states with multithreading
   //void toJson();
   // write the state vectors to json
   //1st
    void writeSimpTabOneEntry(const int &k, const int &q);
    void writeSimpTabAllEntries();
    //2nd
    void writeThetaDTabOneEntry(const int&k,const int&q);
    void writeThetaDTabAllEntries();
    //3rd
    void writeThetaTotTabOneEntry(const int&k, const int &q);
    void writeThetaTotTabAllEntries();
    //4th
    void writeThetaGTabOneEntry(const int&k,const int&q);
    void writeThetaGTabAllEntries();

    std::complex<double> Jkab(const int&k, const int &a, const int &b);//integrand
    double jumpDecision(const double& incr);
    double jump(const double& incr,const double& avg);

    //5th

    void writeBetaOneEntry(const int &k,const int&q);
    void writeBetaAllEntries();
    void jumpAvg();
    void writeBeta00AllEntries();
    //6th
    void writeW();
    //7th
    void plotW();

    //8th
    double xi(const int &k, const int &q);
    //9th
   void writel(const int&q);
    //10th
    void writeRateFuncs();

    //11th
    void plotRateFunction();
    //12th
    void maxW();
    //13th
    void printThetaGTab();
    void printBeta00tab();
    void printThetaDTab();
    void runCalc();




public:

    /*
     * Dictionary to hold all solutions for each k, for each small time step mdt
     * */
    CONSTS CON;
  // std::map<int,std::vector<Eigen::Vector2cd>> solutionAll;//dict, k, [vec0, vec1,...,vec_{QR}]
    std::vector<std::vector<Eigen::Vector2cd>> solutionsAll;
    std::vector<double>tG;//geometric phase at each qds


    std::complex<double> simpsonD(const int &k, const int &a);
   std::vector<std::vector<std::complex<double>>>simpTab;//k=0,1,...,N-1;a=0,1,...,Q-1
   std::vector<std::vector<std::complex<double>>>thetaDTab;//k=0,1,...,N-1;q=0,1,...,Q;
   std::vector<std::vector<std::complex<double>>> thetaTotTab;//k=0,1,...,N-1;q=0,1,...,Q;
   std::vector<std::vector<std::complex<double>>>thetaGTab;//k=0,1,...,N-1;q=0,1,...,Q;
    std::vector<std::vector<double>> beta00;
   std::vector<std::vector<double>> beta;//q=0,1,...,Q;k=0,1,...,N-1;

   std::vector<double>W;//q=0,1,...,Q;
   std::vector<double> rateFunction;//q=0,1,...,Q

};


#endif //KERRKITAEV_SOLVER_HPP
