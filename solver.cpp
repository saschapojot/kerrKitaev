//
// Created by tusa on 16/9/20.
//
#include "solver.hpp"

solver::solver() {
    //init simpTab
    std::vector<std::complex<double>>v0(this->CON.Q,std::complex<double>(0,0));
    for (const int &k : this->CON.kIndAll){
        this->simpTab.push_back(v0);
    }

    //init thetaTab
    std::vector<std::complex<double>>v1(this->CON.Q+1,std::complex<double>(0,0));
    for(const int&k:this->CON.kIndAll){
        this->thetaDTab.push_back(v1);
    }




}

double solver::b0(const int &k) {
    return this->CON.d0 * std::sin((double) k * this->CON.dk);
}

double solver::c0(const int &k) {
    return -this->CON.t0 * std::cos((double) k * this->CON.dk);

}

///
/// \param k
/// \return init vector
Eigen::Vector2cd solver::initVec(const int &k) {

    Eigen::Vector2cd rst;
    double b0Val = this->b0(k);
    double c0Val = this->c0(k);

    double denom2 = 2 * b0Val * b0Val + 2 * c0Val * c0Val + 2 * c0Val * std::sqrt(b0Val * b0Val + c0Val * c0Val);
    if (fabs(denom2) >= this->CON.tol) {
        double denom = std::sqrt(denom2);
        std::complex<double> denomZ{denom, 0};

        std::complex<double> numer1{0, -b0Val};
        std::complex<double> numer2{c0Val + std::sqrt(c0Val * c0Val + b0Val * b0Val), 0};
        rst[0] = numer1 / denomZ;
        rst[1] = numer2 / denomZ;
        return rst;


    }
    std::complex<double> c0 = 1j;

    rst[0] = 0.0 * c0;
    rst[1] = 1.0 + 0.0 * c0;
    return rst;

}

///
/// \param k
/// \param con1
/// \return Linear part of the Hamiltonian after quench
Eigen::Matrix2cd solver::H0(const int &k) {
    Eigen::Matrix2cd rst;
    std::complex<double> idc = 1.0 + 0.0 * 1j;

    rst(0, 0) = (-this->CON.mu1 - this->CON.t1 * std::cos((double) k * this->CON.dk)) * idc;

    std::complex<double> tmp01{0, this->CON.d1 * std::sin(this->CON.dk * (double) k)};
    rst(0, 1) = tmp01;

    rst(1, 0) = -tmp01;
    rst(1, 1) = (this->CON.t1 * std::cos((double) k * this->CON.dk)) * idc;

    return rst;


}

Eigen::Matrix2cd solver::expH0(const int &k) {
    Eigen::Matrix2cd h0Val = this->H0(k);
    Eigen::Matrix2cd z0 = -1.0j * this->CON.dt / 2.0 * h0Val;
    return z0.exp();

}

///
/// \param k
/// \param vecStart
/// \return state vector after 1 step of 2nd order Strang splitting S2
Eigen::Vector2cd solver::S2(const int &k, const Eigen::Vector2cd &vecStart) {
//step 1
    Eigen::Matrix2cd exph0Val = this->expH0(k);
    Eigen::Vector2cd vec1 = exph0Val * vecStart;

//step 2
    std::complex<double> vkm = vec1[0];
    std::complex<double> wkm = vec1[1];
    //calculate etakm
    std::complex<double> tmp1 = -1.0j * this->CON.lmd * std::pow(std::abs(vkm), 2) * this->CON.dt;
    std::complex<double> etakm = vkm * std::exp(tmp1);
    // calculate zetakm
    std::complex<double> tmp2 = -1.0j * this->CON.lmd * std::pow(std::abs(wkm), 2) * this->CON.dt;
    std::complex<double> zetakm = wkm * std::exp(tmp2);

    //step 3
    Eigen::Vector2cd vec2;
    vec2 << etakm, zetakm;
    return exph0Val * vec2;


}

/// calculate the state vectors starting with the kth momentum
/// \param k
void solver::calulateVec(const int &k) {


    //initialization
    Eigen::Vector2cd veck0 = this->initVec(k);

    this->solutionAll[k].push_back(veck0);
    //time step number: 0,1,...,Q*R-1
    for (int m = 0; m < this->CON.Q * this->CON.R; m++) {
        auto vecCurr = this->solutionAll[k].back();
        auto vecNext = this->S2(k, vecCurr);
        this->solutionAll[k].push_back(vecNext);

    }


}

void solver::writeAllVects() {
    std::vector<std::thread> thrds;
    for (const int &k : this->CON.kIndAll) {
        thrds.emplace_back(&solver::calulateVec, this, k);
    }
    for (auto &th:thrds) {
        th.join();
    }

}

std::complex<double> solver::Jkab(const int &k, const int &a, const int &b) {
    Eigen::Vector2cd vec = this->solutionAll[k][a * this->CON.R + b];
    std::complex<double> yk = vec[0];
    std::complex<double> zk = vec[1];

    Eigen::Matrix2cd H1;
    H1(0, 0) = -this->CON.mu1 - this->CON.t1 * std::cos(this->CON.dk * (double) k) +
               this->CON.lmd * std::pow(std::abs(yk), 2);

    H1(0, 1) = std::complex<double>(0, this->CON.d1 * std::sin(this->CON.dk * (double) k));
    H1(1, 0) = -H1(0, 1);
    H1(1, 1) = this->CON.t1 * std::cos(this->CON.dk * (double) k) + this->CON.lmd * std::pow(std::abs(zk), 2);

    std::complex<double> rst;
    rst = vec.adjoint() * H1 * vec;
    rst /= (std::pow(std::abs(yk), 2) + std::pow(std::abs(zk), 2));
    return rst;
}

///
/// \param k
/// \return simpson integration
std::complex<double> solver::simpsonD(const int &k, const int &a) {
    //a=0,1,...,Q-1
    std::complex<double> evenSum, oddSum;

    //compute odd sums
    for (int b = 1; b < this->CON.R; b++) {
        oddSum += this->Jkab(k, a, b);
    }
    //compute even sums
    for (int b = 2; b < this->CON.R; b++) {
        evenSum += this->Jkab(k, a, b);
    }

    std::complex<double> rst = this->Jkab(k, a, 0) + 4.0 * oddSum + 2.0 * evenSum + this->Jkab(k, a, this->CON.R);
    rst*=this->CON.dt/3.0;
    return rst;
}

void solver::writeSimpTabOneEntry(const int &k, const int &a) {
    this->simpTab[k][a]=this->simpsonD(k,a);

}

void solver::writeSimpTabAllEntries() {
    //k=0,1,...,N;
    //a=0,1,...,Q-1;
    std::vector<std::thread> thrdsAll;
    for(const int&k :this->CON.kIndAll){
        for(int a=0;a<this->CON.Q;a++){
            thrdsAll.emplace_back(&solver::writeSimpTabOneEntry,this,k,a);
        }
    }
    for(auto &th:thrdsAll){
        th.join();
    }


}
void solver::writeThetaDTabOneEntry(const int&k, const int&q) {
   //q=1,2,...,Q


}