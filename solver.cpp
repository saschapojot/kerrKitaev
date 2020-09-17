//
// Created by tusa on 16/9/20.
//
#include "solver.hpp"

solver::solver() {
    //init solutionsAll
    std::vector<Eigen::Vector2cd> vecs0(this->CON.Q*this->CON.R+1);
    for(const auto&k:this->CON.kIndAll) {
        this->solutionsAll.push_back(vecs0);
    }
    //init simpTab
    std::vector<std::complex<double>> v0(this->CON.Q, std::complex<double>(0, 0));
    for (const int &k : this->CON.kIndAll) {
        this->simpTab.push_back(v0);
    }

    //init thetaDTab
    std::vector<std::complex<double>> v1(this->CON.Q + 1, std::complex<double>(0, 0));
    for (const int &k:this->CON.kIndAll) {
        this->thetaDTab.push_back(v1);
    }

    //init thetaTotTab
    std::vector<std::complex<double>> v2(this->CON.Q+1, std::complex<double>(0,0));
    for(const auto&k:this->CON.kIndAll){
        this->thetaTotTab.push_back(v2);
    }

    //init thetaGTab
    std::vector<std::complex<double>>v3(this->CON.Q+1,std::complex<double>(0,0));
    for(const auto&k:this->CON.kIndAll){
        this->thetaGTab.push_back(v3);
    }

    //init beta
    std::vector<double> b0(this->CON.N,0.0);
    for(int q=0; q<this->CON.Q+1;q++){
        this->beta.push_back(b0);
    }

    //init W
    this->W=std::vector<double>(this->CON.Q+1,0.0);


}

double solver::b0(const int &k) {
    return this->CON.d0 * std::sin((double) k * this->CON.dk);
}

double solver::c0(const int &k) {
    return -this->CON.t0 * std::cos((double) k * this->CON.dk)-this->CON.mu0/2.0;

}

///
/// \param k
/// \return init vector
Eigen::Vector2cd solver::initVec(const int &k) {

    Eigen::Vector2cd rst;
    double b0Val = this->b0(k);
    double c0Val = this->c0(k);

    double denom2 = 2 * b0Val * b0Val + 2 * c0Val * c0Val + 2 * c0Val * std::sqrt(b0Val * b0Val + c0Val * c0Val);
    if (std::abs(denom2) >= this->CON.tol) {
        double denom = std::sqrt(denom2);
        std::complex<double> denomZ{denom, 0};

        std::complex<double> numer1{0, -b0Val};
        std::complex<double> numer2{c0Val + std::sqrt(c0Val * c0Val + b0Val * b0Val), 0};
        rst[0] = numer1 / denomZ;
        rst[1] = numer2 / denomZ;
        return rst;


    }
    else {
        std::complex<double> c0 = 1j;

        rst[0] = 0.0 * c0;
        rst[1] = 1.0 + 0.0 * c0;
        return rst;
    }

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

    rst(1, 0) = std::conj(tmp01);
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

    this->solutionsAll[k][0]=veck0;
    //time step number: 0,1,...,Q*R-1
    for (int m = 0; m < this->CON.Q * this->CON.R; m++) {
        auto vecCurr = this->solutionsAll[k][m];
        auto vecNext = this->S2(k, vecCurr);
        this->solutionsAll[k][m+1]=vecNext;

        //solutionAll[k] has Q*R vectors

    }


}


//void solver::writeAllVects() {
//    std::vector<std::thread> thrds;
//    for (const int &k : this->CON.kIndAll) {
//        thrds.emplace_back(&solver::calulateVec, this, k);
//    }
//    for (auto &th:thrds) {
//        th.join();
//    }
//
//}


void solver::writeAllVects() {
    int numPtr = -1;
    do {
        std::vector<std::thread> thrdsAll;
        int numPtrTmp = numPtr;
        for (int ki = numPtrTmp + 1; ki < numPtrTmp + this->CON.threadNum + 1 && ki < this->CON.N + 1; ki++) {
            thrdsAll.emplace_back(&solver::calulateVec, this, ki);
            numPtr += 1;
        }
        for (auto &th:thrdsAll) {
            th.join();
        }
    } while (numPtr < this->CON.N);
}

std::complex<double> solver::Jkab(const int &k, const int &a, const int &b) {
    //a=0,1,...,Q-1;b=0,1,...,R;

    Eigen::Vector2cd vec = this->solutionsAll[k][a * this->CON.R + b];
    std::complex<double> yk = vec[0];
    std::complex<double> zk = vec[1];

    Eigen::Matrix2cd H1;
    H1(0, 0) = -this->CON.mu1 - this->CON.t1 * std::cos(this->CON.dk * (double) k) +
               this->CON.lmd * std::pow(std::abs(yk), 2);

    H1(0, 1) = std::complex<double>(0, this->CON.d1 * std::sin(this->CON.dk * (double) k));
    H1(1, 0) = std::conj(H1(0, 1));
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
    //integration over [ads, (a+1)ds]
    std::complex<double> evenSum, oddSum;

    //compute odd sums
    for (int b = 1; b < this->CON.R; b += 2) {
        oddSum += this->Jkab(k, a, b);
    }
    //compute even sums
    for (int b = 2; b < this->CON.R; b += 2) {
        evenSum += this->Jkab(k, a, b);
    }

    std::complex<double> rst = this->Jkab(k, a, 0) + 4.0 * oddSum + 2.0 * evenSum + this->Jkab(k, a, this->CON.R);
    rst *= this->CON.dt / 3.0;
    return rst;
}

void solver::writeSimpTabOneEntry(const int &k, const int &a) {
    this->simpTab[k][a] = this->simpsonD(k, a);

}

//void solver::writeSimpTabAllEntries() {
//    //k=0,1,...,N;
//    //a=0,1,...,Q-1;
//    std::vector<std::thread> thrdsAll;
//    for (const int &k :this->CON.kIndAll) {
//        for (int a = 0; a < this->CON.Q; a++) {
//            thrdsAll.emplace_back(&solver::writeSimpTabOneEntry, this, k, a);
//        }
//    }
//    for (auto &th:thrdsAll) {
//        th.join();
//    }
//
//
//}
void solver::writeSimpTabAllEntries() {
    //k=0,1,...,N;
    //a=0,1,...,Q-1;
    for(int a=0;a<this->CON.Q;a++) {
        int numPtr = -1;
        do {
            std::vector<std::thread> thrds;
            int numPtrTmp = numPtr;
            for (int ki = numPtrTmp + 1; ki < numPtrTmp + this->CON.threadNum + 1 && ki < this->CON.N + 1; ki++) {
                thrds.emplace_back(&solver::writeSimpTabOneEntry, this, ki, a);
                numPtr += 1;
            }
            for (auto &th:thrds) {
                th.join();
            }
        } while (numPtr < this->CON.N);
    }




}
void solver::writeThetaDTabOneEntry(const int &k, const int &q) {
    //q=1,2,...,Q
    std::complex<double> sumSimpTabTmp{0, 0};
    for (int a = 0; a < q; a++) {
        sumSimpTabTmp += this->simpTab[k][a];
    }
    sumSimpTabTmp *= -1.0;

    Eigen::Vector2cd vecqR = this->solutionsAll[k][q * this->CON.R];
    Eigen::Vector2cd vec0 = this->solutionsAll[k][0];


    std::complex<double> tmp2N = (vecqR.adjoint() * vecqR);
    std::complex<double> tmp2D = vec0.adjoint() * vec0;
    std::complex<double> term2 = std::complex<double>(0, 0.5) * std::log(tmp2N / tmp2D);
    this->thetaDTab[k][q] = sumSimpTabTmp + term2;


}

//void solver::writeThetaDTabAllEntries() {
//    //for q=0;
//    for (const auto &k:this->CON.kIndAll) {
//        this->thetaDTab[k][0] = 0;
//    }
//    std::vector<std::thread> thrdsAll;
//    //for q=1,2,...,Q
//
//    for (const auto &k:this->CON.kIndAll) {
//        for (int q = 1; q < this->CON.Q + 1; q++) {
//            thrdsAll.emplace_back(&solver::writeThetaDTabOneEntry, this, k, q);
//        }
//    }
//    for (auto &th:thrdsAll) {
//        th.join();
//    }
//}


void solver::writeThetaDTabAllEntries() {
    //for q=0;
    for (const auto &k:this->CON.kIndAll) {
        this->thetaDTab[k][0] = 0;
    }

    for (int q = 1; q < this->CON.Q + 1; q++) {
        int numPtr = -1;
        do {
            std::vector<std::thread> thrds;
            int numPtrTmp = numPtr;
            for (int ki = numPtrTmp + 1; ki < numPtrTmp + this->CON.threadNum + 1 && ki < this->CON.N + 1; ki++) {
                thrds.emplace_back(&solver::writeThetaDTabOneEntry, this, ki, q);
                numPtr += 1;
            }
            for (auto &th:thrds) {
                th.join();
            }
        } while (numPtr < this->CON.N);
    }

}
void solver::writeThetaTotTabOneEntry(const int &k, const int &q) {
    //k=0,1,...,N;
    //q=0,1,...,Q
    Eigen::Vector2cd vecqR=this->solutionsAll[k][q*this->CON.R];
    Eigen::Vector2cd vec0=this->solutionsAll[k][0];

    std::complex<double> tmp=vec0.adjoint()*vecqR;

   std::complex<double> rst=std::complex<double>(0,-1)*std::log(tmp/std::abs(tmp));
   this->thetaTotTab[k][q]=rst;


}

//void solver::writeThetaTotTabAllEntries() {
//    //k=0,1,..,N;
//    //q=0,1,...,Q;
//    std::vector<std::thread>thrdsAll;
//    for(const auto&k:this->CON.kIndAll){
//        for(int q=0;q<this->CON.Q+1;q++){
//            thrdsAll.emplace_back(&solver::writeThetaTotTabOneEntry,this,k,q);
//        }
//    }
//
//    for(auto &th:thrdsAll){
//        th.join();
//    }
//}

void solver::writeThetaTotTabAllEntries() {
    //k=0,1,...,N;
    //q=0,1,...,Q;
    for (int q = 0; q < this->CON.Q + 1; q++) {
        int numPtr = -1;
        do {
            std::vector<std::thread> thrds;
            int numPtrTmp = numPtr;
            for (int ki = numPtrTmp + 1; ki < numPtrTmp + this->CON.threadNum + 1 && ki < this->CON.N + 1; ki++) {
                thrds.emplace_back(&solver::writeThetaTotTabOneEntry, this, ki, q);
                numPtr += 1;
            }
            for (auto &th:thrds) {
                th.join();
            }
        } while (numPtr < this->CON.N);
    }
}

void solver::writeThetaGTabOneEntry(const int &k, const int &q) {
    //k=0,1,...,N;
    //q=0,1,...,Q;
    std::complex<double> diffTmp = this->thetaTotTab[k][q] - this->thetaDTab[k][q];

    this->thetaGTab[k][q] = diffTmp;
}

//void solver::writeThetaGTabAllEntries() {
//    //k=0,1,...,N;
//    //q=0,1,...,Q;
//    std::vector<std::thread> thrdsAll;
//    for (const auto &k:this->CON.kIndAll) {
//        for (int q = 0; q < this->CON.Q + 1; q++) {
//            thrdsAll.emplace_back(&solver::writeThetaGTabOneEntry, this, k, q);
//        }
//    }
//
//    for (auto &th:thrdsAll) {
//        th.join();
//    }
//
//}
void solver::writeThetaGTabAllEntries() {
    //k=0,1,...,N;
    //q=0,1,...,Q;
    for(const auto &k:this->CON.kIndAll){
        for(int q=0;q<this->CON.Q+1;q++){
            this->writeThetaGTabOneEntry(k,q);

        }
    }
}
double solver::jumpDecision(const double &incr) {
    double tmp = incr / M_PI;

    if (tmp >= this->CON.cutOff) {
        return incr - 2 * M_PI;
    } else if (tmp <= -this->CON.cutOff) {
        return incr + 2 * M_PI;
    } else {
        return incr;
    }
}

void solver::writeBetaOneEntry(const int &k, const int &q) {
    //q=0,1,...,Q;
    //k=0,1,...,N-1;
    double incr=(this->thetaGTab[k+1][q]-this->thetaGTab[k][q]).real();
    double bVal=this->jumpDecision(incr);
    this->beta[q][k]=bVal;

}
//void solver::writeBetaAllEntries() {
//    //q=0,1,...,Q;
//    //k=0,1,...,N-1;
//    std::vector<std::thread> thrdsAll;
//    for (int q = 0; q < this->CON.Q + 1; q++) {
//        for (int k = 0; k < this->CON.N; k++) {
//
//            thrdsAll.emplace_back(&solver::writeBetaOneEntry, this, k, q);
//        }
//    }
//    for (auto &th: thrdsAll) {
//        th.join();
//    }
//
//
//}
void solver::writeBetaAllEntries() {
    //q=0,1,...,Q;
    //k=0,1,...,N-1;
    for(int q=0;q<this->CON.Q+1;q++){
        for(int k=0;k<this->CON.N;k++){
            this->writeBetaOneEntry(k,q);


        }
    }

}
void solver::writeW() {
    //q=0,1,...,Q;
    for (int q = 0; q < this->CON.Q + 1; q++) {
        double wTmp = 0.0;
        for (int k = 0; k < this->CON.N / 2 ; k++) {
            wTmp += this->beta[q][k];
        }
        wTmp /= 2 * M_PI;
        W[q] = wTmp;
    }


}

void solver::plotW() {
    //q=0,1,..,Q;
    std::vector<double>Qtime;
    for(int q=0;q<this->CON.Q+1;q++){
        Qtime.push_back((double)q*this->CON.ds);
    }
    plt::plot(Qtime,this->W);
    plt::save("wn.png");
}

void solver::runCalc() {
    //0th
    this->writeAllVects();
    //1st
    this->writeSimpTabAllEntries();
    //2nd
    this->writeThetaDTabAllEntries();
    //3rd
    this->writeThetaTotTabAllEntries();
    //4th
    this->writeThetaGTabAllEntries();
    //5th
    this->writeBetaAllEntries();
    //6th
    this->writeW();
    //7th
    this->plotW();


}