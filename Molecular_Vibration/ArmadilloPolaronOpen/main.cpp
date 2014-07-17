#include <iostream>
#include <string>
#include <sstream>
#include "armaInclude/armadillo"
#include <fstream>
#include <time.h>
#include <math.h>
using namespace std;
using namespace arma;

//# lcn, lattice configuration number, which is a vec(200),
//# ham, Hamiltonian number, which is a {200,200} mat(200,200).
//# t0, jumping probability between two electrons for anhilation operators and creation operators.
//# te, bk terms, describes ...
//# alpha, electron-lattice coupling constants,
//# k, lattice elastic constant,
//# m, mass of molecues,
//# N, number of lattice,
//# A, lattice distance before dimerization,
//# Ity, Light beam intensity,
//# EEoE, Eigen Energy of Electron, which is {1,200}, vec(200)
//# EVoE, Eigen Vector of Electron, which is {200,200}, mat(200,200), for every energy value, the corresponding eigen vector is {1,200},
//# bcd, bonding charging density, a term is the sum of adjacent EVoE on lattice and times population of electrons on each lattice,
//# F, force on lattice, derived by Feynman-Hellman theorem,
//# v, velocity of lattice,
//# bcdSum, a number used in Feynman-Hellman calculation,
//# Pula, population on lattice,
//# Hlcn, History of lcn,

double delta(int m, int n){
    double result = 0;
    if(m == n)
        result = 1;
    return result;
}


template<int N> string concatenatestr(const string (&a)[N]){
    string s = "";
    for(int i = 0; i < N; ++i)
        s += a[i];
    return s;
}


class SSHH{
public:
    // constructor to initialize necessary data, more importantly input data,
    // need a function to output data into .dat
    // destructors to delete dynamcially allocated memory by declaring to delete its pointer.
    SSHH(const char* dataPlace, int iN, double it0, double ite, double ialpha, double ik, double im, double iA, double iIty):\
        N(iN),t0(it0),te(ite),alpha(ialpha),k(ik),m(im),A(iA),Ity(iIty),bcdSum(0),c1(0),c2(0){

        // initialize lcn,
        vec lcnInit(N,fill::zeros);
        lcn = lcnInit;
        ifstream lcnData(dataPlace,ifstream::in);
        for(int i = 0; i < N; ++i)
            lcnData >> lcn(i);
        // initialize ham,
        mat hamInit(N,N,fill::zeros);
        ham = hamInit;
        // initialize EEoE, EVoE, bcd, F, v;
        EEoE = lcnInit;
        //        BEEoE = lcnInit;
        bcd = lcnInit;
        F = lcnInit;
        v = lcnInit;
        Pula = lcnInit;
        for(int i = 0; i < N/2; ++i)
            Pula(i) = 2;
        //        Pula(N/2 - 1) = 1;// Polaron
        EVoE = hamInit;
        //        BEVoE = hamInit;
//        vec vtempInit(N-2,fill::zeros);
        //        BEEoE = vtempInit;
        BEEoE = lcnInit;
//        mat tempInit(N-2, N-2, fill::zeros);
//        mat tempInit2(N,N,fill::ones);
        B = hamInit;
        BEVoE = hamInit;
        //        B = tempInit;
        //        BEVoE = tempInit;
        //        B = tempInit2;
        //        BEVoE = tempInit2;
        mat tlcn(1000,N,fill::zeros);
        Hlcn = tlcn;
//        HBEVoE1 = tlcn;

    }

    // getEigen:
    // the innermost level is the lattice level of vector elements
    // the outer level is the energy level of eigenvectors
    SSHH& getEigen(){
        for(int d1 = 0; d1 < N; ++d1)
            for(int d2 = 0; d2 < N; ++d2){
                if( d1 + 1 == d2 )
                    ham(d1,d2) = -t0-pow(-1.0,d2 + 1)*( alpha * ( lcn(d2) + lcn(d1) ) + te );
                else if( d1 - 1 == d2)
                    ham(d1,d2) = -t0-pow(-1.0,d1 + 1)*( alpha * ( lcn(d2) + lcn(d1) ) + te );
            }

        eig_sym(EEoE,EVoE,ham,"dc");
        // the same column, the same eigenvalue,
        // row stands for energy
        // row stands for lattice,
        //        EVoE = EVoE.t();
        return *this;
    }

    // get bond charging density
    SSHH& getBCD(){

        bcd.zeros();
        for(int j = 0; j < N - 1; ++j)
            for(int i = 0; i < N; ++i){
                bcd(j) += Pula(i) * EVoE(j,i) * EVoE(j + 1,i);
            }
        for(int j = 0; j < N - 1; ++j)
            bcdSum += bcd(j);
        bcdSum *= alpha*2/N;

        return *this;
    }

    // FH Calculation, it has a time counting factor as an argument, time starts from 0,
    SSHH& getFH(int time){

        F(0) = 2*alpha*bcd(0) - k*(lcn(0)+lcn(1)) - bcdSum;
        // 0.98 is for display,
        v(0) = 0.98*v(0) + F(0) / m;
        lcn(0) += v(0);
        Hlcn(time,0) = lcn(0);
        for(int i = 1; i < N - 1; ++i){
            F(i) = 2*pow(-1.0, i)*alpha*(bcd(i) - bcd(i - 1)) - k*(lcn(i - 1) + lcn(i + 1) + 2*lcn(i));
            v(i) = 0.98*v(i) + F(i) / m;
            lcn(i) += v(i);
            lcn(i) = (lcn(i - 1) + lcn(i) + lcn(i + 1)) / 3.0;
            Hlcn(time, i) = lcn(i);
        }
        F(N - 1) = 2*alpha*bcd(N - 2) - k*(lcn(N - 1) + lcn(N - 2)) - bcdSum;
        v(N - 1) = 0.98*v(N - 1) + F(N - 1) / m;
        lcn(N - 1) += v(N - 1);
        Hlcn(time,N - 1) = lcn(N - 1);

        return *this;
    }

    SSHH& getTrans(){
        double rd = 0; // temp;
        for(int i = 0; i < N; ++i)
            rd += EVoE(i, N/2) * (i + 1) * EVoE(i,N/2 - 1);
        rd *= A;
        double r = rd * rd * 0.00376 * pow( EEoE(N/2) - EEoE(N/2 - 1) ,3);
        Pula(N/2) = Pula(N/2) - r * Pula(N/2) * 0.000001 + Pula(N/2-1) * rd * rd * Ity - Pula(N/2) * rd * rd * Ity;
        Pula(N/2 - 1) = 2 - Pula(N/2);
        //        Pula(N/2 - 1) = 1 - Pula(N/2);

        //        cout << Pula(N/2 - 1) << ", " << Pula(N/2) << endl;
        return *this;
    }


    /*=====new function======*/

    SSHH& getBEigen(){
//    SSHH& getBEigen(int time){
        // IMPORTANT!!! IF THE TASK IS CREATED NEW.
        //        Pula(N/2) = 0.99;
        //        Pula(N/2-1) = 1.0;

        double temp = 0;

        //        for(int m0 = 1; m0 < N - 1; ++m0)
        //            for(int n0 = 1; n0 < N - 1; ++n0){
        //                temp = 0;
        //                B(m0-1, n0-1) = k * (delta(m0, n0 - 1) + 2*delta(m0, n0) + delta(m0, n0 + 1));
        //                for(int i0 = 0; i0 < N/2+1; ++i0)
        //                    for(int j0 = 0; j0 < N; ++j0){
        //                        if( i0 != j0 ){
        //                            c1 = EVoE(m0, i0)*( EVoE(m0 + 1, j0) - EVoE(m0 - 1, j0) ) + EVoE(m0, j0)*( EVoE(m0 + 1, i0) - EVoE(m0 - 1, i0) );
        //                            c2 = EVoE(n0, i0)*( EVoE(n0 + 1, j0) - EVoE(n0 - 1, j0) ) + EVoE(n0, j0)*( EVoE(n0 + 1, i0) - EVoE(n0 - 1, i0) );
        //                            temp += (Pula(i0) * c1 * c2) / (EEoE(i0)-EEoE(j0));
        //                        }
        //                    }
        //                temp *= 2*alpha*alpha*pow(-1.0, m0 + n0);
        //                B(m0-1,n0-1) += temp;
        //            }


// 198 MODEL
//        for(int m0 = 1; m0 < N - 1; ++m0)
//            for(int n0 = 1; n0 < N - 1; ++n0){
//                temp = 0;
//                B(m0-1, n0-1) = k * ( ( 2*delta(m0, n0) + delta(m0, n0 + 1)) * (1 - delta(n0, N - 1)) + (delta(m0,n0) + delta(m0,n0 - 1)) * (1 - delta(n0, 0)));
//                for(int i0 = 0; i0 < N/2+1; ++i0)
//                    for(int j0 = 0; j0 < N; ++j0){
//                        if( i0 != j0 ){

//                            c1 = (1 - delta(m0, N - 1)) * (EVoE(m0 + 1, i0)*EVoE(m0, j0) + EVoE(m0, i0) * EVoE(m0 + 1, j0)) - (1 - delta(m0, 0)) * (EVoE(m0 , i0)*EVoE(m0 - 1, j0) + EVoE(m0 - 1, i0) * EVoE(m0 , j0));
//                            c2 = (1 - delta(n0, N - 1)) * (EVoE(n0 + 1, i0)*EVoE(n0, j0) + EVoE(n0, i0) * EVoE(n0 + 1, j0)) - (1 - delta(n0, 0)) * (EVoE(n0 , i0)*EVoE(n0 - 1, j0) + EVoE(n0 - 1, i0) * EVoE(n0 , j0));
//                            temp += (Pula(i0) * c1 * c2) / (EEoE(i0)-EEoE(j0));
//                        }
//                    }
//                temp *= 2*alpha*alpha*pow(-1.0, m0 + n0);
//                B(m0-1,n0-1) += temp;
//            }
// 200 MODEL

        // make m0 {0, N - 1}
        // make n0 {0, N - 1}
        for(int m0 = 1; m0 < N - 1; ++m0)
            for(int n0 = 1; n0 < N - 1; ++n0){
                temp = 0;
                B(m0, n0) = k * ( ( delta(m0, n0) + delta(m0, n0 + 1)) * (1 - delta(n0, N - 1)) + (delta(m0,n0) + delta(m0,n0 - 1)) * (1 - delta(n0, 0)));
                for(int i0 = 0; i0 < N/2+1; ++i0)
                    for(int j0 = 0; j0 < N; ++j0){
                        if( i0 != j0 ){

                            c1 = (1 - delta(m0, N - 1)) * (EVoE(m0 + 1, i0)*EVoE(m0, j0) + EVoE(m0, i0) * EVoE(m0 + 1, j0)) - (1 - delta(m0, 0)) * (EVoE(m0 , i0)*EVoE(m0 - 1, j0) + EVoE(m0 - 1, i0) * EVoE(m0 , j0));
                            c2 = (1 - delta(n0, N - 1)) * (EVoE(n0 + 1, i0)*EVoE(n0, j0) + EVoE(n0, i0) * EVoE(n0 + 1, j0)) - (1 - delta(n0, 0)) * (EVoE(n0 , i0)*EVoE(n0 - 1, j0) + EVoE(n0 - 1, i0) * EVoE(n0 , j0));
                            temp += (Pula(i0) * c1 * c2) / (EEoE(i0)-EEoE(j0));
                        }
                    }
                temp *= 2*alpha*alpha*pow(-1.0, m0 + n0);
                B(m0,n0) += temp;
            }

        // m0 = 0, n0 -> {1, N - 2}
        for(int n0 = 1; n0 < N - 1; ++n0){
            temp = 0;
            B(0, n0) = k * ( ( delta(0, n0) + delta(0, n0 + 1)) * (1 - delta(n0, N - 1)) + (delta(0,n0) + delta(0,n0 - 1)) * (1 - delta(n0, 0)));
            for(int i0 = 0; i0 < N/2+1; ++i0)
                for(int j0 = 0; j0 < N; ++j0){
                    if( i0 != j0 ){
                        c1 = (EVoE(0 + 1, i0)*EVoE(0, j0) + EVoE(0, i0) * EVoE(0 + 1, j0));
                        c2 = (EVoE(n0 + 1, i0)*EVoE(n0, j0) + EVoE(n0, i0) * EVoE(n0 + 1, j0)) - (EVoE(n0 , i0)*EVoE(n0 - 1, j0) + EVoE(n0 - 1, i0) * EVoE(n0 , j0));
                        temp += (Pula(i0) * c1 * c2) / (EEoE(i0)-EEoE(j0));
                    }
                }
            temp *= 2*alpha*alpha*pow(-1.0, 0 + n0);
            B(0,n0) += temp;
        }

        // m0 = N - 1, n0 -> {1, N - 2}
        for(int n0 = 1; n0 < N - 1; ++n0){
            temp = 0;
            B(N-1, n0) = k * ( ( delta(N-1, n0) + delta(N-1, n0 + 1)) * (1 - delta(n0, N - 1)) + (delta(N-1,n0) + delta(N-1,n0 - 1)) * (1 - delta(n0, 0)));
            for(int i0 = 0; i0 < N/2+1; ++i0)
                for(int j0 = 0; j0 < N; ++j0){
                    if( i0 != j0 ){
                        c1 = - (EVoE(N-1 , i0)*EVoE(N-1 - 1, j0) + EVoE(N-1 - 1, i0) * EVoE(N-1 , j0));
                        c2 = (EVoE(n0 + 1, i0)*EVoE(n0, j0) + EVoE(n0, i0) * EVoE(n0 + 1, j0)) - (EVoE(n0 , i0)*EVoE(n0 - 1, j0) + EVoE(n0 - 1, i0) * EVoE(n0 , j0));
                        temp += (Pula(i0) * c1 * c2) / (EEoE(i0)-EEoE(j0));
                    }
                }
            temp *= 2*alpha*alpha*pow(-1.0, N-1 + n0);
            B(N-1,n0) += temp;
        }

        // n0 = 0, m0 -> {1, N - 2}
        for(int m0 = 1; m0 < N - 1; ++m0){
            temp = 0;
            B(m0, 0) = k * ( ( delta(m0, 0) + delta(m0, 0 + 1)) * (1 - delta(0, N - 1)) + (delta(m0,0) + delta(m0,0 - 1)) * (1 - delta(0, 0)));
            for(int i0 = 0; i0 < N/2+1; ++i0)
                for(int j0 = 0; j0 < N; ++j0){
                    if( i0 != j0 ){
                        c1 = (EVoE(m0 + 1, i0)*EVoE(m0, j0) + EVoE(m0, i0) * EVoE(m0 + 1, j0)) - (EVoE(m0 , i0)*EVoE(m0 - 1, j0) + EVoE(m0 - 1, i0) * EVoE(m0 , j0));
                        c2 = (EVoE(0 + 1, i0)*EVoE(0, j0) + EVoE(0, i0) * EVoE(0 + 1, j0));
                        temp += (Pula(i0) * c1 * c2) / (EEoE(i0)-EEoE(j0));
                    }
                }
            temp *= 2*alpha*alpha*pow(-1.0, m0 + 0);
            B(m0,0) += temp;
        }

        // n0 = N - 1, m0 -> {1, N - 2}
        for(int m0 = 1; m0 < N - 1; ++m0){
            temp = 0;
            B(m0, N - 1) = k * ( ( delta(m0, N - 1) + delta(m0, N - 1 + 1)) * (1 - delta(N - 1, N - 1)) + (delta(m0,N - 1) + delta(m0,N - 1 - 1)) * (1 - delta(N - 1, 0)));
            for(int i0 = 0; i0 < N/2+1; ++i0)
                for(int j0 = 0; j0 < N; ++j0){
                    if( i0 != j0 ){
                        c1 = (EVoE(m0 + 1, i0)*EVoE(m0, j0) + EVoE(m0, i0) * EVoE(m0 + 1, j0)) - (EVoE(m0 , i0)*EVoE(m0 - 1, j0) + EVoE(m0 - 1, i0) * EVoE(m0 , j0));
                        c2 = - (EVoE(N - 1 , i0)*EVoE(N - 1 - 1, j0) + EVoE(N - 1 - 1, i0) * EVoE(N - 1 , j0));
                        temp += (Pula(i0) * c1 * c2) / (EEoE(i0)-EEoE(j0));
                    }
                }
            temp *= 2*alpha*alpha*pow(-1.0, m0 + N - 1);
            B(m0,N - 1) += temp;
        }

        // m0 -> 0, n0 -> 0,
        temp = 0;
        B(0, 0) = k * ( ( delta(0, 0) + delta(0, 0 + 1)) * (1 - delta(0, N - 1)) + (delta(0,0) + delta(0,0 - 1)) * (1 - delta(0, 0)));
        for(int i0 = 0; i0 < N/2+1; ++i0)
            for(int j0 = 0; j0 < N; ++j0){
                if( i0 != j0 ){
                    c1 = (EVoE(0 + 1, i0)*EVoE(0, j0) + EVoE(0, i0) * EVoE(0 + 1, j0));
                    c2 = (EVoE(0 + 1, i0)*EVoE(0, j0) + EVoE(0, i0) * EVoE(0 + 1, j0));
                    temp += (Pula(i0) * c1 * c2) / (EEoE(i0)-EEoE(j0));
                }
            }
        temp *= 2*alpha*alpha*pow(-1.0, 0 + 0);
        B(0,0) += temp;

        // m0 -> 0, n0 -> N - 1
        temp = 0;
        B(0, N - 1) = k * ( ( delta(0, N - 1) + delta(0, N - 1 + 1)) * (1 - delta(N - 1, N - 1)) + (delta(0,N - 1) + delta(0,N - 1 - 1)) * (1 - delta(N - 1, 0)));
        for(int i0 = 0; i0 < N/2+1; ++i0)
            for(int j0 = 0; j0 < N; ++j0){
                if( i0 != j0 ){
                    c1 = (EVoE(0 + 1, i0)*EVoE(0, j0) + EVoE(0, i0) * EVoE(0 + 1, j0));
                    c2 = - (EVoE(N - 1 , i0)*EVoE(N - 1 - 1, j0) + EVoE(N - 1 - 1, i0) * EVoE(N - 1 , j0));
                    temp += (Pula(i0) * c1 * c2) / (EEoE(i0)-EEoE(j0));
                }
            }
        temp *= 2*alpha*alpha*pow(-1.0, 0 + N - 1);
        B(0,N - 1) += temp;

        // n0 -> 0, m0 -> N - 1
        temp = 0;
        B(N - 1, 0) = k * ( ( delta(N - 1, 0) + delta(N - 1, 0 + 1)) * (1 - delta(0, N - 1)) + (delta(N - 1,0) + delta(N - 1,0 - 1)) * (1 - delta(0, 0)));
        for(int i0 = 0; i0 < N/2+1; ++i0)
            for(int j0 = 0; j0 < N; ++j0){
                if( i0 != j0 ){
                    c1 =  - (EVoE(N - 1 , i0)*EVoE(N - 1 - 1, j0) + EVoE(N - 1 - 1, i0) * EVoE(N - 1 , j0));
                    c2 = (EVoE(0 + 1, i0)*EVoE(0, j0) + EVoE(0, i0) * EVoE(0 + 1, j0));
                    temp += (Pula(i0) * c1 * c2) / (EEoE(i0)-EEoE(j0));
                }
            }
        temp *= 2*alpha*alpha*pow(-1.0, N - 1 + 0);
        B(N - 1,0) += temp;

        // n0 -> N - 1, m0 -> N - 1
        temp = 0;
        B(N - 1, N - 1) = k * ( ( delta(N - 1, N - 1) + delta(N - 1, N - 1 + 1)) * (1 - delta(N - 1, N - 1)) + (delta(N - 1,N - 1) + delta(N - 1,N - 1 - 1)) * (1 - delta(N - 1, 0)));
        for(int i0 = 0; i0 < N/2+1; ++i0)
            for(int j0 = 0; j0 < N; ++j0){
                if( i0 != j0 ){
                    c1 =  - (EVoE(N - 1 , i0)*EVoE(N - 1 - 1, j0) + EVoE(N - 1 - 1, i0) * EVoE(N - 1 , j0));
                    c2 =  - (EVoE(N - 1 , i0)*EVoE(N - 1 - 1, j0) + EVoE(N - 1 - 1, i0) * EVoE(N - 1 , j0));
                    temp += (Pula(i0) * c1 * c2) / (EEoE(i0)-EEoE(j0));
                }
            }
        temp *= 2*alpha*alpha*pow(-1.0, N - 1 + N - 1);
        B(N - 1,N - 1) += temp;


// solve it
        eig_sym(BEEoE, BEVoE, B, "dc");

        // record HBEVoE1
//        for(int i = 0; i < N; ++i)
//            HBEVoE1(time,i) = BEVoE(i,1);

        return *this;
    }

    /*========================*/


    // output data
    // as an example, the output function only out puts lcn;
//    SSHH& OutPutData(string lcnOutput, string EVoEOutput, string EEoEOutput){
    SSHH& OutPutData(const char * lcnOutput, const char * EVoEOutput, const char * EEoEOutput/*, const char * HlcnOutput*/){

        // lcn output
        ofstream lcnOut(lcnOutput, ofstream::app);
        for(int i = 0; i < N; ++i){
            lcnOut << lcn(i) << " ";
        }
        // EVoE output
        ofstream EVoEOut(EVoEOutput, ofstream::app);
        for(int i = 0; i < N; ++i)
            for(int j = 0; j < N; ++j){
                if(j == N - 1)
                    EVoEOut << EVoE(i,j) << "\n";
                else
                    EVoEOut << EVoE(i,j) << " ";
            }

        // lcn output
        ofstream EEoEOut(EEoEOutput, ofstream::app);
        for(int i = 0; i < N; ++i){
            EEoEOut << EEoE(i) << " ";
        }

        // Hlcn output
//        ofstream HlcnOut(HlcnOutput, ofstream::app);
//        for(int i = 0; i < 1000; ++i)
//            for(int j = 0; j < N; ++j){
//                if(j == N - 1)
//                    HlcnOut << Hlcn(i,j) << "\n";
//                else
//                    HlcnOut << Hlcn(i,j) << " ";
//            }

        return *this;
    }

//    SSHH& OutPutData2(string EVoEOutput, string EEoEOutput){
    SSHH& OutPutData2(const char * EVoEOutput, const char * EEoEOutput){


        // EVoE output
        // MODEL 198
//        ofstream EVoEOut(EVoEOutput, ofstream::app);
//        for(int i = 0; i < N-2; ++i)
//            for(int j = 0; j < N-2; ++j){
//                if(j == N - 3)
//                    EVoEOut << BEVoE(i,j) << "\n";
//                else
//                    EVoEOut << BEVoE(i,j) << " ";
//            }

        // MODEL 200
        ofstream EVoEOut(EVoEOutput, ofstream::app);
        for(int i = 0; i < N-0; ++i)
            for(int j = 0; j < N-0; ++j){
                if(j == N - 1)
                    EVoEOut << BEVoE(i,j) << "\n";
                else
                    EVoEOut << BEVoE(i,j) << " ";
            }

        // lcn output
        // MODEL 200
        ofstream EEoEOut(EEoEOutput, ofstream::app);
        for(int i = 0; i < N-0; ++i){
            EEoEOut << BEEoE(i) << " ";
        }

        return *this;
    }

//    SSHH& OutPutData3(const char* HBEVoE1Output){
//        // HBEVoE1 output
//        ofstream HBEVoE1Out(HBEVoE1Output, ofstream::app);
//        for(int i = 0; i < 1000; ++i)
//            for(int j = 0; j < N; ++j){
//                if(j == N - 1)
//                    HBEVoE1Out << HBEVoE1(i,j) << "\n";
//                else
//                    HBEVoE1Out << HBEVoE1(i,j) << " ";
//            }

//        return *this;
//    }

protected:
    // ask questions and decide what type to use for each data member?
    int N;
    vec lcn;
    mat ham;
    double t0, te, alpha, k, m, A, bcdSum, Ity;
    vec Pula;
    vec EEoE;
    mat EVoE;
    vec bcd; // where should be N - 1;
    vec F;
    vec v;
    mat Hlcn;
    /*=========add new========*/
    mat B;
    vec BEEoE;
    mat BEVoE;
//    mat HBEVoE1;
    double c1, c2;
};


int main()
{
    clock_t s0, s1;
    s0 = clock();

    SSHH task_Cong("D:\\DataNumerics\\phi.dat", \
                   200, 2.5, 0.06, 4.1, 21, 1345.49, 1.22, 0.000091 * 20);
    for(int i = 0; i < 500; ++i)
        task_Cong.getEigen().getBCD().getFH(i);
    task_Cong.OutPutData("D:\\DataNumerics\\elcn.dat", \
                         "D:\\DataNumerics\\eEVoE.dat", "D:\\DataNumerics\\eEEoE.dat"/*, "D:\\DataNumerics\\eHlcn.dat"*/);

    SSHH task_Cong2("D:\\DataNumerics\\elcn.dat", \
                    200, 2.5, 0.06, 4.1, 21, 1345.49, 1.22, 0.000091 * 20);
    //    for(int i = 0; i < 5; ++i)
    //        task_Cong2.getEigen().getBCD().getFH(i).getTrans();
    //    task_Cong2.OutPutData("D:\\DataNumerics\\elcn2.dat", \
    //                          "D:\\DataNumerics\\eEVoE2.dat", "D:\\DataNumerics\\eEEoE2.dat", "D:\\DataNumerics\\eHlcn2.dat");

//    for(int i = 0; i < 2; ++i)
//        task_Cong2.getEigen().getBCD().getFH(i).getTrans().getBEigen(i);
//    task_Cong2.OutPutData("D:\\DataNumerics\\elcn2.dat", \
//                                      "D:\\DataNumerics\\eEVoE2.dat", "D:\\DataNumerics\\eEEoE2.dat", "D:\\DataNumerics\\eHlcn2.dat").OutPutData2("D:\\DataNumerics\\eBEVoE.dat", "D:\\DataNumerics\\eBEEoE.dat").OutPutData3("D:\\DataNumerics\\eHBEVoE1.dat");



    // ================ address ==================
    // 6 names to be prepared
    string temp;
//    string ads_elcn [3] = {"D:\\DataNumerics\\res",temp,"s\\elcn2.dat"};
//    string ads_EVoE [3] = {"D:\\DataNumerics\\res", temp, "s\\eEVoE2.dat"};
//    string ads_EEoE [3] = {"D:\\DataNumerics\\res", temp, "s\\eEEoE2.dat"};
//    string ads_eBEVoE [3] = {"D:\\DataNumerics\\res", temp, "s\\eBEVoE.dat"};
//    string ads_eBEEoE [3] = {"D:\\DataNumerics\\res", temp, "s\\eBEEoE.dat"};

    string ads_elcn [3] = {"D:\\DataNumerics2\\res",temp,"s\\elcn2.dat"};
    string ads_EVoE [3] = {"D:\\DataNumerics2\\res", temp, "s\\eEVoE2.dat"};
    string ads_EEoE [3] = {"D:\\DataNumerics2\\res", temp, "s\\eEEoE2.dat"};
    string ads_eBEVoE [3] = {"D:\\DataNumerics2\\res", temp, "s\\eBEVoE.dat"};
    string ads_eBEEoE [3] = {"D:\\DataNumerics2\\res", temp, "s\\eBEEoE.dat"};


//    const char* thisadd;
//    thisadd = concatenatestr(ads_elcn).c_str();
//    cout << thisadd << endl;
//    cout << concatenatestr(ads_elcn).c_str() << endl;

    char address [5][100];
    string btemp;
//    string address [5];

//
    for(int i = 0; i < 0; ++i)
        task_Cong2.getEigen().getBCD().getFH(i).getTrans();
    for(int i = 0; i < 31; ++i){
        task_Cong2.getEigen().getBCD().getFH(i).getTrans();
        // mod some number n means to record value every n times,,
        if( (i + 1) % 1 == 0){
            stringstream ss;
            ss << i + 1;
            ss >> temp;
//            cout << temp << endl;
            ads_elcn[1] = temp;
            ads_EVoE[1] = temp;
            ads_EEoE[1] = temp;
            ads_eBEVoE[1] = temp;
            ads_eBEEoE[1] = temp;
            btemp = concatenatestr(ads_elcn);
            strcpy(address[0],btemp.c_str());
            cout << address[0] << endl;
            btemp = concatenatestr(ads_EVoE);
            strcpy(address[1],btemp.c_str());
            cout << address[1] << endl;
            btemp = concatenatestr(ads_EEoE);
            strcpy(address[2],btemp.c_str());
            cout << address[2] << endl;
            btemp = concatenatestr(ads_eBEVoE);
            strcpy(address[3],btemp.c_str());
            cout << address[3] << endl;
            btemp = concatenatestr(ads_eBEEoE);
            strcpy(address[4],btemp.c_str());
            cout << address[4] << endl;
            task_Cong2.getBEigen().OutPutData(address[0], address[1], address[2]).OutPutData2(address[3], address[4]);
        }
    }

    // ===========================================

//    for(int i = 0; i < 5; ++i)
//        task_Cong2.getEigen().getBCD().getFH(i).getTrans();
//    task_Cong2.getBEigen().OutPutData("D:\\DataNumerics\\res1250\\elcn2.dat", \
//                                      "D:\\DataNumerics\\res1250\\eEVoE2.dat", "D:\\DataNumerics\\res1250\\eEEoE2.dat"/*, "D:\\DataNumerics\\eHlcn2.dat"*/).OutPutData2("D:\\DataNumerics\\res1250\\eBEVoE.dat", "D:\\DataNumerics\\res1250\\eBEEoE.dat");
//    task_Cong2.getBEigen().OutPutData("D:\\DataNumerics\\elcn2.dat", \
//                                      "D:\\DataNumerics\\eEVoE2.dat", "D:\\DataNumerics\\eEEoE2.dat"/*, "D:\\DataNumerics\\eHlcn2.dat"*/).OutPutData2("D:\\DataNumerics\\eBEVoE.dat", "D:\\DataNumerics\\eBEEoE.dat");
//        SSHH task_Cong3("D:\\DataNumerics\\elcn2.dat", \
//                        200, 2.5, 0.06, 4.1, 21, 1345.49, 1.22, 0.000091 * 20);
//        task_Cong3.getEigen().getBEigen();
//        task_Cong3.OutPutData2("D:\\DataNumerics\\eBEVoE.dat", "D:\\DataNumerics\\eBEEoE.dat");


    s1 = clock();
    cout << "Whole Process: " << double(s1 - s0) / CLOCKS_PER_SEC << " s." << endl;
    return 0;
}
