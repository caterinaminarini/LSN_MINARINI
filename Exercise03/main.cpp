#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "../ParallelRandomNumberGenerator/random.h"

using namespace std;

void DataBlocking(double*, double*, double*, int, int);
double error (double, double, int);

int main (int argc, char *argv[]){
        
    Random rnd;
    int seed[4];
    int p1, p2;
    ifstream Primes("../ParallelRandomNumberGenerator/Primes");
    if (Primes.is_open()){
       Primes >> p1 >> p2 ;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();
    
    ifstream input("../ParallelRandomNumberGenerator/seed.in");
    string property;
    if (input.is_open()){
       while ( !input.eof() ){
          input >> property;
          if( property == "RANDOMSEED" ){
             input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
             rnd.SetRandom(seed,p1,p2);
          }
       }
       input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;
    
    
    double S0=100.; // asset price at t=0
    double T=1.;  // delivery time
    double K=100.; // strike price
    double r=0.1;  // risk-free interest rate
    double sigma=0.25;  // volatility
    
    // --- PRIMA PARTE ---
    
    int M=100000;
    int N=100;
    
    double* C = new double [M];
    double* sum_prog_C = new double [N];
    double* err_prog_C = new double [N];
    
    double* P = new double [M];
    double* sum_prog_P = new double [N];
    double* err_prog_P = new double [N];
    

    for (int i=0; i<M; i++){
        double S = S0*exp((r-0.5*pow(sigma, 2.))*T+sigma*rnd.Gauss(0., T));
        C[i]=double(exp(-r*T))*double(max(0., S-K));
        P[i]=double(exp(-r*T))*double(max(0., K-S));
    }
    
    DataBlocking(sum_prog_C, err_prog_C, C, M, N);
    DataBlocking(sum_prog_P, err_prog_P, P, M, N);
    
    ofstream uscita;
    uscita.open("results.txt");
    uscita << "#C   " << "    errC  "  << "P   " << "    errP  " << endl;
    for (int i=0; i<N; i++){
        uscita << sum_prog_C[i] << "  " << err_prog_C[i] << "  " << sum_prog_P[i] << "  " << err_prog_P[i] << endl;
    }
    uscita.close();
    
    // --- SECONDA PARTE ---
    
    
    for (int j=0; j<M; j++){
        double s[100];
        for (int i=0; i<100; i++){
            if (i==0) s[i]=S0;
            else {
                s[i]=s[i-1]*exp((r-0.5*pow(sigma, 2))*(T/100.)+sigma*rnd.Gauss(0., 1.)*sqrt(T/100.));
            }
        }
        C[j]=double(exp(-r*T))*double(max(0., s[99]-K));
        P[j]=double(exp(-r*T))*double(max(0., K-s[99]));
    }

    DataBlocking(sum_prog_C, err_prog_C, C, M, N);
    DataBlocking(sum_prog_P, err_prog_P, P, M, N);
    
    uscita.open("results2.txt");
    uscita << "#C   " << "    errC  "  << "P   " << "    errP  " << endl;
    for (int i=0; i<N; i++){
        uscita << sum_prog_C[i] << "  " << err_prog_C[i] << "  " << sum_prog_P[i] << "  " << err_prog_P[i] << endl;
    }
    uscita.close();
    
    rnd.SaveSeed();
    return 0;
}

// ---------

void DataBlocking(double *sum_prog, double* err_prog, double* f, int M, int N){
    for (int i=0; i<N; i++){
            sum_prog[i] = 0.;
            err_prog[i] = 0.;
        }
        
        int L=int(M/N);
        double ave[N];
        double av2[N];
        int k=0;
        double su2_prog[N];

        
        for (int i=0; i<N; i++){
            ave[i]=0.;
            av2[i]=0.;
            su2_prog[i]=0.;
        }
        
        for (int i=0; i<N; i++){
            double sum=0.;
            for (int j=0; j<L; j++){
                k = j + L*i;
                sum += f[k];
            }
            ave[i]=sum/L;
            av2[i]=ave[i]*ave[i];
        }
        
        for (int i=0; i<N; i++){
            for (int j=0; j<i+1; j++){
                sum_prog[i] += ave[j];
                su2_prog[i] += av2[j];
            }
            sum_prog[i] /= double(i+1);
            su2_prog[i] /= double(i+1);
            err_prog[i] = error(sum_prog[i], su2_prog[i], i);
        }
}

double error (double AV, double AV2, int n){
    if (n==0){
        return 0;
    }
    else{
        return pow(((AV2 - AV*AV)/n), 0.5);
    }
}



