#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "../../ParallelRandomNumberGenerator/random.h"

using namespace std;

void DataBlocking(double*, double*, double*, int, int);
double error (double, double, int);

int main (int argc, char *argv[]){
    
    // --- PRIMA PARTE ---
    
    Random rnd;
    int seed[4];
    int p1, p2;
    ifstream Primes("../../ParallelRandomNumberGenerator/Primes");
    if (Primes.is_open()){
       Primes >> p1 >> p2 ;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();
    
    ifstream input("../../ParallelRandomNumberGenerator/seed.in");
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
    
    int M=10000;
    int N=100;
    double* f = new double[M];
    double* sum_prog_f = new double [N];
    double* err_prog_f = new double [N];
    
    for (int i=0; i<M; i++){
        f[i]=(M_PI/2.)*cos((M_PI/2.)*rnd.Rannyu());
    }
    
    DataBlocking(sum_prog_f, err_prog_f, f, M, N);
    
    // --- SECONDA PARTE ---
    
    // Devo generare i punti come p(x)=...
    // Uso reject
    
    double pmax=1./(1. - pow(M_PI, 2.)/24. + pow(M_PI, 4.)/1920.);
    double* g = new double[M];
    double* sum_prog_g = new double [N];
    double* err_prog_g = new double [N];
    double p[M];
    
    for (int i=0; i<M; i++){
        double s = 0.;
        double r = 0.;
        do {
            s = rnd.Rannyu();
            r = rnd.Rannyu();
        } while (s > (1. - pow(M_PI*r, 2.)/8. + pow(M_PI*r, 4.)/384.));
        p[i] = (1. - pow(M_PI*r, 2.)/8. + pow(M_PI*r, 4.)/384.)*pmax;
        g[i] = (M_PI/2.)*cos((M_PI/2.)*r)/p[i];
    }
    
    DataBlocking(sum_prog_g, err_prog_g, g, M, N);
    
    ofstream uscita("results02.1.txt");
    uscita << "#<I>1   " << "    err<I>  " << "<I>2   " << "    err<I>2" << endl;
    for (int i=0; i<N; i++){
        uscita << sum_prog_f[i] << "  " << err_prog_f[i] << "  " << sum_prog_g[i] << "  " << err_prog_g[i] << endl;
    }
    
    delete [] f;
    delete [] sum_prog_f;
    delete [] err_prog_f;
    
    delete [] g;
    delete [] sum_prog_g;
    delete [] err_prog_g;
    
    
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


