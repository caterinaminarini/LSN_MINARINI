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
    
    int M=100000;
    int N=100;
    double* r = new double[M];
    double* sum_prog_r = new double[N];
    double* err_prog_r = new double[N];
    
    
    for (int i=0; i<M; i++){
        r[i]=rnd.Rannyu();
    }
    
    DataBlocking(sum_prog_r, err_prog_r, r, M, N);

    // --- SECONDA PARTE ---
    
    double* f = new double[M];
    double* sum_prog_f = new double[N];
    double* err_prog_f = new double[N];
    
    for (int i=0; i<M; i++){
        f[i]=pow((r[i]-0.5), 2.);
    }
    
    DataBlocking(sum_prog_f, err_prog_f, f, M, N);
    
    // --- TERZA PARTE ---
    
    int M1=100;
    int N1=10000;
    double n[N1];
    double chi2[100];
    int eventi=0;
    
    for (int i=0; i<N1; i++){
        n[i]=0.;
    }
    
    for (int i=0; i<M1; i++){
        chi2[i]=0.;
    }
    
    for (int w=0; w<100; w++){
        double sum=0.;
        for (int s=0; s<100; s++){
            eventi=0;
            for (int j=0; j<N1; j++){
                n[j]=rnd.Rannyu();
                if (n[j]>double(s/100.) && n[j]<double((s+1.)/100.)){
                    eventi++;
                }
            }
            sum+=pow((eventi-double(N1/M1)), 2.);
        }
        chi2[w]=sum/double(N1/M1);
    }
    
    ofstream uscita("results01.1.txt");
    uscita << "#<r>   " << "    err<r>  " << "  <sigma2>    " << "  err<sigma2> " << "  chi2    " << endl;
    for (int i=0; i<N; i++){
        uscita << sum_prog_r[i] << "  " << err_prog_r[i] << "   " << sum_prog_f[i] << "  " << err_prog_f[i] << "    " << chi2[i] << endl;
    }
    
    
    delete [] r;
    delete [] sum_prog_r;
    delete [] err_prog_r;
    delete [] f;
    delete [] sum_prog_f;
    delete [] err_prog_f;
    
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
