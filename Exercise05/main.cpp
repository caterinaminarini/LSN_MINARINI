#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "../ParallelRandomNumberGenerator/random.h"

using namespace std;

double Psi1 (double, double, double);
double Psi2 (double, double, double);
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

    
    // (!!!) ho scritto in unità di a0
    
    
    int M = 1e06;
    int N = 100;
    double* x = new double [M];
    double* y = new double [M];
    double* z = new double [M];
    int accept=0.;
    
    // sto partendo dal raggio a0
    
    
    // ---1S UNIFORM---
    ofstream config;
    config.open("coordinates1.txt");
    config << "#x   " << "  y   " << "  z   " << endl;
    for (int i=0; i<M; i++){
        if (i==0){
            x[0] = 1.;
            y[0] = 1.;
            z[0] = 1.;
        }
        else {
            double xnew = x[i-1] + rnd.Rannyu(-1., 1.);
            double ynew = y[i-1] + rnd.Rannyu(-1., 1.);
            double znew = z[i-1] + rnd.Rannyu(-1., 1.);
            
            double A = min(1., Psi1(xnew, ynew, znew)/Psi1(x[i-1], y[i-1], z[i-1]));
            double s = rnd.Rannyu();
            if (s<=A) {
                x[i] = xnew;
                y[i] = ynew;
                z[i] = znew;
                accept++;
            }
            else{
                x[i] = x[i-1];
                y[i] = y[i-1];
                z[i] = z[i-1];
            }
        }
        config << x[i] << "  " << y[i] << "  " << z[i] << endl;
    }
    config.close();
    cout << "Accettazione1=" << double(accept)/double(M) << endl;
    
    double* f = new double[M];
    double* sum_prog_f = new double[N];
    double* err_prog_f = new double[N];
    
    for (int i=0; i<M; i++){
        f[i] = pow(x[i]*x[i]+y[i]*y[i]+z[i]*z[i], 0.5);
    }
    DataBlocking(sum_prog_f, err_prog_f, f, M, N);
    
    ofstream uscita("results.txt");
    uscita << "#R   " << "    errR  " << endl;
    for (int i=0; i<N; i++){
        uscita << sum_prog_f[i] << "  " << err_prog_f[i] << endl;
    }
    uscita.close();
    
    // ---2P UNIFORM---
    config.open("coordinates2.txt");
    config << "#x   " << "  y   " << "  z   " << endl;
    accept=0;
    for (int i=0; i<M; i++){
        if (i==0){
            x[0] = 4.;
            y[0] = 4.;
            z[0] = 4.;
        }
        else {
            double xnew = x[i-1] + rnd.Rannyu(-3., 3.);
            double ynew = y[i-1] + rnd.Rannyu(-3., 3.);
            double znew = z[i-1] + rnd.Rannyu(-3., 3.);
            
            double A = min(1., Psi2(xnew, ynew, znew)/Psi2(x[i-1], y[i-1], z[i-1]));
            double s = rnd.Rannyu();
            if (s<=A) {
                x[i] = xnew;
                y[i] = ynew;
                z[i] = znew;
                accept++;
            }
            else{
                x[i] = x[i-1];
                y[i] = y[i-1];
                z[i] = z[i-1];
            }
        }
        config << x[i] << "  " << y[i] << "  " << z[i] << endl;
    }
    config.close();
    cout << "Accettazione2=" << double(accept)/double(M) << endl;
        
    for (int i=0; i<M; i++){
        f[i] = pow(x[i]*x[i]+y[i]*y[i]+z[i]*z[i], 0.5);
    }
    DataBlocking(sum_prog_f, err_prog_f, f, M, N);
    
    uscita.open("results2.txt");
    uscita << "#R   " << "    errR  " << endl;
    for (int i=0; i<N; i++){
        uscita << sum_prog_f[i] << "  " << err_prog_f[i] << endl;
    }
    uscita.close();
    
    // ---1S GAUSS---
    config.open("coordinatesGauss.txt");
    config << "#x   " << "  y   " << "  z   " << endl;
    accept=0;
    for (int i=0; i<M; i++){
        if (i==0){
            x[0] = 1.;
            y[0] = 1.;
            z[0] = 1.;
        }
        else {
            double xnew = x[i-1] + rnd.Gauss(0., 0.8);
            double ynew = y[i-1] + rnd.Gauss(0., 0.8);
            double znew = z[i-1] + rnd.Gauss(0., 0.8);
            
            double A = min(1., Psi1(xnew, ynew, znew)/Psi1(x[i-1], y[i-1], z[i-1]));
            double s = rnd.Rannyu();
            if (s<=A) {
                x[i] = xnew;
                y[i] = ynew;
                z[i] = znew;
                accept++;
            }
            else{
                x[i] = x[i-1];
                y[i] = y[i-1];
                z[i] = z[i-1];
            }
        }
        config << x[i] << "  " << y[i] << "  " << z[i] << endl;
    }
    config.close();
    cout << "Accettazione1G=" << double(accept)/double(M) << endl;
    
    for (int i=0; i<M; i++){
        f[i] = pow(x[i]*x[i]+y[i]*y[i]+z[i]*z[i], 0.5);
    }
    DataBlocking(sum_prog_f, err_prog_f, f, M, N);
    
    uscita.open("resultsGauss.txt");
    uscita << "#R   " << "    errR  " << endl;
    for (int i=0; i<N; i++){
        uscita << sum_prog_f[i] << "  " << err_prog_f[i] << endl;
    }
    uscita.close();
    
    // ---2P GAUSS---
    config.open("coordinatesGauss2.txt");
    config << "#x   " << "  y   " << "  z   " << endl;
    accept=0;
    for (int i=0; i<M; i++){
        if (i==0){
            x[0] = 4.;
            y[0] = 4.;
            z[0] = 4.;
        }
        else {
            double xnew = x[i-1] + rnd.Gauss(0., 1.3);
            double ynew = y[i-1] + rnd.Gauss(0., 1.3);
            double znew = z[i-1] + rnd.Gauss(0., 1.3);
            
            double A = min(1., Psi2(xnew, ynew, znew)/Psi2(x[i-1], y[i-1], z[i-1]));
            double s = rnd.Rannyu();
            if (s<=A) {
                x[i] = xnew;
                y[i] = ynew;
                z[i] = znew;
                accept++;
            }
            else{
                x[i] = x[i-1];
                y[i] = y[i-1];
                z[i] = z[i-1];
            }
        }
        config << x[i] << "  " << y[i] << "  " << z[i] << endl;
    }
    config.close();
    cout << "Accettazione2G=" << double(accept)/double(M) << endl;
    
    for (int i=0; i<M; i++){
        f[i] = pow(x[i]*x[i]+y[i]*y[i]+z[i]*z[i], 0.5);
    }
    DataBlocking(sum_prog_f, err_prog_f, f, M, N);
    
    uscita.open("resultsGauss2.txt");
    uscita << "#R   " << "    errR  " << endl;
    for (int i=0; i<N; i++){
        uscita << sum_prog_f[i] << "  " << err_prog_f[i] << endl;
    }
    uscita.close();
    
    delete [] f;
    delete [] x;
    delete [] y;
    delete [] z;
    delete [] sum_prog_f;
    delete [] err_prog_f;
    
    rnd.SaveSeed();
    return 0;
}

double Psi1 (double x, double y, double z){
    double P = exp(-pow(x*x+y*y+z*z, 0.5));
    return P*P;
}

double Psi2 (double x, double y, double z){
    double P = z*exp(-pow(x*x+y*y+z*z, 0.5)/2.);
    return P*P;
}

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


