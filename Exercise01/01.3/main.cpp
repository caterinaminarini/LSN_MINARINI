#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "../../ParallelRandomNumberGenerator/random.h"

using namespace std;

void DataBlocking(double*, double*, double*, int, int);
double error (double, double, int);

int main (int argc, char *argv[]){
        
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
    
    int m=100000;
    int Nthr=10000; //trials
    int Nbl=100; //blocks
    double L=0.9; //needle's lenght
    double d=1.; //distance
    double* pi = new double [m];
    double* sum_prog_pi = new double [Nbl];
    double* err_prog_pi = new double [Nbl];
    
    for (int j=0; j<m; j++){
        double P=0.;
        double Nhit=0.;
        for (int k=0; k<Nthr; k++){
            double c=rnd.Rannyu(0., d);
            double x, y, theta;
            do{
                x=rnd.Rannyu(-1., 1.);
                y=rnd.Rannyu(0., 1.);
                theta=acos(x/pow((x*x+y*y),0.5));
            }while(x*x+y*y>=1.);
            double y_up=c+(L/2.)*sin(theta);
            double y_down=c-(L/2.)*sin(theta);
            if (y_up>=d || y_down<=0.){
                Nhit++;
            }
        }
        P=double(Nhit)/double(Nthr);
        pi[j]=double((2.*L)/(P*d));
    }
    
    DataBlocking(sum_prog_pi, err_prog_pi, pi, m, Nbl);
    
    ofstream uscita("results01.3.txt");
    uscita << "#PI   " << "    err_PI" << endl;
    for (int i=0; i<Nbl; i++){
        uscita << sum_prog_pi[i] << "  " << err_prog_pi[i] << endl;
    }
    
    delete [] pi;
    delete [] sum_prog_pi;
    delete [] err_prog_pi;
    
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

