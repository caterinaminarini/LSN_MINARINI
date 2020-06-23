#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "../../ParallelRandomNumberGenerator/random.h"

using namespace std;

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

    int M=10000;
    double S1[M];
    double S2[M];
    double S10[M];
    double S100[M];

    // --- UNIFORM ---
    
    for (int i=0; i<M; i++) {
        S1[i]=0.;
        S1[i] = rnd.Rannyu();
    }
    
    for (int i=0; i<M; i++){
        S2[i]=0.;
        double sum=0.;
        for (int j=0; j<2;j++){
            double y = rnd.Rannyu();
            sum+=y;
        }
        S2[i]=sum/double(2.);
    }
    
    for (int i=0; i<M; i++){
        S10[i]=0.;
        double sum=0.;
        for (int j=0; j<10;j++){
            double y = rnd.Rannyu();
            sum+=y;
        }
        S10[i]=sum/double(10.);
    }
    
    for (int i=0; i<M; i++){
        S100[i]=0.;
        double sum=0.;
        for (int j=0; j<100;j++){
            double y = rnd.Rannyu();
            sum+=y;
        }
        S100[i]=sum/double(100.);
    }
    
    ofstream uscita("uniform.txt");
    uscita << "#S1   " << "    S2    " << "    S10    " << "    S100" << endl;
    for (int i=0; i<M; i++){
        uscita << S1[i] << "   " << S2[i] << "   " << S10[i] << "   " << S100[i] << endl;
    }
    
    // --- EXPONENTIAL ---
    
    double l=1.;
    
    for (int i=0; i<M; i++) {
        S1[i]=0.;
        S1[i] = rnd.Exp(l);
    }
    
    for (int i=0; i<M; i++){
        S2[i]=0.;
        double sum=0.;
        for (int j=0; j<2;j++){
            double y = rnd.Exp(l);
            sum+=y;
        }
        S2[i]=sum/double(2.);
    }
    
    for (int i=0; i<M; i++){
        S10[i]=0.;
        double sum=0.;
        for (int j=0; j<10;j++){
            double y = rnd.Exp(l);
            sum+=y;
        }
        S10[i]=sum/double(10.);
    }
    
    for (int i=0; i<M; i++){
        S100[i]=0.;
        double sum=0.;
        for (int j=0; j<100;j++){
            double y = rnd.Exp(l);
            sum+=y;
        }
        S100[i]=sum/double(100.);
    }
    
    ofstream uscita2("exponential.txt");
    uscita2 << "#S1   " << "    S2    " << "    S10    " << "    S100" << endl;
    for (int i=0; i<M; i++){
        uscita2 << S1[i] << "   " << S2[i] << "   " << S10[i] << "   " << S100[i] << endl;
    }
    
    // --- CAUCHY ---
    
    double g=1.;
    double mu=0.;
    
    for (int i=0; i<M; i++) {
        S1[i]=0.;
        S1[i] = rnd.Cauchy(g, mu);
    }
    
    for (int i=0; i<M; i++){
        S2[i]=0.;
        double sum=0.;
        for (int j=0; j<2;j++){
            double y = rnd.Cauchy(g, mu);
            sum+=y;
        }
        S2[i]=sum/double(2.);
    }
    
    for (int i=0; i<M; i++){
        S10[i]=0.;
        double sum=0.;
        for (int j=0; j<10;j++){
            double y = rnd.Cauchy(g, mu);
            sum+=y;
        }
        S10[i]=sum/double(10.);
    }
    
    for (int i=0; i<M; i++){
        S100[i]=0.;
        double sum=0.;
        for (int j=0; j<100;j++){
            double y = rnd.Cauchy(g, mu);
            sum+=y;
        }
        S100[i]=sum/double(100.);
    }
    
    ofstream uscita3("cauchy.txt");
    uscita3 << "#S1   " << "    S2    " << "    S10    " << "    S100" << endl;
    for (int i=0; i<M; i++){
        uscita3 << S1[i] << "   " << S2[i] << "   " << S10[i] << "   " << S100[i] << endl;
    }
    
    rnd.SaveSeed();
    return 0;
}
