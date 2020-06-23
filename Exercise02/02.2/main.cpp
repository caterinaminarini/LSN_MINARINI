#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "../../ParallelRandomNumberGenerator/random.h"

using namespace std;

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
    
    
    int M=10000; // throws
    int N=100; // blocks
    int L = int(M/N);
    int steps=100;  // RW lenght
    int dim=3;  // dimensions
    double r[dim];
    double R[steps][M];
    double ave[steps][N];
    double av2[steps][N];
    double sum_prog[steps][N];
    double su2_prog[steps][N];
    double err_prog[steps];
    double sum[N];
    
    // --- RW 3D LATTICE ---
    
    for (int s=0; s<steps; s++){
        err_prog[s]=0.;
        for (int i=0; i<N; i++){
            ave[s][i]=0.;
            av2[s][i]=0.;
            sum_prog[s][i]=0.;
            su2_prog[s][i]=0.;
        }
    }
    
    // devo fare per ogni passo la media a blocchi e poi plottarla in funzione del passo
    for (int i=0; i<M; i++){
        r[0]=0.;
        r[1]=0.;
        r[2]=0.;
        for (int j=0; j<steps; j++){
            R[j][i]=0.;
            int dir = int(rnd.Rannyu(0, 3));
            double x = rnd.Rannyu();
            if (x>0.5) r[dir] += 1.;
            else r[dir] += -1.;
            R[j][i] = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
        }
    }
    int k=0;
    for (int s=0; s<steps; s++){
        for (int i=0; i<N; i++){
            sum[s]=0.;
            for (int j=0; j<L; j++){
                k = j +i*L;
                sum[s] += R[s][k];
            }
            ave[s][i] = sum[s]/L;
            av2[s][i] = ave[s][i]*ave[s][i];
        }
    }
    
    for (int s=0; s<steps; s++){
        for (int i=0; i<N; i++){
            for (int j=0; j<i+1; j++){
                sum_prog[s][i] += ave[s][j];
                su2_prog[s][i] += av2[s][j];
            }
            sum_prog[s][i] /= double(i+1);
            su2_prog[s][i] /= double(i+1);
        }
        err_prog[s] = error(sum_prog[s][N-1], su2_prog[s][N-1], N-1);
    }
    
    ofstream uscita("results02.2.txt");
    uscita << "#step   " << "    R cubo    " << "    sigmaR cubo" << endl;
    for (int i=0; i<steps; i++){
        double a = i+1;
        double b = pow(sum_prog[i][N-1], 0.5);
        double c = (1./(2.*b))*err_prog[i];
        uscita << a << "       " << b << "       " << c << endl;
    }
    
    
    // --- RW 3D CONTINUUM ---
    
    for (int i=0; i<M; i++){
        r[0]=0.;
        r[1]=0.;
        r[2]=0.;
        for (int j=0; j<steps; j++){
            R[j][i]=0.;
            double theta = acos(1.-2.*rnd.Rannyu());
            double phi = 2.*M_PI*rnd.Rannyu();
            double x = rnd.Rannyu();
            if (x>0.5) {
                r[0] += 1.*sin(theta)*cos(phi);
                r[1] += 1.*sin(theta)*sin(phi);
                r[2] += 1.*cos(theta);
            }
            else {
                r[0] += -1.*sin(theta)*cos(phi);
                r[1] += -1.*sin(theta)*sin(phi);
                r[2] += -1.*cos(theta);
            }
            R[j][i] = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
        }
    }
    
    for (int s=0; s<steps; s++){
        err_prog[s]=0.;
        for (int i=0; i<N; i++){
            ave[s][i]=0.;
            av2[s][i]=0.;
            sum_prog[s][i]=0.;
            su2_prog[s][i]=0.;
        }
    }
    
    k=0;
    for (int s=0; s<steps; s++){
        for (int i=0; i<N; i++){
            sum[s]=0.;
            for (int j=0; j<L; j++){
                k = j +i*L;
                sum[s] += R[s][k];
            }
            ave[s][i] = sum[s]/L;
            av2[s][i] = ave[s][i]*ave[s][i];
        }
    }
    
    for (int s=0; s<steps; s++){
        for (int i=0; i<N; i++){
            for (int j=0; j<i+1; j++){
                sum_prog[s][i] += ave[s][j];
                su2_prog[s][i] += av2[s][j];
            }
            sum_prog[s][i] /= double(i+1);
            su2_prog[s][i] /= double(i+1);
        }
        err_prog[s] = error(sum_prog[s][N-1], su2_prog[s][N-1], N-1);
    }
    
    ofstream uscita1("results02.2_cont.txt");
    uscita1 << "#step   " << "    R cont    " << "    sigmaR cont" << endl;
    for (int i=0; i<steps; i++){
        double a = i+1;
        double b = pow(sum_prog[i][N-1], 0.5);
        double c = (1./(2.*b))*err_prog[i];
        uscita1 << a << "       " << b << "       " << c << endl;
    }
        
    rnd.SaveSeed();
    return 0;
}

double error (double AV, double AV2, int n){
    if (n==0){
        return 0;
    }
    else{
        return pow(((AV2 - AV*AV)/n), 0.5);
    }
}
