#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "../ParallelRandomNumberGenerator/random.h"
#include "main.h"

using namespace std;

int main (int argc, char *argv[]){
        
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

    
    // --- 08.1 ---  --- 08.2 ---
    
    Input();
    equilibria_mu_sigma(); // Ho usato il simulated annealing
    
    cout << "mu = " << mu << endl;
    cout << "sigma = " << sigma << endl;
    for(int iblk=1; iblk <= nblk; ++iblk) {
        Reset(iblk);
        for(int istep=1; istep <= nstep; ++istep)
        {
            Move(mu, sigma);
            Fill_istogramma();
            h = H_Psi(x, mu, sigma);
            Accumulate();
        }
        //cout << "Blocco " << iblk << endl;
        //cout << "Accettazione = " << (double)accepted/(double)attempted << endl;
        //cout << "----------" << endl;
        Averages(iblk);
    }
    
    ofstream out;
    out.open("istogramma.dat");
        for (int i=0; i<nbins; i++){
            // la normalizzazione dell'istogramma sarebbe sum(conta[i] * (x_{i+1} - x_{i})), però la differenza fra le x corrisponde a binsize
            // che è costante per tutte le i e la somma dei conta[i] sono i conteggi totali che corrispondono a norm
            out << bin[i] + bin_size/2. << setw(20) << conta[i]/(double)(norm*bin_size) << endl;
        }
    out.close();
    
    delete [] bin;
    delete [] conta;
    
    rnd.SaveSeed();
    return 0;
}



double Psi (double x, double mu, double sigma){
    return exp(-(x-mu)*(x-mu)/(2.*sigma*sigma)) + exp(-(x+mu)*(x+mu)/(2.*sigma*sigma));
}

double Psi_Quadro (double x, double mu, double sigma){
    return Psi(x, mu, sigma)*Psi(x, mu, sigma);
}

double Psi_der2 (double x, double mu, double sigma){
    double S1 = (exp(-(x-mu)*(x-mu)/(2*sigma*sigma)))/(sigma*sigma);
    double S2 = (exp(-(x+mu)*(x+mu)/(2*sigma*sigma)))/(sigma*sigma);
    double S3 = (-1. + (x-mu)*(x-mu)/(sigma*sigma));
    double S4 = (-1. + (x+mu)*(x+mu)/(sigma*sigma));
    return S1*S3 + S2*S4;
}

double H_Psi (double x, double mu, double sigma){
    double A = (-0.5)*Psi_der2(x, mu, sigma) + (pow(x, 4.) - (5./2.)*pow(x, 2.))*(double)Psi(x, mu, sigma);
    return A/Psi(x, mu, sigma);
}

void Input(){
    ifstream ReadInput;
    ReadInput.open("input.dat");
    
    ReadInput >> mu;
    ReadInput >> sigma;
    ReadInput >> x;
    ReadInput >> delta;
    ReadInput >> nblk;
    ReadInput >> nstep;
    
    ReadInput.close();
    
    nbins = 1000;
    sup = 5.;
    inf = -5.;
    bin_size = (sup-inf)/(double)nbins;
    norm =0;
    
    bin = new double [nbins];
    conta = new double [nbins];
    
    for (int i=0; i<nbins; i++){
        bin[i] = inf + i*bin_size;
        conta[i] = 0.;
    }
    
}

void Fill_istogramma(){
    for (int i = 0; i<nbins-1; i++){
        if (bin[i]<=x && bin[i+1]>x){
            conta[i] += 1;
            norm ++;
        }
    }
    if (x>=bin[nbins-1] && x<sup) {
        conta[nbins-1] += 1;
        norm ++;
    }
}


void equilibria_mu_sigma() {
    
    double H, H_new;
    double temp = 25.;
    double energy_new, energy_old, mu_new, sigma_new;
    double beta;
    
    for (int j=0; j<100; j++){
        accepted=0;
        attempted=0;
        beta = 1./temp;
        
        for (int k=0; k<100; k++){
            
            mu_new = mu + rnd.Rannyu(-0.2, 0.2);
            do{
                sigma_new = sigma + rnd.Rannyu(-0.2, 0.2);
            }while(sigma_new<=0.1);
            
            // faccio prima andare un po' la funzione Move affinchè mi campioni i numeri giusti
            for (int p=0; p<10000; p++){
                Move(mu, sigma);
            }
            
            accepted=0;
            attempted=0;
            
            H=0.;
            for (int istep=0; istep<nstep; istep++){
                Move(mu, sigma);
                H += H_Psi(x, mu, sigma);
            }
            energy_old = H/nstep;
            
            // faccio prima andare un po' la funzione Move affinchè mi campioni i numeri giusti
            for (int p=0; p<10000; p++){
                Move(mu, sigma);
            }
            
            accepted=0;
            attempted=0;
            
            H_new=0.;
            for (int istep=0; istep<nstep; istep++){
                Move(mu_new, sigma_new);
                H_new += H_Psi(x, mu_new, sigma_new);
            }
            energy_new = H_new/nstep;
            
            double s = rnd.Rannyu(0.,1.);
            double p = min(1., exp(-(energy_new-energy_old)*(double)(beta)));
            if (s<=p){
                mu = mu_new;
                sigma = sigma_new;
                accepted++;
            }
            attempted ++;
        }
        //cout << "temp = " << temp << endl;
        //cout << "Accettazione = " << (double)accepted/(double)attempted << endl;
        //cout << "----------" << endl;
        temp = temp/1.1;
    }
}

void Move(double mu, double sigma){
    xtest = x + rnd.Rannyu(-delta, delta);
    double A = min(1., Psi_Quadro(xtest, mu, sigma)/Psi_Quadro(x, mu, sigma));
    double s = rnd.Rannyu(0,1);
    if (s<=A) {
        x = xtest;
        accepted++;
    }
    attempted++;
}

void Reset(int iblk) {
    
    if(iblk == 1) {
        H_ave = 0;
        H_av2 = 0;
    }
    
    blk_av = 0;
    blk_norm = 0;
    attempted = 0;
    accepted = 0;
}

void Accumulate(void) {
    blk_av = blk_av + h;
    blk_norm = blk_norm + 1.0;
}

void Averages(int iblk) {
    
    ofstream H, Have;
    H.open ("H_mean.dat", ios::app);
    Have.open("Have.dat", ios::app);
    
    stima_h = blk_av/blk_norm;
    H_ave += stima_h;
    H_av2 += stima_h*stima_h;
    err_h = Error(H_ave, H_av2, iblk);
    H << iblk << setw(20) << H_ave/(double)iblk << setw(20) << err_h << endl;
    
    if (iblk == nblk){
        Have << mu << setw(20) << sigma << setw(20) << H_ave/(double)iblk << setw(20) << err_h << endl;
    }
    
    Have.close();
    H.close();
    
}

double Error(double sum, double sum2, int iblk) {
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}
