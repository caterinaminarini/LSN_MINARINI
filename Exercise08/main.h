#ifndef __main__
#define __main__

Random rnd;
int nblk, nstep, nbins;
double x, xtest, delta, mu, sigma, h;
double sup, inf, bin_size;
int norm;
double* bin;
double* conta;

double H_ave, H_av2, blk_av, blk_norm, stima_h, err_h;
double accepted, attempted;


// funzioni
void Input(void);
double Psi(double, double, double);
double Psi_Quadro(double, double, double);
double Psi_der2(double, double, double);
double H_Psi(double, double, double);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(double, double);
double Error(double,double,int);
void equilibria_mu_sigma(void);
void Fill_istogramma();
#endif
