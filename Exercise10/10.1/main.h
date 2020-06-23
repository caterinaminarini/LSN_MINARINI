#ifndef __main_h__
#define __main_h__
#include "classi.h"
#include "../../ParallelRandomNumberGenerator/random.h"


// Dati del problema
Random rnd;
int n_città=32;
int forma, n_step, n_temp;
double R, temp, beta, accepted, attempted;
double *x = new double [n_città];
double *y = new double [n_città];
Individuo individuo;

void Input();

// mutation operators
void pair_permutation (Individuo&); //permutazione singola
void Try_pair_permutation (Individuo&);
void shift (Individuo&); //shift di un segmento alla fine
void Try_shift (Individuo&);
void permutation (Individuo&); //permuto due segmenti
void Try_permutation (Individuo&);
void inversion (Individuo&); //inverto l'ordine di un segmento
void Try_inversion (Individuo&);
void swap (Individuo&, int, int); //scambia elementi di dati indici del vettore geni (classe individuo)
void check_permutation (Individuo&, Individuo&);

// Stampa
void Stampa_migliore(int, Individuo&);
void Stampa_percorso (Individuo&, double *, double *);

#endif

