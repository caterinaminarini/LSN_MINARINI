#ifndef __main_h__
#define __main_h__
#include "classi.h"
#include "../ParallelRandomNumberGenerator/random.h"


// Dati del problema
Random rnd;
int n_città=32;
int dim_population=100;
int n_generazioni, forma;
double R;
double *x = new double [n_città];
double *y = new double [n_città];
vector <Individuo> popolazione (dim_population);
vector <Individuo> popolazione_appo (dim_population);
Individuo elitario;
vector <double> probability (popolazione.size());

double p_pair_permutation, p_shift, p_permutation, p_inversion, p_crossover, p_elitario;

// mutation and crossover operators
double mutazione;
void pair_permutation (Individuo&); //permutazione singola
void shift (Individuo&); //shift di un segmento alla fine
void permutation (Individuo&); //permuto due segmenti
void inversion (Individuo&); //inverto l'ordine di un segmento
void swap (Individuo&, int, int); //scambia elementi di dati indici del vettore geni (classe individuo)

void Calcola_Normalizzazione(vector <Individuo>&);
int selection (vector <Individuo>&);
void crossover (Individuo&, Individuo&, vector <Individuo>&);


// ordinamento della popolazione e percorso
bool myordine (Individuo&, Individuo&);
bool ordine_indici (Gene&, Gene&);
void Ordina (vector <Individuo>&);
void Ordina_Geni (vector <Gene>&);

void Input();
void Inizializza_popolazione (vector <Individuo>&);
void Compute_length_popolazione (vector <Individuo>&);
void Inserisci_Elitario (vector <Individuo>&);
void Imposta_Elitario (vector <Individuo>&);
void check_popolazione (vector <Individuo>&);
void check_permutation (Individuo&, Individuo&);

// Stampa
void Stampa_migliore(int, Individuo&);
void Stampa_media (int, double);
double Media_migliore (vector <Individuo>&);
void Stampa_percorso (Individuo&, double *, double *);


#endif
