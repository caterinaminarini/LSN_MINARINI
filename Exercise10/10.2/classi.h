#ifndef __classi_h__
#define __classi_h__
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <ctime>
#include <iomanip>

using namespace std;

class Individuo {

public:
    Individuo();
    ~Individuo();
    
    void SetParameters(int); // Set dim and rnd numbers (percorso casuale)
    void SetGene(int i, int valore) {_geni[i] = valore;}
    void check();
    
    int GetGene(int i) {return _geni[i];}
    int GetDim() {return _dim;}
    double Getlength() {return _length;}
    double Computelength(double*, double*);
    
    vector<int>::iterator Get_start(); //ho bisogno che il tipo di ritorno sia un iteratore per usare rotate
    vector<int>::iterator Get_end();
   
private:
    vector <int> _geni;
    int _dim;
    double _length;

};


class Gene {
    
public:
    Gene();
    ~Gene();
    
    int Getcitta() {return _citta;}
    int GetIndice() {return _indice;}
    
    void Setcitta(int citta) {_citta = citta;}
    void SetIndice(int indice) {_indice = indice;}
    
    
private:
    int _citta;
    int _indice;
};


#endif
