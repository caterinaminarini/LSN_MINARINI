#include "classi.h"

Individuo:: Individuo(){
    _dim = 0;
    _length = 0.;
}

Individuo:: ~Individuo(){}

void Individuo:: SetParameters(int dim)
{
    _dim = dim;
    _geni.resize(_dim);
    for (int i=0; i<_dim; i++){
        _geni[i] = i;
    }
    random_shuffle(_geni.begin()+1, _geni.end());
}

void Individuo:: check(){
    int check=0;
    for (int i=0; i<_dim; i++){
        for (int j=0; j<_dim; j++){
            if (_geni[i]==_geni[j] && i!=j) check++;
        }
    }
    if (check!=0) cout << "Non sta rispettando i vincoli!" << endl; //se check = 0, allora rispetta i costrain
}


double Individuo:: Computelength (double *x, double *y){ // Calcola la lunghezza del percorso
    _length = 0.;
    for (int i=0; i<_dim-1; i++){
        _length += pow((x[_geni[i]] - x[_geni[i+1]])*(x[_geni[i]] - x[_geni[i+1]]) + (y[_geni[i]] - y[_geni[i+1]])*(y[_geni[i]] - y[_geni[i+1]]), 0.5);
    }
    _length += pow((x[_geni[_dim-1]] - x[_geni[0]])*(x[_geni[_dim-1]] - x[_geni[0]]) + (y[_geni[_dim-1]] - y[_geni[0]])*(y[_geni[_dim-1]] - y[_geni[0]]), 0.5);
    
    return _length;
}

vector<int>::iterator Individuo::Get_start(){
    return _geni.begin();
}
vector<int>::iterator Individuo::Get_end(){
    return _geni.end();
}
