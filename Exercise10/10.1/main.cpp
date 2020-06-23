#include "main.h"

int main(int argc, char *argv[]){
    
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
    
    //----------------------------------------------------
    
    Input();

    if (forma==1) // Genero città distribuite sulla circonferenza
    {
        for (int i=0; i<n_città; i++){
            double angolo = rnd.Rannyu(0., 2.*M_PI);
            x[i] = R*cos(angolo);
            y[i] = R*sin(angolo);
            //cout << x[i] << setw(20) << y[i] << endl;
        }
    }
    else if (forma==0) // Genero città distribuite nel quadrato
    {
        for (int i=0; i<n_città; i++){
            x[i] = rnd.Rannyu(0., R);
            y[i] = rnd.Rannyu(0., R);
            //cout << x[i] << setw(20) << y[i] << endl;
        }
    }
    
    individuo.SetParameters(n_città);
    
    for (int j=0; j<n_temp; j++){
        
        accepted = 0.;
        attempted = 0.;
        beta = 1./temp;
        
        // Genero un nuovo individuo tramite mutazioni
        for (int i=0; i<n_step; i++){
            Try_pair_permutation(individuo);
            Try_shift(individuo);
            Try_permutation(individuo);
            Try_inversion(individuo);
        }
        
        individuo.Computelength(x,y);
        //cout << "rate = " << (double)(accepted/attempted) << endl;
        Stampa_migliore(j, individuo);
        temp = temp/1.1;
    }
     
    Stampa_percorso(individuo, x, y); //così sto stampando in funzione della temperatura
    //cout << "lunghezza = " << individuo.Getlength() << endl;
    delete [] x;
    delete [] y;
    
    rnd.SaveSeed();
    return 0;
}


// Genero il nuovo individuo tramite delle mutazioni
void Input(){
    ifstream ReadInput;
    ReadInput.open("input.dat");
    ReadInput >> forma;
    ReadInput >> R;
    ReadInput >> temp;
    ReadInput >> n_step;
    ReadInput >> n_temp;
    ReadInput.close();
}

void pair_permutation (Individuo& A){
    Individuo appo;
    appo=A;
    //cout << "sto facendo pair_permutation" << endl;
    int s = (int)(rnd.Rannyu(1, A.GetDim()));
    int t = (int)(rnd.Rannyu(1, A.GetDim()));
    swap(A, s, t);
    check_permutation (appo, A);
}

void Try_pair_permutation (Individuo& A){
    Individuo appo;
    appo=A;
    pair_permutation(appo);
    double s = rnd.Rannyu(0.,1.);
    double p = min(1., exp(-(appo.Computelength(x, y)-A.Computelength(x, y))*(double)(beta)));
    if (s<=p){
        A = appo;
        accepted++;
    }
    attempted ++;
}

void shift (Individuo& A){
    Individuo appo;
    appo=A;
    //cout << "sto facendo shift" << endl;
    int position = (int)(rnd.Rannyu(1, A.GetDim()));
    rotate (A.Get_start() + position, A.Get_start() + position, A.Get_end());
    check_permutation (appo, A);
}

void Try_shift (Individuo& A){
    Individuo appo;
    appo=A;
    shift(appo);
    double s = rnd.Rannyu(0.,1.);
    double p = min(1., exp(-(appo.Computelength(x, y)-A.Computelength(x, y))*(double)(beta)));
    if (s<=p){
        A = appo;
        accepted++;
    }
    attempted ++;
}

void permutation (Individuo& A){
    Individuo appo;
    appo=A;
    //cout << "sto facendo permutation" << endl;
    int position1 = (int)(rnd.Rannyu(1, A.GetDim()/2));
    int position2 = (int)(rnd.Rannyu(A.GetDim()/2, A.GetDim()));
    int lunghezza = (int)(rnd.Rannyu(1, A.GetDim()/4));

    if (lunghezza<(position2-A.GetDim())){
        for (int i=0; i<lunghezza; i++){
            swap(A, position1 + i, position2 + i);
        }
    }
    check_permutation (appo, A);
}

void Try_permutation (Individuo& A){
    Individuo appo;
    appo=A;
    permutation(appo);
    double s = rnd.Rannyu(0.,1.);
    double p = min(1., exp(-(appo.Computelength(x, y)-A.Computelength(x, y))*(double)(beta)));
    if (s<=p){
        A = appo;
        accepted++;
    }
    attempted ++;
}

void inversion (Individuo& A){
    Individuo appo;
    appo=A;
    //cout << "sto facendo inversione" << endl;

    int estremo1 = (int)(rnd.Rannyu(1, A.GetDim()));
    int estremo2 = (int)(rnd.Rannyu(1, A.GetDim()));
    
    if (estremo1<estremo2){
        int lunghezza = estremo2 - estremo1;
        for (int i=0; i<lunghezza/2; i++){
            swap(A, estremo1 + i, estremo2 - i);
        }
    }
    else if (estremo1>estremo2){
        int lunghezza = estremo1 - estremo2;
        for (int i=0; i<lunghezza/2; i++){
            swap(A, estremo2 + i, estremo1 - i);
        }
    }
    check_permutation (appo, A);
}

void Try_inversion (Individuo& A){
    Individuo appo;
    appo=A;
    inversion(appo);
    double s = rnd.Rannyu(0.,1.);
    double p = min(1., exp(-(appo.Computelength(x, y)-A.Computelength(x, y))*(double)(beta)));
    if (s<=p){
        A = appo;
        accepted++;
    }
    attempted ++;
}

void swap (Individuo& A, int t, int z){
    int appo1 = A.GetGene(t);
    int appo2 = A.GetGene(z);
    A.SetGene(t, appo2);
    A.SetGene(z, appo1);
}

void check_permutation (Individuo& a, Individuo& b){
    if (not is_permutation(a.Get_start(), a.Get_end(), b.Get_start())){
        cout << "non è una permutazione" << endl;
        exit(EXIT_SUCCESS);
    }
}

void Stampa_percorso (Individuo& I, double *x, double *y){
    if (forma==1){
        ofstream out ("output_cerchio_best_path.dat",ios::app);
        for (int i=0; i<n_città; i++){
            out << x[I.GetGene(i)] << setw(12) << y[I.GetGene(i)] << endl;
        }
        out.close();
    }
    else if (forma==0){
        ofstream out ("output_quadrato_best_path.dat",ios::app);
        for (int i=0; i<n_città; i++){
            out << x[I.GetGene(i)] << setw(12) << y[I.GetGene(i)] << endl;
        }
        out.close();
    }
}

void Stampa_migliore (int step, Individuo& I){
    if (forma==1){
        ofstream out ("output_cerchio_best_length.dat",ios::app);
        out << step << setw(12) << I.Getlength() << endl;
        out.close();
    }
    else if (forma==0){
        ofstream out ("output_quadrato_best_length.dat",ios::app);
        out << step << setw(12) << I.Getlength() << endl;
        out.close();
    }
}
