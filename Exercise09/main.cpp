#include "main.h"

int main (int argc, char *argv[]){

    int seed[4];
    int p1, p2;
    ifstream Primes("../ParallelRandomNumberGenerator/Primes");
    if (Primes.is_open()){
        Primes >> p1 >> p2;
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
    
    // --- PROBLEMA DEL COMMESSO VIAGGIATORE ---
    
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
    
    Inizializza_popolazione(popolazione);
    check_popolazione(popolazione);
    Inizializza_popolazione(popolazione_appo);
    check_popolazione(popolazione_appo);
    elitario.SetParameters(n_città);
    elitario.check();
    elitario.Computelength(x, y);
    
    Individuo Figlio1;
    Figlio1.SetParameters(n_città);
    Figlio1.check();
    
    Individuo Figlio2;
    Figlio2.SetParameters(n_città);
    Figlio2.check();
    
    Compute_length_popolazione(popolazione);
    Ordina(popolazione);
    Imposta_Elitario(popolazione);
    
    for (int i=0; i<n_generazioni; i++){
        popolazione_appo = popolazione;
        
        // faccio crossover cambiando i due membri peggiori della generazione
        crossover(Figlio1, Figlio2, popolazione_appo);
        popolazione[dim_population-1] = Figlio1;
        popolazione[dim_population-2] = Figlio2;
    
        for (int k=0; k<dim_population; k++){
            mutazione=rnd.Rannyu(0,1);
            int individuo=(int)(rnd.Rannyu(0, dim_population));
            if (mutazione>0.1 && mutazione<0.1+p_pair_permutation) {pair_permutation(popolazione[individuo]);}
            if (mutazione>0.3 && mutazione<0.3+p_shift) {shift(popolazione[individuo]);}
            if (mutazione>0.5 && mutazione<0.5+p_permutation) {permutation(popolazione[individuo]);}
            if (mutazione>0.7 && mutazione<0.7+p_inversion) {inversion(popolazione[individuo]);}
        }
       
        check_popolazione(popolazione);
        Compute_length_popolazione(popolazione);
        Ordina(popolazione);
        Imposta_Elitario(popolazione);

        double estrazione_elitario = rnd.Rannyu();
        if (estrazione_elitario<=p_elitario){
            //cout << "inserisco elitario" << endl;
            Inserisci_Elitario(popolazione);
            Compute_length_popolazione(popolazione);
            Ordina(popolazione);
        }
        
        Stampa_migliore(i, popolazione[0]);
        double Best_mean = Media_migliore(popolazione);
        Stampa_media(i, Best_mean);
    }
        
    Stampa_percorso (popolazione[0], x, y);
    
    delete [] x;
    delete [] y;
    
    rnd.SaveSeed();
    return 0;
}


// ------- FUNZIONI -------

void Input(){
    ifstream ReadInput;
    ReadInput.open("input.dat");
    
    ReadInput >> forma;
    ReadInput >> R;
    ReadInput >> n_generazioni;
    ReadInput >> p_pair_permutation;
    ReadInput >> p_shift;
    ReadInput >> p_permutation;
    ReadInput >> p_inversion;
    ReadInput >> p_crossover;
    ReadInput >> p_elitario;
    
    ReadInput.close();
    
    
    
}

void Inizializza_popolazione (vector <Individuo>& popolazione){
    for (int i=0; i<popolazione.size(); i++){
        if (i==0) {
            popolazione[i].SetParameters(n_città);
            popolazione[i].check();
        }
        else {
            popolazione[i].SetParameters(n_città);
            popolazione[i].check();
            check_permutation(popolazione[i-1], popolazione[i]);
        }
    }
}

void Compute_length_popolazione (vector <Individuo>& popolazione){
    for (int i=0; i<popolazione.size(); i++){
        popolazione[i].Computelength(x, y);
    }
}

void check_popolazione (vector <Individuo>& popolazione){
    for (int i=0; i<popolazione.size(); i++) {
        popolazione[i].check();
    }
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

void shift (Individuo& A){
    Individuo appo;
    appo=A;
    //cout << "sto facendo shift" << endl;
    int position = (int)(rnd.Rannyu(1, A.GetDim()));
    rotate (A.Get_start() + position, A.Get_start() + position, A.Get_end());
    check_permutation (appo, A);
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

void swap (Individuo& A, int t, int z){
    int appo1 = A.GetGene(t);
    int appo2 = A.GetGene(z);
    A.SetGene(t, appo2);
    A.SetGene(z, appo1);
}

bool ordine_length (Individuo& i, Individuo& j){
    return (i.Getlength()<j.Getlength());
}

bool ordine_indici (Gene& i, Gene& j){
    return (i.GetIndice()<j.GetIndice());
}

void Ordina (vector <Individuo>& popolazione){
    sort(popolazione.begin(), popolazione.end(), ordine_length);
}

void Ordina_Geni (vector <Gene>& X){
    sort(X.begin(), X.end(), ordine_indici);
}

void check_permutation (Individuo& a, Individuo& b){
    if (not is_permutation(a.Get_start(), a.Get_end(), b.Get_start())){
        cout << "non è una permutazione" << endl;
        exit(EXIT_SUCCESS);
    }
}

/*int selection (vector <Individuo>& popolazione){
    
    // Se prendo p grandi ho più probabilità di ottenere indici piccoli, quindi lunghezze brevi siccome il vettore è ordinato
    
    double p = 3;
    int j = int(popolazione.size()*pow(rnd.Rannyu(), p));
    
    return j;
}*/


int selection (vector <Individuo>& popolazione){
    double N = 0.; //normalizzazione
    for (int i=0; i<popolazione.size(); i++){
        N += 1./(double)(popolazione[i].Getlength());
    }
    
    //probabilità
    probability[0] = (1./popolazione[0].Getlength())/(double)N;
    for (int i=1; i<popolazione.size(); i++){
        probability[i] = (1./(double)(popolazione[i].Getlength()))/(double)N;
    }
    
    //selezione
    int selected = 0;
    double p_estratto = rnd.Rannyu(0., 1.);
    if (p_estratto<=probability[0]) selected = 0;
    for (int i=1; i<popolazione.size(); i++){
        if (p_estratto>probability[i-1] && p_estratto<=probability[i]) selected = i;
    }
    
    return selected;
}

void crossover (Individuo& A, Individuo& B, vector <Individuo>& popolazione){
    //seleziono un padre e una madre
    int padre = selection(popolazione);
    int madre = selection(popolazione);
    
    // Estraggo numero, se è minore della probabilità allora non faccio crossing-over e i due nuovi individui sono il padre e la madre che ho selezionato
    double p_estratto = rnd.Rannyu(0, 1);
    
    if (p_estratto>p_crossover){
        for (int j=0; j<n_città; j++){
            A.SetGene(j, popolazione[padre].GetGene(j));
            B.SetGene(j, popolazione[madre].GetGene(j));
        }
    }
    
    else if (p_estratto<=p_crossover){
        //cout << "sto facendo crossover" << endl;

        int taglio = (int)(rnd.Rannyu(1, n_città));
        //cout << taglio << endl;
        for (int j=0; j<taglio; j++){
            A.SetGene(j, popolazione[padre].GetGene(j));
            B.SetGene(j, popolazione[madre].GetGene(j));
        }
        
        vector <Gene> geni_padre (n_città-taglio);
        vector <Gene> geni_madre (n_città-taglio);
        
        int t = 0;
        for (int j=taglio; j<n_città; j++){
            geni_padre[t].SetCittà(popolazione[padre].GetGene(j));
            geni_madre[t].SetCittà(popolazione[madre].GetGene(j));
            for (int i=0; i<n_città; i++){
                if (popolazione[padre].GetGene(j)==popolazione[madre].GetGene(i)){
                    geni_padre[t].SetIndice(i);
                }
                if (popolazione[madre].GetGene(j)==popolazione[padre].GetGene(i)){
                    geni_madre[t].SetIndice(i);
                }
            }
            t++;
        }
        
        Ordina_Geni(geni_padre);
        Ordina_Geni(geni_madre);
        
        t=0;
        for (int j=taglio; j<n_città; j++){
            A.SetGene(j, geni_padre[t].GetCittà());
            B.SetGene(j, geni_madre[t].GetCittà());
            t++;
        }
    }
}

void Imposta_Elitario (vector <Individuo>& popolazione){
    if(popolazione[0].Getlength()<elitario.Getlength()){
        elitario = popolazione[0];
    }
}

void Inserisci_Elitario (vector <Individuo>& popolazione){
    popolazione[popolazione.size()-1] = elitario;
}

void Stampa_migliore (int generazione, Individuo& I){
    if (forma==1){
        ofstream out ("output_cerchio_best_length.dat",ios::app);
        out << generazione << setw(12) << I.Getlength() << endl;
        out.close();
    }
    else if (forma==0){
        ofstream out ("output_quadrato_best_length.dat",ios::app);
        out << generazione << setw(12) << I.Getlength() << endl;
        out.close();
    }
}

double Media_migliore (vector <Individuo>& popolazione){
    double media = 0.;
    for (int i=0; i<popolazione.size()/2; i++){
        media += popolazione[i].Getlength();
    }
    media = media/(popolazione.size()/2);
    return media;
}

void Stampa_media (int generazione, double media){
    if (forma==1){
        ofstream out ("output_cerchio_mean_length.dat",ios::app);
        out << generazione << setw(12) << media << endl;
        out.close();
    }
    else if (forma==0){
        ofstream out ("output_quadrato_mean_length.dat",ios::app);
        out << generazione << setw(12) << media << endl;
        out.close();
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

