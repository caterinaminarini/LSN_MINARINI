#include "main.h"
#include "mpi.h"

int main (int argc, char *argv[]){

    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status stat;
    
    int seed[4];
    int p1[4], p2[4];
    ifstream Primes("Primes");
    if (Primes.is_open()){
        for (int i=0; i<4; i++){ Primes >> p1[i] >> p2[i];}
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();
    
    ifstream input("seed.in");
    string property;
    if (input.is_open()){
       while ( !input.eof() ){
          input >> property;
          if( property == "RANDOMSEED" ){
             input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
              for (int i=0; i<4; i++){
		if (rank==i) rnd.SetRandom(seed,p1[i], p2[i]);
		}
          }
       }
       input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;
    
    // --- PROBLEMA DEL COMMESSO VIAGGIATORE ---
    
    Input();
    
    // faccio leggere le citta solo da un nodo e poi tramite MPI_Bcast lo faccio leggere a tutti gli altri perchè ogni nodo avendo un seme del generatore diverso sennò mi farebbe leggere citta diverse per ogni nodo
    
    if (rank==0){
        for (int i=0; i<n_citta; i++){
            x[i] = rnd.Rannyu(0., R);
            y[i] = rnd.Rannyu(0., R);
        }
    }
    
    MPI_Bcast(&x[0], n_citta, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&y[0], n_citta, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    Inizializza_popolazione(popolazione);
    check_popolazione(popolazione);
    Inizializza_popolazione(popolazione_appo);
    check_popolazione(popolazione_appo);
    elitario.SetParameters(n_citta);
    elitario.check();
    elitario.Computelength(x, y);
    
    Individuo Figlio1;
    Figlio1.SetParameters(n_citta);
    Figlio1.check();
    
    Individuo Figlio2;
    Figlio2.SetParameters(n_citta);
    Figlio2.check();

    Compute_length_popolazione(popolazione);
    Ordina(popolazione);
    Imposta_Elitario(popolazione);
    
    // variabili per MPI
    int N_migr = 30;
    int itag=1;

        
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
            Inserisci_Elitario(popolazione);
            Compute_length_popolazione(popolazione);
            Ordina(popolazione);
        }
        
        if (i%N_migr==0){
		vector <int> best_a (n_citta);
		vector <int> best_b (n_citta);
		int a, b;
		
            //genero 2 rank casuali. Come per le coordinate, li genero da un rank essendo che i generatori di numeri casuali sono diversi
	    //poi con Bcast trasferisco a tutti gli altri rank
            
		if (rank==0){
			do{
				a = (int)rnd.Rannyu(0, size);
				b = (int)rnd.Rannyu(0, size);
			}while(a==b);
		}

		MPI_Bcast(&a, 1, MPI_INT, 0, MPI_COMM_WORLD);
   		MPI_Bcast(&b, 1, MPI_INT, 0, MPI_COMM_WORLD);

		//cout << a << " " << b << endl;		

                if (rank==a) {
			for(int i=0; i<n_citta; i++) best_a[i]=popolazione[0].GetGene(i);
		}
		if (rank==b) {
			for(int i=0; i<n_citta; i++) best_b[i]=popolazione[0].GetGene(i);
		}
		
		
		//il rank a manda a b
		if (rank==a){
			//cout << "sono nell'if del rank a" << endl;
                	MPI_Send(&best_a[0], n_citta, MPI_INT, b, itag, MPI_COMM_WORLD);
                	MPI_Recv(&best_b[0], n_citta, MPI_INT, b, itag, MPI_COMM_WORLD, &stat);
			for(int i=0; i<n_citta; i++) popolazione[0].SetGene(i, best_b[i]);
            	}

		//il rank b manda a a
            if (rank==b){
		//cout << "sono nell'if del rank b" << endl;
                MPI_Send(&best_b[0], n_citta, MPI_INT, a, itag, MPI_COMM_WORLD);
                MPI_Recv(&best_a[0], n_citta, MPI_INT, a, itag, MPI_COMM_WORLD, &stat);
		for(int i=0; i<n_citta; i++) popolazione[0].SetGene(i, best_a[i]);
	    }


        }
        
        Ordina(popolazione);
        Stampa_migliore(i, popolazione[0], rank);
    }
   

    Stampa_percorso(popolazione[0], x, y, rank);
    
    
    delete [] x;
    delete [] y;
    
    rnd.SaveSeed();
    MPI_Finalize();
    return 0;
}


// ------- FUNZIONI -------

void Input(){
    ifstream ReadInput;
    ReadInput.open("input.dat");
    
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
    for (int i=0; i<(int)popolazione.size(); i++){
        if (i==0) {
            popolazione[i].SetParameters(n_citta);
            popolazione[i].check();
        }
        else {
            popolazione[i].SetParameters(n_citta);
            popolazione[i].check();
            check_permutation(popolazione[i-1], popolazione[i]);
        }
    }
}

void Compute_length_popolazione (vector <Individuo>& popolazione){
    for (int i=0; i<(int)popolazione.size(); i++){
        popolazione[i].Computelength(x, y);
    }
}

void check_popolazione (vector <Individuo>& popolazione){
    for (int i=0; i<(int)popolazione.size(); i++) {
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

bool ordine_length (Individuo i, Individuo j){
    return (i.Getlength()<j.Getlength());
}

bool ordine_indici (Gene i, Gene j){
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


int selection (vector <Individuo>& popolazione){
    double N = 0.; //normalizzazione
    for (int i=0; i<(int)popolazione.size(); i++){
        N += 1./(double)(popolazione[i].Getlength());
    }
    
    //probabilità
    probability[0] = (1./popolazione[0].Getlength())/(double)N;
    for (int i=1; i<(int)popolazione.size(); i++){
        probability[i] = (1./(double)(popolazione[i].Getlength()))/(double)N;
    }
    
    //selezione
    int selected = 0;
    double p_estratto = rnd.Rannyu(0., 1.);
    if (p_estratto<=probability[0]) selected = 0;
    for (int i=1; i<(int)popolazione.size(); i++){
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
        for (int j=0; j<n_citta; j++){
            A.SetGene(j, popolazione[padre].GetGene(j));
            B.SetGene(j, popolazione[madre].GetGene(j));
        }
    }
    
    else if (p_estratto<=p_crossover){
        //cout << "sto facendo crossover" << endl;

        int taglio = (int)(rnd.Rannyu(1, n_citta));
        //cout << taglio << endl;
        for (int j=0; j<taglio; j++){
            A.SetGene(j, popolazione[padre].GetGene(j));
            B.SetGene(j, popolazione[madre].GetGene(j));
        }
        
        vector <Gene> geni_padre (n_citta-taglio);
        vector <Gene> geni_madre (n_citta-taglio);
        
        int t = 0;
        for (int j=taglio; j<n_citta; j++){
            geni_padre[t].Setcitta(popolazione[padre].GetGene(j));
            geni_madre[t].Setcitta(popolazione[madre].GetGene(j));
            for (int i=0; i<n_citta; i++){
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
        for (int j=taglio; j<n_citta; j++){
            A.SetGene(j, geni_padre[t].Getcitta());
            B.SetGene(j, geni_madre[t].Getcitta());
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

void Stampa_migliore (int generazione, Individuo& I, int rank){
	string k;
	if (rank==0) k = "rank0";
	if (rank==1) k = "rank1";
	if (rank==2) k = "rank2";
	if (rank==3) k = "rank3";
		
	ofstream out ("output_best_length_"+k+".dat",ios::app);
	out << generazione << setw(12) << I.Getlength() << endl;
	out.close();
}

void Stampa_percorso (Individuo& I, double *x, double *y, int rank){
	string k;
	if (rank==0) k = "rank0";
	if (rank==1) k = "rank1";
	if (rank==2) k = "rank2";
	if (rank==3) k = "rank3";
	
	ofstream out ("output_best_path_"+k+".dat",ios::app);			
    	for (int i=0; i<n_citta; i++){
       		out << x[I.GetGene(i)] << setw(12) << y[I.GetGene(i)] << endl;
    	}
    	out.close();
}



