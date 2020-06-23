In ogni cartella dell'esercizio è presente un file Python in cui ho scritto le relazioni e i risultati degli esercizi e qualora ci fossero più esercizi o i codici forniti dal professore sono presenti anche le cartelle dei singoli, altrimenti è tutto racchiuso in una soltanto.

Ho creato un Makefile per ogni esercizio, in cui si compila facendo "make" e si esegue facendo "./" con il nome dell'eseguibile.
Ho lasciato fuori agli esercizi la cartella del generatore random e ho fatto leggere il generatore random dal Makefile e dai singoli file direttamente da questa cartella. L'esercizio 10.2, che ho eseguito da remoto sul computer del laboratorio, è l'unico che contiene al suo interno tutti i file della cartella del generatore random.

Qualora fossero presente dei file "input.dat" è necessario settarli prima di eseguire. Nell'esercitazione 4 e 7 sono presenti tre file di input a seconda della fase da simulare. È necessario cambiare all'interno del main il punto in cui legge il file di input nelle funzione "Input" commentando quelli che non si intende simulare.

Nell'esercitazione 4 ho salvato in cartelle "solid", "liquid", "gas" i file con i risultati sulle grandezze medie e istantanee e fuori da queste ho lasciato i file in cui ci sono i dati con cui ho studiato l'equilibrazione del sistema.

Nell'esercitazione 6 ho salvato i dati in cui ho studiato l'equilibrazione del sistema e la media all'ultimo blocco delle grandezze. Ci sono due file in cui ho usato Metropolis con h=0 e h=0.02, e due file in cui ho usato Gibbs con h=0  e h=0.02.

Nell'esercitazione 7 ho salvato in cartelle "solid", "liquid", "gas" i file con i risultati sulle grandezze medie e i 500000 valori istantanei che ho usato per l'analisi iniziale. Se si rimandasse la simulazione, il codice farebbe salvare i nuovi valori simultanei fuori da queste cartelle, in modo da mantenere quelli precedenti per l'analisi, mentre aggiornerebbe quelli con le medie e i blocchi.

Nell'esercitazione 8 all'interno del codice del Quantum Monte Carlo ho aggiunto la cartella PIGS e PIMC in cui ho spostato copiandoli i risultati ottenuti dalle simulazioni differenti variando i parametri.

Nelle esercitazioni 9 e 10 il main e le funzioni che ho implementato sono contenuti in "main.cpp" (la dichiarazione delle variabili in "main.h"), mentre le classi che ho utilizzato le ho scritte in "classi.cpp" e "classi.h".

Per le ultime due esercitazioni (11 e 12) ho avuto difficoltà a implementare reti con più layer su Python poiché sul mio computer si bloccava e così le ho svolte da Colab. I tre esercizi dell'esercitazione 11 li ho svolti su tre file di Python, mentre per l'esercizio 12, siccome ho lavorato da Colab, ho dovuto usare dei comandi per poter leggere le immagini delle mie cifre da Google Drive. Ho salvato nella cartella 12 le immagini, ma lo script di Python se eseguito tenta di leggerle da Google Drive.

Nei restanti esercizi dovrebbe essere abbastanza intuitivo leggere i file e negli script di Python ho comunque descritto come ho prodotto i dati.