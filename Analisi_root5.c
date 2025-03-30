{
#include <iostream>
#include<sstream> 
#include <memory>
#include <cmath>
#include <vector>
#include <TGraph.h>
#include <math.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>

#include "TList.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TRandom3.h"
#include "TTree.h"
//using namespace std;
//void Calibrazione() {

//  _RAW_285      cosmici   10 V
//  _RAW_286      cosmici   50 V
//  _RAW_287 150 uA 45 kV   10 V
//  _RAW_289   5 uA 45 kV   10 V
//  _RAW_291  10 uA 45 kV   10 V    BAD
//  _RAW_292  15 uA 45 kV   10 V
//  _RAW_293  20 uA 45 kV   10 V
//  _RAW_294  20 uA 45 kV   10 V    BAD
//  _RAW_295  30 uA 45 kV   10 V
//  _RAW_296  40 uA 45 kV   10 V
//  _RAW_298  75 uA 45 kV   10 V
//  _RAW-299 150 uA 45 kV   10 V
//  _RAW_300 150 uA 45 kV   50 V
//  _RAW_301      cosmici   50 V

TF1 *doublegaus = new TF1("doublegaus","gaus(0)+gaus(3)", -1., 1.);
doublegaus->SetParameter(0,20.);
doublegaus->SetParameter(1,15.);
doublegaus->SetParameter(2,2.0);
doublegaus->SetParameter(3,20.);
doublegaus->SetParameter(4,80.);
doublegaus->SetParameter(5,3.0);

TF1 *langaus = new TF1("langaus","landau(0)+gaus(2)", -1., 1.);
langaus->SetParameter(0,20.);
langaus->SetParameter(1,6.);
langaus->SetParameter(2,50.0);
langaus->SetParameter(3,20.);
langaus->SetParameter(4,2.);

TF1 *langaus1 = new TF1("langaus","landau(0)*gaus(2)", -1., 1.);
langaus->SetParameter(0,20.);
langaus->SetParameter(1,6.);
langaus->SetParameter(2,50.0);
langaus->SetParameter(3,20.);
langaus->SetParameter(4,2.);

// TFile è un file di tipo ROOT

// File dove scrivo istogrammi e grafici
   TFile *f = new TFile("C:/root_v5.34.38/macros/New_DAQ_foot_3172_results_1000.root", "RECREATE");  
// File che leggo
   TFile *file = new TFile("C:/root_v5.34.38/macros/foot_3172.root");

// dichiaro un indicatore per un file con etichetta "outfile"
   FILE *outfile;
// dichiaro una variaile carattere con 700 caratteri    
   char fileout[700] = "C:/root_v5.34.38/macros/output.txt";  //file txt che genero di output
// estraggo un oggetto "TTree" chiamato "raw_events" dal file ROOT precedentemente estratto 
  TTree * raw_events = (TTree*)file->Get("raw_events"); //cerca un oggetto con il nome "raw_events" 
// questo albero contiene i dati grezzi degli eventi che verranno analizzati successivamente
  Float_t bias = 50.0;      // sensor bias [V] : bias applicato al sensore
  Float_t Angle = 74.5;		// source angle [°] : angolo della sorgente rispetto al sensore
//
//Inizializzo i vettori dove mettere le letture ADC
//  std::vector<unsigned short> *raw_event = new std::vector<unsigned short>();
  vector<unsigned short> *raw_event = new vector<unsigned short>(); //?da capire meglio a cosa serve crearne uno nuovo
  //vector rappresenta un contenitore dinamico per meorizzare una sequenza di elementi. 
  //unsigned short dice il tipo di elementi che accetta il vettore ,ovvero numeri interi senza segno a 16 bit.
  //*raw_event è un puntatore che fa riferimento al vettore
  //new è un operatore per allocare un nuovo vettore  

  //Associo ogni TBranch a un vettore
  raw_events->SetBranchAddress("RAW Event", &raw_event); //associo il ramo (branch) chiamato "RAW Event" del TTree "raw_events" al puntatore "raw_event".
  //SetBranchAddress è una funzione di ROOT che collega un ramo di un albero a una variabile o a un puntatore
  //Il ramo "RAW Event" contiene i dati grezzi degli eventi , che verranno memorizzati nel vettore raw_event di tipo vector<unsigned short>

  Int_t nentries;
  nentries = raw_events->GetEntries(); //Funzione di ROOT che restituisce il numero totale di enrtrate/eventi contenute nell'albero raw_events
  //Ogni "entrata" corrisponde a un evento registrato nel file ROOT, che può contenere dati come letture ADC,segnali o altre informazioni sperimentali
  Int_t n_pedestal = 500; //numero di eventi da utilizzare per il calcolo del piedistallo
  Int_t n_events = nentries; //numero totale di eventi da analizzare
  n_events = 1000; 
  printf(" n pedestal %d n events %d \n",n_pedestal,n_events); 
  Int_t event_offset = 0; 
  Int_t ped_event_offset = 0;
  Int_t n_chip = 10; //numero di chip
  Int_t n_chip_chan = 64; //numero di canali per chip
  Int_t n_strip_sensor = n_chip * n_chip_chan; //numero totale di strip del sensore 
  Double_t strip_signal[n_strip_sensor]; //array per memorizzare il segnale di ogni strip
  // Double_t è un tipo di dato definito da ROOT che rappresenta un numero in virgola mobile a doppia precisione(64 bit)
  //strip_signal è il nome dell'array dichiarato
  Int_t event_to_display = 535;// evento da visualizzare     
  Double_t threshold = 10.; // soglia per definire inizio di un cluster 
  Double_t common_thre = 10.; //soglia per definire il common mode
  Int_t common_mode_choice = 0; //scelta del common mode
  Int_t n_max_in_event = 60; //numero massimo di cluster in un evento
  Int_t max_cluster_width = 640; //numero massimo di strip in un cluster
  Int_t good_event[n_events]; // array per memorizzare gli eventi buoni
  Double_t local_max[n_events][n_max_in_event][6]; //array tridimensionale per memorizzare i minimi locali
  //6 è il numero di parametri che voglio memorizzare per ogni cluster:
  //1: inizio del cluster (indice della strip iniziale)
  //2: fine del cluster (indice della strip finale)
  //3: massimo del cluster (segnale massimo rilevato nel cluster)
  //4: posizione del massimo (indice della strip con il massimo segnale)
  //5: larghezza del cluster (numero di strip nel cluster)
  //6: integrale del cluster (somma dei segnali delle strip nel cluster)
  
  Double_t cluster_low_threshold = 10.;  //soglia minima per definire la appartenenza di una strip ad un cluster
  Int_t i,j,k,j_index,l; 
  // i indice per iterare sugli eventi;
  // j indice per iterare sulle strip;
  // k indice per iterare sui chip ; 
  // j_index indice globale di una strip, calcolato combinando k e j;
  // l indice generico per iterare su cluster o strip in un intervallo specifico.

  
  char histo_title[3000]; // array di caratteri per memorizzare il titolo degli istogrammi 
  char histo_base[50]; //array di caratteri per memorizzare il suffisso del titolo degli istogrammi
  //Definisco una serie di istogrammi ROOT
  TH1F *hRMS = new TH1F("RMS","RMS",100,0.,10.);
  TH1F *hMean = new TH1F("Mean","Mean",75,100.,400.);
  TH1F *hcm100 = new TH1F("CM100","CM100",20,-0.5,0.5); //Istogramma per rappresentare il common mode(rumore comune) con un intervallo ristretto di 100 eventi
  // "CM100" è il nome dell'istogramma,"CM100" è anche il titolo dell'istogramma, 20 è il numero di bin e [-0.5,0.5] è l'intervallo dell'asse x
  // Scopo: analizzare il rumore comune in un intervallo ristretto di eventi per identificare eventuali anomalie o rumori indesiderati.
  TH1F *hcm_chip[n_chip]; // array di istogrammi per rappresentare il common mode calcolato per ogni chip
  // Scopo: analizzare il rumore comune separatamente per ciascun chip.
  TH1F *hcm_prof_chip[n_chip]; // array di istogrammi(uno per ogni chip) per rappresentare il profilo del common mode in funzione del numero di eventi
  // Scopo: analizzare l'andamento del rumore comune nel tempo per ciascun chip.
  TH1F *common_ADC_hist[2]; // array di istogrammi per rappresentare il common mode calcolato per ogni ADC
  // Scopo: analizzare il rumore comune separatamente per ciascun ADC.
  TH1F *common_chip[n_chip]; // array di istogrammi per rappresentare il common mode calcolato per ogni chip
//// Scopo: ?RIDONDANZA con hcm_chip?
  TH1F *strip_value[n_strip_sensor]; // array di istogrammi per rappresentare il segnale di ogni strip
  // Scopo: analizzare il segnale di ogni strip separatamente.
  TH1F *common_ADC_1_distr = new TH1F("Common noise ADC 1","CM ADC 1 distribution",400,-25.,25.); // istogramma per rappresentare il common mode calcolato per il primo ADC
  // Scopo: analizzare il rumore comune separatamente per il primo ADC.
  TH1F *common_ADC_2_distr = new TH1F("Common noise ADC 2","CM ADC 2 distribution",400,-25.,25.); // istogramma per rappresentare il common mode calcolato per il secondo ADC
  // Scopo: analizzare il rumore comune separatamente per il secondo ADC.

  // dichiaro array di istogrammi per il piedistallo
  for(k=0; k<n_chip; k++) // sto ciclando sui 10 chip
  // quindi sto ciclando su tutto il sensore
   { TH1F *hcm_chip[k] = new TH1F("CM x chip k","CM x chip k",100,-20,20); // e per ogni chip creo un istogramma per il common mode
    TH1F *hcm_prof_chip[k] = new TH1F("CM profile x chip k","CM profile x chip k",n_events,0,float(n_events));// e per ogni chip creo un profilo del common mode in funzione del numero di eventi
    }
  // Scopo: analizzare il rumore comune separatamente per ciascun chip e nel tempo.
  TH1F *n_strip_common_distr = new TH1F("CM","CM",70,0.,70.); 
  // istogramma per rappresentare il numero di strip utilizzate per calcolare il common mode con 70 bin e intervallo da 0 a 70  
  TH1F *n_total_strip_common_distr = new TH1F("total strip CM","total strip CM",700,0.,700.);
  // istogramma per rappresentare il numero totale di strip utilizzate per calcolare il common mode con 700 bin e intervallo da 0 a 700
  TH1F *hcm = new TH1F("h4","CM",100,-20,20);
  // istogramma per rappresentare il common mode con 100 bin e intervallo da -20 a 20
  TH1F *signal_event = new TH1F("signal_event","signal_event",n_events,0,float(n_events));
  // istogramma per rappresentare il segnale dell'evento con n_events bin e intervallo da 0 a n_events
  TH1F *rms_signal_distr = new TH1F("RMS signal distribution","Event RMS signal distribution",1000,0.,200.);
  // istogramma per rappresentare la distribuzione del segnale RMS degli eventi con 1000 bin e intervallo da 0 a 200
  TH1F *common_tot_rms_distr = new TH1F("RMS of common mode distribution","Event RMS signal distribution",1000,0.,100.);
  // istogramma per rappresentare la distribuzione del rumore comune RMS degli eventi con 1000 bin e intervallo da 0 a 100
  TH1F *strip_total_distr = new TH1F("Strip total distribution","Strip total distribution",2500,-50.,200.);
  // istogramma per rappresentare la distribuzione totale delle strip con 2500 bin e intervallo da -50 a 200
  TH1F *strip_occupancy = new TH1F("Strip occupancy","Max strip occupancy",n_strip_sensor,0.,float(n_strip_sensor)); 
  TH1F *strip_occupancy_max = new TH1F("Max strip occupancy","Strip occupancy",n_strip_sensor,0.,float(n_strip_sensor));
  TH1F *n_cluster_profile = new TH1F("Number of clusters vs event number","Number of clusters vs event number",n_events,0.,float(n_events));
  TH1F *n_cluster_distr = new TH1F("Distribuzioone numero di cluster","Number of clusters distributionr",100,0.,100.);
  TH1F *n_strip_over_threshold_distr = new TH1F("N strip over threshold","Number of strip over threshold distributionr",640,0.,640.);
  TH1F *cluster_width_distr = new TH1F("Cluster width distribution","Cluster width distribution",64,0.,64.);
  TH1F *cluster_strip_signal_distr = new TH1F("Cluster strip signal distribution","Cluster strip signal distribution",20000,0.,200.);
  TH1F *cluster_integral_signal_distr = new TH1F("Cluster width distribution","Cluster width distribution",2000,0.,500.);
  TProfile *hcm_prof = new TProfile("CM profile","CM profile",n_events,0,float(n_events),-1000.,1000.);
  TH1F *hfit = new TH1F("h5","par1fit",10000,-1,1);
  TH1F *hchi = new TH1F("h6","chi_quadro",100,-10,10);
  TH2 *common_tot_vs_rms = new TH2F("Total common noise vs its rms","",400,-20.,20.,400,0.,20.);
  //TCanvas * c = new TCanvas("c", "c", 1000, 1000); // crea un canvas per visualizzare gli istogrammi

  // dichiaro array di istogrammi per il piedistallo
  TH1 *h[n_strip_sensor]; // array di istogrammi per il piedistallo
  // h è un array di puntatori a istogrammi
  // ogni elemento dell'array h è un puntatore a un istogramma di tipo TH1
  // TH1 è una classe di ROOT che rappresenta un istogramma unidimensionale
  TH1 *h_prof[n_strip_sensor];
  for (Int_t k=0; k<n_strip_sensor; k++) // sto ciclando su tutto il sensore
   { h[k] = new TH1F("histo","piedistallo",350,50.,400.);  
    h_prof[k] = new TH1F("histo","piedistallo profile",n_pedestal,0.,float(n_pedestal)); 
    strip_value[k] = new TH1F("strip signal","strip signal",800,-20.,20.);
   }

  // CALCOLO PIEDISTALLO
  // cicli dove riempio 640 istogrammi e 640 grafici uno per ogni strip
  unsigned short adc[n_strip_sensor]; // array per memorizzare i valori ADC
  // adc è un array di tipo unsigned short che memorizza i valori letti da ogni strip
  for (Int_t i=ped_event_offset;i<n_pedestal;i++) // ciclo sugli eventi per il calcolo del piedistallo
   { raw_events->GetEntry(i); // leggo l'evento i dal file ROOT
    for (Int_t j=0; j<n_strip_sensor; j++) // ciclo su ogni strip
     {
      adc[j] = (*raw_event)[j]; // leggo il valore ADC per ogni strip
      h[j]->Fill(adc[j]); // riempio l'istogramma h[j] con il valore ADC
      h_prof[j]->Fill(float(i),adc[j]); // riempio il profilo h_prof[j] con l'evento e il valore ADC
     }
    }
  // salvo il valor medio del piedistallo per ogni strip
  Float_t Mean[n_strip_sensor];// array per memorizzare la media del piedistallo
  // Mean è un array di tipo Float_t che memorizza la media del piedistallo per ogni strip
  // Float_t è un tipo di dato definito da ROOT che rappresenta un numero in virgola mobile a singola precisione(32 bit)
  Float_t Noise[n_strip_sensor];
  // array per memorizzare il rumore del piedistallo
  // Noise è un array di tipo Float_t che memorizza il rumore del piedistallo per ogni strip
  for (j=0; j<n_strip_sensor; j++) 
   {
    Mean[j]=h[j]->GetMean(); // calcolo la media del piedistallo per ogni strip
    // GetMean() è una funzione di ROOT che restituisce la media dei valori contenuti nell'istogramma
    // h[j]->GetMean() calcola la media dei valori contenuti nell'istogramma h[j]
    // e la memorizza nell'array Mean[j]
    // h[j] è un puntatore all'istogramma j-esimo  
    Noise[j]=h[j]->GetRMS(); // calcolo il rumore del piedistallo per ogni strip 
    // h[j]->GetRMS() calcola la deviazione standard dei valori contenuti nell'istogramma h[j]
    // e la memorizza nell'array Noise[j]
    } 
//
// definizioni variabili di singolo evento
   Float_t signal[n_strip_sensor]; // array per memorizzare il segnale di ogni strip
  // signal è un array di tipo Float_t che memorizza il segnale per ogni strip
  // Float_t è un tipo di dato definito da ROOT che rappresenta un numero in virgola mobile a singola precisione(32 bit)

// Calcolo common noise
// grandezze calcolate:  
// 0: Common mode per tutto il sensore (1)
// 1: Common mode per ogni ADC   (2)
// 2: Common mode per ogni chip  (10)

  Int_t ev_to_display = 0;      // evento da visualizzare
  Float_t common_tot;           // array dove inserire il common mode unico [0]
  Float_t common_tot2;          // array dove inserire il common mode unico [0]
  Int_t n_total_strip_common;	 	// numero di strip usate per il common mode unico
  Float_t common_ADC[2];				// array dove inserire il common mode per ogni ADC [1]
  Int_t n_ADC_strip_common[2];  // numero di strip usate per il common mode per ogni ADC
  Float_t common[n_chip];				// array dove inserire il common mode per ogni chip [2]
  Int_t n_strip_common[n_chip]; // numero di strip usate per il common mode per ogni chip
 
 //
  for ( Int_t i=event_offset; i<n_events; i++)
   { //printf(" dentro loop eventi \n");
    raw_events->GetEntry(i);
	  good_event[i] = 1;
    common_tot = 0.;
    common_tot2 = 0.;
	  n_total_strip_common = n_strip_sensor;// numero totale di strip usate per il common mode unico
    common_ADC[0] = 0.;// common mode per il primo ADC
	  n_ADC_strip_common[0] = n_chip_chan*5;// numero di strip usate per il common mode per il primo ADC 
    common_ADC[1] = 0.;// common mode per il secondo ADC
	  n_ADC_strip_common[1] = n_chip_chan*5;// numero di strip usate per il common mode per il secondo ADC 
     
// creo un array di caratteri per il titolo dell'istogramma
// event display histogram  
	  strcpy(histo_title,"event_display_"); // strcpy è una funzione che copia una stringa in un'altra
    stringstream histo_suffix; // creo un oggetto stringstream per convertire l'intero i in una stringa
    histo_suffix << i; // inserisco l'evento i nella stringa
	  histo_suffix >> histo_base; // converti l'intero i in una stringa
	  strcat(histo_title,histo_base); // concatena la stringa "event_display_" con il numero dell'evento i
    TH1F *ev_display = new TH1F(histo_title,histo_title,n_strip_sensor,0.,float(n_strip_sensor)); // creo un istogramma 1D per visualizzare il segnale di ogni strip
    TH1F *ev_raw_display = new TH1F(histo_title,histo_title,n_strip_sensor,0.,float(n_strip_sensor)); // creo un istogramma 1D per visualizzare il segnale grezzo di ogni strip
//	printf(" evento %d prima chip \n",i);

    for( k=0; k<n_chip; k++)
     { common[k] = 0;// inizializzo il common mode per ogni chip
	     n_strip_common[k] = n_chip_chan; 
       for (j=0; j<n_chip_chan; j++)
	      { j_index = j + k*n_chip_chan;
//	printf(" i %d chip %d strip %d strip ID %d \n",i,k,j,j_index);
          ev_raw_display->Fill(float(j_index),signal[j_index]);
		      signal[j_index] = ((*raw_event)[j_index]);
          ev_raw_display->Fill(float(j_index),signal[j_index]);
		      signal[j_index] = signal[j_index]-Mean[j_index];
          if(signal[j_index] < common_thre) 
		       { common_tot = common_tot+signal[j_index];
		         common_tot2 = common_tot2+signal[j_index]**2;
//		  printf(" i %d j_index %d signal %f common %f common2 %f \n",i,j_index,signal[j_index],common_tot,common_tot2);
		         common[k] = common[k] + signal[j_index];
		         if(k < 5 ) { common_ADC[0] = common_ADC[0] + signal[j_index]; }
		         else { common_ADC[1] = common_ADC[1] + signal[j_index]; }
		        }
		      else
		       { n_strip_common[k] = n_strip_common[k] -1;
		         n_total_strip_common = n_total_strip_common -1;
		      if(k < 5 )   { n_ADC_strip_common[0] = n_ADC_strip_common[0] -1;}
		      else {n_ADC_strip_common[1] = n_ADC_strip_common[0] -1;}
            }
        } // fine ciclo strip

	     if(n_strip_common[k] > 0) 
	      { if(n_strip_common[k] < (n_chip_chan/2)) 
	         { common[k] = common[k] /n_strip_common[k]; 
		         good_event[i] = 0;
		        }
	        else
		       { common[k] = common[k] /n_strip_common[k];
 	           hcm_chip[k]->Fill(common[k]);
	           hcm_prof_chip[k]->Fill(float(i),common[k]);
	          }
	       } 
	    }// fine ciclo chip

	  common_ADC[0] = common_ADC[0]/n_ADC_strip_common[0];
	  common_ADC[1] = common_ADC[1]/n_ADC_strip_common[1];
	  common_tot = common_tot/n_total_strip_common;
	  common_tot_rms = sqrt(common_tot2/n_total_strip_common - common_tot**2);
    n_strip_common_distr->Fill(float(n_total_strip_common));
    hcm->Fill(common_tot);
    hcm_prof->Fill(float(i),common_tot);
    common_tot_rms_distr->Fill(common_tot_rms);
	  common_tot_vs_rms->Fill(common_tot,common_tot_rms);
    
	  printf(" i %d common %f common_rms %f n_strip %d \n",i,common_tot,common_tot_rms,n_total_strip_common);
    //i = indice dell'evento
    //common_tot = valore del common mode calcolato per l'evento i
    //common_tot_rms = valore del rumore comune calcolato per l'evento i
    //n_total_strip_common = numero totale di strip usate per calcolare il common mode

	  common_ADC_1_distr->Fill(common_ADC[0]); // riempio l'istogramma con il valore del common mode per il primo ADC
	  common_ADC_2_distr->Fill(common_ADC[1]); // riempio l'istogramma con il valore del common mode per il secondo ADC
//
// start event processing
//
    if(good_event > 0)// se l'evento è buono 
	   { if(common_tot_rms < 6.) //se ha una certa qualità, altrimenti lo buttiamo
	    { Int_t i_local_max = 0;  
	      Int_t in_cluster = 0;  	// flag per indicare che stiamo processando un cluster
     // Tutti azzeramenti dei parametri del cluster: 
	      local_max[i][i_local_max][0] = 0.;							//  start cluster 
	      local_max[i][i_local_max][1] = 0.;    					//  end cluster 
	      local_max[i][i_local_max][2] = 0.;							//  cluster maximum
	      local_max[i][i_local_max][3] = 0.;							//  position of max in cluster
	      local_max[i][i_local_max][4] = 0.;    					//  cluster width 
	      local_max[i][i_local_max][5] = 0.;  						//  integral of cluster signal
	      Int_t n_strip_over_threshold = 0;  //  number of strip over threshold
	      Double_t ev_signal = 0.;									//  Calcolo media e RMS signal per identificare bad events
	      Double_t ev_signal2 = 0.;									//  Calcolo media e RMS signal per identificare bad events
	      Double_t rms_signal = 0.;									//  Calcolo media e RMS signal per identificare bad events

	      for(k=0; k<max_cluster_width; k++) { strip_signal[k] = 0.; } // azzero il segnale delle strip 
        for( k=0; k<n_chip; k++)// ciclo sui chip 
         { for (j=0; j<n_chip_chan; j++)// ciclo sui canali
	          { j_index = j + k*n_chip_chan; //j_index va da 0 a 640, quindi sto scorrendo su tutte le strip del sensore
		          signal[j_index] = (signal[j_index]-common[k])*good_event[i]; //rimuovo il common mode dal segnale della strip filtrando eventi cattivi azzerando il segnale se l'evento è cattivo
			        strip_value[j_index]->Fill(signal[j_index]);
		          ev_display->Fill(float(j_index),signal[j_index]);//riempio istogramma con indice di strip e valore del segnale sulla strip 
              //Scopo : rappresentare segnale di strip ed ogni bin corrisponde ad una strip e il valore del bin è il segnale della strip
		          ev_signal = ev_signal + signal[j_index]; //accumula la somma totale dei segnali di tutte le strip per l'evento corrente
              //Scopo : calcolare media del segnale
		          ev_signal2 = ev_signal2 + signal[j_index]**2;	 //accumula la somma dei quadrati	dei segnali di tutte le strip per l'evento corrente 
              //Scopo : calcolare deviazione standard del segnale
	  	        strip_total_distr->Fill(signal[j_index]); //riempio istogramma per rappresentare la distribuzione totale dei segnali delle strip

           // printf(" evento %d k %d j %d indxex %d prima massimo relativo \n",i,k,j,j_index);
           // i = indice dell'evento
           // k = indice del chip
           // j = indice della strip
           // j_index = indice globale della strip (k*64+j) 

// look for relative maximum
              if(signal[j_index] > cluster_low_threshold) //inizio il cluster 
		           { if(in_cluster < 1) {local_max[i][i_local_max][0] = float(j_index);} // start cluster
			           if(in_cluster < max_cluster_width) //sempre vero in quanto in_cluster < 640
                  { if(signal[j_index] > cluster_low_threshold || signal[j_index] < -cluster_low_threshold)// verifico se  il valore del segnale di una strip è sopra la soglia positiva o sotto la soglia negativa 
                    {n_strip_over_threshold = n_strip_over_threshold + 1;}//incremento contatore di strip sopra soglia
		                in_cluster = in_cluster +1; //incremento contatore di strip che fanno parte del cluster
                    // in_cluster = numero di strip che fanno parte del cluster
			              printf(" evento %d chip %d canale %d strip %d dentro loop massimo relativo \n",i,k,j,j_index);
                    // i = indice dell'evento
                    // k = indice del chip
                    // j = indice della strip sul canale del chip
                    // j_index = indice globale della strip (k*64+j)
                    
		                if(signal[j_index] > local_max[i][i_local_max][2]) //se il segnale di strip è maggiore dell'attuale massimo locale
		                 { local_max[i][i_local_max][3] = float(j_index); //allora memorizzo la posizione del massimo 
		                   local_max[i][i_local_max][2] = signal[j_index]; //allora memorizzo il nuovo massimo locale
                       //quindi io ho spostato il valore del massimo e la posizione del massimo.

                    // printf(" inside cluster j_index %d signal %5.2f local_max %5.2f position %5.0f \n",j_index,signal[j_index],local_max[i][i_local_max][2],local_max[i][i_local_max][3]);
		                  }
		                   local_max[i][i_local_max][1] = float(j_index);   // end cluster 
		                   local_max[i][i_local_max][4] = local_max[i][i_local_max][1]-local_max[i][i_local_max][0]+1;   // cluster width = end - start +1
		                   local_max[i][i_local_max][5] = local_max[i][i_local_max][5]+signal[j_index];   // integral of cluster signal = sum of signals in cluster
		                   cluster_strip_signal_distr->Fill(signal[j_index]); //costruisco un istogramma per rappresentare la distribuzione del segnale delle strip del cluster
	                     strip_occupancy->Fill(float(j_index)); //riempio l'istogramma strip_occupancy con il valore del j_index

                    // printf(" evento %d strip %d # cluster %d signal %5.1f pos %5.1f local max %5.1f cluster width %5.0f \n",i,j_index,i_local_max,signal[j_index],local_max[i][i_local_max][2],local_max[i][i_local_max][3],in_cluster);
			             }
	 	            }
		          else //altrimenti non inizio il cluster 
		           { if(in_cluster > 0 ) 
		              { in_cluster = 0;      		  // dichiaro chiuso quel cluster 
		                Int_t l1 = 0;
		                strip_occupancy_max->Fill(local_max[i][i_local_max][3]); //riempio l'istogramma strip_occupancy_max con il valore del massimo locale

                 // fprintf(outfile," %d  %d %d %d ",i,i_local_max+1,int(local_max[i][i_local_max][0])-2,int(local_max[i][i_local_max][1])+2);

                 // Specie di salvaguardia:
                    if(local_max[i][i_local_max][0] < 2) { Int_t lower_cluster = 0;} // se lo start < 2, cioè vale 0 o 1, allora il lower_cluster lo faccio partire da 0
                    else {Int_t lower_cluster = int(local_max[i][i_local_max][0])-2;} // altrimenti prendo lo start meno 2
                    if(local_max[i][i_local_max][1] > n_strip_sensor-2) { Int_t upper_cluster = 640;} // se l'end > 638, allora lo faccio finire a 640
                    else {Int_t upper_cluster = int(local_max[i][i_local_max][1])+2;} // altrimenti prendo l'end più 2
			              for(l=lower_cluster; l<upper_cluster; l++) 	
			               { // fprintf(outfile,"%.2f ",signal[l]);	  
			                strip_signal[l1] = signal[l];
			                l1 = l1 + 1;
		   	              }
                       // fprintf(outfile,"\n");
                     
			                ev_number = i;
			                n_cluster = i_local_max;
			                start_cluster = lower_cluster;
			                stop_cluster = upper_cluster;
			                cluster_width = upper_cluster - lower_cluster +1;
			                common_noise = common_tot;
			                common_noise_rms = common_tot_rms;
                   // t_ridotto.Fill();
		                  for(l1=0; l1<cluster_width; l1++) { strip_signal[l1] = 0.; } 
                   // local_max[i][i_local_max][3] = float(j_index-1);
			             // Int_t test_index = 0;
			             // printf(" %d  %d %d %d ",i,i_local_max+1,int(local_max[i][i_local_max][0])-2,int(local_max[i][i_local_max][1])+2);
	                 // printf(" evento %d j_index %d numero cluster in frame %d start %5.0f stop %5.0f cluster max pos %5.0f max value %5.1f \n",i,j_index,i_local_max+1,local_max[i][i_local_max][0],local_max[i][i_local_max][1],local_max[i][i_local_max][3],local_max[i][i_local_max][2]); 
		               // for(test_index=(int(local_max[i][i_local_max][0])-2); test_index<(int(local_max[i][i_local_max][1])+3); test_index++)
			             // { printf(" %d  %5.1f \n",test_index,signal[test_index]);}
		                  i_local_max = i_local_max + 1;    // incremento contatore per il prossimo cluster
		               }
		            }
		        }
	        }
	    }
 	   if(i_local_max > 0 && i_local_max < 8)
	    { cluster_width_distr->Fill(local_max[i][i_local_max-1][4]);
	      cluster_integral_signal_distr->Fill(local_max[i][i_local_max-1][5]);
        if(i_local_max > 5.) { printf( "   evento troppi cluster %d n strip over threshold %d \n",i,n_strip_over_threshold);}
	     }
	   if(good_event > 0) 
	    { ev_signal = ev_signal/n_strip_sensor;
	      rms_signal = sqrt(ev_signal2/n_strip_sensor - ev_signal**2);
	      rms_signal_distr->Fill(rms_signal);
	      n_cluster_profile->Fill(float(i),float(i_local_max));   // profilo # cluster nel frame
	      n_cluster_distr->Fill(float(i_local_max));   // distribuzione # cluster nel frame
	      n_strip_over_threshold_distr->Fill(float(n_strip_over_threshold));
	     }
  // if(i_local_max > 4 ) { ev_display->Write();}
     f.cd();
	   if(i < 100)
	    { ev_display->Write();
	      ev_raw_display->Write();
	     }
	   delete ev_display;
	   delete ev_raw_display;
	   printf(" wevento %d \n",i);
	   file.cd();
	    }

  // if((i - i/10) < 1) {cout << " evento  " << i << " processato " << endl;}
    }

  // entro nel TFile dove voglio scrivere gli istogrammi e i grafici
  // f_tree->cd();
  printf(" fine data processing \n");
  file->Close();
  printf(" close input data file \n"); 
  f->cd();
  common_tot_vs_rms->Write();
  hcm->Write();
  hcm_prof->Write();
  common_tot_rms_distr->Write();
//  common_tot_vs_rms->Write();
//  common_ADC_1_distr->Write();
//  common_ADC_2_distr->Write();
  printf(" prima scrive strip signal \n");
  for(i=0; i<n_strip_sensor; i++) { strip_value[i]->Write();}
  printf(" dopo strip signal \n");
  for(k=0; k<n_chip; k++) 
   { hcm_chip[k]->Write();
    hcm_prof_chip[k]->Write();
   }
  strip_total_distr->Write();
  n_cluster_profile->Write();
  cluster_width_distr->Write();
  cluster_strip_signal_distr->Write();
  cluster_integral_signal_distr->Write();
  n_cluster_distr->Write();
  strip_occupancy->Write();
  strip_occupancy_max->Write();
  f->Close();
//  f1->Close();
//   
// chiusura file ridotto
//
//  fclose(outfile);*/
}
