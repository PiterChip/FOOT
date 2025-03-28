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


// File dove scrivo istogrammi e grafici
  TFile *f = new TFile("C:/root_v5.34.38/macros/PROVA.root", "RECREATE");
// File che leggo
  TFile *file = new TFile("C:/root_v5.34.38/macros/foot_3172.root");
  FILE *outfile; 
  char fileout[700] = "C:/root_v5.34.38/macros/nnn.txt";  //file txt che genero di output
  outfile = fopen(fileout,"a");
  TTree * raw_events = (TTree*)file->Get("raw_events");
  Float_t bias = 50.0;      // sensor bias [V]
  Float_t Angle = 74.5;		// source angle [°]

//Inizializzo i vettori dove mettere le letture ADC
//  std::vector<unsigned short> *raw_event = new std::vector<unsigned short>();
  vector<unsigned short> *raw_event = new vector<unsigned short>();
  //Associo ogni TBranch a un vettore
  raw_events->SetBranchAddress("RAW Event", &raw_event);

  Int_t nentries;
  nentries = raw_events->GetEntries();
  Int_t n_pedestal = 500;
  Int_t n_events = nentries;
  n_events = 999;
  printf(" n pedestal %d n events %d nentries %d \n",n_pedestal,n_events, nentries);
  Int_t event_offset = 0;
  Int_t ped_event_offset = 0;
  Int_t n_chip = 10;
  Int_t n_chip_chan = 64;
  Int_t n_strip_sensor = n_chip * n_chip_chan;
  Double_t strip_signal[n_strip_sensor];
  Int_t event_to_display = 535;    
  Double_t threshold = 10.;   // soglia per definire inizio di un cluster
  Double_t common_thre = 10.;
  Int_t common_mode_choice = 0;
  Int_t n_max_in_event = 10;
  Int_t max_cluster_width = 640;
  Int_t good_event[n_events];
  Double_t local_max[n_events][n_max_in_event][6];
  Double_t cluster_low_threshold = 3.;            //soglia minima per definire la appartenenza di una strip ad un cluster
  Int_t i,j,k,j_index,l;
  char histo_title[3000];
  char histo_base[50];
  
  TH1F *hRMS = new TH1F("RMS","RMS",100,0.,10.);
  TH1F *hMean = new TH1F("Mean","Mean",75,100.,400.);
  TH1F *hcm100 = new TH1F("CM100","CM100",20,-0.5,0.5);
  TH1F *hcm_chip[n_chip];
  TH1F *hcm_prof_chip[n_chip];
  TH1F *common_ADC_hist[2];
  TH1F *common_chip[n_chip];
  TH1F *strip_value[n_strip_sensor];
  TH1F *common_ADC_1_distr = new TH1F("Common noise ADC 1","CM ADC 1 distribution",400,-25.,25.);
  TH1F *common_ADC_2_distr = new TH1F("Common noise ADC 2","CM ADC 2 distribution",400,-25.,25.);
  for(k=0; k<n_chip; k++) 
  { TH1F *hcm_chip[k] = new TH1F("CM x chip k","CM x chip k",100,-20,20);
    TH1F *hcm_prof_chip[k] = new TH1F("CM profile x chip k","CM profile x chip k",n_events,0,float(n_events));
  }
  TH1F *n_strip_common_distr = new TH1F("CM","CM",70,0.,70.);
  TH1F *n_total_strip_common_distr = new TH1F("total strip CM","total strip CM",700,0.,700.);
  TH1F *hcm = new TH1F("h4","CM",100,-20,20);
  TH1F *signal_event = new TH1F("signal_event","signal_event",n_events,0,float(n_events));
  TH1F *rms_signal_distr = new TH1F("RMS signal distribution","Event RMS signal distribution",1000,0.,200.);
  TH1F *common_tot_rms_distr = new TH1F("RMS of common mode distribution","Event RMS signal distribution",1000,0.,100.);
  TH1F *strip_total_distr = new TH1F("Strip total distribution","Strip total distribution",2500,-50.,200.);
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
//  TCanvas * c = new TCanvas("c", "c", 1000, 1000);

  // dichiaro array di istogrammi per il piedistallo
  TH1 *h[n_strip_sensor];
  TH1 *h_prof[n_strip_sensor];
  for (Int_t k=0; k<n_strip_sensor; k++)
  { h[k] = new TH1F("histo","piedistallo",350,50.,400.);
    h_prof[k] = new TH1F("histo","piedistallo profile",n_pedestal,0.,float(n_pedestal));
    strip_value[k] = new TH1F("strip signal","strip signal",800,-20.,20.);
  }

  // CALCOLO PIEDISTALLO
  // cicli dove riempo 640 istogrammi e 640 grafici uno per ogni strip
  unsigned short adc[n_strip_sensor];
  for (Int_t i=ped_event_offset;i<n_pedestal;i++)
  { raw_events->GetEntry(i);
    for (Int_t j=0; j<n_strip_sensor; j++){
      adc[j] = (*raw_event)[j];
      h[j]->Fill(adc[j]);
      h_prof[j]->Fill(float(i),adc[j]);
    }
  }
  // salvo il valor medio del piedistallo per ogni strip
  Float_t Mean[n_strip_sensor];
  Float_t Noise[n_strip_sensor];
  for (j=0; j<n_strip_sensor; j++){
    Mean[j]=h[j]->GetMean();
    Noise[j]=h[j]->GetRMS();
  }
//
//  definizioni variabili di singolo evento

  Float_t signal[n_strip_sensor];

// Calcolo common noise
// grandezze calcolate:  
// 0: Common mode per tutto il sensore (1)
// 1: Common mode per ogni ADC   (2)
// 2: Common mode per ogni chip  (10)

  Int_t ev_to_display = 0;
  Float_t common_tot;                  //  array dove inserire il common mode unico [0]
  Float_t common_tot2;                 //  array dove inserire il common mode unico [0]
  Int_t n_total_strip_common;	 	     // numero di strip usate per il common mode unico
  Float_t common_ADC[2];				 //  array dove inserire il common mode per ogni ADC [1]
  Int_t n_ADC_strip_common[2];         // numero di strip usate per il common mode per ogni ADC
  Float_t common[n_chip];				 //  array dove inserire il common mode per ogni chip [2]
  Int_t n_strip_common[n_chip];        // numero di strip usate per il common mode per ogni chip
 
 //
  for ( Int_t i=event_offset; i<n_events; i++)
  { //printf(" dentro loop eventi \n");
    raw_events->GetEntry(i);
	  good_event[i] = 1;
    common_tot = 0.;
    common_tot2 = 0.;
	  n_total_strip_common = n_strip_sensor;
    common_ADC[0] = 0.;
	  n_ADC_strip_common[0] = n_chip_chan*5;
    common_ADC[1] = 0.;
	  n_ADC_strip_common[1] = n_chip_chan*5;
    // event display histogram  
	  strcpy(histo_title,"event_display_");
    stringstream histo_suffix;
    histo_suffix << i;
	  histo_suffix >> histo_base; 
	  strcat(histo_title,histo_base);
    TH1F *ev_display = new TH1F(histo_title,histo_title,n_strip_sensor,0.,float(n_strip_sensor));
    TH1F *ev_raw_display = new TH1F(histo_title,histo_title,n_strip_sensor,0.,float(n_strip_sensor));
    //	printf(" evento %d prima chip \n",i);
    for( k=0; k<n_chip; k++)
    { common[k] = 0;
	    n_strip_common[k] = n_chip_chan;
      for (j=0; j<n_chip_chan; j++)
	    { j_index = j + k*n_chip_chan;
        //	    printf(" i %d chip %d strip %d strip ID %d \n",i,k,j,j_index);
        ev_raw_display->Fill(float(j_index),signal[j_index]);
		    signal[j_index] = ((*raw_event)[j_index]);
        ev_raw_display->Fill(float(j_index),signal[j_index]);
		    signal[j_index] = signal[j_index]-Mean[j_index];
        if(signal[j_index] < common_thre) 
		    { common_tot = common_tot+signal[j_index];
		      common_tot2 = common_tot2+signal[j_index]**2;
          //		   printf(" i %d j_index %d signal %f common %f common2 %f \n",i,j_index,signal[j_index],common_tot,common_tot2);
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
      }
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
	  }
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
	  common_ADC_1_distr->Fill(common_ADC[0]);
	  common_ADC_2_distr->Fill(common_ADC[1]);
//
//  start event processing
//
    if(good_event[i] > 0) 
	  { if(common_tot_rms < 6.)   // rms del common mode <6 (impostazione per grarantirci qualità) altrimenti l'evento lo buttiamo
	    { Int_t i_local_max = 0;
	      Int_t in_cluster = 0;  	// flag per indicare che stiamo processando un cluster 
	      local_max[i][i_local_max][0] = 0.;							//  start cluster 
	      local_max[i][i_local_max][1] = 0.;    						//  end cluster 
	      local_max[i][i_local_max][2] = 0.;							//  cluster maximum
	      local_max[i][i_local_max][3] = 0.;							//  position of max in cluster
	      local_max[i][i_local_max][4] = 0.;    						//  cluster width 
	      local_max[i][i_local_max][5] = 0.;  						//  integral of cluster signal
	      Int_t n_strip_over_threshold = 0;
	      Double_t ev_signal = 0.;									//  Calcolo media e RMS signal per identificare bad events
	      Double_t ev_signal2 = 0.;									//  Calcolo media e RMS signal per identificare bad events
	      Double_t rms_signal = 0.;									//  Calcolo media e RMS signal per identificare bad events
	      
        for(k=0; k<max_cluster_width; k++) { strip_signal[k] = 0.; } //sto ciclando su? 
        for( k=0; k<n_chip; k++) //sto ciclando sui 10 chip
        { for (j=0; j<n_chip_chan; j++) //sto ciclando sulle 64 strip di ogni chip
          //quindi sto ciclando su tutto il sensore
	        { j_index = j + k*n_chip_chan;
		        signal[j_index] = (signal[j_index]-common[k])*good_event[i];
			      strip_value[j_index]->Fill(signal[j_index]);
		        ev_display->Fill(float(j_index),signal[j_index]);
		        ev_signal = ev_signal + signal[j_index];
		        ev_signal2 = ev_signal2 + signal[j_index]**2;		
	  	      strip_total_distr->Fill(signal[j_index]);
            //printf(" evento %d k %d j %d indxex %d prima massimo relativo \n",i,k,j,j_index);
            // look for relative maximum
           if(signal[j_index] > cluster_low_threshold)    
		       { if(in_cluster < 1) {local_max[i][i_local_max][0] = float(j_index);   } // start cluster
			       if(in_cluster < max_cluster_width)  
		         { if(signal[j_index] > cluster_low_threshold) {n_strip_over_threshold = n_strip_over_threshold + 1;}
		           in_cluster = in_cluster +1;
			         // printf(" evento %d k %d j %d indxex %d  signal %f dentro loop massimo relativo \n",i,k,j,j_index,signal[j_index]);
		           if(signal[j_index] > local_max[i][i_local_max][2])  //se il segnale su questa microstrip è > del massimo locale 
		           { local_max[i][i_local_max][3] = float(j_index);    //allora il massimo diventa questo che ho trovato adesso 
		              local_max[i][i_local_max][2] = signal[j_index];  //e la posizione del massimo diventa questa che ho trovato adesso
                  //quindi io ho spostato il valore del massimo e la posizione del massimo.
		  	          // printf(" inside cluster j_index %d signal %5.2f local_max %5.2f position %5.0f \n",j_index,signal[j_index],local_max[i][i_local_max][2],local_max[i][i_local_max][3]);
		           }
		           local_max[i][i_local_max][1] = float(j_index);   // current end of cluster,
               // cioè ho trovato un'altra strip che fa parte del cluster e questa per il momento 
               // è l'ultima del cluster però devo andare a vedere anche quelle dopo. 
               // Quindi questa è una cosa che mi dice che ogni volta che sarò maggiore io aumenterò questa,
               // alla fine mi dirà: arrivo all'ultima poi sto sotto la soglia e quindi questa diventa l'ultima del cluster che potrebbe anche essere la massima.
               // Il massimo può coincidere con la prima o con l'ultima strip del cluster, io il massimo devo sapere dove è ma anche dove siano l'inizio e la fine di un cluster;
               // quindi ogni strip che aggiungo devo controllare se per caso sia il massimo e la aggiungo alla fine del cluster.

		           local_max[i][i_local_max][4] = local_max[i][i_local_max][1]-local_max[i][i_local_max][0]+1;   // cluster width : aumento la larghezza del cluster di 1
		           local_max[i][i_local_max][5] = local_max[i][i_local_max][5]+signal[j_index];   // integral of cluster signal : faccio l'integrale di tutto il segnale che sta sul cluster
		           
               cluster_strip_signal_distr->Fill(signal[j_index]); //riempio degli istogrammi
               
	             strip_occupancy->Fill(float(j_index)); //riempio degli istogrammi
               
		           //printf(" evento %d strip %d # cluster %d signal %5.1f pos %5.1f local max %5.1f  %5.0f \n",i,j_index,i_local_max,signal[j_index],local_max[i][i_local_max][2],local_max[i][i_local_max][3],in_cluster);
			       }
	 	       }
		       else 
		       { if(in_cluster > 0 ) 
		         { in_cluster = 0;      		  // dichiaro chiuso quel cluster 
		           Int_t l1 = 0;
		           strip_occupancy_max->Fill(local_max[i][i_local_max][3]);
               //if(local_max[i][i_local_max][0] < 3) { Int_t lower_cluster = 0;} else {Int_t lower_cluster = int(local_max[i][i_local_max][0])-3;}
               //if(local_max[i][i_local_max][1] > n_strip_sensor-3) { Int_t upper_cluster = 640;} else {Int_t upper_cluster = int(local_max[i][i_local_max][1])+3;}
               if((local_max[i][i_local_max][0] > 3 ) && (local_max[i][i_local_max][1] < (n_strip_sensor-3)) ) 
               { lower_cluster = int(local_max[i][i_local_max][0])-3;
                 upper_cluster = int(local_max[i][i_local_max][0])+3;
                 printf("lower_cluster %d upper_cluster %d \n",lower_cluster,upper_cluster);
                 fprintf(outfile," %d  %d %d %d ",i,i_local_max+1,int(local_max[i][i_local_max][0])-3,int(local_max[i][i_local_max][1])+3);
                 for(l=lower_cluster-1; l<upper_cluster+1; l++) 	
			           { fprintf(outfile,"%.2f ",signal[l]);	  
			             strip_signal[l1] = signal[l];
			             l1 = l1 + 1;
		   	         }
			           fprintf(outfile,"\n");
                 printf("fines crittura cluster \n");
                 cluster_width = upper_cluster - lower_cluster;
		             for(l1=0; l1<cluster_width; l1++) { strip_signal[l1] = 0.; } 
		             i_local_max = i_local_max + 1;    // incremento contatore per il prossimo cluster
		           }
		         }
		       }
           if(i>996 ) {printf("fine blocco j_index %d j %d k %d signal %f \n",j_index, j, k, signal[l]);}
	       }
	     }
       printf("fine blocco 1 \n");
 	     if(i_local_max > 0 && i_local_max < 8)
	     { cluster_width_distr->Fill(local_max[i][i_local_max-1][4]);
	       cluster_integral_signal_distr->Fill(local_max[i][i_local_max-1][5]);
         if(i_local_max > 5.) { printf( "   evento troppi cluster %d n strip over threshold %d \n",i,n_strip_over_threshold);}
	     }
	     if(good_event[i] > 0) 
	     { ev_signal = ev_signal/n_strip_sensor;
	       rms_signal = sqrt(ev_signal2/n_strip_sensor - ev_signal**2);
	       rms_signal_distr->Fill(rms_signal);
	       n_cluster_profile->Fill(float(i),float(i_local_max));   // profilo # cluster nel frame
	       n_cluster_distr->Fill(float(i_local_max));   // distribuzione # cluster nel frame
	       n_strip_over_threshold_distr->Fill(float(n_strip_over_threshold));
	     }
       //	  if(i_local_max > 4 ) { ev_display->Write();}
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
     //
     //	if((i - i/10) < 1) {cout << " evento  " << i << " processato " << endl;}
   }
 }
 //entro nel TFile dove voglio scrivere gli istogrammi e i grafici  
 //  f_tree->cd();
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
 fclose(outfile);
}



