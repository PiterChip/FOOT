{
//   example of macro to read data from an ascii file and
//   create a root file with an histogram and an ntuple.

#include <iostream>
#include <sstream> 
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


gROOT->Reset();
gStyle->SetOptStat(1111111);
gStyle->SetOptFit(111);

#include "Riostream.h"
#define max_row 1536
#define max_col 2048
#define max_sensors  4

TF1 *doublegaus = new TF1("doublegaus","gaus(0)+gaus(3)", -1., 1.);
doublegaus->SetParameter(0,20.);
doublegaus->SetParameter(1,60.);
doublegaus->SetParameter(2,3.0);
doublegaus->SetParameter(3,20.);
doublegaus->SetParameter(4,80.);
doublegaus->SetParameter(5,3.0);

  Int_t event_offset = 0;  
  Int_t ped_event_offset = 0;
  Int_t n_chip = 10;
  Int_t n_chip_chan = 64;
  Int_t n_strip_sensor = n_chip * n_chip_chan;
  Double_t strip_signal[n_strip_sensor];
  Int_t event_to_display = 535;    
  Double_t threshold = 10.;
  Double_t common_thre = 10.;
  Int_t common_mode_choice = 0;
  Int_t n_max_in_event = 60;
  Int_t max_cluster_width = 640;

// File che leggo
  Float_t bias ;      // sensor bias [V]
  Float_t Angle ;		// source angle [�]
  Int_t ev_number;
  Int_t n_cluster;
  Int_t start_cluster;
  Int_t stop_cluster;
  Int_t cluster_width;
  Int_t n_cluster;
  Double_t common_noise ;
  Double_t common_noise_rms ;

  Int_t n_chip = 10;
  Int_t n_chip_chan = 64;
  Int_t n_strip_sensor = n_chip * n_chip_chan;
  Double_t strip_threshold = 5.;
  Float_t strip_bad[n_strip_sensor];
  for(i=0; i<n_strip_sensor; i++) { strip_bad[i] = 1.;}
//  strip_bad[206] = 0.;
//
//  Histogram definition
//
  Double_t sum[60];
  TH1F *sum_x_strip_signal_inner_distr[60];
  TH1F *strip_signal_vs_position_distr[60];
  TH1F *strip_signal_distr = new TH1F("Strip signal distribution","Strip signal distribution",2500,-50.,200.);
  TH1F *strip_signal_distr_1 = new TH1F("Strip signal distribution","Strip signal distribution",2500,-50.,200.);
  TH1F *strip_signal_extreme_distr_1 = new TH1F("Strip signal distribution","Strip signal distribution",2500,-50.,200.);
  TH1F *strip_signal_extreme_distr_2 = new TH1F("Strip signal distribution","Strip signal distribution",2500,-50.,200.);
  TH1F *strip_signal_inner_distr = new TH1F("inner Strip signal distribution","Strip signal distribution",2500,-50.,200.);
  for(i=0; i<60; i++)
  { TH1F *sum_x_strip_signal_inner_distr[i] = new TH1F("Inner sum 1 Strip signal distribution","Strip signal distribution",(100+20*i),1.,(101.+20.*i));
    sum[i] = 0.;
    TH1F *strip_signal_vs_position_distr[i] = new TH1F("Strip signal vs distribution","Strip signal vs strip position distribution",1000,10.,110.);
  }
  TH1F *cluster_width_distr = new TH1F("Cluster width distribution","Cluster width distribution",100,0.,100.);
  TH1F *strip_occupancy = new TH1F("Strip Occupancy","Strip Occupancy",n_strip_sensor,0.,float(n_strip_sensor));
  
// input file
  char conf_file[100], fileio[100], fileout[100], pippo[100];
  strcpy(conf_file,"C:/root_v5.34.38/macros/reduced_PPFM001B_PL_tutti_i_txt.txt");
  printf(" input file \n",conf_file);
  ifstream ino;
  ino.open(conf_file);
  int n_events = 126097;
  int event_number, numero_cluster, start_cluster, stop_cluster, strip_pos;
  float strip_value[640];

  for(int event=0; event<n_events; event++)
  { ino >> event_number;
    ino >> numero_cluster;
  	ino >> start_cluster;
	  ino >> stop_cluster;
	for(i=0; i<640; i++) {strip_value[i] = 0.;}
//	printf(" event number %d \n",event_number);
  int cluster_width_1;
	cluster_width_1 = (stop_cluster - start_cluster +1);
	cluster_width_distr.Fill(cluster_width_1);
  for(int strip_pos=start_cluster; strip_pos<(stop_cluster+1); strip_pos++)
    { ino >> strip_value[strip_pos];
      strip_signal_distr->Fill(strip_value[strip_pos]);
      if(cluster_width_1 > 	9 ) { strip_signal_distr_1->Fill(strip_value[strip_pos]);}
    }
    strip_signal_extreme_distr_1->Fill(strip_value[start_cluster]);
    strip_signal_extreme_distr_2->Fill(strip_value[stop_cluster]);	
  }
/*

  Double_t scale_factor = cluster_width_distr->GetEntries();
  scale_factor = 1./scale_factor;
  printf(" scale_factor %f \n",scale_factor);
  TH1F *cluster_width_norm_distr = (TH1F*)cluster_width_distr->Clone();
  cluster_width_norm_distr->Scale(scale_factor);*/
  
} // chiudi  programma

