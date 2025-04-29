#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TH1F.h>
#include <TF1.h>
#include <TRandom.h>
#include <TLegend.h>
#include <TString.h>
#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include "TFrame.h"
#include "TBenchmark.h"
#include "TFile.h"
#include "TROOT.h"
#include "TError.h"
#include "TInterpreter.h"
#include "TSystem.h"                             
#include "TPaveText.h"
#include "algorithm"
#include <TSpline.h>
using namespace std;


void Analisi_Strip() {

    std::ifstream infile("C:/root_v5.34.38/macros/reduced_PPFM001B_PL_tutti_i_txt.txt");  
    if (!infile.is_open()) 
    {
        std::cerr << "Errore nell'aprire il file!" << std::endl;
        return;
    }

    double adc = 0; // accumulatore di segnale sul cluster
    double width = 0; // serve per il cluster width

    // Istogrammi bidimensionali
    TH2 *ADC_vs_ClusterWidth = new TH2F("ADC vs  Width", "ADC vs WIDTH", 3500., 0., 3500., 500., 0., 500.);
    ADC_vs_ClusterWidth->GetXaxis()->SetTitle("Cluster Width [um]");
    ADC_vs_ClusterWidth->GetYaxis()->SetTitle("ADC [...]");

    TProfile *profile = new TProfile("ADC vs Width Profile ", "ADC vs WIDTH Profile", 3500., 0., 3500. , 0., 500.);
    profile->GetXaxis()->SetTitle("Cluster Width [um]");
    profile->GetYaxis()->SetTitle("ADC [...]");


    //std::vector<double> signal, width;   // salvo signal e width     
    std::string line;     
    int j = 0;

    // Leggi i dati
    while (std::getline(infile, line)) // mentre legge una riga dal file di input "infile" e la memorizza nella stringa "line"
    {
        std::stringstream ss(line); // oggetto "ss" inizializzato con il contenuto della riga appena letta
        if (j>0)  //non salta nessuna riga
        {               
            double evento, max, start_strip, end_strip;
            double nums[40]; //Array per memorizzare i segnali sulle strip del cluster letti sulla riga
            ss >> evento >> max >> start_strip >> end_strip;

            //Leggo i tot numeri nella riga
            for (int i = 0; i < 40; i++) 
            {
                ss >> nums[i];
            }  

            // Sommo i valori maggiori di 10
            adc = 0; // Reset 
            for (int i = 0; i < 40; i++) 
            {
              if (nums[i] > 10) 
              {
                 adc += nums[i]; //riempio l'accumulatore di segnale
              }
            }
            
            // Calcola la larghezza del cluster: 
            width = ((end_strip - start_strip) -5) * 150;  //+1 perchè voglio start ed end compresi nel cluster e -6 perchè non voglio i bordi ; *150 serve per tener conto della larghezza di una strip espressa in micrometri

            // Riempimento degli istogrammi
            ADC_vs_ClusterWidth->Fill(width, adc);
            profile->Fill(width, adc);                    
        }        
         
        j++;
    }

    infile.close();

    // Salva gli istogrammi in un file .root e .pdf
    TFile *outputFile = new TFile("C:/root_v5.34.38/macros/Grafici/reduced_tutti.root", "RECREATE");
    ADC_vs_ClusterWidth->Write();
    profile->Write();
    outputFile->Close();
    
    TCanvas* c1 = new TCanvas("c1", "ADC-Width Graph");

    /*TGraph* gr = new TGraph(spessore.size(), &spessore[0], &sum[0]);
    gr->SetLineColor(600); gr->SetMarkerStyle(20); gr->SetMarkerColor(600);

    TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->AddEntry(gr, " bias voltage: 50 V");
   

    gr->SetTitle("ADC vs Spessore;Spessore [um];ADC");
    gr->Draw("ALP");   //se non vuoi la linea, AP
    leg->Draw("Same");*/

    
    //ADC_vs_ClusterWidth->Draw("COLZ");
    

    //Salvo come files .root  e .pdf il plot 
    c1->SaveAs("C:/root_v5.34.38/macros/Grafici/reduced_tutti.pdf");
    

 }