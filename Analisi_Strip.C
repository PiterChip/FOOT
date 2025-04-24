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

    std::ifstream infile("C:/root_v5.34.38/macros/reduced_PPFM001B_PL_007.txt");  //nel file testo cambiare , con . !!!!!!!

    if (!infile.is_open()) {
        std::cerr << "Errore nell'aprire il file!" << std::endl;
        return;
    }

    double x = 0; //serve per la somma (accumulatore)
    std::vector<double> sum, spessore;   //ci salvo x e x*150 um
    std::string line;
    

    int j = 0; 
    // Leggi i dati
    while (std::getline(infile, line)) {
        std::stringstream ss(line);
        if (j >= 0) {                                            //non salta nessuna riga   
            double evento, max, start_strip, end_strip, num1, num2, num3, num4, num5;
            ss >> evento >> max >> start_strip >> end_strip >> num1 >> num2 >> num3 >> num4 >> num5;
            if (num1 > 10) {
                x = x + num1;
            }
            if (num2 > 10) {
                x = x + num2;
            }
            if (num3 > 10) {
                x = x + num3;
            }
            if (num4 > 10) {
                x = x + num4;
            }
            if (num5 > 10) {
                x = x + num5;
            }
            else {
                x = x + 0;
            }

        }
        sum.push_back(x);
        spessore.push_back(x * 150);  //tutti gli elementi in micron
        x = 0;

        j++;
    }

    infile.close();


    //FACCIO GRAFICO ADC VS SPESSORE
    TCanvas* c1 = new TCanvas("c1", "ADC-Spessore Graph");
    TGraph* gr = new TGraph(spessore.size(), &spessore[0], &sum[0]);
    gr->SetLineColor(600); gr->SetMarkerStyle(20); gr->SetMarkerColor(600);

    TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->AddEntry(gr, "50 V");
   

    gr->SetTitle("ADC vs Spessore;Spessore [um];ADC");
    gr->Draw("ALP");   //se non vuoi la linea, AP
    leg->Draw("Same");


    //Salvo come files .root il plot senza piedistallo
    c1->SaveAs("C:/root_v5.34.38/macros.root");
    c1->SaveAs("C:/root_v5.34.38/macros.pdf","RECREATE");


}