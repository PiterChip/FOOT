
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

double calculateArea(const std::vector<double>& x, const std::vector<double>& y, double xmin, double xmax) 
{
    double area = 0.0;
    for (size_t i = 0; i < x.size() - 1; ++i) 
    {
        if (x[i] >= xmin && x[i + 1] <= xmax) 
        {
            double dx = x[i + 1] - x[i];
            double avg_y = (y[i] + y[i + 1]) / 2.0;
            area += dx * avg_y;
        }
    }
    return area;
}


void Analisi_Profilo_Efficienza_di_Carica() 
{ 
   // Percorso del file
   std::string filePath = "C:/root_v5.34.38/macros/Correzione_eta_.txt";

   // Apri il file
   std::ifstream infile(filePath.c_str());
   if (!infile.is_open()) 
    {
      std::cerr << "Errore nell'aprire il file: " << filePath << std::endl;
      return;
    }

   // Ignora la prima riga
   std::string line;
   std::getline(infile, line);

   // Vettori per i dati
   std::vector<double> x_values;
   std::vector<double> y_values_with_correction;

   int rowCount = 0;

   // Leggi il file riga per riga
   while (std::getline(infile, line)) 
    {
      rowCount++;
      if (rowCount > 12000) break;

      std::stringstream ss(line);
      double col1, col2, col3;
      char delimiter;

      // Parsing della riga
      ss >> col1 >> delimiter >> col2 >> delimiter >> col3;
      if (ss.fail())
       {
          std::cerr << "Errore nel parsing della riga: " << line << std::endl;
          continue;
       }

      // Calcola la posizione e verifica la condizione
      double position = col1 * 150;
      if (fabs(col2 - 4.25) < 1e-6) 
       {
          x_values.push_back(position);
          y_values_with_correction.push_back(col3);
       }
    }

    infile.close();

    // Controlla che i vettori non siano vuoti
    if (x_values.empty() || y_values_with_correction.empty()) 
    {
      std::cerr << "Errore: i vettori sono vuoti!" << std::endl;
      return;
    }

    // Trova il massimo valore di y_values_with_correction
    double max_y = *std::max_element(y_values_with_correction.begin(), y_values_with_correction.end());

    // Crea un vettore con valori costanti pari al massimo
    std::vector<double> y_values_constant(x_values.size(), max_y);


    // Crea il primo grafico
    TGraph *graph_with_correction = new TGraph(x_values.size(), &x_values[0], &y_values_with_correction[0]);
    graph_with_correction->SetTitle("Profilo di Efficienza di Raccolta di Carica (Con Correzione)");
    graph_with_correction->GetXaxis()->SetTitle("Coordinata Spaziale [um]");
    graph_with_correction->GetYaxis()->SetTitle("Efficienza di Carica");
    graph_with_correction->SetLineColor(kBlue);
    graph_with_correction->SetMarkerStyle(21);

    // Crea il secondo grafico
    TGraph *graph_constant = new TGraph(x_values.size(), &x_values[0], &y_values_constant[0]);
    graph_constant->SetLineColor(kRed);
    graph_constant->SetLineStyle(2); // Linea tratteggiata
    graph_constant->SetMarkerStyle(22);


    // Disegna i grafici sullo stesso canvas
    TCanvas *canvas = new TCanvas("canvas", "Grafico di Efficienza", 800, 600);
    graph_with_correction->Draw("ALP"); // Disegna il primo grafico
    graph_constant->Draw("LP SAME");   // Disegna il secondo grafico sopra il primo

    // Aggiungi una legenda
    TLegend *legend = new TLegend(0.7, 0.8, 0.9, 0.9);
    legend->AddEntry(graph_with_correction, "Efficienza Corretta", "l");
    legend->AddEntry(graph_constant, "Massimo Costante", "l");
    legend->Draw();

    // Calcola le aree sottese ai grafici
    double xmin = 1.875;
    double xmax = 148.125;
    double area_with_correction = calculateArea(x_values, y_values_with_correction, xmin, xmax);
    double area_constant = calculateArea(x_values, y_values_constant, xmin, xmax);

    // Calcola il rapporto
    double ratio = area_constant / area_with_correction;

    // Stampa i risultati
    std::cout << "Area sotto il grafico con correzione: " << area_with_correction << std::endl;
    std::cout << "Area sotto il grafico costante: " << area_constant << std::endl;
    std::cout << "Rapporto (costante / correzione): " << ratio << std::endl;

    // Salva il grafico
    TFile *outputFile = new TFile("c:/root_v5.34.38/macros/Profili_Efficienza_Grafici.root", "RECREATE");
    graph_with_correction->Write("Graph_With_Correction");
    graph_constant->Write("Graph_Constant");
    canvas->Write();
    outputFile->Close();






}
    
  
     
 


