#include <TFile.h>
#include <TH1F.h>
#include <TEfficiency.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <unordered_set>

#include "utils/definitions.h"

void eff_pur_plots(bool withElectrons=true) {
  // Open ROOT file and retrieve histograms
  TFile *file = TFile::Open((Utils::Definitions::stackRoot + "/output/PrCheckerPlots.root").c_str());
  if (!file || file->IsZombie()) {
      std::cerr << "Error: Could not open ROOT file." << std::endl;
      return;
  }

  // set up canvases for forward and backward region, and types of plots
  std::vector<std::string> types = {"Pt", "Eta", "docaz"};
  std::vector<std::string> x_axes = {"pT (MeV)", "eta", "docaz (mm)"};
  std::string base_string = "velo_kalman_validator/VeloTracks_eta";
  TCanvas *canvas = new TCanvas("c", "Efficiencies", 1000, 800);
  canvas->Divide(2, 3);
  int c_idx = 1;


  for (int isForward = 0; isForward < 2; isForward++) {
    for (int typeIndex = 0; typeIndex < types.size(); typeIndex++) {
      std::string type = types[typeIndex];
      std::string base_string_region = (isForward ? base_string + "25" : base_string + "-5-2");
      
      TH1F *nonElectronsReconstructible =
        (TH1F*)file->Get((base_string_region + "_notElectrons_" + type + "_reconstructible").c_str());
      TH1F *electronsReconstructible =
        (TH1F*)file->Get((base_string_region + "_electrons_" + type + "_reconstructible").c_str());
      TH1F *nonElectronsReconstructed =
        (TH1F*)file->Get((base_string_region + "_notElectrons_" + type + "_reconstructed").c_str());
      TH1F *electronsReconstructed =
        (TH1F*)file->Get((base_string_region + "_electrons_" + type + "_reconstructed").c_str());

      // Check if histograms are properly loaded
      if (!nonElectronsReconstructible) {
        std::cerr << "Error: Hist "
                  << (base_string_region + "_notElectrons_" + type + "_reconstructible").c_str()
                  << "not found in root file ../../output/PrCheckerPlots.root" << std::endl;
        file->Close();
        return;
      } else if (!electronsReconstructible) {
        std::cerr << "Error: Hist "
                  << (base_string_region + "_electrons_" + type + "_reconstructible").c_str()
                  << "not found in root file ../../output/PrCheckerPlots.root" << std::endl;
        file->Close();
        return;
      } else if (!nonElectronsReconstructed) {
        std::cerr << "Error: Hist "
                  << (base_string_region + "_notElectrons_" + type + "_reconstructed").c_str()
                  << "not found in root file ../../output/PrCheckerPlots.root" << std::endl;
        file->Close();
        return;
      } else if (!electronsReconstructed) {
        std::cerr << "Error: Hist "
                  << (base_string_region + "_electrons_" + type + "_reconstructed").c_str()
                  << "not found in root file ../../output/PrCheckerPlots.root" << std::endl;
        file->Close();
        return;
      }

      // Add the two numerators and denominators
      TH1F *reconstructedTotal = (TH1F*)nonElectronsReconstructed->Clone("reconstructedTotal");
      if (withElectrons) reconstructedTotal->Add(electronsReconstructed);
      
      TH1F *reconstructibleTotal = (TH1F*)nonElectronsReconstructible->Clone("reconstructibleTotal");
      if (withElectrons) reconstructibleTotal->Add(electronsReconstructible);

      // Create TEfficiency object
      if (TEfficiency::CheckConsistency(*reconstructedTotal, *reconstructibleTotal)) {
        TEfficiency *eff = new TEfficiency(*reconstructedTotal, *reconstructibleTotal);
        eff->SetLineColor(kMagenta);
        eff->SetMarkerColor(kMagenta);
        // Draw the efficiency plot
        canvas->cd(c_idx);
        std::string region = (isForward ? "forward" : "backward");
        std::string title = "Efficiency in " + region + " region wrt " + type;
        eff->SetTitle((title + ";" + x_axes[typeIndex] + ";Efficiency").c_str());
        eff->Draw("AP");

        gPad->Update();
        // for pT and eta, cut at 0.8 to have the same
        // TODO: Change to c_idx < 5 for electrons too when better?
        if (withElectrons) {
          if (c_idx < 3 || c_idx == 4) {  // top two are pT
            eff->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.8, 1.01);
            // graph->GetYaxis()->SetNdivisions(10, "N");  // steps of 0.02
          }
        } else {
          if (c_idx < 5) {  // top two are pT
            eff->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.8, 1.01);
            // graph->GetYaxis()->SetNdivisions(10, "N");  // steps of 0.02
          }
        }
        c_idx += 2;  // move to next row
      } else {
        std::cerr << "Histograms are not consistent for TEfficiency calculation." << std::endl;
      }
    }  // types
    c_idx = 2;  // move to second column
  }  // isForward
  canvas->SaveAs((Utils::Definitions::analysisRoot + "/output/efficiency_plots.pdf").c_str());
  // Clean up
  file->Close();
  delete file;
  delete canvas;
}

void allen_post_processing(bool withElectrons=true) {
  gROOT->SetBatch();
  eff_pur_plots(withElectrons);
  get_clone_rates();
}

int main() {
    allen_post_processing();
    return 0;
}