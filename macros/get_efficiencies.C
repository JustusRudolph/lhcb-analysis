#include <TFile.h>
#include <TH1F.h>
#include <TEfficiency.h>
#include <TCanvas.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <unordered_set>

#include "utils/definitions.h"
#include "utils/basic_functions.h"

void plot_efficiencies(unsigned nEvents, unsigned max_scatter, unsigned max_dt,
                       bool isForward=true, bool withElectrons=true) {
  gROOT->SetBatch();  // so stuff isn't autoplotted
  // Open ROOT file and retrieve histograms
  std::string region = (isForward ? "forward" : "backward");
  TString suffix = Utils::Functions::get_suffix(nEvents, max_scatter, max_dt);
  TString input_suffix = suffix + ".root";
  TFile *file = TFile::Open(
    TString(Utils::Definitions::stackRoot + "/output/PrCheckerPlots") + input_suffix);
  if (!file || file->IsZombie()) {
      std::cerr << "Error: Could not open ROOT file." << std::endl;
      return;
  }
  // get ghost rates, these are calculated separately
  TFile* clone_histo_file = TFile::Open(
    TString(Utils::Definitions::analysisRoot + "/hists/clones/mc_hists") + input_suffix);
  TProfile* ghostRates = (TProfile*) clone_histo_file->Get("ghost_rates");


  // set up canvases for forward and backward region, and types of plots
  std::vector<std::string> types = {"Pt", "Eta", "docaz"};
  std::vector<std::string> x_axes = {"p_{T} (MeV)", "#eta", "DOCA_{z} (mm)"};
  std::string base_string = "velo_kalman_validator/VeloTracks_eta";
  TCanvas* canvas = new TCanvas("canvas", "Forward Region Efficiencies", 1000, 800);
  canvas->Divide(2, 2);
  int c_idx = 1;

  for (int typeIndex = 0; typeIndex < types.size(); typeIndex++, c_idx++) {
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
      std::string title = "Efficiency vs " + x_axes[typeIndex];
      eff->SetTitle((title + ";" + x_axes[typeIndex] + ";Efficiency").c_str());
      eff->Draw("AP");

      gPad->Update();
      // for pT and eta, cut y at 0.8 to have the same range
      if (type == "Eta") {
        if (isForward) {
          eff->GetPaintedGraph()->GetXaxis()->SetRangeUser(2., 5.);
        } else {
          eff->GetPaintedGraph()->GetXaxis()->SetRangeUser(-5., -2.);
        }
         eff->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.8, 1.01);
      } else if (type == "docaz") {
        eff->GetPaintedGraph()->GetYaxis()->SetRangeUser(0, 1.01);
        eff->GetPaintedGraph()->GetXaxis()->SetRangeUser(0, 7);
      } else if (type == "Pt") {
        eff->GetPaintedGraph()->GetXaxis()->SetRangeUser(0., 5000);
        eff->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.8, 1.01);
      }
    } else {
      std::cerr << "Histograms are not consistent for TEfficiency calculation." << std::endl;
    }
  }  // types
  canvas->cd(4);
  // Draw ghost rates
  gStyle->SetOptStat(0); // remove the info box
  ghostRates->SetTitle("Ghost Rates;#eta;Ghost Rate");
  ghostRates->SetLineColor(kMagenta);
  ghostRates->GetYaxis()->SetRangeUser(0., 0.1);
  if (isForward)
    ghostRates->GetXaxis()->SetRangeUser(1., 5.5);
  else
    ghostRates->GetXaxis()->SetRangeUser(-5.5, -1.);
  ghostRates->Draw();
  canvas->SaveAs(TString((Utils::Definitions::analysisRoot + "/output/efficiency/" +
                          region + "_efficiency_plots").c_str()) + suffix + ".pdf");
  // Clean up
  file->Close();
  delete file;
  delete canvas;
}

void get_efficiencies(unsigned nEvents=5000, unsigned max_scatter=80000,
                      unsigned max_dt=0, bool withElectrons=false) {
  // plot both front and backwards regions
  plot_efficiencies(nEvents, max_scatter, max_dt, true, withElectrons);
  plot_efficiencies(nEvents, max_scatter, max_dt, false, withElectrons);
}