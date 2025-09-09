#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TString.h>
#include <iostream>
#include <algorithm>
#include <vector>

#include "utils/definitions.h"
#include "utils/basic_functions.h"

void overflow_underflow_analysis(const char* filename = "fake_clusters.root", 
                                const char* treename = "FakeClusters", 
                                const char* branchname = "t",
                                float underflow_cutoff = 0.0,
                                float overflow_cutoff = 4.0) {
  
  gROOT->SetBatch();  // so stuff isn't autoplotted
  gStyle->SetOptStat(0); // remove the info box
  
  // Open the ROOT file
  TFile* file = TFile::Open((Utils::Definitions::stackRoot + filename).c_str());
  if (!file || file->IsZombie()) {
    std::cerr << "Error: Could not open ROOT file: " << filename << std::endl;
    return;
  }

  // Get the tree
  TTree* tree = (TTree*) file->Get(treename);
  if (!tree) {
    std::cerr << "Error: Could not find tree: " << treename << std::endl;
    file->Close();
    delete file;
    return;
  }

  // Set up branch reading
  float branchValue;
  tree->SetBranchAddress(branchname, &branchValue);

  // First pass: find min and max values
  Long64_t nEntries = tree->GetEntries();
  std::vector<float> values;
  values.reserve(nEntries);
  
  std::cout << "Reading " << nEntries << " entries to find min/max values..." << std::endl;
  
  for (Long64_t i = 0; i < nEntries; i++) {
    tree->GetEntry(i);
    values.push_back(branchValue);
  }

  // Find min and max
  auto minmax = std::minmax_element(values.begin(), values.end());
  float minVal = *minmax.first;
  float maxVal = *minmax.second;
  
  std::cout << "Branch statistics:" << std::endl;
  std::cout << "  Minimum value: " << minVal << std::endl;
  std::cout << "  Maximum value: " << maxVal << std::endl;
  std::cout << "  Number of entries: " << nEntries << std::endl;

  // Create histograms for overflow and underflow
  // Underflow: everything below underflow_cutoff
  TH1D* h_underflow = new TH1D("h_underflow", 
                               Form("Underflow distribution (%s < %.1f);%s;Counts", branchname, underflow_cutoff, branchname),
                               100, minVal, underflow_cutoff);
  
  // Overflow: everything above overflow_cutoff
  TH1D* h_overflow = new TH1D("h_overflow", 
                              Form("Overflow distribution (%s > %.1f);%s;Counts", branchname, overflow_cutoff, branchname),
                              100, overflow_cutoff, maxVal);

  // Second pass: fill histograms
  std::cout << "Filling histograms..." << std::endl;
  
  int underflowCount = 0;
  int overflowCount = 0;
  
  for (float value : values) {
    if (value < underflow_cutoff) {
      h_underflow->Fill(value);
      underflowCount++;
    } else if (value > overflow_cutoff) {
      h_overflow->Fill(value);
      overflowCount++;
    }
  }

  std::cout << "Histogram statistics:" << std::endl;
  std::cout << "  Underflow entries (< " << underflow_cutoff << "): " << underflowCount << " (" 
            << (100.0 * underflowCount / nEntries) << "%)" << std::endl;
  std::cout << "  Overflow entries (> " << overflow_cutoff << "): " << overflowCount << " (" 
            << (100.0 * overflowCount / nEntries) << "%)" << std::endl;
  std::cout << "  Normal range entries (" << underflow_cutoff << "-" << overflow_cutoff << "): " 
            << (nEntries - underflowCount - overflowCount) << " (" 
            << (100.0 * (nEntries - underflowCount - overflowCount) / nEntries) << "%)" << std::endl;

  // Create canvas and draw histograms
  TCanvas* canvas = new TCanvas("canvas", "Overflow and Underflow Analysis", 1200, 600);
  canvas->Divide(2, 1);

  // Draw underflow histogram
  canvas->cd(1);
  h_underflow->SetLineColor(kRed);
  h_underflow->SetFillColor(kRed);
  h_underflow->SetFillStyle(3001);
  h_underflow->Draw();
  gPad->SetLogy();

  // Draw overflow histogram
  canvas->cd(2);
  h_overflow->SetLineColor(kBlue);
  h_overflow->SetFillColor(kBlue);
  h_overflow->SetFillStyle(3001);
  h_overflow->Draw();
  gPad->SetLogy();

  // Save the canvas
  TString outputPath = (Utils::Definitions::analysisRoot + "output/4d_tracking/overflow_underflow_analysis.pdf").c_str();
  canvas->SaveAs(outputPath);
  std::cout << "Canvas saved to: " << outputPath << std::endl;

  // Also save histograms to ROOT file
  TString histOutputPath = (Utils::Definitions::analysisRoot + "hists/4d_tracking/overflow_underflow_histograms.root").c_str();
  TFile* histFile = TFile::Open(histOutputPath, "RECREATE");
  h_underflow->Write();
  h_overflow->Write();
  histFile->Close();
  delete histFile;
  std::cout << "Histograms saved to: " << histOutputPath << std::endl;

  // Clean up
  delete canvas;
  delete h_underflow;
  delete h_overflow;
  file->Close();
  delete file;
} 