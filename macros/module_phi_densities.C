#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TMath.h>

#include <cstdlib>
#include <iostream>
#include <string>

void module_phi_densities(unsigned minMod, unsigned maxMod, unsigned step=2,
                          unsigned runNoCheck=-1, unsigned evNoCheck=-1) {
  gROOT->SetBatch();  // don't show hists while running
  std::string stackRoot = std::getenv("STACK_ROOT");
  std::string analysisRoot = std::getenv("ANALYSIS_ROOT");

  TFile *file = TFile::Open((stackRoot + "/fake_clusters.root").c_str());
  if (!file || file->IsZombie()) {
      std::cerr << "Error: Could not open ROOT file." << std::endl;
      return;
  }
  // Create Histograms for each module
  std::vector<TH1D*> moduleDensitiesByPhi;
  for (unsigned i = minMod; i <= maxMod; i += step) {
    std::string histName = "moduleDensityByPhi_" + std::to_string(i);
    moduleDensitiesByPhi.push_back(
      new TH1D(histName.c_str(),
               ("Occupancy in module " + std::to_string(i) + "by #phi;#phi;Counts").c_str(),
               100, -3.15, 3.15));
  }
  // get data from tree
  unsigned lhcbid, evNo, runNo;
  float x, y;
  TTree* fakeClusTree = (TTree*) file->Get("FakeClusters");
  if (!fakeClusTree) {
      std::cerr << "TTree not found!" << std::endl;
      return;
  }
  fakeClusTree->SetBranchAddress("runNo", &runNo);
  fakeClusTree->SetBranchAddress("evNo", &evNo);
  fakeClusTree->SetBranchAddress("id", &lhcbid);
  fakeClusTree->SetBranchAddress("x", &x);
  fakeClusTree->SetBranchAddress("y", &y);

  unsigned nEntries = fakeClusTree->GetEntries();

  for (unsigned i = 0; i < nEntries; i++) {
    fakeClusTree->GetEntry(i);
    // only take from specific run and event if specified
    if (runNoCheck >=0 && runNo != runNoCheck) continue;
    if (evNoCheck >=0 && evNo != evNoCheck) continue;

    unsigned tvMod = (lhcbid >> 12) & 0x3F;
    if (tvMod < minMod || tvMod > maxMod) continue;
    unsigned modIdx = (tvMod - minMod) / step;
    float phi = TMath::ATan2(y, x);
    moduleDensitiesByPhi[modIdx]->Fill(phi);

    // print status
    if (i && (i % 1000000 == 0))
      printf("Finished with %d %% of fake clusters.\n", (int) (100. * (float) i / nEntries));
  }
  // get maximum value for scaling
  unsigned maxVal = 0;
  for (unsigned i = 0; i < moduleDensitiesByPhi.size(); i++) {
    if (moduleDensitiesByPhi[i]->GetMaximum() > maxVal)
      maxVal = moduleDensitiesByPhi[i]->GetMaximum();
  }
  // set maximum value & general title for scaling only on first histogram
  moduleDensitiesByPhi[0]->SetMaximum(maxVal * 1.1);
  moduleDensitiesByPhi[0]->SetTitle("Module Occupancy by #phi");

  // now plot the histograms on the same canvas with different colours
  TCanvas* canvas = new TCanvas("c", "Module Occupancy by #phi", 1000, 800);
  gStyle->SetOptStat(0);  // remove the info box

  // make vector of different colours, cap at 8 to not bloat the plot
  std::vector<int> colours = {kBlack, kRed, kBlue, kGreen, kMagenta, kCyan, kOrange, kViolet};
  TLegend* legend = new TLegend(0.4, 0.6, 0.6, 0.8);

  for (unsigned i = 0; i < moduleDensitiesByPhi.size(); i++) {
    moduleDensitiesByPhi[i]->SetLineColor(colours[i % colours.size()]);
    // draw histograms, adapt size if plot gets larger
    if (i == 0) moduleDensitiesByPhi[i]->Draw();
    else moduleDensitiesByPhi[i]->Draw("SAME");
    legend->AddEntry(moduleDensitiesByPhi[i],
                     ("Module " + std::to_string(minMod + i * step)).c_str(), "l");
  }
  legend->Draw();
  canvas->SetTitle("Module Occupancy by #phi");

  // save the canvas
  canvas->SaveAs((analysisRoot + "/output/module_phi_densities.pdf").c_str());
  // clean up
  file->Close();
  delete file;
  delete canvas;
}