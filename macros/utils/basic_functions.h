#pragma once

#include <TH1D.h>
#include <TString.h>

#include "definitions.h"

namespace Utils::Functions {
  
  void printBinEdges(TH1D* hist) {
    int nBins = hist->GetNbinsX();
    for (int i = 1; i <= nBins; ++i) {
      std::cout << "Bin " << i << ": " << hist->GetBinLowEdge(i) << " to "
                << hist->GetBinLowEdge(i) + hist->GetBinWidth(i) << std::endl;
    }
  }

  bool inEtaAcceptance(float eta) {
    float absEta = std::abs(eta);
    return (absEta >= 2. && absEta <= 5.);
  }

  void disableProfileMarkers(TProfile* p) {
      p->SetMarkerStyle(0);
      p->SetMarkerSize(0);
      p->SetMarkerColor(0);
  }

  TH1D* profileToHist(TProfile* prof, const char* name) {
      int nbins = prof->GetNbinsX();
      double xmin = prof->GetXaxis()->GetXmin();
      double xmax = prof->GetXaxis()->GetXmax();

      TH1D* h = new TH1D(name, prof->GetTitle(), nbins, xmin, xmax);

      for (int i = 1; i <= nbins; ++i) {
          h->SetBinContent(i, prof->GetBinContent(i));
          h->SetBinError(i, prof->GetBinError(i));
      }

      return h;
  }

  /*
   * Function to generate a suffix for the data file based on the number of events,
    * maximum deflection (in nanometres), and maximum time (in picoseconds).
   */
  inline TString get_suffix(unsigned nEvents, unsigned max_scatter, unsigned max_dt) {
      return Form("_%uev_%unm_%ups", nEvents, max_scatter, max_dt);
  }

  std::vector<double> get_module_to_z_vector() {
  // Get module to z conversion histogram
  TFile *moduleFile = TFile::Open(
    (Utils::Definitions::analysisRoot + "hists/module_mc_info.root").c_str());
  if (!moduleFile || moduleFile->IsZombie()) {
      std::cerr << "Error: Could not open ROOT file." << std::endl;
      return {};
  }
  TH1D* moduleToZ = (TH1D*) moduleFile->Get("module_to_z");
  // TODO maybe at some point: Why does the below not work?
  // first entry in GetArray is underflow (last is overflow)
  // std::vector<double> moduleToZVec(moduleToZ->GetArray() + 1, 
  //                                  moduleToZ->GetArray() + 1 + moduleToZ->GetNbinsX());
  std::vector<double> moduleToZVec(moduleToZ->GetNbinsX());
  for (int i = 0; i < moduleToZ->GetNbinsX(); i++) {
    moduleToZVec[i] = moduleToZ->GetBinContent(i+1);
  }
  assert(moduleToZVec.size() == 64);
  // print for testing
  for (unsigned i = 0; i < moduleToZVec.size(); i++) {
    std::cout << "Module " << i << ": z=" << moduleToZVec[i] << "mm\n";
  }
  // clean up
  moduleFile->Close();
  delete moduleFile;

  return moduleToZVec;
  }
} // namespace Utils::Functions