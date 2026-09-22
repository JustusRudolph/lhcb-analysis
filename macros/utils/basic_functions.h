#pragma once

#include <TH1D.h>
#include <TString.h>

#include <bitset>
#include <vector>

#include "definitions.h"
#include "Hit.h"

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

  std::vector<double> get_module_to_z_vector(bool print=false) {
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
  if (print) {
    for (unsigned i = 0; i < moduleToZVec.size(); i++) {
      std::cout << "Module " << i << ": z=" << moduleToZVec[i] << "mm\n";
    }
  }
  // clean up
  moduleFile->Close();
  delete moduleFile;

  return moduleToZVec;
  }

  /*
   * Return only the first hit found on each module, keeping the order of the input.
   * The module ID is the 6-bit field at bit 12 of the LHCbID, so hits sharing a module
   * also sit on the same z-plane. Two such hits in a row make the dz of an
   * extrapolation zero, which is what we want to keep out of the forwarding dts.
   */
  inline std::vector<Hit::BaseHit> get_unique_module_hits(const std::vector<Hit::BaseHit>& hits) {
    std::bitset<Utils::Definitions::kModules> seen;
    std::vector<Hit::BaseHit> unique_hits;
    unique_hits.reserve(hits.size());
    for (const Hit::BaseHit& hit : hits) {
      unsigned module = (hit.id >> 12) & 0x3F;  // 64 modules, i.e. 3F
      if (seen.test(module)) continue;
      seen.set(module);
      unique_hits.push_back(hit);
    }
    return unique_hits;
  }
} // namespace Utils::Functions