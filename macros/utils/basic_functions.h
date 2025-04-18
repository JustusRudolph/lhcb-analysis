#pragma once

#include <TH1D.h>

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
} // namespace Utils::Functions