#pragma once

#include <TCanvas.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TProfile.h>
#include <TString.h>

#include <iostream>
#include <string>
#include <vector>

#include "definitions.h"

/*
 * Drawing of the efficiency, ghost rate and clone rate canvas. Shared between
 * plot_efficiencies.C (one dataset) and compare_efficiencies.C (several of them), so that the
 * layout and styling only have to be changed in one place.
 */
namespace Utils::EfficiencyPlots {

  // has to stay in sync with the names written by get_eff_clone_ghosts.C
  const std::vector<std::string> kEfficiencyTypes = {"Pt", "Eta", "docaz"};
  const std::vector<int> kColors =
    {kMagenta, kBlue, kRed, kGreen + 2, kOrange + 7, kBlack, kCyan + 2};
  const std::vector<int> kMarkers = {20, 21, 22, 23, 33, 34, 29};

  /*
   * Open one of the files written by get_eff_clone_ghosts.C. histName is everything after
   * mc_hists_, i.e. the data suffix without its leading underscore.
   */
  inline TFile* openHistFile(const TString& histName) {
    TString path =
      TString((Utils::Definitions::analysisRoot + "/hists/eff_clone_ghosts/mc_hists_").c_str()) +
      histName + ".root";
    TFile* file = TFile::Open(path);
    if (!file || file->IsZombie()) {
      std::cerr << "Error: Could not open " << path << std::endl;
      return nullptr;
    }
    return file;
  }

  // ranges are per variable, they are only known once the pad has been painted
  inline void setEfficiencyRanges(TEfficiency* eff, const std::string& type,
                                  bool isForward, const float minEff) {
    TGraphAsymmErrors* graph = eff->GetPaintedGraph();
    if (!graph) return;  // only painted after the pad has been updated
    if (type == "Eta") {
      if (isForward) graph->GetXaxis()->SetRangeUser(2., 5.);
      else graph->GetXaxis()->SetRangeUser(-5., -2.);
      graph->GetYaxis()->SetRangeUser(minEff, 1.05);
    } else if (type == "docaz") {
      graph->GetXaxis()->SetRangeUser(0., 7.);
      graph->GetYaxis()->SetRangeUser(minEff, 1.05);
    } else if (type == "Pt") {
      graph->GetXaxis()->SetRangeUser(0., 5000.);
      graph->GetYaxis()->SetRangeUser(minEff, 1.05);
    }
  }

  /*
   * Overlay one variable of all given files into the current pad. With a single file this is
   * simply that file's efficiency, with the legend still telling which dataset it is.
   */
  inline void drawEfficiency(const std::vector<TFile*>& files,
                             const std::vector<TString>& labels,
                             const float minEff,
                             const std::string& type, bool isForward) {
    std::string region = (isForward ? "forward" : "backward");
    std::string histName = "efficiency_" + region + "_" + type;
    TLegend* legend = new TLegend(0.6, 0.15, 0.88, 0.35);
    legend->SetBorderSize(0);
    TEfficiency* firstDrawn = nullptr;
    for (unsigned i = 0; i < files.size(); i++) {
      TEfficiency* eff = (TEfficiency*) files[i]->Get(histName.c_str());
      if (!eff) {
        std::cerr << "Warning: " << histName << " not found in " << files[i]->GetName() << "\n";
        continue;
      }
      // same name in every file, rename so the pad keeps them apart
      eff->SetName(Form("%s_%u", histName.c_str(), i));
      int color = kColors[i % kColors.size()];
      eff->SetLineColor(color);
      eff->SetMarkerColor(color);
      eff->SetMarkerStyle(kMarkers[i % kMarkers.size()]);
      eff->SetMarkerSize(0.7);
      eff->Draw(firstDrawn ? "SAME P" : "AP");
      if (labels.size() > i)
        legend->AddEntry(eff, labels[i], "lp");
      if (!firstDrawn) firstDrawn = eff;
    }
    if (!firstDrawn) return;  // nothing in this pad
    gPad->Update();  // painted graph of the first one carries the axes
    setEfficiencyRanges(firstDrawn, type, isForward, minEff);
    if (labels.size() > 0) legend->Draw();
    gPad->Update();
  }

  /*
   * Ghost and clone rates are TProfiles rather than efficiencies, so they need their own
   * overlay. The eta profiles hold both regions, so their x range is what picks the region out.
   */
  inline void drawRate(const std::vector<TFile*>& files, const std::vector<TString>& labels,
                       const std::string& histName, const std::string& title,
                       bool isVsEta, bool isForward, double yMax) {
    TLegend* legend = new TLegend(0.6, 0.65, 0.88, 0.85);
    legend->SetBorderSize(0);
    bool anyDrawn = false;
    for (unsigned i = 0; i < files.size(); i++) {
      TProfile* rate = (TProfile*) files[i]->Get(histName.c_str());
      if (!rate) {
        std::cerr << "Warning: " << histName << " not found in " << files[i]->GetName() << "\n";
        continue;
      }
      rate->SetName(Form("%s_%u", histName.c_str(), i));
      int color = kColors[i % kColors.size()];
      rate->SetLineColor(color);
      rate->SetMarkerColor(color);
      rate->SetMarkerStyle(kMarkers[i % kMarkers.size()]);
      rate->SetMarkerSize(0.7);
      rate->SetTitle(title.c_str());
      rate->GetYaxis()->SetRangeUser(0., yMax);
      if (isVsEta) {
        if (isForward) rate->GetXaxis()->SetRangeUser(1., 5.5);
        else rate->GetXaxis()->SetRangeUser(-5.5, -1.);
      } else {
        rate->GetXaxis()->SetRangeUser(0., 5000.);
      }
      rate->Draw(anyDrawn ? "SAME" : "");
      
      if (labels.size() > i)
        legend->AddEntry(rate, labels[i], "lp");
      anyDrawn = true;
    }
    if (anyDrawn && labels.size() > 0) legend->Draw();
  }

  /*
   * The whole canvas for one region: efficiencies on the top row, ghost and clone rates on the
   * bottom one.
   */
  inline void drawRegion(const std::vector<TFile*>& files,
                         const std::vector<TString>& labels,
                         const std::vector<float>& minEffs,
                         bool isForward, TCanvas* canvas) {
    canvas->Divide(3, 2);
    for (unsigned typeIndex = 0; typeIndex < kEfficiencyTypes.size(); typeIndex++) {
      canvas->cd(typeIndex + 1);
      drawEfficiency(files, labels, minEffs[typeIndex],
                     kEfficiencyTypes[typeIndex], isForward);
    }
    canvas->cd(4);
    drawRate(files, labels, "ghost_rates", "Ghost Rates;#eta;Ghost Rate", true, isForward, 0.1);
    canvas->cd(5);
    drawRate(files, labels, "clone_rate", "Clone Rates;#eta;Clone Rate", true, isForward, 0.1);
    canvas->cd(6);
    drawRate(files, labels, "clone_rate_pt", "Clone Rates;p_{T} (MeV);Clone Rate",
             false, isForward, 0.1);
  }

  /*
   * Write one PDF per region, named <outputBase>_forward.pdf and <outputBase>_backward.pdf.
   */
  inline void drawAndSave(const std::vector<TFile*>& files,
                          const std::vector<TString>& labels,
                          const std::vector<float>& minEffs,
                          const TString& outputBase) {
    for (bool isForward : {true, false}) {
      TString region = (isForward ? "forward" : "backward");
      TCanvas* canvas =
        new TCanvas("canvas_" + region, region + " Region Efficiencies", 1500, 800);
      drawRegion(files, labels, minEffs, isForward, canvas);
      canvas->SaveAs(outputBase + "_" + region + ".pdf");
      delete canvas;
    }
    std::cout << "Wrote " << outputBase << "_{forward,backward}.pdf" << std::endl;
  }

}  // namespace Utils::EfficiencyPlots
