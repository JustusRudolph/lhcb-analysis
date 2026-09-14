#pragma once

#include <TCanvas.h>
#include <TEfficiency.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TProfile.h>
#include <TString.h>

#include <algorithm>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
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

  // the x range that is drawn, per variable and region
  inline std::pair<double, double> efficiencyXRange(const std::string& type, bool isForward) {
    if (type == "Eta") return isForward ? std::make_pair(2., 5.) : std::make_pair(-5., -2.);
    if (type == "docaz") return std::make_pair(0., 7.);
    if (type == "Pt") return std::make_pair(0., 5000.);
    return std::make_pair(0., 0.);
  }

  /*
   * The y range that fits everything drawn in the pad: 0.05 below the lowest and 0.05 above the
   * highest efficiency of any of the datasets, with the top at least 1.01 so that a perfect
   * efficiency stays visible. Only bins inside the drawn x range and with something
   * reconstructible in them count, empty bins would otherwise pull the bottom down to zero.
   */
  inline std::pair<double, double> efficiencyYRange(const std::vector<TEfficiency*>& effs,
                                                    double xMin, double xMax) {
    double lowest = 1., highest = 0.;
    bool anyFilled = false;
    for (TEfficiency* eff : effs) {
      const TH1* total = eff->GetTotalHistogram();
      for (int bin = 1; bin <= total->GetNbinsX(); bin++) {
        if (total->GetBinContent(bin) == 0) continue;  // nothing reconstructible here
        double x = total->GetBinCenter(bin);
        if (x < xMin || x > xMax) continue;  // outside of what is drawn
        double value = eff->GetEfficiency(bin);
        lowest = std::min(lowest, value);
        highest = std::max(highest, value);
        anyFilled = true;
      }
    }
    if (!anyFilled) return std::make_pair(0., 1.01);
    return std::make_pair(lowest - 0.05, std::max(highest + 0.05, 1.01));
  }

  // the x range that is drawn for the rates, the eta profiles hold both regions
  inline std::pair<double, double> rateXRange(bool isVsEta, bool isForward) {
    if (!isVsEta) return std::make_pair(0., 5000.);
    return isForward ? std::make_pair(1., 5.5) : std::make_pair(-5.5, -1.);
  }

  /*
   * Same idea as efficiencyYRange, but the bottom is kept at -0.01 at the lowest so that rates
   * sitting at zero do not push the axis far into the negatives. Bins without entries are
   * skipped, a TProfile reports those as zero.
   */
  inline std::pair<double, double> rateYRange(const std::vector<TProfile*>& rates,
                                              double xMin, double xMax) {
    double lowest = std::numeric_limits<double>::max();
    double highest = std::numeric_limits<double>::lowest();
    bool anyFilled = false;
    for (TProfile* rate : rates) {
      for (int bin = 1; bin <= rate->GetNbinsX(); bin++) {
        if (rate->GetBinEntries(bin) == 0) continue;  // nothing filled this bin
        double x = rate->GetBinCenter(bin);
        if (x < xMin || x > xMax) continue;  // outside of what is drawn
        double value = rate->GetBinContent(bin);
        lowest = std::min(lowest, value);
        highest = std::max(highest, value);
        anyFilled = true;
      }
    }
    if (!anyFilled) return std::make_pair(-0.01, 0.1);
    return std::make_pair(std::max(lowest - 0.05, -0.01), highest + 0.05);
  }

  // ranges are per variable, they are only known once the pad has been painted
  inline void setEfficiencyRanges(TEfficiency* eff, const std::string& type, bool isForward,
                                  const std::pair<double, double>& yRange) {
    TGraphAsymmErrors* graph = eff->GetPaintedGraph();
    if (!graph) return;  // only painted after the pad has been updated
    std::pair<double, double> xRange = efficiencyXRange(type, isForward);
    graph->GetXaxis()->SetRangeUser(xRange.first, xRange.second);
    graph->GetYaxis()->SetRangeUser(yRange.first, yRange.second);
  }

  /*
   * Overlay one variable of all given files into the current pad. With a single file this is
   * simply that file's efficiency, with the legend still telling which dataset it is.
   */
  inline void drawEfficiency(const std::vector<TFile*>& files,
                             const std::vector<TString>& labels,
                             const std::string& type, bool isForward) {
    std::string region = (isForward ? "forward" : "backward");
    std::string histName = "efficiency_" + region + "_" + type;
    TLegend* legend = new TLegend(0.6, 0.15, 0.88, 0.35);
    legend->SetBorderSize(0);
    TEfficiency* firstDrawn = nullptr;
    std::vector<TEfficiency*> drawn{};  // all of them decide the y range
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
      drawn.push_back(eff);
      if (!firstDrawn) firstDrawn = eff;
    }
    if (!firstDrawn) return;  // nothing in this pad
    gPad->Update();  // painted graph of the first one carries the axes
    std::pair<double, double> xRange = efficiencyXRange(type, isForward);
    setEfficiencyRanges(firstDrawn, type, isForward,
                        efficiencyYRange(drawn, xRange.first, xRange.second));
    if (labels.size() > 0) legend->Draw();
    gPad->Update();
  }

  /*
   * Ghost and clone rates are TProfiles rather than efficiencies, so they need their own
   * overlay. The eta profiles hold both regions, so their x range is what picks the region out.
   */
  inline void drawRate(const std::vector<TFile*>& files, const std::vector<TString>& labels,
                       const std::string& histName, const std::string& title,
                       bool isVsEta, bool isForward) {
    TLegend* legend = new TLegend(0.6, 0.65, 0.88, 0.85);
    legend->SetBorderSize(0);
    // fetch and style first, the y range needs all of them before anything is drawn
    std::vector<TProfile*> rates{};
    std::vector<unsigned> fileIndices{};
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
      rates.push_back(rate);
      fileIndices.push_back(i);
    }
    if (rates.empty()) return;  // nothing in this pad
    std::pair<double, double> xRange = rateXRange(isVsEta, isForward);
    std::pair<double, double> yRange = rateYRange(rates, xRange.first, xRange.second);
    for (unsigned i = 0; i < rates.size(); i++) {
      rates[i]->GetXaxis()->SetRangeUser(xRange.first, xRange.second);
      rates[i]->GetYaxis()->SetRangeUser(yRange.first, yRange.second);
      rates[i]->Draw(i ? "SAME" : "");
      if (labels.size() > fileIndices[i])
        legend->AddEntry(rates[i], labels[fileIndices[i]], "lp");
    }
    if (labels.size() > 0) legend->Draw();
  }

  /*
   * The whole canvas for one region: efficiencies on the top row, ghost and clone rates on the
   * bottom one.
   */
  inline void drawRegion(const std::vector<TFile*>& files,
                         const std::vector<TString>& labels,
                         bool isForward, TCanvas* canvas) {
    canvas->Divide(3, 2);
    for (unsigned typeIndex = 0; typeIndex < kEfficiencyTypes.size(); typeIndex++) {
      canvas->cd(typeIndex + 1);
      drawEfficiency(files, labels, kEfficiencyTypes[typeIndex], isForward);
    }
    canvas->cd(4);
    drawRate(files, labels, "ghost_rates", "Ghost Rates;#eta;Ghost Rate", true, isForward);
    canvas->cd(5);
    drawRate(files, labels, "clone_rate", "Clone Rates;#eta;Clone Rate", true, isForward);
    canvas->cd(6);
    drawRate(files, labels, "clone_rate_pt", "Clone Rates;p_{T} (MeV);Clone Rate",
             false, isForward);
  }

  /*
   * Write one PDF per region, named <outputBase>_forward.pdf and <outputBase>_backward.pdf.
   */
  inline void drawAndSave(const std::vector<TFile*>& files,
                          const std::vector<TString>& labels,
                          const TString& outputBase) {
    for (bool isForward : {true, false}) {
      TString region = (isForward ? "forward" : "backward");
      TCanvas* canvas =
        new TCanvas("canvas_" + region, region + " Region Efficiencies", 1500, 800);
      drawRegion(files, labels, isForward, canvas);
      canvas->SaveAs(outputBase + "_" + region + ".pdf");
      delete canvas;
    }
    std::cout << "Wrote " << outputBase << "_{forward,backward}.pdf" << std::endl;
  }

}  // namespace Utils::EfficiencyPlots
