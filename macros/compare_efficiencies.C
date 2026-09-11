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

#include "utils/definitions.h"
#include "utils/basic_functions.h"

// has to stay in sync with the names written by get_eff_clone_ghosts.C
const std::vector<std::string> kEfficiencyTypes = {"Pt", "Eta", "docaz"};
const std::vector<int> kCompareColors =
  {kMagenta, kBlue, kRed, kGreen + 2, kOrange + 7, kBlack, kCyan + 2};
const std::vector<int> kCompareMarkers = {20, 21, 22, 23, 33, 34, 29};

// ranges are per variable, they are only known once the pad has been painted
void set_efficiency_ranges(TEfficiency* eff, const std::string& type, bool isForward) {
  TGraphAsymmErrors* graph = eff->GetPaintedGraph();
  if (!graph) return;  // only painted after the pad has been updated
  if (type == "Eta") {
    if (isForward) graph->GetXaxis()->SetRangeUser(2., 5.);
    else graph->GetXaxis()->SetRangeUser(-5., -2.);
    graph->GetYaxis()->SetRangeUser(0.8, 1.01);
  } else if (type == "docaz") {
    graph->GetXaxis()->SetRangeUser(0., 7.);
    graph->GetYaxis()->SetRangeUser(0., 1.01);
  } else if (type == "Pt") {
    graph->GetXaxis()->SetRangeUser(0., 5000.);
    graph->GetYaxis()->SetRangeUser(0.8, 1.01);
  }
}

/*
 * Overlay one variable of all given files into the current pad.
 */
void draw_efficiency_comparison(const std::vector<TFile*>& files,
                                const std::vector<TString>& labels,
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
    int color = kCompareColors[i % kCompareColors.size()];
    eff->SetLineColor(color);
    eff->SetMarkerColor(color);
    eff->SetMarkerStyle(kCompareMarkers[i % kCompareMarkers.size()]);
    eff->SetMarkerSize(0.7);
    eff->Draw(firstDrawn ? "SAME P" : "AP");
    legend->AddEntry(eff, labels[i], "lp");
    if (!firstDrawn) firstDrawn = eff;
  }
  if (!firstDrawn) return;  // nothing in this pad
  gPad->Update();  // painted graph of the first one carries the axes
  set_efficiency_ranges(firstDrawn, type, isForward);
  legend->Draw();
  gPad->Update();
}

/*
 * Ghost and clone rates are TProfiles rather than efficiencies, so they need their own overlay.
 * The eta profiles hold both regions, so their x range is what picks the region out.
 */
void draw_rate_comparison(const std::vector<TFile*>& files, const std::vector<TString>& labels,
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
    int color = kCompareColors[i % kCompareColors.size()];
    rate->SetLineColor(color);
    rate->SetMarkerColor(color);
    rate->SetMarkerStyle(kCompareMarkers[i % kCompareMarkers.size()]);
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
    legend->AddEntry(rate, labels[i], "lp");
    anyDrawn = true;
  }
  if (anyDrawn) legend->Draw();
}

/*
 * Compare the efficiencies, ghost and clone rates of two or more datasets written by
 * get_eff_clone_ghosts.C. histNames are file names (without .root) inside hists/eff_clone_ghosts,
 * labels
 * are what ends up in the legend and default to the file names. One PDF is written per region.
 *
 * e.g. root -l 'compare_efficiencies.C({"mc_hists_5000ev_20nm_0ps",
 *                                       "mc_hists_5000ev_80000nm_0ps"}, {"20nm", "80000nm"})'
 */
void compare_efficiencies(std::vector<TString> histNames, std::vector<TString> labels={},
                          TString outName="") {
  gROOT->SetBatch();  // so stuff isn't autoplotted
  gStyle->SetOptStat(0);  // remove the info box
  if (histNames.size() < 2) {
    std::cerr << "Error: give at least two hist names to compare." << std::endl;
    return;
  }
  TString histPath = (Utils::Definitions::analysisRoot + "/hists/eff_clone_ghosts/mc_hists_").c_str();
  std::vector<TFile*> files{};
  std::vector<TString> usedLabels{};
  for (unsigned i = 0; i < histNames.size(); i++) {
    TFile* file = TFile::Open(histPath + histNames[i] + ".root");
    if (!file || file->IsZombie()) {
      std::cerr << "Error: Could not open " << histPath + histNames[i] + ".root" << std::endl;
      continue;
    }
    files.push_back(file);
    // fall back to the file name if no (or not enough) labels were given
    usedLabels.push_back(i < labels.size() ? labels[i] : histNames[i]);
  }
  if (files.size() < 2) {
    std::cerr << "Error: less than two of the given hists could be opened." << std::endl;
    return;
  }
  if (outName.IsNull()) {
    outName = histNames[0];
    for (unsigned i = 1; i < histNames.size(); i++) outName += "_vs_" + histNames[i];
  }
  TString outputBase =
    TString((Utils::Definitions::analysisRoot + "/output/efficiency/").c_str()) + outName;

  for (bool isForward : {true, false}) {
    TString region = (isForward ? "forward" : "backward");
    TCanvas* canvas = new TCanvas("canvas_" + region, region + " Region Efficiencies", 1500, 800);
    canvas->Divide(3, 2);
    // top row: efficiencies, bottom row: ghost and clone rates
    for (unsigned typeIndex = 0; typeIndex < kEfficiencyTypes.size(); typeIndex++) {
      canvas->cd(typeIndex + 1);
      draw_efficiency_comparison(files, usedLabels, kEfficiencyTypes[typeIndex], isForward);
    }
    canvas->cd(4);
    draw_rate_comparison(files, usedLabels, "ghost_rates", "Ghost Rates;#eta;Ghost Rate",
                         true, isForward, 0.1);
    canvas->cd(5);
    draw_rate_comparison(files, usedLabels, "clone_rate", "Clone Rates;#eta;Clone Rate",
                         true, isForward, 0.1);
    canvas->cd(6);
    draw_rate_comparison(files, usedLabels, "clone_rate_pt",
                         "Clone Rates;p_{T} (MeV);Clone Rate", false, isForward, 0.1);
    canvas->SaveAs(outputBase + "_" + region + ".pdf");
    delete canvas;
  }
  std::cout << "Wrote comparison to " << outputBase << "_{forward,backward}.pdf" << std::endl;

  for (TFile* file : files) {
    file->Close();
    delete file;
  }
}
