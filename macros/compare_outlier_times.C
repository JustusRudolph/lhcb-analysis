#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TString.h>
#include <TSystem.h>
#include <iostream>
#include <string>
#include <vector>

#include "utils/definitions.h"
#include "utils/basic_functions.h"

// colours of the four curves in the seeding and forwarding pads, as in plot_outlier_times.C
const std::vector<int> kSeedingColors = {kBlue, kRed, kBlue, kRed};
const std::vector<int> kForwardingColors = {kBlack, kRed, kBlue, kGreen};

/*
 * Shape ratio of one histogram between two files. Both are normalised to unit integral first,
 * so the ratio shows where the distribution moved rather than how much statistics each run had.
 * Returns nullptr if the histogram is missing or empty in either file.
 */
template <typename T>
T* make_ratio(TFile* numerFile, TFile* denomFile, const std::string& name, unsigned index) {
  T* numer = (T*) numerFile->Get(name.c_str());
  T* denom = (T*) denomFile->Get(name.c_str());
  if (!numer || !denom) {
    std::cerr << "Warning: " << name << " missing in one of the files, skipping.\n";
    return nullptr;
  }
  if (numer->Integral() == 0 || denom->Integral() == 0) {
    std::cerr << "Warning: " << name << " is empty in one of the files, skipping.\n";
    return nullptr;
  }
  T* ratio = (T*) numer->Clone(Form("ratio_%s_%u", name.c_str(), index));
  ratio->SetDirectory(nullptr);
  ratio->Scale(1.0 / numer->Integral());
  T* denomNormalised = (T*) denom->Clone(Form("denom_%s_%u", name.c_str(), index));
  denomNormalised->SetDirectory(nullptr);
  denomNormalised->Scale(1.0 / denom->Integral());
  ratio->Divide(denomNormalised);
  delete denomNormalised;
  return ratio;
}

/*
 * Overlay the ratios of a group of 1D histograms into the current pad, with a line at unity.
 */
void draw_ratio_group(TFile* numerFile, TFile* denomFile, const std::vector<std::string>& names,
                      const std::vector<std::string>& entryLabels, const std::vector<int>& colors,
                      const std::string& title, double xMin, double xMax,
                      unsigned nSolidCurves) {
  TLegend* legend = new TLegend(0.6, 0.65, 0.88, 0.85);
  legend->SetBorderSize(0);
  TH1D* firstDrawn = nullptr;
  for (unsigned i = 0; i < names.size(); i++) {
    TH1D* ratio = make_ratio<TH1D>(numerFile, denomFile, names[i], i);
    if (!ratio) continue;
    // curves past nSolidCurves reuse a colour, so they are drawn transparent as before
    if (i < nSolidCurves) ratio->SetLineColor(colors[i % colors.size()]);
    else ratio->SetLineColorAlpha(colors[i % colors.size()], 0.35);
    ratio->Draw(firstDrawn ? "HIST SAME" : "HIST");
    legend->AddEntry(ratio, entryLabels[i].c_str(), "l");
    if (!firstDrawn) firstDrawn = ratio;
  }
  if (!firstDrawn) return;
  firstDrawn->SetTitle(title.c_str());
  firstDrawn->GetYaxis()->SetTitle("Ratio");
  firstDrawn->GetYaxis()->SetRangeUser(0., 2.);
  firstDrawn->GetXaxis()->SetRangeUser(xMin, xMax);
  // unity line to read the deviation against
  TLine* unity = new TLine(xMin, 1., xMax, 1.);
  unity->SetLineStyle(2);
  unity->SetLineColor(kGray + 2);
  unity->Draw();
  legend->Draw();
}

/*
 * Draw the ratio of one 2D histogram. A ratio has no reason to be log scaled, unlike the
 * distributions themselves.
 */
void draw_ratio_2d(TFile* numerFile, TFile* denomFile, const std::string& name,
                   const std::string& title) {
  TH2D* ratio = make_ratio<TH2D>(numerFile, denomFile, name, 0);
  if (!ratio) return;
  ratio->SetTitle(title.c_str());
  ratio->GetZaxis()->SetRangeUser(0., 2.);
  ratio->Draw("COLZ");
}

/*
 * Compare two datasets written by get_outlier_times.C by taking the ratio of every figure
 * plot_outlier_times.C draws. histNames are file names (without .root) inside
 * hists/4d_tracking, labels default to those names. Every ratio is first over second.
 *
 * e.g. root -l 'compare_outlier_times.C({"outlier_times_20nm", "outlier_times_80000nm"},
 *                                       {"20nm", "80000nm"})'
 */
void compare_outlier_times(std::vector<TString> histNames, std::vector<TString> labels={},
                           TString outName="") {
  gROOT->SetBatch();  // so stuff isn't autoplotted
  gStyle->SetOptStat(0); // remove the info box for the plots
  if (histNames.size() != 2) {
    std::cerr << "Error: give exactly two hist names, the ratio is the first over the second."
              << std::endl;
    return;
  }
  TString histPath = (Utils::Definitions::analysisRoot + "hists/4d_tracking/").c_str();
  std::vector<TFile*> files{};
  for (unsigned i = 0; i < histNames.size(); i++) {
    TFile* file = TFile::Open(histPath + histNames[i] + ".root");
    if (!file || file->IsZombie()) {
      std::cerr << "Error: Could not open " << histPath + histNames[i] + ".root" << std::endl;
      return;
    }
    files.push_back(file);
  }
  // fall back to the file names if no (or not enough) labels were given
  std::vector<TString> usedLabels{};
  for (unsigned i = 0; i < histNames.size(); i++)
    usedLabels.push_back(i < labels.size() ? labels[i] : histNames[i]);
  std::string ratioOf =
    " (" + std::string(usedLabels[0].Data()) + " / " + std::string(usedLabels[1].Data()) + ")";

  TCanvas* canvas = new TCanvas("canvas", "Outlier Times Ratios", 1200, 1200);
  canvas->Divide(2, 3);

  canvas->cd(1);
  draw_ratio_group(files[0], files[1],
                   {"h_dt0_forward", "h_dt2_forward", "h_dt0_backward", "h_dt2_backward"},
                   {"#Deltat_{0}^{f}", "#Deltat_{2}^{f}", "#Deltat_{0}^{b}", "#Deltat_{2}^{b}"},
                   kSeedingColors, "Seeding time scatter ratio" + ratioOf + ";#Deltat (ns)",
                   -0.3, 0.3, 2);  // forward curves solid, backward ones transparent

  canvas->cd(2);
  draw_ratio_group(files[0], files[1],
                   {"h_dt_forwarding_forward", "h_dt_forwarding_forward_h5",
                    "h_dt_forwarding_forward_h10", "h_dt_forwarding_forward_h15"},
                   {"All forwarded", "Hit index 5", "Hit index 10", "Hit index 15"},
                   kForwardingColors,
                   "#Deltat forwarding ratio, forward" + ratioOf + ";#Deltat_{forward} (ns)",
                   -0.5, 0.5, 4);

  canvas->cd(3);
  draw_ratio_2d(files[0], files[1], "h_dt_forwarding_vs_moduleID",
                "#Deltat vs module ID ratio" + ratioOf + ";Module ID;#Deltat (ns)");

  canvas->cd(4);
  draw_ratio_group(files[0], files[1],
                   {"h_dt_forwarding_backward", "h_dt_forwarding_backward_h5",
                    "h_dt_forwarding_backward_h10", "h_dt_forwarding_backward_h15"},
                   {"All forwarded", "Hit index 5", "Hit index 10", "Hit index 15"},
                   kForwardingColors,
                   "#Deltat forwarding ratio, backward" + ratioOf + ";#Deltat_{forward} (ns)",
                   -0.5, 0.5, 4);

  canvas->cd(5);
  draw_ratio_2d(files[0], files[1], "h_nthHits_vs_moduleID",
                "Hit number vs module ID ratio" + ratioOf + ";Module ID;Hit number in track");

  canvas->cd(6);
  draw_ratio_2d(files[0], files[1], "h_dt_forwarding_vs_nthHit",
                "#Deltat vs hit number ratio" + ratioOf + ";Nth hit;#Deltat (ns)");

  if (outName.IsNull()) outName = "outlier_times_" + histNames[0] + "_vs_" + histNames[1];
  TString outputDir = (Utils::Definitions::analysisRoot + "output/4d_tracking").c_str();
  gSystem->mkdir(outputDir, true);
  TString outputPath = outputDir + "/" + outName + ".pdf";
  canvas->SaveAs(outputPath);
  std::cout << "Canvas saved to: " << outputPath << std::endl;

  delete canvas;
  for (TFile* file : files) {
    file->Close();
    delete file;
  }
}
