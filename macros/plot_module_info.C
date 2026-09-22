#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TString.h>
#include <TSystem.h>

#include <iostream>

#include "utils/definitions.h"
#include "utils/basic_functions.h"

/*
 * Draw the hit timestamp distribution and the mean timestamp per module written by
 * get_module_info.C. Both come from the same fake cluster sample, the first as a raw
 * distribution and the second profiled against module ID.
 */
void plot_module_info(bool logY=true) {
  gROOT->SetBatch();  // so stuff isn't autoplotted
  gStyle->SetOptStat(0);  // remove the info box
  TString histPath = (Utils::Definitions::analysisRoot + "hists/module_mc_info.root").c_str();
  TFile* file = TFile::Open(histPath);
  if (!file || file->IsZombie()) {
    std::cerr << "Error: Could not open " << histPath << std::endl;
    return;
  }
  TH1D* hitTimes = (TH1D*) file->Get("hit_times");
  TProfile* moduleTimes = (TProfile*) file->Get("module_times");
  if (!hitTimes || !moduleTimes) {
    std::cerr << "Error: hit_times or module_times not found in " << histPath
              << ", re-run get_module_info.C" << std::endl;
    return;
  }
  TString outputDir = (Utils::Definitions::analysisRoot + "output/hits").c_str();
  gSystem->mkdir(outputDir, true);

  TCanvas* c_times = new TCanvas("c_times", "Hit Timestamps", 800, 600);
  hitTimes->SetLineColor(kBlue);
  hitTimes->Draw("HIST");
  // the distribution has a long tail, so a log y is usually the readable choice
  if (logY) gPad->SetLogy();
  c_times->SaveAs(outputDir + "/hit_times.pdf");

  TCanvas* c_module_times = new TCanvas("c_module_times", "Mean Time per Module", 800, 600);
  moduleTimes->SetLineColor(kBlue);
  moduleTimes->SetMarkerColor(kBlue);
  moduleTimes->SetMarkerStyle(20);
  moduleTimes->SetMarkerSize(0.7);
  moduleTimes->Draw("E1");
  c_module_times->SaveAs(outputDir + "/module_times.pdf");

  std::cout << "Wrote " << outputDir << "/hit_times.pdf and module_times.pdf" << std::endl;
  delete c_times;
  delete c_module_times;
  file->Close();
  delete file;
}
