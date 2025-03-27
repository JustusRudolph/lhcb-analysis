#include <TFile.h>
#include <TH1F.h>
#include <TEfficiency.h>
#include <TCanvas.h>

#include <cstdlib>
#include <iostream>
#include <string>

void plot_clone_histos() {
  gROOT->SetBatch();  // so stuff isn't autoplotted
  std::string analysisRoot = std::getenv("ANALYSIS_ROOT");
  TFile* file = TFile::Open((analysisRoot + "/hists/clones/mc_hists.root").c_str());

  TProfile* hDuplicateIDRates = (TProfile*) file->Get("duplicate_match_id");
  TProfile* hUniqueIDRates = (TProfile*) file->Get("unique_id_match_rate");
  TProfile* hLongestMatchedTrackRate = (TProfile*) file->Get("longest_matched_track_fraction");
  TProfile* hLongestMatchedTrackRateClones = (TProfile*) file->Get("longest_matched_track_fraction_1_clones");
  TProfile* hLongestMatchedTrackRate5Clones = (TProfile*) file->Get("longest_matched_track_fraction_5_clones");
  TProfile* hLongestTrackTagged = (TProfile*) file->Get("longest_match_tag_rate");
  TProfile* hNMatches = (TProfile*) file->Get("number_of_matches_of_MC_track");

  // Now plot
  TCanvas* canvas = new TCanvas("c", "MC Clone Information", 1000, 1200);
  canvas->Divide(2, 3);

  gStyle->SetOptStat(0);  // remove the info box
  canvas->cd(1);
  hDuplicateIDRates->Draw();
  canvas->cd(2);
  hUniqueIDRates->Draw();

  canvas->cd(3);
  // set line colours and plot all histograms
  hLongestMatchedTrackRate->SetLineColor(kBlack);
  hLongestMatchedTrackRate->Draw();

  hLongestMatchedTrackRateClones->SetLineColor(kMagenta);
  hLongestMatchedTrackRateClones->Draw("SAME");

  hLongestMatchedTrackRate5Clones->SetLineColor(kBlue);
  hLongestMatchedTrackRate5Clones->Draw("SAME");
  // Add legend
  TLegend* legend = new TLegend(0.4, 0.4, 0.2, 0.2);
  legend->AddEntry(hLongestMatchedTrackRate, "No clones", "l");
  legend->AddEntry(hLongestMatchedTrackRateClones, "1 #leq N_{clones} #leq 5", "l");
  legend->AddEntry(hLongestMatchedTrackRate5Clones, "N_{clones} #geq 5", "l");
  legend->Draw();

  canvas->cd(4);
  hLongestTrackTagged->Draw();
  canvas->cd(5);
  hNMatches->Draw();

  canvas->SaveAs((analysisRoot + "/output/mc_clone_plots.pdf").c_str());

  // clean up
  file->Close();
  delete file;
  delete canvas;
  delete legend;
}