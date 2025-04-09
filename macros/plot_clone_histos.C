#include <TFile.h>
#include <TH1F.h>
#include <TEfficiency.h>
#include <TCanvas.h>

#include <cstdlib>
#include <iostream>
#include <string>

void disableProfileMarkers(TProfile* p) {
    p->SetMarkerStyle(0);
    p->SetMarkerSize(0);
    p->SetMarkerColor(0);
}

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
  TProfile* ghostRates = (TProfile*) file->Get("ghost_rates");
  TProfile* seedingClonesMCByEta = (TProfile*) file->Get("clone_types_mc_by_eta");
  TProfile* tripletClonesMCByEta = (TProfile*) file->Get("triplet_clones_mc_by_eta");
  TProfile* moduleOverlapClonesMCByEta = (TProfile*) file->Get("module_overlap_clones_mc_by_eta");
  TProfile* otherClonesMCByEta = (TProfile*) file->Get("other_clones_mc_by_eta");
  TProfile* seedingClonesRecoByEta = (TProfile*) file->Get("clone_types_reco_by_eta");
  TProfile* tripletClonesRecoByEta = (TProfile*) file->Get("triplet_clones_reco_by_eta");
  TProfile* moduleOverlapClonesRecoByEta = (TProfile*) file->Get("module_overlap_clones_reco_by_eta");
  TProfile* otherClonesRecoByEta = (TProfile*) file->Get("other_clones_reco_by_eta");

  // Now plot
  TCanvas* canvas = new TCanvas("c", "MC Clone Information", 1200, 1800);
  canvas->Divide(2, 4);

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
  TLegend* legend = new TLegend(0.2, 0.2, 0.4, 0.4);
  legend->AddEntry(hLongestMatchedTrackRate, "No clones", "l");
  legend->AddEntry(hLongestMatchedTrackRateClones, "1 #leq N_{clones} #leq 5", "l");
  legend->AddEntry(hLongestMatchedTrackRate5Clones, "N_{clones} > 5", "l");
  legend->Draw();

  canvas->cd(4);
  hLongestTrackTagged->Draw();
  canvas->cd(5);
  hNMatches->Draw();
  canvas->cd(6);
  ghostRates->SetLineColor(kMagenta);
  ghostRates->GetYaxis()->SetRangeUser(0., 0.2);
  ghostRates->Draw();

  canvas->cd(7);
  seedingClonesMCByEta->SetTitle("MC particle clone rate by #eta");
  seedingClonesMCByEta->SetLineColor(kRed);
  disableProfileMarkers(seedingClonesMCByEta);
  seedingClonesMCByEta->Draw();
  tripletClonesMCByEta->SetLineColor(kBlue);
  disableProfileMarkers(tripletClonesMCByEta);
  tripletClonesMCByEta->Draw("SAME");
  moduleOverlapClonesMCByEta->SetLineColor(kGreen);
  disableProfileMarkers(moduleOverlapClonesMCByEta);
  moduleOverlapClonesMCByEta->Draw("SAME");
  otherClonesMCByEta->SetLineColor(kBlack);
  disableProfileMarkers(otherClonesMCByEta);
  otherClonesMCByEta->Draw("SAME");
  // Add legend
  TLegend* cloneLegend = new TLegend(0.2, 0.3, 0.4, 0.5);
  cloneLegend->AddEntry(seedingClonesMCByEta, "Seeding", "l");
  cloneLegend->AddEntry(tripletClonesMCByEta, "Triplet", "l");
  cloneLegend->AddEntry(moduleOverlapClonesMCByEta, "Module overlap", "l");
  cloneLegend->AddEntry(otherClonesMCByEta, "Other", "l");
  cloneLegend->Draw();

  canvas->cd(8);
  seedingClonesRecoByEta->SetTitle("Reconstructed track clone rates by #eta");
  seedingClonesRecoByEta->SetLineColor(kRed);
  disableProfileMarkers(seedingClonesRecoByEta);
  seedingClonesRecoByEta->Draw();
  tripletClonesRecoByEta->SetLineColor(kBlue);
  disableProfileMarkers(tripletClonesRecoByEta);
  tripletClonesRecoByEta->Draw("SAME");
  moduleOverlapClonesRecoByEta->SetLineColor(kGreen);
  disableProfileMarkers(moduleOverlapClonesRecoByEta);
  moduleOverlapClonesRecoByEta->Draw("SAME");
  otherClonesRecoByEta->SetLineColor(kBlack);
  disableProfileMarkers(otherClonesRecoByEta);
  otherClonesRecoByEta->Draw("SAME");
  // Add legend
  TLegend* recoCloneLegend = new TLegend(0.2, 0.3, 0.4, 0.5);
  recoCloneLegend->AddEntry(seedingClonesRecoByEta, "Seeding", "l");
  recoCloneLegend->AddEntry(tripletClonesRecoByEta, "Triplet", "l");
  recoCloneLegend->AddEntry(moduleOverlapClonesRecoByEta, "Module overlap", "l");
  recoCloneLegend->AddEntry(otherClonesRecoByEta, "Other", "l");
  recoCloneLegend->Draw();

  canvas->SaveAs((analysisRoot + "/output/mc_clone_plots.pdf").c_str());

  // clean up
  file->Close();
  delete file;
  delete canvas;
  delete legend;
}