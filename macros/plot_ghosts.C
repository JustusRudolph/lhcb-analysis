#include <TFile.h>
#include <TH1F.h>
#include <TString.h>

#include "utils/definitions.h"
#include "utils/basic_functions.h"
#include "utils/Hit.h"

void plot_ghosts(unsigned nEvents=5000, unsigned max_scatter=80000, unsigned max_dt=0) {
  gROOT->SetBatch();  // so stuff isn't autoplotted
  gStyle->SetOptStat(0); // remove the info box
  gStyle->SetAxisMaxDigits(3);
  TString suffix = Utils::Functions::get_suffix(nEvents, max_scatter, max_dt);
  TString input_suffix = suffix + ".root";
  TFile* file = TFile::Open(
    TString((Utils::Definitions::analysisRoot + "hists/ghosts/basic_hists").c_str()) +
    input_suffix);
  if (!file || file->IsZombie()) {
    std::cerr << "Error: Could not open ROOT file." << std::endl;
    return;
  }

  TH1D* h_dt_0_matched = (TH1D*) file->Get("h_dt_0_matched");
  TH1D* h_dt_0_fake = (TH1D*) file->Get("h_dt_0_fake");
  TH1D* h_dt_1_matched = (TH1D*) file->Get("h_dt_1_matched");
  TH1D* h_dt_1_fake = (TH1D*) file->Get("h_dt_1_fake");
  TH1D* h_dt_seeds_matched = (TH1D*) file->Get("h_dt_seeds_matched");
  TH1D* h_dt_seeds_fake = (TH1D*) file->Get("h_dt_seeds_fake");
  TH1D* h_dt_seeds_total_matched = (TH1D*) file->Get("h_dt_seeds_total_matched");
  TH1D* h_dt_seeds_total_fake = (TH1D*) file->Get("h_dt_seeds_total_fake");
  TH1D* h_dt_seeds_total_scaled_matched = (TH1D*) file->Get("h_dt_seeds_total_scaled_matched");
  TH1D* h_dt_seeds_total_scaled_fake = (TH1D*) file->Get("h_dt_seeds_total_scaled_fake");
  TH1D* h_chi_sq_matched = (TH1D*) file->Get("h_chi_sq_matched");
  TH1D* h_chi_sq_fake = (TH1D*) file->Get("h_chi_sq_fake");
  TH1D* h_dr_sq_dz_sq_matched = (TH1D*) file->Get("h_dr_sq_dz_sq_matched");
  TH1D* h_dr_sq_dz_sq_fake = (TH1D*) file->Get("h_dr_sq_dz_dq_fake");
  TH1D* h_dr_sq_dz_matched = (TH1D*) file->Get("h_dr_sq_dz_matched");
  TH1D* h_dr_sq_dz_fake = (TH1D*) file->Get("h_dr_sq_dz_fake");
  TH1D* h_eta_matched = (TH1D*) file->Get("h_eta_matched");
  TH1D* h_eta_fake = (TH1D*) file->Get("h_eta_fake");
  TH1D* h_size_matched = (TH1D*) file->Get("h_size_matched");
  TH1D* h_size_fake = (TH1D*) file->Get("h_size_fake");
  TH1D* h_split_index_matched = (TH1D*) file->Get("h_split_index_matched");
  TH1D* h_split_index_fake = (TH1D*) file->Get("h_split_index_fake");
  TH1D* h_split_index_matched_2 = (TH1D*) file->Get("h_split_index_matched_2");
  TH1D* h_split_index_fake_2 = (TH1D*) file->Get("h_split_index_fake_2");
  TH1D* h_split_index_matched_3 = (TH1D*) file->Get("h_split_index_matched_3");
  TH1D* h_split_index_fake_3 = (TH1D*) file->Get("h_split_index_fake_3");
  TH1D* h_split_index_fake_4 = (TH1D*) file->Get("h_split_index_fake_4");
  TH1D* h_split_index_fake_5 = (TH1D*) file->Get("h_split_index_fake_5");
  TH1D* h_number_mcs_fake = (TH1D*) file->Get("h_number_mcs_fake");
  TH1D* h_number_mcs_matched = (TH1D*) file->Get("h_number_mcs_matched");
  TH1D* h_number_mc_splits_matched = (TH1D*) file->Get("h_number_mc_splits_matched");
  TH1D* h_number_mc_splits_fake = (TH1D*) file->Get("h_number_mc_splits_fake");
  TProfile* p_recod_hits_fake_by_size = (TProfile*) file->Get("p_recod_hits_fake_by_size");
  TH1D* h_dt_splits_matched = (TH1D*) file->Get("h_dt_splits_matched");
  TH1D* h_dt_splits_fake = (TH1D*) file->Get("h_dt_splits_fake");
  TH1D* h_fraction_mc_sizes_2 = (TH1D*) file->Get("h_fraction_mc_sizes_2");
  TH1D* h_fraction_mc_sizes_3 = (TH1D*) file->Get("h_fraction_mc_sizes_3");
  TH1D* h_fraction_mc_sizes_4 = (TH1D*) file->Get("h_fraction_mc_sizes_4");

  printf("Number of entries in fake 2,3,4,5 splits: %f, %f, %f, %f\n",
         h_split_index_fake_2->GetEntries(), h_split_index_fake_3->GetEntries(),
         h_split_index_fake_4->GetEntries(), h_split_index_fake_5->GetEntries());
  printf("Number of entries in matched 2,3 splits: %f, %f\n",
         h_split_index_matched_2->GetEntries(), h_split_index_matched_3->GetEntries()); 

  // on same axes we want (1) dt, (2) dt_seeds, (3) chi_sq, (4) dr_sq_dz_sq, (5) dr_sq_dz & (6) eta
  TCanvas* canvas = new TCanvas("canvas", "Ghosts", 800, 1200);
  canvas->Divide(2, 3);

  canvas->cd(1);  // --------------- (0,0)
  h_dt_0_matched->SetLineColor(kRed);
  h_dt_0_matched->Scale(1. / h_dt_0_matched->GetEntries());
  h_dt_0_fake->SetLineColor(kBlue);
  h_dt_0_fake->Scale(1. / h_dt_0_fake->GetEntries());
  h_dt_1_matched->SetLineColor(kRed);
  h_dt_1_matched->Scale(1. / h_dt_1_matched->GetEntries());
  h_dt_1_fake->SetLineColor(kBlue);
  h_dt_1_fake->Scale(1. / h_dt_1_fake->GetEntries());
  h_dt_1_matched->SetLineStyle(2); // dashed
  h_dt_1_fake->SetLineStyle(2); // dashed
  gPad->SetLogy();
  // draw
  h_dt_0_matched->Draw();
  h_dt_0_fake->Draw("SAME");
  h_dt_1_matched->Draw("SAME");
  h_dt_1_fake->Draw("SAME");
  h_dt_0_matched->SetTitle("Time scatter of subsequent hits;dt (ps);Density");
  // legend
  TLegend* dtLegend = new TLegend(0.4, 0.15, 0.6, 0.35);
  dtLegend->AddEntry(h_dt_0_matched, "Matched (0th order)", "l");
  dtLegend->AddEntry(h_dt_0_fake, "Fake (0th order)", "l");
  dtLegend->AddEntry(h_dt_1_matched, "Matched (1st order)", "l");
  dtLegend->AddEntry(h_dt_1_fake, "Fake (1st order)", "l");
  dtLegend->Draw();

  canvas->cd(2);  // -------------- (0,1)
  h_dt_seeds_matched->SetLineColor(kRed);
  h_dt_seeds_matched->Scale(1. / h_dt_seeds_matched->GetEntries());
  h_dt_seeds_fake->SetLineColor(kBlue);
  h_dt_seeds_fake->Scale(1. / h_dt_seeds_fake->GetEntries());
  h_dt_seeds_matched->Draw();
  h_dt_seeds_fake->Draw("SAME");
  h_dt_seeds_matched->SetTitle("Time scatter between 0th & 2nd hits in seeds;dt (ps);Density");
  gPad->SetLogy();
  // legend
  TLegend* basicFakeMatchedLegend = new TLegend(0.4, 0.15, 0.6, 0.35);
  basicFakeMatchedLegend->AddEntry(h_dt_seeds_matched, "Matched", "l");
  basicFakeMatchedLegend->AddEntry(h_dt_seeds_fake, "Fake", "l");
  basicFakeMatchedLegend->Draw();

  canvas->cd(3);  // -------------- (1,0)
  h_chi_sq_matched->SetLineColor(kRed);
  h_chi_sq_matched->Scale(1. / h_chi_sq_matched->GetEntries());
  h_chi_sq_fake->SetLineColor(kBlue);
  h_chi_sq_fake->Scale(1. / h_chi_sq_fake->GetEntries());
  h_chi_sq_matched->Draw();
  h_chi_sq_fake->Draw("SAME");
  h_chi_sq_matched->SetTitle("#chi^{2} of tracks;#chi^{2};Density");
  basicFakeMatchedLegend->Draw(); // same legend, just fake & matched
  gPad->SetLogy();

  canvas->cd(4);  // -------------- (1,1)
  // Time scatter of seeds in more detail
  h_dt_seeds_total_matched->SetLineColor(kRed);
  h_dt_seeds_total_matched->Scale(1. / h_dt_seeds_total_matched->GetEntries());
  h_dt_seeds_total_fake->SetLineColor(kBlue);
  h_dt_seeds_total_fake->Scale(1. / h_dt_seeds_total_fake->GetEntries());
  h_dt_seeds_total_scaled_matched->SetLineColor(kOrange);
  h_dt_seeds_total_scaled_matched->Scale(1. / h_dt_seeds_total_scaled_matched->GetEntries());
  h_dt_seeds_total_scaled_fake->SetLineColor(kCyan);
  h_dt_seeds_total_scaled_fake->Scale(1. / h_dt_seeds_total_scaled_fake->GetEntries());
  h_dt_seeds_total_matched->SetTitle("Sum of square #Delta t^{2} of 0th to 1st & 2nd hits in seeds;#Delta t^{2} (ps^{2});Density");
  h_dt_seeds_total_matched->Draw();
  h_dt_seeds_total_fake->Draw("SAME");
  h_dt_seeds_total_scaled_matched->Draw("SAME");
  h_dt_seeds_total_scaled_fake->Draw("SAME");
  gPad->SetLogy();
  // legend
  TLegend* dtSeedsLegend = new TLegend(0.5, 0.6, 0.8, 0.85);
  dtSeedsLegend->AddEntry(h_dt_seeds_total_matched, "Matched (unscaled)", "l");
  dtSeedsLegend->AddEntry(h_dt_seeds_total_fake, "Fake (unscaled)", "l");
  dtSeedsLegend->AddEntry(h_dt_seeds_total_scaled_matched, "Matched (scaled)", "l");
  dtSeedsLegend->AddEntry(h_dt_seeds_total_scaled_fake, "Fake (scaled)", "l");
  dtSeedsLegend->Draw();

  canvas->cd(5);  // -------------- (2,0)
  // The following is dr^2 wrt z^2. If you want wrt z, comment and uncomment the block below it
  h_dr_sq_dz_sq_matched->SetLineColor(kRed);
  h_dr_sq_dz_sq_matched->Scale(1. / h_dr_sq_dz_sq_matched->GetEntries());
  h_dr_sq_dz_sq_fake->SetLineColor(kBlue);
  h_dr_sq_dz_sq_fake->Scale(1. / h_dr_sq_dz_sq_fake->GetEntries());
  h_dr_sq_dz_sq_matched->Draw();
  h_dr_sq_dz_sq_fake->Draw("SAME");
  h_dr_sq_dz_sq_matched->SetTitle("Deflection of tracks (wrt z^{2});#Delta r^{2}/z^{2} (unitless);Density");
  basicFakeMatchedLegend->Draw(); // same legend, just fake & matched
  gPad->SetLogy();

  // h_dr_sq_dz_matched->SetLineColor(kRed);
  // h_dr_sq_dz_matched->Scale(1. / h_dr_sq_dz_matched->GetEntries());
  // h_dr_sq_dz_fake->SetLineColor(kBlue);
  // h_dr_sq_dz_fake->Scale(1. / h_dr_sq_dz_fake->GetEntries());
  // h_dr_sq_dz_matched->Draw();
  // h_dr_sq_dz_fake->Draw("SAME");
  // h_dr_sq_dz_matched->SetTitle("Deflection of tracks;#Delta r^{2}/z (mm);Density");
  // basicFakeMatchedLegend->Draw(); // same legend, just fake & matched
  // gPad->SetLogy();

  canvas->cd(6);  // -------------- (2,1)
  TH1D* h_eta_total = (TH1D*) h_eta_fake->Clone("h_eta_total");
  h_eta_total->Add(h_eta_matched);
  h_eta_fake->Divide(h_eta_total);
  gPad->SetLogy();
  h_eta_fake->Draw();
  h_eta_fake->SetTitle("Ghost rate by #eta;#eta;Ghost rate");

  // Save the canvas
  TString outPrefix = (Utils::Definitions::analysisRoot + "output/ghosts/basic_distributions").c_str();
  canvas->SaveAs(outPrefix + input_suffix + ".pdf");

  // Clean up canvas
  delete canvas;

  // Now also plot the size of the tracks and the MC hit usage
  TCanvas* canvas2 = new TCanvas("canvas2", "Ghosts", 800, 1200);
  canvas2->Divide(2, 3);
  canvas2->cd(1);  // -------------- (0,0)
  h_size_fake->Scale(1. / h_size_fake->GetEntries());
  h_size_fake->SetLineColor(kBlue);
  h_size_fake->Draw();
  h_size_matched->SetLineColor(kRed);
  h_size_matched->Scale(1. / h_size_matched->GetEntries());
  h_size_matched->Draw("SAME");
  h_size_fake->SetTitle("Track size distribution (density);N_{Hits};density");
  basicFakeMatchedLegend->Draw(); // same legend, just fake & matched

  canvas2->cd(2);  // -------------- (0,1)
  h_number_mcs_matched->SetLineColor(kRed);
  h_number_mcs_matched->Scale(1. / h_number_mcs_matched->GetEntries());
  h_number_mcs_matched->GetYaxis()->SetRangeUser(0, 1);
  h_number_mcs_matched->Draw();
  h_number_mcs_matched->SetTitle("Density of number of MC particles & splits for tracks;N;Density");
  h_number_mcs_fake->SetLineColor(kBlue);
  h_number_mcs_fake->Scale(1. / h_number_mcs_fake->GetEntries());
  h_number_mcs_fake->Draw("SAME");
  h_number_mc_splits_matched->SetLineColor(kOrange);
  h_number_mc_splits_matched->Scale(1. / h_number_mcs_matched->GetEntries());
  h_number_mc_splits_matched->Draw("SAME");
  h_number_mc_splits_fake->SetLineColor(kCyan);
  h_number_mc_splits_fake->Scale(1. / h_number_mcs_fake->GetEntries());
  h_number_mc_splits_fake->Draw("SAME");
  // legend
  TLegend* mcNLegend = new TLegend(0.45, 0.6, 0.8, 0.8);
  mcNLegend->AddEntry(h_number_mcs_matched, "N_{MCs} (Match)", "l");
  mcNLegend->AddEntry(h_number_mcs_fake, "N_{MCs} (Fake)", "l");
  mcNLegend->AddEntry(h_number_mc_splits_matched, "N_{splits} + 1 (Match)", "l");
  mcNLegend->AddEntry(h_number_mc_splits_fake, "N_{splits} + 1 (Fake)", "l");
  mcNLegend->Draw();


  canvas2->cd(3);  // -------------- (1,0)
  h_fraction_mc_sizes_2->SetLineColor(kRed);
  h_fraction_mc_sizes_2->Scale(1. / h_fraction_mc_sizes_2->GetEntries());
  h_fraction_mc_sizes_3->SetLineColor(kBlue);
  h_fraction_mc_sizes_3->Scale(1. / h_fraction_mc_sizes_3->GetEntries());
  h_fraction_mc_sizes_4->SetLineColor(kGreen);
  h_fraction_mc_sizes_4->Scale(1. / h_fraction_mc_sizes_4->GetEntries());
  h_fraction_mc_sizes_2->Draw();
  h_fraction_mc_sizes_3->Draw("SAME");
  h_fraction_mc_sizes_4->Draw("SAME");
  h_fraction_mc_sizes_2->SetTitle("Distribution of MC fractions of track (normalised);n_{hits}^{MC} / n_{hits}^{track};Fraction Density");
  h_fraction_mc_sizes_2->GetXaxis()->SetRangeUser(0, 0.7);
  // need legend for this
  TLegend* fracLegend = new TLegend(0.2, 0.6, 0.4, 0.8);
  fracLegend->AddEntry(h_fraction_mc_sizes_2, "2 MCs", "l");
  fracLegend->AddEntry(h_fraction_mc_sizes_3, "3 MCs", "l");
  fracLegend->AddEntry(h_fraction_mc_sizes_4, "4 MCs", "l");
  fracLegend->Draw();

  canvas2->cd(4);  // -------------- (1,1)
  // Where the splits occur
  h_split_index_matched->SetLineColor(kRed);
  h_split_index_matched->Scale(1. / h_split_index_matched->GetEntries());
  h_split_index_fake->SetLineColor(kBlue);
  h_split_index_fake->Scale(1. / h_split_index_fake->GetEntries());
  h_split_index_matched->Draw();
  h_split_index_fake->Draw("SAME");
  h_split_index_matched->SetTitle("Distribution of where tracks split;index in track;Density");
  h_split_index_matched->GetYaxis()->SetRangeUser(0, 0.35);
  basicFakeMatchedLegend->Draw();

  canvas2->cd(5);  // -------------- (2,0)
  // time scatter at splits (not including first step: Hit 0->1)
  h_dt_splits_matched->SetLineColor(kRed);
  h_dt_splits_matched->Scale(1. / h_dt_splits_matched->GetEntries());
  h_dt_splits_fake->SetLineColor(kBlue);
  h_dt_splits_fake->Scale(1. / h_dt_splits_fake->GetEntries());
  h_dt_splits_matched->Draw();
  h_dt_splits_fake->Draw("SAME");
  h_dt_splits_matched->SetTitle("Time scatter of split hits (0th order);dt (ps);Density");
  basicFakeMatchedLegend->Draw();

  canvas2->cd(6);  // -------------- (2,1)
  p_recod_hits_fake_by_size->SetLineColor(kBlue);
  p_recod_hits_fake_by_size->SetTitle("Fraction of hits in fake tracks belonging to non-recod MCs;N_{Hits};N_{reco}/N_{Hits}");
  p_recod_hits_fake_by_size->Draw();

  // Save the canvas
  TString outPrefix2 = (Utils::Definitions::analysisRoot + "output/ghosts/size_distributions").c_str();
  canvas2->SaveAs(outPrefix2 + input_suffix + ".pdf");
  // Clean up canvas
  delete canvas2;

  // Now a canvas specifically for split info
  TCanvas* canvas3 = new TCanvas("canvas3", "Ghosts", 800, 800);
  canvas3->Divide(2, 2);
  canvas3->cd(1);  // -------------- (0,0)
  // matched track split indices
  h_split_index_matched_2->SetLineColor(kRed);
  h_split_index_matched_2->Scale(1. / h_split_index_matched_2->GetEntries());
  h_split_index_matched_3->SetLineColor(kBlue);
  printf("Entries: %f\n", h_split_index_matched_3->GetEntries());
  h_split_index_matched_3->Scale(1. / h_split_index_matched_3->GetEntries());
  h_split_index_matched_2->Draw();
  h_split_index_matched_3->Draw("SAME");
  h_split_index_matched_2->SetTitle("Distribution of index of split for matched tracks;index in track;Density");
  h_split_index_matched_2->GetYaxis()->SetRangeUser(0, 0.35);
  // legend
  TLegend* splitMatchedLegend = new TLegend(0.5, 0.6, 0.7, 0.8);
  splitMatchedLegend->AddEntry(h_split_index_matched_2, "1 Split", "l");
  splitMatchedLegend->AddEntry(h_split_index_matched_3, "2 Splits", "l");
  splitMatchedLegend->Draw();

  canvas3->cd(2);  // -------------- (0,1)
  // fake track split indices
  h_split_index_fake_2->SetLineColor(kRed);
  h_split_index_fake_2->Scale(1. / h_split_index_fake_2->GetEntries());
  h_split_index_fake_3->SetLineColor(kBlue);
  h_split_index_fake_3->Scale(1. / h_split_index_fake_3->GetEntries());
  h_split_index_fake_4->SetLineColor(kGreen);
  h_split_index_fake_4->Scale(1. / h_split_index_fake_4->GetEntries());
  h_split_index_fake_5->SetLineColor(kMagenta);
  h_split_index_fake_5->Scale(1. / h_split_index_fake_5->GetEntries());
  h_split_index_fake_2->Draw();
  h_split_index_fake_3->Draw("SAME");
  h_split_index_fake_4->Draw("SAME");
  h_split_index_fake_5->Draw("SAME");
  h_split_index_fake_2->SetTitle("Distribution of index of split for fake tracks;index in track;Density");
  h_split_index_fake_2->GetYaxis()->SetRangeUser(0, 0.35);
  // legend
  TLegend* splitFakeLegend = new TLegend(0.5, 0.6, 0.7, 0.8);
  splitFakeLegend->AddEntry(h_split_index_fake_2, "1 Split", "l");
  splitFakeLegend->AddEntry(h_split_index_fake_3, "2 Splits", "l");
  splitFakeLegend->AddEntry(h_split_index_fake_4, "3 Splits", "l");
  splitFakeLegend->AddEntry(h_split_index_fake_5, "4 Splits", "l");
  splitFakeLegend->Draw();

  canvas3->cd(3);  // -------------- (1,0)
  // empty for now

  canvas3->cd(4);  // -------------- (1,1)
  // empty for now

  // Save the canvas
  TString outPrefix3 = (Utils::Definitions::analysisRoot + "output/ghosts/split_distributions").c_str();
  canvas3->SaveAs(outPrefix3 + input_suffix + ".pdf");
  // Clean up canvas
  delete canvas3;

  delete file;
}