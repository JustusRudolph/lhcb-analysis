#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TProfile.h>
#include <TString.h>
#include <TSystem.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "utils/definitions.h"
#include "utils/basic_functions.h"

/*
 * Times are quoted in the legends in ps, rounded to two significant figures,
 * e.g. 0.1234 ns becomes 120 ps.
 */
int to_ps_two_sig_figs(double t_ns) {
  double t_ps = t_ns * 1000.;
  if (t_ps == 0.) return 0;
  double magnitude = std::pow(10., std::floor(std::log10(std::abs(t_ps))) - 1.);
  return (int) std::lround(std::lround(t_ps / magnitude) * magnitude);
}

// every histogram the canvas needs, in the order they are fetched and checked
const std::vector<std::string> kOutlierTimeHists = {
  "h_dt0_forward", "h_dt2_forward", "h_dt0_backward", "h_dt2_backward",
  "h_dt_forwarding_forward", "h_dt_forwarding_forward_h5",
  "h_dt_forwarding_forward_h10", "h_dt_forwarding_forward_h15",
  "h_dt_forwarding_backward", "h_dt_forwarding_backward_h5",
  "h_dt_forwarding_backward_h10", "h_dt_forwarding_backward_h15",
  "h_dt_forwarding_vs_moduleID", "h_nthHits_vs_moduleID", "h_dt_forwarding_vs_nthHit",
  "p_dt_forwarding_vs_moduleID", "p_inv_beta_vs_moduleID", "h_inv_beta_vs_moduleID",
  "p_dt_vs_nth_hit_forward", "p_dt_vs_nth_hit_backward"};

/*
 * Draw the canvas that check_outlier_times.C used to draw directly, from the histograms
 * written by get_outlier_times.C. The dataset is picked the same way as there, histName
 * overrides that if the hists were written under a different name.
 */
void plot_outlier_times(unsigned nEvents=5000, unsigned max_scatter=80000, unsigned max_dt=0,
                        TString mc_file_suffix="", TString histName="") {
  gROOT->SetBatch();  // so stuff isn't autoplotted
  gStyle->SetOptStat(0); // remove the info box for the plots
  TString suffix;
  if (histName.IsNull()) {
    if (mc_file_suffix.IsNull()) {
      suffix = Utils::Functions::get_suffix(nEvents, max_scatter, max_dt);
    } else {
      suffix = Form("_%uev_%s", nEvents, mc_file_suffix.Data());
    }
    histName = TString("outlier_times_hists") + suffix;
  }
  TString histPath =
    TString((Utils::Definitions::analysisRoot + "hists/4d_tracking/").c_str()) + histName + ".root";
  TFile* file = TFile::Open(histPath);
  if (!file || file->IsZombie()) {
    std::cerr << "Error: Could not open " << histPath << std::endl;
    return;
  }
  TH1D* h_dt0_forward = (TH1D*) file->Get("h_dt0_forward");
  TH1D* h_dt2_forward = (TH1D*) file->Get("h_dt2_forward");
  TH1D* h_dt0_backward = (TH1D*) file->Get("h_dt0_backward");
  TH1D* h_dt2_backward = (TH1D*) file->Get("h_dt2_backward");
  TH1D* h_dt_forwarding_forward = (TH1D*) file->Get("h_dt_forwarding_forward");
  TH1D* h_dt_forwarding_forward_h5 = (TH1D*) file->Get("h_dt_forwarding_forward_h5");
  TH1D* h_dt_forwarding_forward_h10 = (TH1D*) file->Get("h_dt_forwarding_forward_h10");
  TH1D* h_dt_forwarding_forward_h15 = (TH1D*) file->Get("h_dt_forwarding_forward_h15");
  TH1D* h_dt_forwarding_backward = (TH1D*) file->Get("h_dt_forwarding_backward");
  TH1D* h_dt_forwarding_backward_h5 = (TH1D*) file->Get("h_dt_forwarding_backward_h5");
  TH1D* h_dt_forwarding_backward_h10 = (TH1D*) file->Get("h_dt_forwarding_backward_h10");
  TH1D* h_dt_forwarding_backward_h15 = (TH1D*) file->Get("h_dt_forwarding_backward_h15");
  TH2D* h_dt_forwarding_vs_moduleID = (TH2D*) file->Get("h_dt_forwarding_vs_moduleID");
  TH2D* h_nthHits_vs_moduleID = (TH2D*) file->Get("h_nthHits_vs_moduleID");
  TH2D* h_dt_forwarding_vs_nthHit = (TH2D*) file->Get("h_dt_forwarding_vs_nthHit");
  TProfile* p_dt_forwarding_vs_moduleID = (TProfile*) file->Get("p_dt_forwarding_vs_moduleID");
  TProfile* p_inv_beta_vs_moduleID = (TProfile*) file->Get("p_inv_beta_vs_moduleID");
  TH2D* h_inv_beta_vs_moduleID = (TH2D*) file->Get("h_inv_beta_vs_moduleID");
  TProfile* p_dt_vs_nth_hit_forward = (TProfile*) file->Get("p_dt_vs_nth_hit_forward");
  TProfile* p_dt_vs_nth_hit_backward = (TProfile*) file->Get("p_dt_vs_nth_hit_backward");
  std::vector<TH1*> allHists = {
    h_dt0_forward, h_dt2_forward, h_dt0_backward, h_dt2_backward,
    h_dt_forwarding_forward, h_dt_forwarding_forward_h5,
    h_dt_forwarding_forward_h10, h_dt_forwarding_forward_h15,
    h_dt_forwarding_backward, h_dt_forwarding_backward_h5,
    h_dt_forwarding_backward_h10, h_dt_forwarding_backward_h15,
    h_dt_forwarding_vs_moduleID, h_nthHits_vs_moduleID, h_dt_forwarding_vs_nthHit,
    p_dt_forwarding_vs_moduleID, p_inv_beta_vs_moduleID, h_inv_beta_vs_moduleID,
    p_dt_vs_nth_hit_forward, p_dt_vs_nth_hit_backward};
  for (unsigned i = 0; i < allHists.size(); i++) {
    if (!allHists[i]) {
      std::cerr << "Error: " << kOutlierTimeHists[i] << " not found in " << histPath << std::endl;
      return;
    }
  }

  TCanvas* canvas = new TCanvas("canvas", "Outlier Times", 2500, 800);
  canvas->Divide(5, 2);
  canvas->cd(1);
  // draw h0 and h2
  h_dt0_forward->SetTitle("Seeding time scatter from t_{0}^{est} and t_{2}^{est} (normalised)");
  h_dt0_forward->SetLineColor(kBlue);
  h_dt2_forward->SetLineColor(kRed);
  // make colours on backward plots same but more transparent
  h_dt0_backward->SetLineColorAlpha(kBlue, 0.25);
  h_dt2_backward->SetLineColorAlpha(kRed, 0.25);
  h_dt0_forward->Scale(1.0 / h_dt0_forward->Integral());
  h_dt2_forward->Scale(1.0 / h_dt2_forward->Integral());
  h_dt0_backward->Scale(1.0 / h_dt0_backward->Integral());
  h_dt2_backward->Scale(1.0 / h_dt2_backward->Integral());
  h_dt0_forward->Draw("HIST");
  h_dt0_forward->GetXaxis()->SetTitle("#Deltat (ns)");
  h_dt2_forward->Draw("HIST SAME");
  h_dt0_backward->Draw("HIST SAME");
  h_dt2_backward->Draw("HIST SAME");
  h_dt0_forward->GetYaxis()->SetRangeUser(0, 0.15);
  h_dt0_forward->GetXaxis()->SetRangeUser(-0.3, 0.3);
  // draw legend
  TLegend* legend_dt_seeding = new TLegend(0.6, 0.6, 0.88, 0.8);
  legend_dt_seeding->SetBorderSize(0);
  legend_dt_seeding->AddEntry(h_dt0_forward, Form("#Deltat_{0}^{f} (#mu, #sigma) = (%d, %d) ps",
    to_ps_two_sig_figs(h_dt0_forward->GetMean()), to_ps_two_sig_figs(h_dt0_forward->GetStdDev()) ), "l");
  legend_dt_seeding->AddEntry(h_dt2_forward, Form("#Deltat_{2}^{f} (#mu, #sigma) = (%d, %d) ps",
    to_ps_two_sig_figs(h_dt2_forward->GetMean()), to_ps_two_sig_figs(h_dt2_forward->GetStdDev()) ), "l");
  legend_dt_seeding->AddEntry(h_dt0_backward, Form("#Deltat_{0}^{b} (#mu, #sigma) = (%d, %d) ps",
    to_ps_two_sig_figs(h_dt0_backward->GetMean()), to_ps_two_sig_figs(h_dt0_backward->GetStdDev()) ), "l");
  legend_dt_seeding->AddEntry(h_dt2_backward, Form("#Deltat_{2}^{b} (#mu, #sigma) = (%d, %d) ps",
    to_ps_two_sig_figs(h_dt2_backward->GetMean()), to_ps_two_sig_figs(h_dt2_backward->GetStdDev()) ), "l");
  legend_dt_seeding->Draw();
  
  canvas->cd(2);
  // plot the 3 histograms for hit indices 5, 10, and 15 on the same canvas
  // use copy of h_dt_forwarding_forward to not modify previous plot (if we use it later again)
  TH1D* h_dt_forwarding_forward_copy = (TH1D*) h_dt_forwarding_forward->Clone();
  h_dt_forwarding_forward_copy->SetLineColor(kBlack);
  h_dt_forwarding_forward_h5->SetLineColor(kRed);
  h_dt_forwarding_forward_h10->SetLineColor(kBlue);
  h_dt_forwarding_forward_h15->SetLineColor(kGreen);
  // normalise the histograms
  h_dt_forwarding_forward_copy->Scale(1.0 / h_dt_forwarding_forward_copy->Integral());
  h_dt_forwarding_forward_h5->Scale(1.0 / h_dt_forwarding_forward_h5->Integral());
  h_dt_forwarding_forward_h10->Scale(1.0 / h_dt_forwarding_forward_h10->Integral());
  h_dt_forwarding_forward_h15->Scale(1.0 / h_dt_forwarding_forward_h15->Integral());
  // draw as histograms not as lines
  h_dt_forwarding_forward_copy->Draw("HIST");
  h_dt_forwarding_forward_h5->Draw("HIST SAME");
  h_dt_forwarding_forward_h10->Draw("HIST SAME");
  h_dt_forwarding_forward_h15->Draw("HIST SAME");
  h_dt_forwarding_forward_copy->SetTitle("#Deltat scatter for forwarded hits (norm, forward)");
  h_dt_forwarding_forward_copy->GetXaxis()->SetTitle("#Deltat_{forward} (ns)");
  h_dt_forwarding_forward_copy->GetYaxis()->SetRangeUser(0, 0.15);
  h_dt_forwarding_forward_copy->GetXaxis()->SetRangeUser(-0.5, 0.5);
  // legend
  TLegend* legend = new TLegend(0.6, 0.6, 0.88, 0.8);
  legend->SetBorderSize(0);
  legend->AddEntry(h_dt_forwarding_forward_copy, Form("All forwarded (#mu, #sigma) = (%d, %d) ps",
    to_ps_two_sig_figs(h_dt_forwarding_forward->GetMean()), to_ps_two_sig_figs(h_dt_forwarding_forward->GetStdDev())), "l");
  legend->AddEntry(h_dt_forwarding_forward_h5, Form("Hit index 5 (#mu, #sigma) = (%d, %d) ps",
    to_ps_two_sig_figs(h_dt_forwarding_forward_h5->GetMean()), to_ps_two_sig_figs(h_dt_forwarding_forward_h5->GetStdDev())), "l");
  legend->AddEntry(h_dt_forwarding_forward_h10, Form("Hit index 10 (#mu, #sigma) = (%d, %d) ps",
    to_ps_two_sig_figs(h_dt_forwarding_forward_h10->GetMean()), to_ps_two_sig_figs(h_dt_forwarding_forward_h10->GetStdDev())), "l");
  legend->AddEntry(h_dt_forwarding_forward_h15, Form("Hit index 15 (#mu, #sigma) = (%d, %d) ps",
    to_ps_two_sig_figs(h_dt_forwarding_forward_h15->GetMean()), to_ps_two_sig_figs(h_dt_forwarding_forward_h15->GetStdDev())), "l");
  legend->Draw();

  canvas->cd(3);
  // plot 2D histogram of dt vs module ID for forwarded hits
  // log scale on z axis
  h_dt_forwarding_vs_moduleID->Draw("COLZ");
  h_dt_forwarding_vs_moduleID->GetYaxis()->SetTitle("#Deltat (ns)");
  h_dt_forwarding_vs_moduleID->SetTitle("#Deltat scatter vs module ID for all forwarded hits");
  gPad->SetLogz();

  canvas->cd(4);
  TH1D* h_dt_forwarding_backward_copy = (TH1D*) h_dt_forwarding_backward->Clone();
  h_dt_forwarding_backward_copy->SetLineColor(kBlack);
  h_dt_forwarding_backward_h5->SetLineColor(kRed);
  h_dt_forwarding_backward_h10->SetLineColor(kBlue);
  h_dt_forwarding_backward_h15->SetLineColor(kGreen);
  // normalise the histograms
  h_dt_forwarding_backward_copy->Scale(1.0 / h_dt_forwarding_backward_copy->Integral());
  h_dt_forwarding_backward_h5->Scale(1.0 / h_dt_forwarding_backward_h5->Integral());
  h_dt_forwarding_backward_h10->Scale(1.0 / h_dt_forwarding_backward_h10->Integral());
  h_dt_forwarding_backward_h15->Scale(1.0 / h_dt_forwarding_backward_h15->Integral());
  // draw as histograms not as lines
  h_dt_forwarding_backward_copy->Draw("HIST");
  h_dt_forwarding_backward_h5->Draw("HIST SAME");
  h_dt_forwarding_backward_h10->Draw("HIST SAME");
  // h_dt_forwarding_backward_h15->Draw("HIST SAME");
  h_dt_forwarding_backward_copy->SetTitle("#Deltat scatter for forwarded hits (norm, backward)");
  h_dt_forwarding_backward_copy->GetXaxis()->SetTitle("#Deltat_{forward} (ns)");
  h_dt_forwarding_backward_copy->GetYaxis()->SetRangeUser(0, 0.15);
  h_dt_forwarding_backward_copy->GetXaxis()->SetRangeUser(-0.5, 0.5);
  // legend

  TLegend* legend_backward = new TLegend(0.55, 0.7, 0.89, 0.8);
  legend_backward->SetBorderSize(0);
  legend_backward->AddEntry(h_dt_forwarding_backward_copy, Form("All forwarded (#mu, #sigma) = (%d, %d) ps",
    to_ps_two_sig_figs(h_dt_forwarding_backward->GetMean()), to_ps_two_sig_figs(h_dt_forwarding_backward->GetStdDev())), "l");
  legend_backward->AddEntry(h_dt_forwarding_backward_h5, Form("Hit index 5 (#mu, #sigma) = (%d, %d) ps",
    to_ps_two_sig_figs(h_dt_forwarding_backward_h5->GetMean()), to_ps_two_sig_figs(h_dt_forwarding_backward_h5->GetStdDev())), "l");
  legend_backward->AddEntry(h_dt_forwarding_backward_h10, Form("Hit index 10 (#mu, #sigma) = (%d, %d) ps",
    to_ps_two_sig_figs(h_dt_forwarding_backward_h10->GetMean()), to_ps_two_sig_figs(h_dt_forwarding_backward_h10->GetStdDev())), "l");
  // legend_backward->AddEntry(h_dt_forwarding_backward_h15, Form("Hit index 15 (#mu, #sigma) = (%d, %d) ps",
  //   to_ps_two_sig_figs(h_dt_forwarding_backward_h15->GetMean()), to_ps_two_sig_figs(h_dt_forwarding_backward_h15->GetStdDev())), "l");
  legend_backward->Draw();

  canvas->cd(5);
  // plot 2D histogram of n hits vs module ID for forwarded hits
  h_nthHits_vs_moduleID->Draw("COLZ");
  h_nthHits_vs_moduleID->GetYaxis()->SetTitle("Hit number in track");
  h_nthHits_vs_moduleID->SetTitle("Hit number vs module ID for forwarded hits");
  gPad->SetLogz();

  // get plot of number of outliers per track & set log scale y axis
  // gPad->SetLogy();
  // h_n_outliers->SetTitle("Number of outliers per track");
  // h_n_outliers->Draw();
  
  canvas->cd(6);
  h_dt_forwarding_vs_nthHit->Draw("COLZ");
  gPad->SetLogz();
  
  // gPad->SetLogy();
  // h_outlier_position->SetTitle("Position of outlier in track");
  // h_outlier_position->Draw();

  canvas->cd(7);
  // mean with the standard deviation as the error bars, per module the hit is extrapolated to
  p_dt_forwarding_vs_moduleID->SetLineColor(kBlue);
  p_dt_forwarding_vs_moduleID->SetMarkerColor(kBlue);
  p_dt_forwarding_vs_moduleID->SetMarkerStyle(20);
  p_dt_forwarding_vs_moduleID->SetMarkerSize(0.7);
  p_dt_forwarding_vs_moduleID->Draw("E1");

  canvas->cd(8);
  // The 2D holds the full c/v distribution per module. Its median is taken rather than its
  // mean: the distribution is a sharp peak near 1 with a tail of genuinely slow particles, and
  // a mean over that tail describes no track in particular. Inverting a median is exact,
  // median(v/c) = 1 / median(c/v), so the inversion introduces no bias of its own.
  TH1D* h_beta_vs_moduleID = new TH1D("h_beta_vs_moduleID",
    "Median speed of forwarded hits vs module ID;Module ID;#beta",
    h_inv_beta_vs_moduleID->GetNbinsX(),
    h_inv_beta_vs_moduleID->GetXaxis()->GetXmin(),
    h_inv_beta_vs_moduleID->GetXaxis()->GetXmax());
  h_beta_vs_moduleID->SetDirectory(nullptr);
  double lowestBeta = 1., highestBeta = 1.;  // start at 1 so the reference line stays in range
  for (int bin = 1; bin <= h_inv_beta_vs_moduleID->GetNbinsX(); bin++) {
    TH1D* invBetaInModule =
      h_inv_beta_vs_moduleID->ProjectionY(Form("inv_beta_module_%d", bin), bin, bin);
    invBetaInModule->SetDirectory(nullptr);
    double nEntries = invBetaInModule->Integral();
    if (nEntries < 2.) { delete invBetaInModule; continue; }
    double quantile = 0.5, median = 0.;
    invBetaInModule->GetQuantiles(1, &median, &quantile);
    if (median <= 0.) { delete invBetaInModule; continue; }
    // the error on a median is about 1.25 times the one on a mean for a gaussian core
    double medianError = 1.2533 * invBetaInModule->GetStdDev() / std::sqrt(nEntries);
    double beta = 1. / median;
    double error = medianError / (median * median);
    h_beta_vs_moduleID->SetBinContent(bin, beta);
    h_beta_vs_moduleID->SetBinError(bin, error);
    lowestBeta = std::min(lowestBeta, beta - error);
    highestBeta = std::max(highestBeta, beta + error);
    delete invBetaInModule;
  }
  h_beta_vs_moduleID->SetLineColor(kBlue);
  h_beta_vs_moduleID->SetMarkerColor(kBlue);
  h_beta_vs_moduleID->SetMarkerStyle(20);
  h_beta_vs_moduleID->SetMarkerSize(0.7);
  // the deviation from 1 is small, so zoom onto the points while keeping 1 in view
  double betaMargin = 0.05 * (highestBeta - lowestBeta);
  if (betaMargin == 0.) betaMargin = 0.001;  // everything at exactly 1
  h_beta_vs_moduleID->GetYaxis()->SetRangeUser(lowestBeta - betaMargin, highestBeta + betaMargin);
  h_beta_vs_moduleID->Draw("E1");
  gPad->Update();
  TLine* speedOfLight = new TLine(-0.5, 1., 63.5, 1.);
  speedOfLight->SetLineStyle(2);
  speedOfLight->SetLineColor(kGray + 2);
  speedOfLight->Draw();

  canvas->cd(9);
  // mean dt against position in the track, the filter's convergence seen directly
  p_dt_vs_nth_hit_forward->SetLineColor(kBlue);
  p_dt_vs_nth_hit_forward->SetMarkerColor(kBlue);
  p_dt_vs_nth_hit_forward->SetMarkerStyle(20);
  p_dt_vs_nth_hit_forward->SetMarkerSize(0.7);
  p_dt_vs_nth_hit_forward->Draw("E1");

  canvas->cd(10);
  p_dt_vs_nth_hit_backward->SetLineColor(kRed);
  p_dt_vs_nth_hit_backward->SetMarkerColor(kRed);
  p_dt_vs_nth_hit_backward->SetMarkerStyle(21);
  p_dt_vs_nth_hit_backward->SetMarkerSize(0.7);
  p_dt_vs_nth_hit_backward->Draw("E1");

  // save canvas
  TString outputDir = (Utils::Definitions::analysisRoot + "output/4d_tracking").c_str();
  gSystem->mkdir(outputDir, true);
  TString outputPath = outputDir + "/outlier_times" + suffix + ".pdf";
  canvas->SaveAs(outputPath);
  std::cout << "Canvas saved to: " << outputPath << std::endl;


  // clean up
  delete canvas;
  file->Close();
  delete file;
}
