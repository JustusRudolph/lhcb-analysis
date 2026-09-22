#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TGaxis.h>
#include <TProfile.h>
#include <TLegend.h>
#include <TString.h>
#include <TTree.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <string>

#include "utils/definitions.h"
#include "utils/basic_functions.h"

// tMin and tMax bound the timestamp histogram, the out of range counts are printed so it is
// clear how much of the sample falls outside them
void get_plot_hit_density(unsigned nEvents=2000, unsigned t_res_ps=50, bool plot_only=false,
                          float tMin=-5., float tMax=15., int nTimeBins=1000) {
  gStyle->SetOptStat(0);  // remove the info box
  TString output_suffix = Form("_%uev_%ups", nEvents, t_res_ps);
  TString root_output_path = Form((Utils::Definitions::analysisRoot + "hists/hits/hit_density%s.root").c_str(), output_suffix.Data());
  TH1D *h_hit_density = nullptr;  // placeholder, filled either from file or from monitoring tree
  // number of velo hits per MC particle, split by direction and by whether electrons are kept
  TH1D *h_hits_forward_all = nullptr;
  TH1D *h_hits_forward_no_electrons = nullptr;
  TH1D *h_hits_backward_all = nullptr;
  TH1D *h_hits_backward_no_electrons = nullptr;
  // every hit timestamp, and the mean timestamp per module
  TH1D *h_hit_times = nullptr;
  TProfile *p_module_times = nullptr;

  if (plot_only) {
    printf("Plotting from file %s...\n", root_output_path.Data());
    TFile *root_file = TFile::Open(root_output_path);
    if (!root_file || root_file->IsZombie()) {
      std::cerr << "Error: Could not open ROOT file at " << root_output_path << std::endl;
      return;
    }
    h_hit_density = (TH1D*) root_file->Get("hit_density");
    h_hits_forward_all = (TH1D*) root_file->Get("hits_per_track_forward_all");
    h_hits_forward_no_electrons = (TH1D*) root_file->Get("hits_per_track_forward_no_electrons");
    h_hits_backward_all = (TH1D*) root_file->Get("hits_per_track_backward_all");
    h_hits_backward_no_electrons = (TH1D*) root_file->Get("hits_per_track_backward_no_electrons");
    // detach before closing, the file owns whatever is read out of it
    if (h_hit_density) h_hit_density->SetDirectory(0);
    if (h_hits_forward_all) h_hits_forward_all->SetDirectory(0);
    if (h_hits_forward_no_electrons) h_hits_forward_no_electrons->SetDirectory(0);
    if (h_hits_backward_all) h_hits_backward_all->SetDirectory(0);
    if (h_hits_backward_no_electrons) h_hits_backward_no_electrons->SetDirectory(0);
    h_hit_times = (TH1D*) root_file->Get("hit_times");
    p_module_times = (TProfile*) root_file->Get("module_times");
    if (h_hit_times) h_hit_times->SetDirectory(0);
    if (p_module_times) p_module_times->SetDirectory(0);
    root_file->Close();
  } else {
    TString input_filepath = Form((Utils::Definitions::stackRoot + "monitoring/monitoring%ups.root").c_str(), t_res_ps);
    printf("Plotting from monitoring file %s...\n", input_filepath.Data());
    TFile *monitoring_file = TFile::Open(input_filepath);
    if (!monitoring_file || monitoring_file->IsZombie()) {
      std::cerr << "Error: Could not open ROOT file." << std::endl;
      return;
    }
    TTree *hit_tree = (TTree*) monitoring_file->Get("decode_tvclusters/monitor_tree_after");
    if (!hit_tree) {
      std::cerr << "Error: monitor_tree_after TTree not found!" << std::endl;
      return;
    }
    // ----------------- BRANCH DATA ----------------
    printf("Setting branch addresses...\n");
    unsigned hit_lhcbid;
    hit_tree->SetBranchAddress("lhcbID", &hit_lhcbid);
    // make empty histogram to fill
    h_hit_density = new TH1D("hit_density", "Hit Density",
                             Utils::Definitions::kModules, 0,
                             Utils::Definitions::kModules);
    h_hit_density->SetDirectory(0);
    printf("Filling histogram...\n");
    for (unsigned i_hit = 0; i_hit < hit_tree->GetEntries(); i_hit++) {
      hit_tree->GetEntry(i_hit);
      h_hit_density->Fill((hit_lhcbid >> 12) & 0x3F);
    }
    printf("h_hit_density->GetEntries() = %f\n", h_hit_density->GetEntries());

    // hits per track come from the MC truth, the monitoring tree has no particle association
    TString mc_filepath = Form((Utils::Definitions::stackRoot +
                                "output/MCData_Checking_%uev_%ups_1x.root").c_str(),
                               nEvents, t_res_ps);
    printf("Getting hits per track from %s...\n", mc_filepath.Data());
    TFile *mc_file = TFile::Open(mc_filepath);
    if (!mc_file || mc_file->IsZombie()) {
      std::cerr << "Error: Could not open " << mc_filepath << std::endl;
      return;
    }
    TTree *mc_track_tree = (TTree*) mc_file->Get("MCTrackData");
    if (!mc_track_tree) {
      std::cerr << "Error: MCTrackData TTree not found!" << std::endl;
      return;
    }
    unsigned n_velo_hits;
    int mc_pid;
    float mc_eta;
    std::vector<unsigned> *mc_lhcbid = nullptr;
    std::vector<float> *mc_t = nullptr;
    mc_track_tree->SetBranchAddress("nHitsVelo", &n_velo_hits);
    mc_track_tree->SetBranchAddress("pid", &mc_pid);
    mc_track_tree->SetBranchAddress("eta", &mc_eta);
    mc_track_tree->SetBranchAddress("lhcbid", &mc_lhcbid);
    mc_track_tree->SetBranchAddress("t", &mc_t);
    h_hit_times = new TH1D("hit_times", "Hit timestamps;t (ns);Hits", nTimeBins, tMin, tMax);
    p_module_times = new TProfile("module_times",
      "Mean hit timestamp by module;module ID;t (ns)",
      Utils::Definitions::kModules, -0.5, Utils::Definitions::kModules - 0.5);
    h_hit_times->SetDirectory(0);
    p_module_times->SetDirectory(0);
    float tMinSeen = 1e9f, tMaxSeen = -1e9f;
    // one bin per hit count, so the bins sit on the integers
    h_hits_forward_all = new TH1D("hits_per_track_forward_all",
      "Hits per track;N_{hits}^{velo};MC particles", 41, -0.5, 40.5);
    h_hits_forward_no_electrons = new TH1D("hits_per_track_forward_no_electrons",
      "Hits per track;N_{hits}^{velo};MC particles", 41, -0.5, 40.5);
    h_hits_backward_all = new TH1D("hits_per_track_backward_all",
      "Hits per track;N_{hits}^{velo};MC particles", 41, -0.5, 40.5);
    h_hits_backward_no_electrons = new TH1D("hits_per_track_backward_no_electrons",
      "Hits per track;N_{hits}^{velo};MC particles", 41, -0.5, 40.5);
    for (TH1D *hist : {h_hits_forward_all, h_hits_forward_no_electrons,
                       h_hits_backward_all, h_hits_backward_no_electrons}) {
      hist->SetDirectory(0);
      hist->Sumw2();
    }
    for (unsigned i_track = 0; i_track < mc_track_tree->GetEntries(); i_track++) {
      mc_track_tree->GetEntry(i_track);
      // every hit of the track, regardless of direction or particle type
      for (unsigned i_hit = 0; i_hit < mc_lhcbid->size(); i_hit++) {
        float hit_t = mc_t->at(i_hit);
        h_hit_times->Fill(hit_t);
        p_module_times->Fill((mc_lhcbid->at(i_hit) >> 12) & 0x3F, hit_t);
        if (hit_t < tMinSeen) tMinSeen = hit_t;
        if (hit_t > tMaxSeen) tMaxSeen = hit_t;
      }
      int mc_pid_abs = mc_pid < 0 ? -mc_pid : mc_pid;
      bool isForward = mc_eta > 0;
      if (isForward) {
        h_hits_forward_all->Fill(n_velo_hits);
        if (mc_pid_abs != 11) h_hits_forward_no_electrons->Fill(n_velo_hits);  // no e+ or e-
      } else {
        h_hits_backward_all->Fill(n_velo_hits);
        if (mc_pid_abs != 11) h_hits_backward_no_electrons->Fill(n_velo_hits);
      }
    }
    printf("MC particles: %.0f forward, %.0f backward\n",
           h_hits_forward_all->GetEntries(), h_hits_backward_all->GetEntries());
    printf("Hit timestamps: seen %.3f .. %.3f ns, histogram covers %.3f .. %.3f ns\n",
           tMinSeen, tMaxSeen, tMin, tMax);
    printf("  %.0f below range, %.0f above range, out of %.0f hits\n",
           h_hit_times->GetBinContent(0), h_hit_times->GetBinContent(nTimeBins + 1),
           h_hit_times->GetEntries());
    mc_file->Close();

    printf("Saving to file %s...\n", root_output_path.Data());
    TFile *root_file = TFile::Open(root_output_path, "RECREATE");
    h_hit_density->Write();
    h_hits_forward_all->Write();
    h_hits_forward_no_electrons->Write();
    h_hits_backward_all->Write();
    h_hits_backward_no_electrons->Write();
    h_hit_times->Write();
    p_module_times->Write();
    printf("h_hit_density->GetEntries() = %f\n", h_hit_density->GetEntries());
    root_file->Close();
    printf("h_hit_density->GetEntries() = %f\n", h_hit_density->GetEntries());
    monitoring_file->Close();
    printf("h_hit_density->GetEntries() = %f\n", h_hit_density->GetEntries());
    printf("Written files and closed them...\n");
  }
  printf("h_hit_density->GetEntries() = %f\n", h_hit_density->GetEntries());
  printf("Making plots...\n");
  TCanvas *c = new TCanvas("c", "Hit Density", 800, 600);
  if (h_hit_density == nullptr) {
    printf("Error: Histogram not found!\n");
    return;
  }
  printf("h_hit_density->GetEntries() = %f\n", h_hit_density->GetEntries());
  // scale density by number of events
  h_hit_density->Scale(1.0 / nEvents);
  printf("Drawing histogram...\n");
  h_hit_density->Draw();
  c->SaveAs(Form((Utils::Definitions::analysisRoot + "output/hits/hit_density%s.pdf").c_str(), output_suffix.Data()));
  printf("Plots made and saved to %s\n", Form((Utils::Definitions::analysisRoot + "output/hits/hit_density%s.pdf").c_str(), output_suffix.Data()));
  delete c;

  // hits per track, drawn as points with error bars rather than filled bars. Forward and
  // backward differ by orders of magnitude in yield, so backward gets its own axis on the
  // right: it is drawn scaled into the forward frame and a TGaxis shows its true values.
  if (h_hits_forward_all == nullptr || h_hits_forward_no_electrons == nullptr ||
      h_hits_backward_all == nullptr || h_hits_backward_no_electrons == nullptr) {
    printf("Error: Hits per track histograms not found!\n");
    return;
  }
  TCanvas *c_hits = new TCanvas("c_hits", "Hits per Track", 800, 600);
  c_hits->SetRightMargin(0.14);  // room for the backward axis
  // per event, as for the hit density above. Sumw2 was set at fill time so the errors scale too
  for (TH1D *hist : {h_hits_forward_all, h_hits_forward_no_electrons,
                     h_hits_backward_all, h_hits_backward_no_electrons}) {
    hist->Scale(1.0 / nEvents);
  }
  // blue is forward and red backward, filled markers keep electrons and open ones drop them
  h_hits_forward_all->SetLineColor(kBlue);
  h_hits_forward_all->SetMarkerColor(kBlue);
  h_hits_forward_all->SetMarkerStyle(20);
  h_hits_forward_no_electrons->SetLineColor(kBlue);
  h_hits_forward_no_electrons->SetMarkerColor(kBlue);
  h_hits_forward_no_electrons->SetMarkerStyle(24);
  h_hits_backward_all->SetLineColor(kRed);
  h_hits_backward_all->SetMarkerColor(kRed);
  h_hits_backward_all->SetMarkerStyle(21);
  h_hits_backward_no_electrons->SetLineColor(kRed);
  h_hits_backward_no_electrons->SetMarkerColor(kRed);
  h_hits_backward_no_electrons->SetMarkerStyle(25);
  for (TH1D *hist : {h_hits_forward_all, h_hits_forward_no_electrons,
                     h_hits_backward_all, h_hits_backward_no_electrons}) {
    hist->SetMarkerSize(0.7);
  }

  double maxForward = std::max(h_hits_forward_all->GetMaximum(),
                               h_hits_forward_no_electrons->GetMaximum()) * 1.15;
  double maxBackward = std::max(h_hits_backward_all->GetMaximum(),
                                h_hits_backward_no_electrons->GetMaximum()) * 1.15;
  h_hits_forward_all->SetTitle("Hits per track;N_{hits}^{velo};Forward MC particles / event");
  h_hits_forward_all->SetMinimum(0.);
  h_hits_forward_all->SetMaximum(maxForward);
  h_hits_forward_all->Draw("E1");
  h_hits_forward_no_electrons->Draw("E1 SAME");
  // clone the backward ones so the originals keep their true values
  double backwardScale = (maxBackward > 0.) ? maxForward / maxBackward : 1.;
  TH1D *h_backward_all_scaled = (TH1D*) h_hits_backward_all->Clone("backward_all_scaled");
  TH1D *h_backward_no_electrons_scaled =
    (TH1D*) h_hits_backward_no_electrons->Clone("backward_no_electrons_scaled");
  h_backward_all_scaled->SetDirectory(0);
  h_backward_no_electrons_scaled->SetDirectory(0);
  h_backward_all_scaled->Scale(backwardScale);
  h_backward_no_electrons_scaled->Scale(backwardScale);
  h_backward_all_scaled->Draw("E1 SAME");
  h_backward_no_electrons_scaled->Draw("E1 SAME");
  gPad->Update();
  // right hand axis carrying the unscaled backward numbers
  TGaxis *backwardAxis = new TGaxis(gPad->GetUxmax(), 0., gPad->GetUxmax(), maxForward,
                                    0., maxBackward, 510, "+L");
  backwardAxis->SetTitle("Backward MC particles / event");
  backwardAxis->SetLineColor(kRed);
  backwardAxis->SetLabelColor(kRed);
  backwardAxis->SetTitleColor(kRed);
  backwardAxis->SetLabelFont(h_hits_forward_all->GetYaxis()->GetLabelFont());
  backwardAxis->SetLabelSize(h_hits_forward_all->GetYaxis()->GetLabelSize());
  backwardAxis->SetTitleFont(h_hits_forward_all->GetYaxis()->GetTitleFont());
  backwardAxis->SetTitleSize(h_hits_forward_all->GetYaxis()->GetTitleSize());
  backwardAxis->Draw();
  TLegend *legend_hits = new TLegend(0.55, 0.68, 0.85, 0.87);
  legend_hits->SetBorderSize(0);
  legend_hits->AddEntry(h_hits_forward_all, "Forward, all particles", "lp");
  legend_hits->AddEntry(h_hits_forward_no_electrons, "Forward, no electrons", "lp");
  legend_hits->AddEntry(h_backward_all_scaled, "Backward, all particles", "lp");
  legend_hits->AddEntry(h_backward_no_electrons_scaled, "Backward, no electrons", "lp");
  legend_hits->Draw();
  TString hits_output_path = Form((Utils::Definitions::analysisRoot +
                                   "output/hits/hits_per_track%s.pdf").c_str(),
                                  output_suffix.Data());
  c_hits->SaveAs(hits_output_path);
  printf("Hits per track plot saved to %s\n", hits_output_path.Data());
  delete c_hits;

  // hit timestamps, and their mean per module
  if (h_hit_times == nullptr || p_module_times == nullptr) {
    printf("Error: Timestamp histograms not found!\n");
    return;
  }
  TCanvas *c_times = new TCanvas("c_times", "Hit Timestamps", 800, 600);
  h_hit_times->SetLineColor(kBlue);
  h_hit_times->Draw("HIST");
  gPad->SetLogy();  // the distribution has a long tail
  TString times_output_path = Form((Utils::Definitions::analysisRoot +
                                    "output/hits/hit_times%s.pdf").c_str(),
                                   output_suffix.Data());
  c_times->SaveAs(times_output_path);
  delete c_times;

  TCanvas *c_module_times = new TCanvas("c_module_times", "Mean Time per Module", 800, 600);
  p_module_times->SetLineColor(kBlue);
  p_module_times->SetMarkerColor(kBlue);
  p_module_times->SetMarkerStyle(20);
  p_module_times->SetMarkerSize(0.7);
  p_module_times->Draw("E1");
  TString module_times_output_path = Form((Utils::Definitions::analysisRoot +
                                           "output/hits/module_times%s.pdf").c_str(),
                                          output_suffix.Data());
  c_module_times->SaveAs(module_times_output_path);
  printf("Timestamp plots saved to %s and %s\n",
         times_output_path.Data(), module_times_output_path.Data());
  delete c_module_times;
}
