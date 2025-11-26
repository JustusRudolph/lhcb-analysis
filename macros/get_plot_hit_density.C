#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <string>

#include "utils/definitions.h"
#include "utils/basic_functions.h"

void get_plot_hit_density(unsigned nEvents=2000, unsigned t_res_ps=50, bool plot_only=false) {
  TString output_suffix = Form("_%uev_%ups", nEvents, t_res_ps);
  TString root_output_path = Form((Utils::Definitions::analysisRoot + "hists/hits/hit_density%s.root").c_str(), output_suffix.Data());
  TH1D *h_hit_density = nullptr;  // placeholder, filled either from file or from monitoring tree

  if (plot_only) {
    printf("Plotting from file %s...\n", root_output_path.Data());
    TFile *root_file = TFile::Open(root_output_path);
    if (!root_file || root_file->IsZombie()) {
      std::cerr << "Error: Could not open ROOT file at " << root_output_path << std::endl;
      return;
    }
    h_hit_density = (TH1D*) root_file->Get("hit_density");
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
    printf("Saving to file %s...\n", root_output_path.Data());
    TFile *root_file = TFile::Open(root_output_path, "RECREATE");
    h_hit_density->Write();
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
}