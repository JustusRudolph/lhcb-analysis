#include <TFile.h>
#include <TMath.h>
#include <TProfile.h>
#include <TSystem.h>
#include <TTree.h>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <string>
#include <tuple>

#include "utils/definitions.h"
#include "utils/allen_functions.h"
#include "utils/basic_functions.h"

/*
 * Fill the seeding/forwarding #Deltat histograms of one dataset and write them to
 * hists/4d_tracking, so that plot_outlier_times.C draws them and compare_outlier_times.C
 * can take ratios between two datasets. outName defaults to the usual data suffix.
 */
// max_dt is in picoseconds and scatter in micrometers
void get_outlier_times(unsigned nEvents=5000, unsigned max_scatter=80000, unsigned max_dt=0,
                       TString mc_file_suffix="", float t_res=0.05, float t_stddevs=3.,
                       float tWack=700.0, TString outName="") {
  TString suffix;
  if (mc_file_suffix.IsNull()) {
    suffix = Utils::Functions::get_suffix(nEvents, max_scatter, max_dt);
  } else {
    suffix = Form("_%uev_%s", nEvents, mc_file_suffix.Data());
  }
  TString input_prefix = (Utils::Definitions::stackRoot + "output/MCData_Checking").c_str();
  TString filepath = input_prefix + suffix + ".root";

  TFile *file = TFile::Open(filepath);
  if (!file || file->IsZombie()) {
      std::cerr << "Error: Could not open " << filepath << std::endl;
      return;
  }
  if (outName.IsNull()) outName = TString("outlier_times_hists") + suffix;
  // get mc track tree
  TTree* mcTrackTree = (TTree*) file->Get("MCTrackData");
  if (!mcTrackTree) {
      std::cerr << "TTree not found!" << std::endl;
      return;
  }

  // define output histograms for dts (0, 2, and forwarding)
  TH1D* h_dt0_forward = new TH1D("h_dt0_forward", "h_dt0_forward;dt_{0} (ns)", 100, -1, 1);
  TH1D* h_dt2_forward = new TH1D("h_dt2_forward", "h_dt2_forward;dt_{2} (ns)", 100, -1, 1);
  TH1D* h_dt_forwarding_forward = new TH1D("h_dt_forwarding_forward", "dt scatter for forwarded hits (outlier times not used in time filtering);dt_{forward} (ns)", 100, -1, 1);
  TH1D* h_dt_forwarding_forward_h5 = new TH1D("h_dt_forwarding_forward_h5", "dt scatter for hit index 5;dt_{forward} (ns)", 100, -1, 1);
  TH1D* h_dt_forwarding_forward_h10 = new TH1D("h_dt_forwarding_forward_h10", "dt scatter for hit index 10;dt_{forward} (ns)", 100, -1, 1);
  TH1D* h_dt_forwarding_forward_h15 = new TH1D("h_dt_forwarding_forward_h15", "dt scatter for hit index 15;dt_{forward} (ns)", 100, -1, 1);
  
  TH1D* h_dt0_backward = new TH1D("h_dt0_backward", "h_dt0_backward;dt_{0} (ns)", 100, -1, 1);
  TH1D* h_dt2_backward = new TH1D("h_dt2_backward", "h_dt2_backward;dt_{2} (ns)", 100, -1, 1);
  TH1D* h_dt_forwarding_backward = new TH1D("h_dt_forwarding_backward", "dt scatter for backwards tracks;dt_{forward} (ns)", 100, -1, 1);
  TH1D* h_dt_forwarding_backward_h5 = new TH1D("h_dt_forwarding_backward_h5", "dt scatter for hit index 5;dt_{forward} (ns)", 100, -1, 1);
  TH1D* h_dt_forwarding_backward_h10 = new TH1D("h_dt_forwarding_backward_h10", "dt scatter for hit index 10;dt_{forward} (ns)", 100, -1, 1);
  TH1D* h_dt_forwarding_backward_h15 = new TH1D("h_dt_forwarding_backward_h15", "dt scatter for hit index 15;dt_{forward} (ns)", 100, -1, 1);

  TH2D* h_dt_forwarding_vs_moduleID = new TH2D("h_dt_forwarding_vs_moduleID", "#Deltat scatter vs module ID for forwarded hits;Module ID;dt_{forward} (ns)", 64, -0.5, 63.5, 100, -1, 1);
  // "s" makes the error bars the standard deviation rather than the error on the mean, so the
  // profile carries both. Filled directly rather than profiled off the 2D above, which would
  // silently drop everything beyond its +-1 ns range.
  TProfile* p_dt_forwarding_vs_moduleID = new TProfile("p_dt_forwarding_vs_moduleID", "#Deltat mean and spread vs module ID of the extrapolated hit;Module ID;dt_{forward} (ns)", 64, -0.5, 63.5, "s");
  TH2D* h_nthHits_vs_moduleID = new TH2D("h_nthHits_vs_moduleID", "Number of hit at module ID for forwarded hits;Module ID;n_{hits}", 64, -0.5, 63.5, 40, -0.5, 39.5);
  TH2D* h_dt_forwarding_vs_nthHit = new TH2D("h_dt_forwarding_vs_nthHit", "#Deltat scatter vs hit number for forwarded hits;Nth hit;dt_{forward} (ns)", 40, -0.5, 39.5, 100, -1, 1);

  TH1D* h_n_outliers = new TH1D("h_n_outliers", "Number of outliers per track;n_{outliers}", 10, 0, 10);
  unsigned max_track_length_plotting = 20;
  TH1D* h_outlier_position = new TH1D("h_outlier_position", "Position of outlier in track;Track index",
                                      max_track_length_plotting, 0, max_track_length_plotting);

  // for now, just check pid of the concerned cluster
  unsigned runNo, evNo, nMatches;
  int pid;
  std::vector<unsigned>* lhcbid = nullptr;
  std::vector<float>* x = nullptr;
  std::vector<float>* y = nullptr;
  std::vector<float>* z = nullptr;
  std::vector<float>* t = nullptr;
  float mcEta;
  
  mcTrackTree->SetBranchAddress("runNo", &runNo);
  mcTrackTree->SetBranchAddress("evNo", &evNo);
  mcTrackTree->SetBranchAddress("nMatches", &nMatches);
  mcTrackTree->SetBranchAddress("lhcbid", &lhcbid);
  mcTrackTree->SetBranchAddress("pid", &pid);
  mcTrackTree->SetBranchAddress("x", &x);
  mcTrackTree->SetBranchAddress("y", &y);
  mcTrackTree->SetBranchAddress("z", &z);
  mcTrackTree->SetBranchAddress("t", &t);
  mcTrackTree->SetBranchAddress("eta", &mcEta);

  float min_forward_abs_dt = 1000000.0f;
  float max_forward_abs_dt = 0.0f;
  int n_not_in_eta_acceptance = 0;
  unsigned n_tracks_hit01_same_module = 0;
  unsigned n_tracks_too_few_unique_modules = 0;
  float max_dt_sq = t_res * t_res * t_stddevs * t_stddevs;
  printf("max_dt_sq: %.5f\n", max_dt_sq);

  printf("Going through %lld MCPs\n", mcTrackTree->GetEntries());
  unsigned nOutliers = 0; // define as number of MCPs with hits where t>3ns
  for (unsigned i_track = 0; i_track < mcTrackTree->GetEntries(); i_track++) {
    unsigned nOutliersForwarding = 0;  // this is equivalent to forwarding outlier detection in Allen
    mcTrackTree->GetEntry(i_track);
    if (lhcbid->size() < 3) continue; // need at least 3 hits to do anything with tracking
    if (!Utils::Functions::inEtaAcceptance(mcEta)) {
      n_not_in_eta_acceptance++;
      continue; // only look in nominal eta acceptance
    }
    // track is forward if hits go in increasing module numbers
    if ( ( (lhcbid->at(0) >> 12) & 0x3F ) == ( (lhcbid->at(1) >> 12) & 0x3F ) ) {
      n_tracks_hit01_same_module++;
      continue;
    }
    // NOTE: THE FOLLOWING IS REDUNDANT BECAUSE MCP HITS ARE SORTED BY MODULE FROM 0 TO 63, NOT IN TIME
    bool mcHitsSortedForward = ( (lhcbid->at(0) >> 12) & 0x3F ) < ( (lhcbid->at(1) >> 12) & 0x3F );

    // get the dts for each seeding/forwarding stage & write to right histogram
    float filtered_t = 0.0f;
    // printf("BEFORE: t0, t1, t2: %.3f, %.3f, %.3f\n", t->at(0), t->at(1), t->at(2));
    // Collect the hits in the order the sweep visits them -- backwards tracks in the
    // stored order, forward tracks from the last hit to the first -- and then keep only
    // the first hit on each module. Two hits on one module sit on the same z-plane, and
    // a zero dz sends the extrapolated dt to infinity.
    std::vector<Hit::BaseHit> track_hits;
    track_hits.reserve(lhcbid->size());
    for (unsigned i_hit = 0; i_hit < lhcbid->size(); i_hit++) {
      unsigned index_to_use = !mcHitsSortedForward ? i_hit : lhcbid->size() - i_hit - 1;
      track_hits.emplace_back(lhcbid->at(index_to_use), x->at(index_to_use),
                              y->at(index_to_use), z->at(index_to_use),
                              t->at(index_to_use));
    }
    track_hits = Utils::Functions::get_unique_module_hits(track_hits);
    unsigned kHits = track_hits.size();
    if (kHits < 3) {  // filtering can leave too few hits to seed with
      n_tracks_too_few_unique_modules++;
      continue;
    }
    Hit::BaseHit h0 = track_hits[0];
    Hit::BaseHit h1 = track_hits[1];
    Hit::BaseHit h2 = track_hits[2];
    // get the dts for the seeding stage
    auto res_seed = Utils::Functions::get_seeding_dts(h0, h1, h2);
    if (std::get<3>(res_seed) > 0) { // positive drhodz is forward track
      h_dt0_forward->Fill(std::get<1>(res_seed));
      h_dt2_forward->Fill(std::get<2>(res_seed));
    } else {
      h_dt0_backward->Fill(std::get<1>(res_seed));
      h_dt2_backward->Fill(std::get<2>(res_seed));
    }
    filtered_t = std::get<0>(res_seed);
    // if (lhcbid->at(2) == 2952835394 && evNo == 31) {
    //   printf("h0: %u, %.3f, %.3f, %.3f, %.3f\n", h0.id, h0.x, h0.y, h0.z, h0.t);
    //   printf("h1: %u, %.3f, %.3f, %.3f, %.3f\n", h1.id, h1.x, h1.y, h1.z, h1.t);
    //   printf("h2: %u, %.3f, %.3f, %.3f, %.3f\n", h2.id, h2.x, h2.y, h2.z, h2.t);
    //   printf("Number of hits in track: %u\n", lhcbid->size());
    // }
    // printf("From BaseHits: t0, t1, t2: %.3f, %.3f, %.3f\n", h0.t, h1.t, h2.t);
    // get the dts for the forwarding stage
    for (unsigned i_hit = 3; i_hit < kHits; i_hit++) {
      h0 = h1;
      h1 = h2;
      h2 = track_hits[i_hit];
      auto res = Utils::Functions::get_next_filtered_t(h0, h1, h2, filtered_t, i_hit);
      float dt_next_hit = std::get<1>(res);
      
      if (std::abs(dt_next_hit) < min_forward_abs_dt) min_forward_abs_dt = std::abs(dt_next_hit);
      if (std::abs(dt_next_hit) > max_forward_abs_dt) max_forward_abs_dt = std::abs(dt_next_hit);

      // check if outlier
      if ( ( dt_next_hit * dt_next_hit ) > ( max_dt_sq * (1. + 1. / i_hit) ) ) {
        // we have an outlier so we just use the filtered time as the estimate
        // recall dt = t_est - h2.t, so t_est = h2.t + dt
        nOutliersForwarding++;
        filtered_t = h2.t + dt_next_hit;
        // fill histogram for position of outlier
        if (i_hit > max_track_length_plotting) h_outlier_position->Fill(max_track_length_plotting);
        else h_outlier_position->Fill(i_hit);

      } else {
        filtered_t = std::get<0>(res);
      }
      filtered_t = std::get<0>(res);
      // fill all forwarding histograms
      if (std::get<3>(res_seed) > 0) {
        h_dt_forwarding_forward->Fill(dt_next_hit);
        if (i_hit == 5) h_dt_forwarding_forward_h5->Fill(dt_next_hit);
        if (i_hit == 10) h_dt_forwarding_forward_h10->Fill(dt_next_hit);
        if (i_hit == 15) h_dt_forwarding_forward_h15->Fill(dt_next_hit);
      } else {
        h_dt_forwarding_backward->Fill(dt_next_hit);
        if (i_hit == 5) h_dt_forwarding_backward_h5->Fill(dt_next_hit);
        if (i_hit == 10) h_dt_forwarding_backward_h10->Fill(dt_next_hit);
        if (i_hit == 15) h_dt_forwarding_backward_h15->Fill(dt_next_hit);
      }
      h_dt_forwarding_vs_moduleID->Fill( (h2.id >> 12) & 0x3F, dt_next_hit );
      p_dt_forwarding_vs_moduleID->Fill( (h2.id >> 12) & 0x3F, dt_next_hit );
      // fill with i+1 because we use the dt of the next hit
      h_nthHits_vs_moduleID->Fill( (h2.id >> 12) & 0x3F, i_hit + 1 );
      h_dt_forwarding_vs_nthHit->Fill( i_hit + 1, dt_next_hit );
    }
    // printf("\n");
    // printf("%u outliers in track %u\n", nOutliersForwarding, i_track);
    h_n_outliers->Fill(nOutliersForwarding); // resets to 0 for each track  

    // Then do some printing for outrageous times
    bool hasOutlier = false;

    bool hasCrazyTime = false;
    for (unsigned i_hit = 0; i_hit < lhcbid->size(); i_hit++) {
      if (t->at(i_hit) > tWack) {
        hasCrazyTime = true;
        break;
      }
      if (t->at(i_hit) > 3 || t->at(i_hit) < -3) {
        hasOutlier = true;
      }
    }
    nOutliers += hasOutlier;
    // ignore if no wacky time in MCP
    if (!hasCrazyTime) continue;
    for (unsigned i_hit = 0; i_hit < lhcbid->size(); i_hit++) {
      unsigned modID = (lhcbid->at(i_hit) >> 12) & 0x3F;
      float phi = TMath::ATan2(y->at(i_hit), x->at(i_hit));
      printf("Run %u, Event %u, MCP %u at eta=%.3f with hit id %u at (x=%.3f, y=%.3f, t=%.3f, phi=%.3f) in Module %u with PID: %d\n",
              runNo, evNo, i_track, mcEta, lhcbid->at(i_hit), x->at(i_hit), y->at(i_hit), t->at(i_hit), phi, modID, pid);
    }
    // printf("\tThis MCP with %u matches has hits in the following %lu modules: {",
    //         nMatches, lhcbid->size());
    // for (unsigned i_hit2 = 0; i_hit2 < lhcbid->size(); i_hit2++) {
    //   printf("%u", (lhcbid->at(i_hit2) >> 12) & 0x3F);
    //   if (i_hit2 != lhcbid->size() - 1) printf(", ");
    // }
    // printf("}\n\n");
  }
  printf("Found %u outlier MCPs from a total of %lld\n", nOutliers, mcTrackTree->GetEntries());
  printf("Min forward abs dt: %.3f ns, Max forward abs dt: %.3f ns\n", min_forward_abs_dt, max_forward_abs_dt);
  printf("Number of MCPs not in eta acceptance: %u\n", n_not_in_eta_acceptance);
  printf("Total number of accepted forwarded hits: %u\n", (unsigned) h_dt_forwarding_forward->GetEntries());
  printf("Number of tracks with first two hits in the same module: %u\n", n_tracks_hit01_same_module);
  printf("Number of tracks left with fewer than three unique modules: %u\n", n_tracks_too_few_unique_modules);
  // save histograms, all of them so that the plotting macros do not need the MC file
  TString histOutputDir = (Utils::Definitions::analysisRoot + "hists/4d_tracking").c_str();
  gSystem->mkdir(histOutputDir, true);  // TFile does not create the directory itself
  TString histOutputPath = histOutputDir + "/" + outName + ".root";
  TFile* histFile = TFile::Open(histOutputPath, "RECREATE");
  h_dt0_forward->Write();
  h_dt2_forward->Write();
  h_dt_forwarding_forward->Write();
  h_dt_forwarding_forward_h5->Write();
  h_dt_forwarding_forward_h10->Write();
  h_dt_forwarding_forward_h15->Write();
  h_dt0_backward->Write();
  h_dt2_backward->Write();
  h_dt_forwarding_backward->Write();
  h_dt_forwarding_backward_h5->Write();
  h_dt_forwarding_backward_h10->Write();
  h_dt_forwarding_backward_h15->Write();
  h_dt_forwarding_vs_moduleID->Write();
  p_dt_forwarding_vs_moduleID->Write();
  h_nthHits_vs_moduleID->Write();
  h_dt_forwarding_vs_nthHit->Write();
  h_n_outliers->Write();
  h_outlier_position->Write();
  std::cout << "Histograms saved to: " << histOutputPath << std::endl;
  // clean up
  histFile->Close();
  delete histFile;
  file->Close();
  delete file;
  return;
}