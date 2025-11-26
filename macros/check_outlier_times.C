#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <string>
#include <tuple>

#include "utils/definitions.h"
#include "utils/allen_functions.h"
#include "utils/basic_functions.h"

void check_outlier_times(std::string mcFilePath="MCData_Checking.root",
                         float t_res=0.05, float t_stddevs=3., float tWack=700.0) {
  gStyle->SetOptStat(0); // remove the info box for the plots

  TFile *file = TFile::Open((Utils::Definitions::stackRoot + mcFilePath).c_str());
  if (!file || file->IsZombie()) {
      std::cerr << "Error: Could not open ROOT file." << std::endl;
      return;
  }
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
    Hit::BaseHit h0, h1, h2;
    unsigned kHits = lhcbid->size();
    if (!mcHitsSortedForward) {  // backwards track, i.e. will sweep in already ordered direction
      h0 = {lhcbid->at(0), x->at(0), y->at(0), z->at(0), t->at(0)};
      h1 = {lhcbid->at(1), x->at(1), y->at(1), z->at(1), t->at(1)};
      h2 = {lhcbid->at(2), x->at(2), y->at(2), z->at(2), t->at(2)};
    } else {  // forward track, i.e. sweep from last to first hit
      h0 = {lhcbid->at(kHits - 1), x->at(kHits - 1), y->at(kHits - 1), z->at(kHits - 1), t->at(kHits - 1)};
      h1 = {lhcbid->at(kHits - 2), x->at(kHits - 2), y->at(kHits - 2), z->at(kHits - 2), t->at(kHits - 2)};
      h2 = {lhcbid->at(kHits - 3), x->at(kHits - 3), y->at(kHits - 3), z->at(kHits - 3), t->at(kHits - 3)};
    }
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
      unsigned index_to_use = !mcHitsSortedForward ? i_hit : kHits - i_hit - 1;
      h2 = Hit::BaseHit(lhcbid->at(index_to_use), x->at(index_to_use),
                        y->at(index_to_use), z->at(index_to_use),
                        t->at(index_to_use));
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
      h_dt_forwarding_vs_moduleID->Fill( (lhcbid->at(index_to_use) >> 12) & 0x3F, dt_next_hit );
      // fill with i+1 because we use the dt of the next hit
      h_nthHits_vs_moduleID->Fill( (lhcbid->at(index_to_use) >> 12) & 0x3F, i_hit + 1 );
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
  // save histograms and make canvas
  TCanvas* canvas = new TCanvas("canvas", "Outlier Times", 1200, 1200);
  canvas->Divide(2, 3);
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
    (int) ( h_dt0_forward->GetMean() * 1000 ), (int) (h_dt0_forward->GetStdDev() * 1000) ), "l");
  legend_dt_seeding->AddEntry(h_dt2_forward, Form("#Deltat_{2}^{f} (#mu, #sigma) = (%d, %d) ps",
    (int) ( h_dt2_forward->GetMean() * 1000 ), (int) (h_dt2_forward->GetStdDev() * 1000) ), "l");
  legend_dt_seeding->AddEntry(h_dt0_backward, Form("#Deltat_{0}^{b} (#mu, #sigma) = (%d, %d) ps",
    (int) ( h_dt0_backward->GetMean() * 1000 ), (int) (h_dt0_backward->GetStdDev() * 1000) ), "l");
  legend_dt_seeding->AddEntry(h_dt2_backward, Form("#Deltat_{2}^{b} (#mu, #sigma) = (%d, %d) ps",
    (int) ( h_dt2_backward->GetMean() * 1000 ), (int) (h_dt2_backward->GetStdDev() * 1000) ), "l");
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
  float mu_all = h_dt_forwarding_forward->GetMean();
  float mu_5 = h_dt_forwarding_forward_h5->GetMean();
  float mu_10 = h_dt_forwarding_forward_h10->GetMean();
  float mu_15 = h_dt_forwarding_forward_h15->GetMean();
  TLegend* legend = new TLegend(0.6, 0.6, 0.88, 0.8);
  legend->SetBorderSize(0);
  legend->AddEntry(h_dt_forwarding_forward_copy, Form("All forwarded (#mu = %d ps)", (int) (mu_all * 1000)), "l");
  legend->AddEntry(h_dt_forwarding_forward_h5, Form("Hit index 5 (#mu = %d ps)", (int) (mu_5 * 1000)), "l");
  legend->AddEntry(h_dt_forwarding_forward_h10, Form("Hit index 10 (#mu = %d ps)", (int) (mu_10 * 1000)), "l");
  legend->AddEntry(h_dt_forwarding_forward_h15, Form("Hit index 15 (#mu = %d ps)", (int) (mu_15 * 1000)), "l");
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
  float mu_all_backward = h_dt_forwarding_backward->GetMean();
  float mu_5_backward = h_dt_forwarding_backward_h5->GetMean();
  float mu_10_backward = h_dt_forwarding_backward_h10->GetMean();
  float mu_15_backward = h_dt_forwarding_backward_h15->GetMean();

  TLegend* legend_backward = new TLegend(0.55, 0.7, 0.89, 0.8);
  legend_backward->SetBorderSize(0);
  legend_backward->AddEntry(h_dt_forwarding_backward_copy, Form("All forwarded (#mu = %d ps)", (int) (mu_all_backward * 1000)), "l");
  legend_backward->AddEntry(h_dt_forwarding_backward_h5, Form("Hit index 5 (#mu = %d ps)", (int) (mu_5_backward * 1000)), "l");
  legend_backward->AddEntry(h_dt_forwarding_backward_h10, Form("Hit index 10 (#mu = %d ps)", (int) (mu_10_backward * 1000)), "l");
  // legend_backward->AddEntry(h_dt_forwarding_backward_h15, Form("Hit index 15 (#mu = %d ps)", (int) (mu_15_backward * 1000)), "l");
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

  // save canvas
  TString outputPath = (Utils::Definitions::analysisRoot + "output/4d_tracking/outlier_times.pdf").c_str();
  canvas->SaveAs(outputPath);
  std::cout << "Canvas saved to: " << outputPath << std::endl;

  // save histograms
  TString histOutputPath = (Utils::Definitions::analysisRoot + "hists/4d_tracking/outlier_times_hists.root").c_str();
  TFile* histFile = TFile::Open(histOutputPath, "RECREATE");
  h_dt0_forward->Write();
  h_dt2_forward->Write();
  h_dt_forwarding_forward->Write();
  h_n_outliers->Write();
  h_outlier_position->Write();
  // clean up
  delete canvas;
  histFile->Close();
  delete histFile;
  file->Close();
  delete file;
  return;
}