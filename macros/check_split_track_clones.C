#include <TFile.h>
#include <TH1F.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <string>

#include "utils/definitions.h"
#include "utils/basic_functions.h"
#include "utils/Hit.h"

void check_split_track_clones(unsigned nEvents=5000, unsigned max_scatter=80,
                              unsigned max_dt=0, unsigned kEventsPerRun=200) {
  gROOT->SetBatch();  // so stuff isn't autoplotted
  TString suffix = Utils::Functions::get_suffix(nEvents, max_scatter, max_dt);
  TString input_suffix = suffix + ".root";
  TString input_prefix = (Utils::Definitions::stackRoot + "output/MCData_Checking").c_str();
  TString filepath = input_prefix + input_suffix;
  TFile* file = TFile::Open(filepath);
  if (!file || file->IsZombie()) {
    std::cerr << "Error: Could not open ROOT file." << std::endl;
    return;
  }
  // Get module to z conversion histogram
  TFile *moduleFile = TFile::Open(
    (Utils::Definitions::analysisRoot + "hists/module_mc_info.root").c_str());
  if (!moduleFile || moduleFile->IsZombie()) {
      std::cerr << "Error: Could not open ROOT file." << std::endl;
      return;
  }
  TH1D* moduleToZ = (TH1D*) moduleFile->Get("module_to_z");
  // Trees that will be relevant
  TTree* mcTrackTree = (TTree*) file->Get("MCTrackData");
  if (!mcTrackTree) {
      std::cerr << "MCTrackData TTree not found!" << std::endl;
      return;
  }
  TTree* recoTrackTree = (TTree*) file->Get("RecoTrackData");
  if (!recoTrackTree) {
      std::cerr << "RecoTrackData TTree not found!" << std::endl;
      return;
  }
  TTree* recoEventTree = (TTree*) file->Get("RecoEventData");
  if (!recoEventTree) {
      std::cerr << "RecoEventData TTree not found!" << std::endl;
      return;
  }
  // ----------- BINS ---------------
  // phi bins
  int nPhiBins = 50;
  std::vector<float> phiBinEdges(nPhiBins + 1);
  float phiMax{0.1}, phiMin{-0.1};  // short range for granularity
  float phiStep = (phiMax - phiMin) / (nPhiBins);
  for (unsigned i = 0; i <= nPhiBins; i++) {
    phiBinEdges[i] = phiMin + phiStep * i;
  }
  // pt bins
  int nPTBins = 50;
  std::vector<float> ptBinEdges(nPTBins + 1);
  float ptMax{5000.}, ptMin{0.};  // only up to TeV
  float ptStep = (ptMax - ptMin) / (nPTBins);
  for (unsigned i = 0; i <= nPTBins; i++) {
    ptBinEdges[i] = ptMin + ptStep * i;
  }
  // dr bins (deflection)
  int nDeflectionBins = 50;
  std::vector<float> deflectionBinEdges(nDeflectionBins + 1);
  float deflectionMax{0.5}, deflectionMin{0.};
  float deflectionStep = (deflectionMax - deflectionMin) / (nDeflectionBins);
  for (unsigned i = 0; i <= nDeflectionBins; i++) {
    deflectionBinEdges[i] = deflectionMin + deflectionStep * i;
  }
  // deflection bins scaled by dz
  std::vector<float> deflectionBinEdgesDZScaled(nDeflectionBins + 1);
  std::vector<float> deflectionBinEdgesDZSQScaled(nDeflectionBins + 1);
  for (unsigned i = 0; i <= nDeflectionBins; i++) {
    // roughly 40mm per module pair but scaled more for larger displaced ones
    deflectionBinEdgesDZScaled[i] = deflectionBinEdges[i] / 100.;
    deflectionBinEdgesDZSQScaled[i] = deflectionBinEdges[i] / (100. * 100.);
  }
  // ----------------- BRANCH DATA ----------------
  unsigned mcTrackRunNo, mcTrackEvNo, nMatches, recoEvOffset, nMCVeloHits;
  int mcMatchIdxForReco;  // to access right reco track with offset
  float recoEta, mcPT;
  std::vector<unsigned>* matchedRecoTrackIndices = nullptr;
  std::vector<unsigned>* recoTrackLHCbIDs = nullptr;
  std::vector<unsigned>* mcTrackLHCbIDs = nullptr;
  std::vector<float>* x_mc = nullptr;
  std::vector<float>* y_mc = nullptr;
  std::vector<float>* x_reco = nullptr;
  std::vector<float>* y_reco = nullptr;
  // Set addresses
  // MC Track Tree
  mcTrackTree->SetBranchAddress("runNo", &mcTrackRunNo);
  mcTrackTree->SetBranchAddress("evNo", &mcTrackEvNo);
  mcTrackTree->SetBranchAddress("nHitsVelo", &nMCVeloHits);
  mcTrackTree->SetBranchAddress("nMatches", &nMatches);
  mcTrackTree->SetBranchAddress("pt", &mcPT);
  mcTrackTree->SetBranchAddress("matchedTracks", &matchedRecoTrackIndices);
  mcTrackTree->SetBranchAddress("lhcbid", &mcTrackLHCbIDs);
  mcTrackTree->SetBranchAddress("x", &x_mc);
  mcTrackTree->SetBranchAddress("y", &y_mc);
  // Reco Trees
  recoEventTree->SetBranchAddress("globalTrackOffset", &recoEvOffset);
  recoTrackTree->SetBranchAddress("lhcbid", &recoTrackLHCbIDs);
  recoTrackTree->SetBranchAddress("mcMatchIdx", &mcMatchIdxForReco);
  recoTrackTree->SetBranchAddress("eta", &recoEta);
  recoTrackTree->SetBranchAddress("x", &x_reco);
  recoTrackTree->SetBranchAddress("y", &y_reco);

  // Here we want to go through specifically the split tracks
  // Create profile for delta phi and 2d hist for pt/phi
  TH1D* h_deltaPhi_splitTrackClone = new TH1D(
    "delta_phi", "#Delta#phi Distribution;#eta;#Delta#phi", nPhiBins, phiBinEdges.data());
  TH1D* h_deltaPhi_splitTrackClone_1Missed = new TH1D(
    "delta_phi_1Missed", "#Delta#phi Distribution;#eta;#Delta#phi", nPhiBins, phiBinEdges.data());
  TH1D* h_deltaPhi_splitTrackClone_2Missed = new TH1D(
    "delta_phi_2Missed", "#Delta#phi Distribution;#eta;#Delta#phi", nPhiBins, phiBinEdges.data());
  TH1D* h_deflection_splitTrackClone = new TH1D(
    "deflection", "Deflection Distribution;#eta;Deflection", nDeflectionBins, deflectionBinEdges.data());
  TH1D* h_deflection_splitTrackClone_1Missed = new TH1D(
    "deflection_1Missed", "Deflection Distribution;#eta;Deflection", nDeflectionBins, deflectionBinEdges.data());
  TH1D* h_deflection_splitTrackClone_2Missed = new TH1D(
    "deflection_2Missed", "Deflection Distribution;#eta;Deflection", nDeflectionBins, deflectionBinEdges.data());
  TH1D* h_deflection_per_z_splitTrackClone = new TH1D(
    "deflection_per_z", "Deflection per z Distribution;#eta;Deflection per z", nDeflectionBins, deflectionBinEdgesDZScaled.data());
  TH1D* h_deflection_per_z_splitTrackClone_1Missed = new TH1D(
    "deflection_per_z_1Missed", "Deflection per z Distribution;#eta;Deflection per z", nDeflectionBins, deflectionBinEdgesDZScaled.data());
  TH1D* h_deflection_per_z_splitTrackClone_2Missed = new TH1D(
    "deflection_per_z_2Missed", "Deflection per z Distribution;#eta;Deflection per z", nDeflectionBins, deflectionBinEdgesDZScaled.data());
  TH1D* h_deflection_per_z_sq_splitTrackClone = new TH1D(
    "deflection_per_z_sq", "Deflection per z^2 Distribution;#eta;Deflection per z^2", nDeflectionBins, deflectionBinEdgesDZSQScaled.data());
  TH1D* h_deflection_per_z_sq_splitTrackClone_1Missed = new TH1D(
    "deflection_per_z_sq_1Missed", "Deflection per z^2 Distribution;#eta;Deflection per z^2", nDeflectionBins, deflectionBinEdgesDZSQScaled.data());
  TH1D* h_deflection_per_z_sq_splitTrackClone_2Missed = new TH1D(
    "deflection_per_z_sq_2Missed", "Deflection per z^2 Distribution;#eta;Deflection per z^2", nDeflectionBins, deflectionBinEdgesDZSQScaled.data());
  TH1D* h_pT_reference = new TH1D(
    "pT_reference", "pT Distribution;#eta;pT", nPTBins, ptBinEdges.data());
  TH1D* h_pT_splitTrackClone = new TH1D(
    "pT", "pT Distribution;#eta;pT", nPTBins, ptBinEdges.data());
  TH1D* h_pT_splitTrackClone_1Missed = new TH1D(
    "pT_1Missed", "pT Distribution;#eta;pT", nPTBins, ptBinEdges.data());
  TH1D* h_pT_splitTrackClone_2Missed = new TH1D(
    "pT_2Missed", "pT Distribution;#eta;pT", nPTBins, ptBinEdges.data());
  TH2D* h_phi_pT_splitTrackClone = new TH2D(
    "delta_phi_pT", "pT vs phi (LO Split Track);#phi;pT", nPhiBins, phiBinEdges.data(), nPTBins, ptBinEdges.data());
  TH2D* h_phi_pT_splitTrackClone_1Missed = new TH2D(
    "delta_phi_pT_1Missed", "pT vs phi (NLO Split Track);#phi;pT", nPhiBins, phiBinEdges.data(), nPTBins, ptBinEdges.data());
  TH2D* h_phi_pT_splitTrackClone_2Missed = new TH2D(
    "delta_phi_pT_2Missed", "pT vs phi (NNLO Split Track);#phi;pT", nPhiBins, phiBinEdges.data(), nPTBins, ptBinEdges.data());

  // Go through MC particles
  unsigned nMCTracks = mcTrackTree->GetEntries();
  printf("Going through %u MC particles now.\n", nMCTracks);
  for (unsigned i_mct = 0; i_mct < nMCTracks; i_mct++) {
    // print progress
    if (i_mct && (i_mct % 1000000 == 0))
      printf("Finished with %d%% of MC Tracks.\n", (int) (100. * (float) i_mct / nMCTracks));
    
    mcTrackTree->GetEntry(i_mct);
    h_pT_reference->Fill(mcPT);  // fill always, not just with clones
    // don't continue if without clones
    if (nMatches < 2 || !matchedRecoTrackIndices) continue;
    unsigned kTotalIDs = 0;
    // To check if MC track to be classified as split
    std::unordered_set<unsigned> mcTrackLHCbIDsSet(mcTrackLHCbIDs->begin(),
                                                   mcTrackLHCbIDs->end());

    Hit::BaseHit t0_penultimateHit, t0_lastHit, t1_firstHit;
     // delta z, x, y, & phi between two tracks (first point on 2nd)
    float dz_to_next_module, dx, dy, dPhi;
    for (unsigned int i_rt = 0; i_rt < nMatches; i_rt++) {
      unsigned matchIdx = matchedRecoTrackIndices->at(i_rt);
      // get relevant reco track index
      unsigned evIdx = (mcTrackRunNo - 1) * kEventsPerRun + (mcTrackEvNo - 1);
      recoEventTree->GetEntry(evIdx);
      unsigned recoTrackIdx = recoEvOffset + matchIdx;
      recoTrackTree->GetEntry(recoTrackIdx);
      if (!recoTrackLHCbIDs) {
        std::cout << "Nothing found in reco track tree for lhcbid at index " << recoTrackIdx << ".\n";
        continue;  // if nothing there, this is weird...
      }
      unsigned nRecoHits = recoTrackLHCbIDs->size();
      kTotalIDs += nRecoHits;
      for (unsigned i_hit_reco = 0; i_hit_reco < nRecoHits; i_hit_reco++) {
        unsigned reco_lhcbid = recoTrackLHCbIDs->at(i_hit_reco);
        mcTrackLHCbIDsSet.erase(reco_lhcbid);
        unsigned moduleNumber = (reco_lhcbid >> 12) & 0x3F;
        float z = moduleToZ->GetBinContent(moduleNumber + 1);  // +1 because ROOT 1 indexed 
        if (i_rt == 0 && i_hit_reco == nRecoHits - 2) {
          t0_penultimateHit = Hit::BaseHit(reco_lhcbid, x_reco->at(i_hit_reco),
                                           y_reco->at(i_hit_reco), z, 0.);
        } else if (i_rt == 0 && i_hit_reco == nRecoHits - 1) {
          // both last hits required to calculate estimated phi
          t0_lastHit = Hit::BaseHit(reco_lhcbid, x_reco->at(i_hit_reco),
                                    y_reco->at(i_hit_reco), z, 0.);
          
          float dz = t0_lastHit.z - t0_penultimateHit.z;
          // using two modules until next hit as the estimate is gotten like that too:
          // module_pair_data[shared::next_module_pair].z[h0_module % 2] - h0.z
          // remember to reset to zero before doing next difference of split track
          dz_to_next_module = 0;
          if (dz < 0 && moduleNumber > 1) {
            // track moving forwards, trace out detector backwards (this is what we shoudl always do)
            // jump by steps of 2 since dealing with module pairs
            dz_to_next_module = moduleToZ->GetBinContent(moduleNumber - 1) - t0_penultimateHit.z;
          } else {
            // detector being traced forwards
            dz_to_next_module = moduleToZ->GetBinContent(moduleNumber + 3) - t0_penultimateHit.z;
          }
          auto estimate = Hit::estimateNextPhi(t0_penultimateHit, t0_lastHit, dz_to_next_module);
          dx = std::get<0>(estimate);
          dy = std::get<1>(estimate);
          dPhi = std::get<2>(estimate);
        } else if (i_rt == 1 && i_hit_reco == 0) {
          t1_firstHit = Hit::BaseHit(reco_lhcbid, x_reco->at(i_hit_reco),
                                     y_reco->at(i_hit_reco), z, 0.);
          // subtract to get deltas once in next track
          dx -= x_reco->at(i_hit_reco);
          dy -= y_reco->at(i_hit_reco);
          dPhi -= TMath::ATan2(y_reco->at(i_hit_reco), x_reco->at(i_hit_reco));
        }
      }
    }
    // all MC hits reconstructed within all reco matched tracks
    bool isSplitTrack = mcTrackLHCbIDsSet.size() == 0 && kTotalIDs == nMCVeloHits &&
                        nMatches > 1;
    // Either one hit missing between, one overlap, or one extra hit (covered by case 2)
    bool isSplitTrack_1Missed = ( (mcTrackLHCbIDsSet.size() == 1 && kTotalIDs == (nMCVeloHits - 1)) ||
                                  (mcTrackLHCbIDsSet.size() == 0 && kTotalIDs == (nMCVeloHits + 1)) ) &&
                                   nMatches > 1;
    // Two hits missing, two overlap (or two extra), or one missing and one extra
    bool isSplitTrack_2Missed = ( (mcTrackLHCbIDsSet.size() == 2 && kTotalIDs == (nMCVeloHits - 2)) ||
                                  (mcTrackLHCbIDsSet.size() == 0 && kTotalIDs == (nMCVeloHits + 2)) ||
                                  (mcTrackLHCbIDsSet.size() == 1 && kTotalIDs == (nMCVeloHits + 1)) ) &&
                                   nMatches > 1;

    // if not split track, continue
    if ( !(isSplitTrack || isSplitTrack_1Missed || isSplitTrack_2Missed) ) continue;
    
    float dr_sq = dx * dx + dy * dy;
    float dr_sq_per_z = dr_sq / abs(dz_to_next_module);
    float dr_sq_per_z_sq = dr_sq / (dz_to_next_module * dz_to_next_module);
    // Fill histograms
    if (isSplitTrack) {
      h_deltaPhi_splitTrackClone->Fill(dPhi);
      h_deflection_splitTrackClone->Fill(dr_sq);
      h_deflection_per_z_splitTrackClone->Fill(dr_sq_per_z);
      h_deflection_per_z_sq_splitTrackClone->Fill(dr_sq_per_z_sq);
      h_pT_splitTrackClone->Fill(mcPT);
      h_phi_pT_splitTrackClone->Fill(dPhi, mcPT);
    } else if (isSplitTrack_1Missed) {
      h_deltaPhi_splitTrackClone_1Missed->Fill(dPhi);
      h_deflection_splitTrackClone_1Missed->Fill(dr_sq);
      h_deflection_per_z_splitTrackClone_1Missed->Fill(dr_sq_per_z);
      h_deflection_per_z_sq_splitTrackClone_1Missed->Fill(dr_sq_per_z_sq);
      h_pT_splitTrackClone_1Missed->Fill(mcPT);
      h_phi_pT_splitTrackClone_1Missed->Fill(dPhi, mcPT);
    } else if (isSplitTrack_2Missed) {
      h_deltaPhi_splitTrackClone_2Missed->Fill(dPhi);
      h_deflection_splitTrackClone_2Missed->Fill(dr_sq);
      h_deflection_per_z_splitTrackClone_2Missed->Fill(dr_sq_per_z);
      h_deflection_per_z_sq_splitTrackClone_2Missed->Fill(dr_sq_per_z_sq);
      h_pT_splitTrackClone_2Missed->Fill(mcPT);
      h_phi_pT_splitTrackClone_2Missed->Fill(dPhi, mcPT);
    }
  }  // MC Particles loop
  // Write histograms and clean up
  TString outPrefix =
    (Utils::Definitions::analysisRoot + "/hists/clones/split_track_hists").c_str();
  TFile* outFile = new TFile(outPrefix + input_suffix, "RECREATE");
  h_deltaPhi_splitTrackClone->Write();
  h_deltaPhi_splitTrackClone_1Missed->Write();
  h_deltaPhi_splitTrackClone_2Missed->Write();
  h_deflection_splitTrackClone->Write();
  h_deflection_splitTrackClone_1Missed->Write();
  h_deflection_splitTrackClone_2Missed->Write();
  h_deflection_per_z_splitTrackClone->Write();
  h_deflection_per_z_splitTrackClone_1Missed->Write();
  h_deflection_per_z_splitTrackClone_2Missed->Write();
  h_deflection_per_z_sq_splitTrackClone->Write();
  h_deflection_per_z_sq_splitTrackClone_1Missed->Write();
  h_deflection_per_z_sq_splitTrackClone_2Missed->Write();
  h_pT_reference->Write();
  h_pT_splitTrackClone->Write();
  h_pT_splitTrackClone_1Missed->Write();
  h_pT_splitTrackClone_2Missed->Write();
  h_phi_pT_splitTrackClone->Write();
  h_phi_pT_splitTrackClone_1Missed->Write();
  h_phi_pT_splitTrackClone_2Missed->Write();

  outFile->Close();
  delete outFile;
  moduleFile->Close();
  delete moduleFile;
  file->Close();
  delete file;
}