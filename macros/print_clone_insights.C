#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <string>

std::string stackRoot = std::getenv("STACK_ROOT");
std::string analysisRoot = std::getenv("ANALYSIS_ROOT");

void print_clone_insights(unsigned kEventsPerRun=200, unsigned nClonesMinimum=8,
                          unsigned nClonesMaximum=20, unsigned nHitsMinimum=3,
                          unsigned nHitsMaximum=8, std::string mcFilePath="MCData_Checking.root") {
  TFile *file = TFile::Open((stackRoot + mcFilePath).c_str());
  if (!file || file->IsZombie()) {
      std::cerr << "Error: Could not open ROOT file." << std::endl;
      return;
  }
  // get all trees
  TTree* mcTrackTree = (TTree*) file->Get("MCTrackData");
  if (!mcTrackTree) {
      std::cerr << "TTree not found!" << std::endl;
      return;
  }
  TTree* mcEventTree = (TTree*) file->Get("MCEventData");
  if (!mcEventTree) {
      std::cerr << "TTree not found!" << std::endl;
      return;
  }
  TTree* recoTrackTree = (TTree*) file->Get("RecoTrackData");
  if (!recoTrackTree) {
      std::cerr << "TTree not found!" << std::endl;
      return;
  }
  TTree* recoEventTree = (TTree*) file->Get("RecoEventData");
  if (!recoEventTree) {
      std::cerr << "TTree not found!" << std::endl;
      return;
  }
  unsigned nMatches, mcTrackEvNo, mcTrackRunNo, recoEvOffset, nMCVeloHits,
           recoTrackEvNo, recoTrackRunNo;
  float doca_t_reco, mcEta, recoEta, recoChi2, pt;
  bool isBackward;
  std::vector<unsigned>* matchedRecoTrackIndices = nullptr;
  std::vector<unsigned>* recoTrackLHCbIDs = nullptr;
  std::vector<unsigned>* mcTrackLHCbIDs = nullptr;
  std::vector<float>* t_mc = nullptr;
  std::vector<float>* x_mc = nullptr;
  std::vector<float>* y_mc = nullptr;
  std::vector<float>* t_reco = nullptr;
  std::vector<float>* x_reco = nullptr;
  std::vector<float>* y_reco = nullptr;

  float t_min = 1, t_max = 0;

  // MC Track Tree
  mcTrackTree->SetBranchAddress("matchedTracks", &matchedRecoTrackIndices);
  mcTrackTree->SetBranchAddress("nMatches", &nMatches);
  mcTrackTree->SetBranchAddress("evNo", &mcTrackEvNo);
  mcTrackTree->SetBranchAddress("runNo", &mcTrackRunNo);
  mcTrackTree->SetBranchAddress("nHitsVelo", &nMCVeloHits);
  mcTrackTree->SetBranchAddress("lhcbid", &mcTrackLHCbIDs);
  mcTrackTree->SetBranchAddress("t", &t_mc);
  mcTrackTree->SetBranchAddress("x", &x_mc);
  mcTrackTree->SetBranchAddress("y", &y_mc);
  mcTrackTree->SetBranchAddress("eta", &mcEta);
  mcTrackTree->SetBranchAddress("pt", &pt);
  // Reco Event Tree
  recoEventTree->SetBranchAddress("globalTrackOffset", &recoEvOffset);
  // Reco Track Tree
  recoTrackTree->SetBranchAddress("evNo", &recoTrackEvNo);
  recoTrackTree->SetBranchAddress("runNo", &recoTrackRunNo);
  recoTrackTree->SetBranchAddress("lhcbid", &recoTrackLHCbIDs);
  recoTrackTree->SetBranchAddress("isBackward", &isBackward);
  recoTrackTree->SetBranchAddress("t_k", &doca_t_reco);
  recoTrackTree->SetBranchAddress("t", &t_reco);
  recoTrackTree->SetBranchAddress("x", &x_reco);
  recoTrackTree->SetBranchAddress("y", &y_reco);
  recoTrackTree->SetBranchAddress("eta", &recoEta);
  recoTrackTree->SetBranchAddress("chi2", &recoChi2);

  unsigned nRecoTrackLongerThanMC = 0;
  // Now let's go through all MC tracks
  unsigned nMCTracks = mcTrackTree->GetEntries();
  for (unsigned i_mct = 0; i_mct < nMCTracks; i_mct++) {
    mcTrackTree->GetEntry(i_mct);
    // first just check the time range
    for (unsigned int i_lhcbid = 0; i_lhcbid < mcTrackLHCbIDs->size(); i_lhcbid++) {
      float t = t_mc->at(i_lhcbid);
      if (t < t_min) t_min = t;
      else if (t > t_max) t_max = t;
    }
    if (nMCVeloHits > nHitsMaximum || nMCVeloHits < nHitsMinimum ||
        nMatches <= nClonesMinimum ||  nMatches > (nClonesMaximum+1)) continue;
  
    printf("Run %u, Ev %u, MC Track %u at eta %.3f & pT %.3f with %u hits {",
           mcTrackRunNo, mcTrackEvNo, i_mct, mcEta, pt, nMCVeloHits);
    for (unsigned int i_lhcbid = 0; i_lhcbid < mcTrackLHCbIDs->size(); i_lhcbid++) {
      unsigned lhcbid = mcTrackLHCbIDs->at(i_lhcbid);
      unsigned moduleNumber = (lhcbid >> 12) & 0x3F;
      float phi = TMath::ATan2(y_mc->at(i_lhcbid), x_mc->at(i_lhcbid));
      printf("%u (id=%u, t=%.3f)", moduleNumber, lhcbid, t_mc->at(i_lhcbid));
      if (i_lhcbid != mcTrackLHCbIDs->size() - 1) printf(", ");
    }
    printf("} has %u clones.\n", nMatches-1);
    for (unsigned matchIdx : *matchedRecoTrackIndices) {
      // check if backward or not
      if (isBackward) printf("\t[B] ");
      else printf("\t[F] ");
      unsigned evIdx = (mcTrackRunNo - 1) * kEventsPerRun + (mcTrackEvNo - 1);
      recoEventTree->GetEntry(evIdx);
      unsigned recoTrackIdx = recoEvOffset + matchIdx;
      recoTrackTree->GetEntry(recoTrackIdx);
      if (!recoTrackLHCbIDs) {
        std::cout << "Nothing found in reco track tree for lhcbid at index " << recoTrackIdx << ".\n";
        continue;  // if nothing there, this is weird...
      }
      unsigned nRecoHits = recoTrackLHCbIDs->size();
      nRecoTrackLongerThanMC += (nRecoHits > nMCVeloHits);

      printf("Run %u, Ev %u, Reco Track %u with %u hits, chi2 %.3f & eta %.3f: {",
             recoTrackRunNo, recoTrackEvNo, matchIdx, nRecoHits, recoChi2, recoEta);
      for (unsigned int i_lhcbid = 0; i_lhcbid < nRecoHits; i_lhcbid++) {
        unsigned lhcbid = recoTrackLHCbIDs->at(i_lhcbid);
        unsigned moduleNumber = (lhcbid >> 12) & 0x3F;
        float phi = TMath::ATan2(y_reco->at(i_lhcbid), x_reco->at(i_lhcbid));
        printf("%u (phi=%.3f, t=%.3f)", moduleNumber, phi, t_reco->at(i_lhcbid));
        if (i_lhcbid != nRecoHits - 1) printf(", ");
      }
      printf("}\n");
    }
    printf("\n");
  }
  printf("Largest and smallest times throughout all events: (%.3f, %.3f)\n", t_min, t_max);
  printf("Number of reco tracks longer than MC tracks: %u\n", nRecoTrackLongerThanMC);

  return;
}