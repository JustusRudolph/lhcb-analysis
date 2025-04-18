#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <string>

#include "utils/definitions.h"
#include "utils/basic_functions.h"
#include "utils/Hit.h"

void print_clone_insights(unsigned kEventsPerRun=200, unsigned nClonesMinimum=8,
                          unsigned nClonesMaximum=20, unsigned nHitsMinimum=3,
                          unsigned nHitsMaximum=8, bool printRecoIDs=false,
                          std::string mcFilePath="output/MCData_Checking5000.root") {
  // first get module to z conversion histogram
  TFile *moduleFile = TFile::Open(
    (Utils::Definitions::analysisRoot + "hists/module_mc_info.root").c_str());
  if (!moduleFile || moduleFile->IsZombie()) {
      std::cerr << "Error: Could not open ROOT file." << std::endl;
      return;
  }
  TH1D* moduleToZ = (TH1D*) moduleFile->Get("module_to_z");
  // print the bin edges
  // for (int i = 1; i <= moduleToZ->GetNbinsX(); ++i) {
  //   std::cout << "Bin " << i << ": " << moduleToZ->GetBinCenter(i)
  //             << " with value " << moduleToZ->GetBinContent(i) << std::endl;
  // }
  TFile *file = TFile::Open(
    (Utils::Definitions::stackRoot + mcFilePath).c_str());
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

    TString mcTrackString =
      Form("Run %u, Ev %u, MC Track %u at eta %.3f & pT %.3f with %u hits {",
           mcTrackRunNo, mcTrackEvNo, i_mct, mcEta, pt, nMCVeloHits);
    for (unsigned int i_lhcbid = 0; i_lhcbid < mcTrackLHCbIDs->size(); i_lhcbid++) {
      unsigned mc_lhcbid = mcTrackLHCbIDs->at(i_lhcbid);
      unsigned moduleNumber = (mc_lhcbid >> 12) & 0x3F;
      float phi = TMath::ATan2(y_mc->at(i_lhcbid), x_mc->at(i_lhcbid));
      if (printRecoIDs)
        mcTrackString += Form("%u (%u)", moduleNumber, mc_lhcbid);
      else
        mcTrackString +=
          Form("%u (phi=%.3f, t=%.3f)", moduleNumber, phi, t_mc->at(i_lhcbid));
      if (i_lhcbid != mcTrackLHCbIDs->size() - 1) mcTrackString += ", ";
    }
    mcTrackString += Form("} has %u clones.\n", nMatches-1);
    // Track clone type to print at the end
    unsigned nTriplets = 0;
    unsigned kTotalIDs = 0;
    std::unordered_set<unsigned> mcTrackLHCbIDsSet(mcTrackLHCbIDs->begin(),
                                                   mcTrackLHCbIDs->end());
    // track the lhcbIDs of the first three hits of all matched tracks to check for seeding clones
    std::unordered_set<unsigned> firstHitLHCbIDs;
    std::unordered_set<unsigned> secondHitLHCbIDs;
    std::unordered_set<unsigned> thirdHitLHCbIDs;
    // keep track of last two hits for extrapolation simulation
    Hit::BaseHit h0, h1;
    for (unsigned matchIdx : *matchedRecoTrackIndices) {
      // check if backward or not (using eta is more reliable)
      if (recoEta < 0) mcTrackString += "\t[B] ";
      else mcTrackString += "\t[F] ";
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
      nTriplets += (nRecoHits == 3);
      kTotalIDs += nRecoHits;

      mcTrackString +=
        Form("Run %u, Ev %u, Reco Track %u with %u hits, chi2 %.3f & eta %.3f: {",
             recoTrackRunNo, recoTrackEvNo, matchIdx, nRecoHits, recoChi2, recoEta);
      for (unsigned int i_hit_reco = 0; i_hit_reco < nRecoHits; i_hit_reco++) {
        unsigned reco_lhcbid = recoTrackLHCbIDs->at(i_hit_reco);
        unsigned moduleNumber = (reco_lhcbid >> 12) & 0x3F;
        float phi = TMath::ATan2(y_reco->at(i_hit_reco), x_reco->at(i_hit_reco));
        mcTrackLHCbIDsSet.erase(reco_lhcbid);
        // print either the reco ID or the phi and time
        if (printRecoIDs)
          mcTrackString += Form("%u (%u", moduleNumber, reco_lhcbid);
        else
          mcTrackString +=
            Form("%u (phi=%.3f, t=%.3f", moduleNumber, phi, t_reco->at(i_hit_reco));
        
        // Also check first three hits to see if we get seeding clones
        if (i_hit_reco == 0) {
          firstHitLHCbIDs.insert(reco_lhcbid);
        } else if (i_hit_reco == 1) {
          secondHitLHCbIDs.insert(reco_lhcbid);
        } else if (i_hit_reco == 2) {
          thirdHitLHCbIDs.insert(reco_lhcbid);
        }
        // predict next phi in case two hits already exist
        float z = moduleToZ->GetBinContent(moduleNumber + 1);  // +1 because ROOT 1 indexed
        if (i_hit_reco == 0) {
          h1 = Hit::BaseHit(reco_lhcbid, x_reco->at(i_hit_reco), y_reco->at(i_hit_reco), z,
                   t_reco->at(i_hit_reco));
        } else {
          h0 = Hit::BaseHit(h1);
          h1 = Hit::BaseHit(reco_lhcbid, x_reco->at(i_hit_reco), y_reco->at(i_hit_reco), z,
                   t_reco->at(i_hit_reco));
        }
        if (i_hit_reco) {  // after i_hit_reco=1 we have two hits
          float dz = h1.z - h0.z;
          float dx = (h1.x - h0.x);
          float dy = (h1.y - h0.y);
          float dxdz = dx / dz;
          float dydz = dy / dz;
          // take difference to h0, that way we get longer interpolation (same as in real alg)
          float dz_to_next_module = moduleToZ->GetBinContent(moduleNumber + 2) - h0.z;
          float predx = dxdz * dz_to_next_module;
          float predy = dydz * dz_to_next_module;
          float x_prediction = h0.x + predx;
          float y_prediction = h0.y + predy;
          float track_extrapolation_phi = TMath::ATan2(y_prediction, x_prediction);
          mcTrackString += Form(", pred_phi=%.3f)", track_extrapolation_phi);
        } else
          mcTrackString += ")";  // close the phi, t paranthesis
        if (i_hit_reco != nRecoHits - 1) mcTrackString += ", ";
      }
      mcTrackString += "}\n";
    }
    // Determine what type of clone (directly copied from get_clone_rates.C)
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
    bool isAnySplitTrack = isSplitTrack || isSplitTrack_1Missed || isSplitTrack_2Missed;
    std::vector<unsigned> firstThreeModules(3, 0);
    bool firstThreeHitSameModules = true;
    for (auto r_lhcbid : firstHitLHCbIDs) {
      unsigned modID = (r_lhcbid >> 12) & 0x3F;
      if (firstThreeModules[0] == 0) firstThreeModules[0] = modID;
      else if (firstThreeModules[0] != modID) {
        firstThreeHitSameModules = false;
        break;
      }
    }
    // only check next two if previously already same module
    if (firstThreeHitSameModules) {
      for (auto r_lhcbid : secondHitLHCbIDs) {
        unsigned modID = (r_lhcbid >> 12) & 0x3F;
        if (firstThreeModules[1] == 0) firstThreeModules[1] = modID;
        else if (firstThreeModules[1] != modID) {
          firstThreeHitSameModules = false;
          break;
        }
      }
    }
    if (firstThreeHitSameModules) {
      for (auto r_lhcbid : thirdHitLHCbIDs) {
        unsigned modID = (r_lhcbid >> 12) & 0x3F;
        if (firstThreeModules[2] == 0) firstThreeModules[2] = modID;
        else if (firstThreeModules[2] != modID) {
          firstThreeHitSameModules = false;
          break;
        }
      }
    }
    bool isSeedingClone{false}, isSeedingClone2{false},
         isSeedingClone3{false}, isSeedingClone4{false};
    // set the right seeding clone. All will have x-Many-y hits & at least two matches
    // only look at those with at most 4-Many-4 setups
    bool canBeSeedingClone = firstHitLHCbIDs.size() < 5 && thirdHitLHCbIDs.size() < 5 &&
                             secondHitLHCbIDs.size() == nMatches && firstThreeHitSameModules &&
                             nMatches > 1 && !isAnySplitTrack;
    if (canBeSeedingClone) {
      if (firstHitLHCbIDs.size() == 4 || thirdHitLHCbIDs.size() == 4) {
        isSeedingClone4 = true;
      } else if (firstHitLHCbIDs.size() == 3 || thirdHitLHCbIDs.size() == 3) {
        isSeedingClone3 = true;
      } else if (firstHitLHCbIDs.size() == 2 || thirdHitLHCbIDs.size() == 2) {
        isSeedingClone2 = true;
        // this last one could just be an else tbf
      } else if (firstHitLHCbIDs.size() == 1 && thirdHitLHCbIDs.size() == 1) {
        isSeedingClone = true;
      }
    }
    // the shape of a seeding clone makes it ineligible for a triplet clone
    // Need at least four hits, otherwise every clone is a triplet clone
    bool isAnySeedingClone = isSeedingClone || isSeedingClone2 ||
                             isSeedingClone3 || isSeedingClone4;
    bool isCloneOfTriplets = nTriplets == nMatches && nMCVeloHits > 3 && nMatches > 1 &&
                             !isAnySeedingClone && !isAnySplitTrack;
    // Often the first one is a bit longer
    bool isCloneOfTripletsPlus = nTriplets == (nMatches - 1) && nMCVeloHits > 3 &&
                                 nMatches > 1 && !isAnySeedingClone && !isAnySplitTrack;
    bool isAnyTripletClone = isCloneOfTriplets || isCloneOfTripletsPlus;                            
    bool containsOddModules = false, containsEvenModules = false;
    // fill odd and even
    for (unsigned i_hit = 0; i_hit < mcTrackLHCbIDs->size(); i_hit++) {
      unsigned modID = (mcTrackLHCbIDs->at(i_hit) >> 12) & 0x3F;
      if (modID % 2 == 0) containsEvenModules = true;
      else containsOddModules = true;
    }
    // can only be tagged as a module overlap clone if not already assigned to something else
    bool isCloneOfModuleOverlap = (nMatches > 1) &&
                                  !(isAnySplitTrack || isAnySeedingClone || isAnyTripletClone) &&
                                  containsEvenModules && containsOddModules;
    bool isOtherClone = !(isAnySplitTrack || isAnySeedingClone ||
                          isAnyTripletClone || isCloneOfModuleOverlap)
                          && nMatches > 1;
    mcTrackString += "\tThis MC Particle has clone type: ";
    if (isAnySplitTrack) {
      mcTrackString +=
        Form("Split track clone");
    } else if (isAnySeedingClone) {
      mcTrackString +=
        Form("Seeding clone, Sizes: %lu, %lu, %lu",
             firstHitLHCbIDs.size(), secondHitLHCbIDs.size(), thirdHitLHCbIDs.size());
    } else if (isCloneOfTriplets) {
      mcTrackString += "Triplet clone";
    } else if (isCloneOfModuleOverlap) {
      mcTrackString += Form("Module overlap clone, Odd: %u, Even %u",
                            containsOddModules, containsEvenModules);
    } else if (isOtherClone) {
      mcTrackString +=
        Form("Other clone, Sizes: %lu, %lu, %lu",
             firstHitLHCbIDs.size(), secondHitLHCbIDs.size(), thirdHitLHCbIDs.size());
    }
    mcTrackString += "\n\n";

    // print the string (change the if here for what you want)
    if (isSplitTrack) {
      std::cout << mcTrackString;
    }
  }
  printf("Largest and smallest times throughout all events: (%.3f, %.3f)\n", t_min, t_max);
  printf("Number of reco tracks longer than MC tracks: %u\n", nRecoTrackLongerThanMC);

  return;
}