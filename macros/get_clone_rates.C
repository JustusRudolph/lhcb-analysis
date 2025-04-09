#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <string>

std::string stackRoot = std::getenv("STACK_ROOT");
std::string analysisRoot = std::getenv("ANALYSIS_ROOT");

void get_clone_rates(unsigned kEventsPerRun=200) {
  TFile *file = TFile::Open((stackRoot + "/MCData_Checking.root").c_str());
  if (!file || file->IsZombie()) {
      std::cerr << "Error: Could not open ROOT file." << std::endl;
      return;
  }
  // 1. Draw the duplication rate of IDs for each MC track by eta
  // first get eta bins
  int nBins = 50;
  std::vector<float> etaBinEdges(nBins+1);
  float etaMax{6.}, etaMin{-5.};
  float etaStep = (etaMax - etaMin) / (nBins);
  for (unsigned i = 0; i <= nBins; i++) {
    etaBinEdges[i] = etaMin + etaStep * i;
    // std::cout << etaBinEdges[i];
    // if (i != (nBins - 1 ))
    //   std::cout << ", ";
  }
  // std::cout << "\n";

  unsigned nReco{0}, nNonReco{0}, nMCTracksWithClones{0}, nMCTracksWith5Clones{0},
           nMCTracksWith10Clones{0};
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
  // define everything that you will need
  unsigned nMatches, mcTrackEvNo, mcTrackRunNo, recoEvOffset, nMCVeloHits;
  int mcMatchIdxForReco;
  float mcEta, recoEta;
  std::vector<unsigned>* matchedRecoTrackIndices = nullptr;
  std::vector<unsigned>* recoTrackLHCbIDs = nullptr;
  std::vector<unsigned>* mcTrackLHCbIDs = nullptr;
  // MC Event Tree
  // MC Track Tree
  mcTrackTree->SetBranchAddress("matchedTracks", &matchedRecoTrackIndices);
  mcTrackTree->SetBranchAddress("nMatches", &nMatches);
  mcTrackTree->SetBranchAddress("evNo", &mcTrackEvNo);
  mcTrackTree->SetBranchAddress("runNo", &mcTrackRunNo);
  mcTrackTree->SetBranchAddress("eta", &mcEta);
  mcTrackTree->SetBranchAddress("nHitsVelo", &nMCVeloHits);
  mcTrackTree->SetBranchAddress("lhcbid", &mcTrackLHCbIDs);
  // Reco Event Tree
  recoEventTree->SetBranchAddress("globalTrackOffset", &recoEvOffset);
  // Reco Track Tree
  recoTrackTree->SetBranchAddress("lhcbid", &recoTrackLHCbIDs);
  recoTrackTree->SetBranchAddress("mcMatchIdx", &mcMatchIdxForReco);
  recoTrackTree->SetBranchAddress("eta", &recoEta);

  unsigned nMCTracks = mcTrackTree->GetEntries();
  for (unsigned i_mct = 0; i_mct < nMCTracks; i_mct++) {
    mcTrackTree->GetEntry(i_mct);
    nReco += (nMatches > 0);
    nNonReco += (nMatches == 0);
    nMCTracksWithClones += (nMatches > 1);
    nMCTracksWith5Clones += (nMatches > 5);
    nMCTracksWith10Clones += (nMatches > 10);
  }
  unsigned nRecoWithoutClones = nReco - nMCTracksWithClones;

  TProfile* hDuplicateIDRates = new TProfile(
    "duplicate_match_id", "Duplicate match ID rate;#eta;DuplicateIDRate",
    nBins, etaBinEdges.data());
  TProfile* hUniqueIDRates = new TProfile(
    "unique_id_match_rate", "Unique ID match rate;#eta;Total ID reconstruction rate",
    nBins, etaBinEdges.data());
  TProfile* hLongestMatchedTrackRate = new TProfile(
    "longest_matched_track_fraction", "Longest matched track fraction;#eta;Longest Track Size (norm)",
    nBins, etaBinEdges.data());
  TProfile* hLongestMatchedTrackRateClones = new TProfile(
    "longest_matched_track_fraction_1_clones",
    "Longest matched track fraction;#eta;Longest Track Size (norm)",
    nBins, etaBinEdges.data());
  TProfile* hLongestMatchedTrackRate5Clones = new TProfile(
    "longest_matched_track_fraction_5_clones",
    "Longest matched track fraction;#eta;Longest Track Size (norm)",
    nBins, etaBinEdges.data());
  // following is whether the longest track was also tagged as the main match
  TH1D* hLongestTrackTagged = new TProfile(
    "longest_match_tag_rate", "Longest match tag rate;#eta;Longest Track Tag Rate",
    nBins, etaBinEdges.data());
  TProfile* hNMatches = new TProfile(
    "number_of_matches_of_MC_track", "Number of Matches of MC Track;#eta;N_{matches}",
    nBins, etaBinEdges.data());
  TProfile* ghostRates = new TProfile(
    "ghost_rates", "Ghost Rate wrt #eta (reconstructed tracks);#eta;Ghost Rate",
    nBins, etaBinEdges.data());
  // clone types by MC and reco
  // by MC describes the rate of an MC track having at least one clone,
  // and by reco describes the rate of it being a clone.
  TProfile* seedingClonesMCByEta = new TProfile(
    "clone_types_mc_by_eta", "Clone Types by #eta;#eta;Clone Rate",
    nBins, etaBinEdges.data());
  TProfile* tripletClonesMCByEta = new TProfile(
    "triplet_clones_mc_by_eta", "Triplet Clones by #eta;#eta;Clone Rate",
    nBins, etaBinEdges.data());
  TProfile* moduleOverlapClonesMCByEta = new TProfile(
    "module_overlap_clones_mc_by_eta", "Module Overlap Clones by #eta;#eta;Clone Rate",
    nBins, etaBinEdges.data());
  TProfile* otherClonesMCByEta = new TProfile(
    "other_clones_mc_by_eta", "Other Clones by #eta;#eta;Clone Rate",
    nBins, etaBinEdges.data());
  TProfile* seedingClonesRecoByEta = new TProfile(
    "clone_types_reco_by_eta", "Clone Types by #eta;#eta;Clone Rate",
    nBins, etaBinEdges.data());
  TProfile* tripletClonesRecoByEta = new TProfile(
    "triplet_clones_reco_by_eta", "Triplet Clones by #eta;#eta;Clone Rate",
    nBins, etaBinEdges.data());
  TProfile* moduleOverlapClonesRecoByEta = new TProfile(
    "module_overlap_clones_reco_by_eta", "Module Overlap Clones by #eta;#eta;Clone Rate",
    nBins, etaBinEdges.data());
  TProfile* otherClonesRecoByEta = new TProfile(
    "other_clones_reco_by_eta", "Other Clones by #eta;#eta;Clone Rate",
    nBins, etaBinEdges.data());

  printf("Going through MC particles now.\n");
  for (unsigned i_mct = 0; i_mct < nMCTracks; i_mct++) {
    mcTrackTree->GetEntry(i_mct);
    if (!nMatches || !matchedRecoTrackIndices) continue;
    std::unordered_set<unsigned> uniqueIDs;
    unsigned longestTrackSize = 0;
    unsigned kTotalIDs = 0;
    unsigned longestTrackTagged = 1;  // represent bool as int for histogram
    // define bools for different types of clones. Set to true initially if there are clones
    bool isCloneOfTriplets = nMatches > 1;

    // track the lhcbIDs of the second hits of all matched tracks to check for seeding clones
    std::unordered_set<unsigned> secondHitLHCbIDs;
    for (unsigned matchIdx : *matchedRecoTrackIndices) {
      unsigned evIdx = (mcTrackRunNo - 1) * kEventsPerRun + (mcTrackEvNo - 1);
      recoEventTree->GetEntry(evIdx);
      unsigned recoTrackIdx = recoEvOffset + matchIdx;
      recoTrackTree->GetEntry(recoTrackIdx);
      if (!recoTrackLHCbIDs) {
        std::cout << "Nothing found in reco track tree for lhcbid at index " << recoTrackIdx << ".\n";
        continue;  // if nothing there, this is weird...
      }
      isCloneOfTriplets = isCloneOfTriplets && (recoTrackLHCbIDs->size() == 3);
      kTotalIDs += recoTrackLHCbIDs->size();
      if (recoTrackLHCbIDs->size() > longestTrackSize) {
        // if longest track size not zero, that means it has already been changed
        // hence a later track is bigger, which shouldn't be
        if (longestTrackSize) longestTrackTagged = 0;
        longestTrackSize = recoTrackLHCbIDs->size();
      }
      for (unsigned i_hit_reco = 0; i_hit_reco < recoTrackLHCbIDs->size(); i_hit_reco++) {
        unsigned reco_lhcbid = recoTrackLHCbIDs->at(i_hit_reco);
        uniqueIDs.insert(reco_lhcbid);
        if (i_hit_reco == 1) {
          secondHitLHCbIDs.insert(reco_lhcbid);
        }
      }
    }
    bool isSeedingClone = nMatches > 1 && secondHitLHCbIDs.size() == nMatches;
    bool containsOddModules = false, containsEvenModules = false;
    // fill odd and even
    for (unsigned i_hit = 0; i_hit < mcTrackLHCbIDs->size(); i_hit++) {
      unsigned modID = (mcTrackLHCbIDs->at(i_hit) >> 12) & 0x3F;
      if (modID % 2 == 0) containsEvenModules = true;
      else containsOddModules = true;
    }
    // can only be tagged as a module overlap clone if not already assigned to something else
    bool isCloneOfModuleOverlap = (nMatches > 1) &&
                                  !(isSeedingClone || isCloneOfTriplets) &&
                                  (containsEvenModules != containsOddModules);

    float idDuplRate = 1. - (float) uniqueIDs.size() / (float) kTotalIDs;
    float idMatchRate = (float) uniqueIDs.size() / (float) nMCVeloHits;
    float longestTrackRate = (float) longestTrackSize / (float) nMCVeloHits;

    hDuplicateIDRates->Fill(mcEta, idDuplRate);
    hUniqueIDRates->Fill(mcEta, idMatchRate);

    if (nMatches > 5) {  // at least 5 clones
      hLongestMatchedTrackRate5Clones->Fill(mcEta, longestTrackRate);
    } else if (nMatches > 1) {  // at least one clone but less than 5
      hLongestMatchedTrackRateClones->Fill(mcEta, longestTrackRate);
    } else {  // will be exactly one match now because already checked not zero
      hLongestMatchedTrackRate->Fill(mcEta, longestTrackRate);
    }
    // fill the clone types only if at least one match, otherwise non-reconstructed
    if (nMatches) {
      bool isOtherClone = !(isSeedingClone || isCloneOfTriplets || isCloneOfModuleOverlap)
                          && nMatches > 1;
      // for MC just fill whether it has clones or not
      seedingClonesMCByEta->Fill(mcEta, isSeedingClone);
      tripletClonesMCByEta->Fill(mcEta, isCloneOfTriplets);
      moduleOverlapClonesMCByEta->Fill(mcEta, isCloneOfModuleOverlap);
      otherClonesMCByEta->Fill(mcEta, isOtherClone);
      // weight the contribution by the number of clones for reco clone plots
      seedingClonesRecoByEta->Fill(mcEta, 10 * isSeedingClone);
      tripletClonesRecoByEta->Fill(mcEta, 10 * isCloneOfTriplets);
      moduleOverlapClonesRecoByEta->Fill(mcEta, 10 * isCloneOfModuleOverlap);
      otherClonesRecoByEta->Fill(mcEta, 10 * isOtherClone);
    }

    hLongestTrackTagged->Fill(mcEta, longestTrackTagged);
    hNMatches->Fill(mcEta, nMatches);

    if (i_mct && (i_mct % 20000 == 0))
      printf("Finished with %d%% of MC Tracks.\n", (int) (100. * (float) i_mct / nMCTracks));
  }
  // now lastly go through all reco tracks and find the ghost rates
  unsigned nRecoTracks = recoTrackTree->GetEntries();
  for (unsigned i_recot = 0; i_recot < nRecoTracks; i_recot++) {
    recoTrackTree->GetEntry(i_recot);
    unsigned wasTagged = mcMatchIdxForReco != -1;
    ghostRates->Fill(recoEta, 1. - wasTagged);
  }
  // write histograms
  TFile* outFile = new TFile((analysisRoot + "/hists/clones/mc_hists.root").c_str(), "RECREATE");
  hDuplicateIDRates->Write();
  hUniqueIDRates->Write();
  hLongestMatchedTrackRate->Write();
  hLongestMatchedTrackRateClones->Write();
  hLongestMatchedTrackRate5Clones->Write();
  hLongestTrackTagged->Write();
  hNMatches->Write();
  ghostRates->Write();
  seedingClonesMCByEta->Write();
  tripletClonesMCByEta->Write();
  moduleOverlapClonesMCByEta->Write();
  otherClonesMCByEta->Write();
  seedingClonesRecoByEta->Write();
  tripletClonesRecoByEta->Write();
  moduleOverlapClonesRecoByEta->Write();
  otherClonesRecoByEta->Write();
  // Clean up outfile
  outFile->Close();
  delete outFile;

  // clean up file
  file->Close();
  delete file;
}