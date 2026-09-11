#include <TEfficiency.h>
#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TString.h>
#include <TSystem.h>
#include <bitset>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <string>

#include "utils/definitions.h"
#include "utils/basic_functions.h"

/*
 * Number of different modules the MC track left a hit in. The LHCb reconstructibility is
 * defined on separate sensors rather than on the raw hit count, so two hits in one module
 * must not count twice.
 */
unsigned count_distinct_modules(const std::vector<unsigned>* lhcbIDs) {
  if (!lhcbIDs) return 0;
  std::bitset<Utils::Definitions::kModules> modules;
  for (unsigned lhcbID : *lhcbIDs) modules.set((lhcbID >> 12) & 0x3F);
  return modules.count();
}

// max_dt is in picoseconds and scatter in micrometers
void get_eff_clone_ghosts(unsigned nEvents=5000, unsigned max_scatter=80000, unsigned max_dt=0,
                          TString mc_file_suffix="", unsigned kEventsPerRun=200) {
  TString input_suffix;
  if (mc_file_suffix.IsNull()) {
    TString suffix = Utils::Functions::get_suffix(nEvents, max_scatter, max_dt);
    input_suffix = suffix + ".root";
  } else {
    input_suffix = Form("_%uev_%s.root", nEvents, mc_file_suffix.Data());
  }
  TString input_prefix = (Utils::Definitions::stackRoot + "output/MCData_Checking").c_str();
  TString filepath = input_prefix + input_suffix;

  TFile *file = TFile::Open(filepath);
  if (!file || file->IsZombie()) {
    std::cerr << "Error: Could not open ROOT file." << std::endl;
    return;
  }
  // 1. Draw the duplication rate of IDs for each MC track by eta
  // first get eta bins
  int nEtaBins = 50;
  std::vector<float> etaBinEdges(nEtaBins+1);
  float etaMax{6.}, etaMin{-5.};
  float etaStep = (etaMax - etaMin) / (nEtaBins);
  for (unsigned i = 0; i <= nEtaBins; i++) {
    etaBinEdges[i] = etaMin + etaStep * i;
    // std::cout << etaBinEdges[i];
    // if (i != (nEtaBins - 1 ))
    //   std::cout << ", ";
  }
  // std::cout << "\n";
  // pt bins, same range as the efficiency plots use
  int nPTBins = 50;
  std::vector<float> ptBinEdges(nPTBins+1);
  float ptMax{5000.}, ptMin{0.};
  float ptStep = (ptMax - ptMin) / (nPTBins);
  for (unsigned i = 0; i <= nPTBins; i++) {
    ptBinEdges[i] = ptMin + ptStep * i;
  }
  // docaz bins for the efficiency
  int nDocaZBins = 50;
  std::vector<float> docaZBinEdges(nDocaZBins+1);
  float docaZMax{10.}, docaZMin{0.};
  float docaZStep = (docaZMax - docaZMin) / (nDocaZBins);
  for (unsigned i = 0; i <= nDocaZBins; i++) {
    docaZBinEdges[i] = docaZMin + docaZStep * i;
  }
  int nHitMax = 15, nHitMin = 3;
  int nHitBins = nHitMax - nHitMin + 1;
  float epsilon = 1e-6;
  std::vector<float> hitBinEdges(nHitBins+1);
  for (unsigned i = 0; i <= nHitBins; i++) {
    hitBinEdges[i] = 3. + i - epsilon;  // move low bin edge left to avoid floating point errors
  }

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
  bool isClone, hasVelo;
  unsigned nMatches, mcTrackEvNo, mcTrackRunNo, recoEvOffset, nMCVeloHits;
  int mcMatchIdxForReco, mcPID, runNoReco, evNoReco;
  float mcEta, recoEta, mcPT, mcDocaZ;
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
  mcTrackTree->SetBranchAddress("pt", &mcPT);
  mcTrackTree->SetBranchAddress("docaz", &mcDocaZ);
  mcTrackTree->SetBranchAddress("hasVelo", &hasVelo);
  mcTrackTree->SetBranchAddress("pid", &mcPID);
  mcTrackTree->SetBranchAddress("nHitsVelo", &nMCVeloHits);
  mcTrackTree->SetBranchAddress("lhcbid", &mcTrackLHCbIDs);
  // Reco Event Tree
  recoEventTree->SetBranchAddress("runNo", &runNoReco);
  recoEventTree->SetBranchAddress("evNo", &evNoReco);
  recoEventTree->SetBranchAddress("globalTrackOffset", &recoEvOffset);
  // Reco Track Tree
  recoTrackTree->SetBranchAddress("lhcbid", &recoTrackLHCbIDs);
  recoTrackTree->SetBranchAddress("mcMatchIdx", &mcMatchIdxForReco);
  recoTrackTree->SetBranchAddress("eta", &recoEta);
  recoTrackTree->SetBranchAddress("isClone", &isClone);


  // generic hit distribution for all MCs sliced up with same hit bins
  TH1D* hMCHitDistribution = new TH1D(
    "mc_hit_distribution", "MC Hit Distribution;N_{hits};Frequency",
    nHitBins, hitBinEdges.data());
  unsigned nHasVelo{0}, nHasAtLeastOneHit{0},
           nHasAtLeastTwoHits{0}, nHasAtLeastThreeHits{0};
  unsigned nMCTracks = mcTrackTree->GetEntries();
  for (unsigned i_mct = 0; i_mct < nMCTracks; i_mct++) {
    mcTrackTree->GetEntry(i_mct);
    nReco += (nMatches > 0);
    nNonReco += (nMatches == 0);
    nMCTracksWithClones += (nMatches > 1);
    nMCTracksWith5Clones += (nMatches > 5);
    nMCTracksWith10Clones += (nMatches > 10);
    nHasVelo += hasVelo;
    nHasAtLeastOneHit += (nMCVeloHits > 0);
    nHasAtLeastTwoHits += (nMCVeloHits > 1);
    nHasAtLeastThreeHits += (nMCVeloHits > 2);
    if (nMCVeloHits < nHitMax) {
      hMCHitDistribution->Fill(nMCVeloHits);
    } else {
      hMCHitDistribution->Fill(nHitMax);
    }
  }
  printf("Out of %u MC tracks, %u have Velo hits, %u, %u, %u have at least 1,2,3 hits.\n",
         nMCTracks, nHasVelo, nHasAtLeastOneHit, nHasAtLeastTwoHits, nHasAtLeastThreeHits);
  unsigned nRecoWithoutClones = nReco - nMCTracksWithClones;

  // create map of run number to run index
  unsigned nRuns = nEvents / kEventsPerRun;
  std::vector<unsigned> runNumberToRunIndex(nRuns, UINT_MAX); // Initialize with invalid index
  // loop over all entries, every 200 of which will be a new run
  for (unsigned i_run = 0; i_run < nRuns; i_run++) {
    recoEventTree->GetEntry(i_run * kEventsPerRun);
    // run numbers not zero indexed by default since Gauss0.sim "broken"
    runNumberToRunIndex[runNoReco - 1] = i_run;  // now zero indexed
  }

  TProfile* hDuplicateIDRates = new TProfile(
    "duplicate_match_id", "Duplicate match ID rate;#eta;DuplicateIDRate",
    nEtaBins, etaBinEdges.data());
  TProfile* hUniqueIDRates = new TProfile(
    "unique_id_match_rate", "Unique ID match rate;#eta;Total ID reconstruction rate",
    nEtaBins, etaBinEdges.data());
  TProfile* hLongestMatchedTrackRate = new TProfile(
    "longest_matched_track_fraction", "Longest matched track fraction;#eta;Longest Track Size (norm)",
    nEtaBins, etaBinEdges.data());
  TProfile* hLongestMatchedTrackRateClones = new TProfile(
    "longest_matched_track_fraction_1_clones",
    "Longest matched track fraction;#eta;Longest Track Size (norm)",
    nEtaBins, etaBinEdges.data());
  TProfile* hLongestMatchedTrackRate5Clones = new TProfile(
    "longest_matched_track_fraction_5_clones",
    "Longest matched track fraction;#eta;Longest Track Size (norm)",
    nEtaBins, etaBinEdges.data());
  // clone distributions
  TH1D* hCloneDistribution = new TH1D(
    "clone_distribution", "Distribution of clones;N_{clones};Frequency", 11, 0, 11);
  TH1D* hCloneDistributionEta35_5 = new TH1D(
    "clone_distribution_eta35_5", "Distribution of clones;N_{clones};Frequency", 11, 0, 11);
  TH1D* hSplitTrackCloneDistribution = new TH1D(
    "split_track_clone_distribution", "Distribution of clones;N_{clones};Frequency", 10, 1, 11);
  TH1D* hSplitTrackClone_1MissedDistribution = new TH1D(
    "split_track_clone_1missed_distribution", "Distribution of clones;N_{clones};Frequency", 10, 1, 11);
  TH1D* hSplitTrackClone_2MissedDistribution = new TH1D(
    "split_track_clone_2missed_distribution", "Distribution of clones;N_{clones};Frequency", 10, 1, 11);
  TH1D* hTripletCloneDistribution = new TH1D(
    "triplet_clone_distribution", "Distribution of clones;N_{clones};Frequency", 10, 1, 11);
  TH1D* hTripletClonePlusDistribution = new TH1D(
    "triplet_clone_plus_distribution", "Distribution of clones;N_{clones};Frequency", 10, 1, 11);
  TH1D* hSeedingCloneDistribution = new TH1D(
    "seeding_clone_distribution", "Distribution of clones;N_{clones};Frequency", 10, 1, 11);
  TH1D* hSeedingClone2Distribution = new TH1D(
    "seeding_clone2_distribution", "Distribution of clones;N_{clones};Frequency", 10, 1, 11);
  TH1D* hSeedingClone3Distribution = new TH1D(
    "seeding_clone3_distribution", "Distribution of clones;N_{clones};Frequency", 10, 1, 11);
  TH1D* hSeedingClone4Distribution = new TH1D(
    "seeding_clone4_distribution", "Distribution of clones;N_{clones};Frequency", 10, 1, 11);
  TH1D* hOverlapCloneDistribution = new TH1D(
    "overlap_clone_distribution", "Distribution of clones;N_{clones};Frequency", 10, 1, 11);
  TH1D* hOtherCloneDistribution = new TH1D(
    "other_clone_distribution", "Distribution of clones;N_{clones};Frequency", 10, 1, 11);
  Utils::Functions::printBinEdges(hSeedingCloneDistribution);
  TProfile* hNMatches = new TProfile(
    "number_of_matches_of_MC_track", "Number of Matches of MC Track;#eta;N_{matches}",
    nEtaBins, etaBinEdges.data());
  TProfile* ghostRates = new TProfile(
    "ghost_rates", "Ghost Rate wrt #eta (reconstructed tracks);#eta;Ghost Rate",
    nEtaBins, etaBinEdges.data());
  TProfile* cloneRate = new TProfile(
    "clone_rate", "Clone Rate;#eta;Clone Rate",
    nEtaBins, etaBinEdges.data());
  TProfile* cloneRatePT = new TProfile(
    "clone_rate_pt", "Clone Rate;p_{T} (MeV);Clone Rate",
    nPTBins, ptBinEdges.data());
  // Efficiency inputs: every reconstructible MC track goes into the denominator, the ones
  // with at least one match into the numerator. Per region and variable, in that order.
  const std::vector<std::string> kEffTypes = {"Pt", "Eta", "docaz"};
  const std::vector<std::string> kEffAxes = {"p_{T} (MeV)", "#eta", "DOCA_{z} (mm)"};
  std::vector<TH1D*> effReconstructible{}, effReconstructed{};
  for (unsigned region = 0; region < 2; region++) {
    std::string regionName = (region == 0 ? "forward" : "backward");
    for (unsigned t = 0; t < kEffTypes.size(); t++) {
      std::string base = "efficiency_" + regionName + "_" + kEffTypes[t];
      std::string title = "Efficiency vs " + kEffAxes[t] + ";" + kEffAxes[t] + ";Efficiency";
      int nBins = (t == 0 ? nPTBins : (t == 1 ? nEtaBins : nDocaZBins));
      const float* edges = (t == 0 ? ptBinEdges.data()
                                   : (t == 1 ? etaBinEdges.data() : docaZBinEdges.data()));
      effReconstructible.push_back(
        new TH1D((base + "_reconstructible").c_str(), title.c_str(), nBins, edges));
      effReconstructed.push_back(
        new TH1D((base + "_reconstructed").c_str(), title.c_str(), nBins, edges));
    }
  }
  // clone types by MC and reco
  // by MC describes the rate of an MC track having at least one clone,
  // and by reco describes the rate of it being a clone.
  TProfile* seedingClonesMCByEta = new TProfile(
    "seeding_clones_mc_by_eta", "Clone Types by #eta;#eta;Clone Rate",
    nEtaBins, etaBinEdges.data());
  TProfile* seedingClones2MCByEta = new TProfile(
    "seeding_clones2_mc_by_eta", "Clone Types by #eta;#eta;Clone Rate",
    nEtaBins, etaBinEdges.data());
  TProfile* seedingClones3MCByEta = new TProfile(
    "seeding_clones3_mc_by_eta", "Clone Types by #eta;#eta;Clone Rate",
    nEtaBins, etaBinEdges.data());
  TProfile* seedingClones4MCByEta = new TProfile(
    "seeding_clones4_mc_by_eta", "Clone Types by #eta;#eta;Clone Rate",
    nEtaBins, etaBinEdges.data());
  TProfile* splitTrackClonesMCByEta = new TProfile(
    "split_track_clones_mc_by_eta", "Split Track Clones by #eta;#eta;Clone Rate",
    nEtaBins, etaBinEdges.data());
  TProfile* splitTrackClones_1MissedMCByEta = new TProfile(
    "split_track_clones_1missed_mc_by_eta", "Split Track Clones by #eta;#eta;Clone Rate",
    nEtaBins, etaBinEdges.data());
  TProfile* splitTrackClones_2MissedMCByEta = new TProfile(
    "split_track_clones_2missed_mc_by_eta", "Split Track Clones by #eta;#eta;Clone Rate",
    nEtaBins, etaBinEdges.data());
  TProfile* tripletClonesMCByEta = new TProfile(
    "triplet_clones_mc_by_eta", "Triplet Clones by #eta;#eta;Clone Rate",
    nEtaBins, etaBinEdges.data());
  TProfile* tripletClonesPlusMCByEta = new TProfile(
    "triplet_clones_plus_mc_by_eta", "Triplet Clones by #eta;#eta;Clone Rate",
    nEtaBins, etaBinEdges.data());
  TProfile* moduleOverlapClonesMCByEta = new TProfile(
    "module_overlap_clones_mc_by_eta", "Module Overlap Clones by #eta;#eta;Clone Rate",
    nEtaBins, etaBinEdges.data());
  TProfile* otherClonesMCByEta = new TProfile(
    "other_clones_mc_by_eta", "Other Clones by #eta;#eta;Clone Rate",
    nEtaBins, etaBinEdges.data());
  TProfile* seedingClonesRecoByEta = new TProfile(
    "seeding_clones_reco_by_eta", "Clone Types by #eta;#eta;Clone Rate",
    nEtaBins, etaBinEdges.data());
  TProfile* seedingClones2RecoByEta = new TProfile(
    "seeding_clones2_reco_by_eta", "Clone Types by #eta;#eta;Clone Rate",
    nEtaBins, etaBinEdges.data());
  TProfile* seedingClones3RecoByEta = new TProfile(
    "seeding_clones3_reco_by_eta", "Clone Types by #eta;#eta;Clone Rate",
    nEtaBins, etaBinEdges.data());
  TProfile* seedingClones4RecoByEta = new TProfile(
    "seeding_clones4_reco_by_eta", "Clone Types by #eta;#eta;Clone Rate",
    nEtaBins, etaBinEdges.data());
  TProfile* splitTrackClonesRecoByEta = new TProfile(
    "split_track_clones_reco_by_eta", "Split Track Clones by #eta;#eta;Clone Rate",
    nEtaBins, etaBinEdges.data());
  TProfile* splitTrackClones_1MissedRecoByEta = new TProfile(
    "split_track_clones_1missed_reco_by_eta", "Split Track Clones by #eta;#eta;Clone Rate",
    nEtaBins, etaBinEdges.data());
  TProfile* splitTrackClones_2MissedRecoByEta = new TProfile(
    "split_track_clones_2missed_reco_by_eta", "Split Track Clones by #eta;#eta;Clone Rate",
    nEtaBins, etaBinEdges.data());
  TProfile* tripletClonesRecoByEta = new TProfile(
    "triplet_clones_reco_by_eta", "Triplet Clones by #eta;#eta;Clone Rate",
    nEtaBins, etaBinEdges.data());
  TProfile* tripletClonesPlusRecoByEta = new TProfile(
    "triplet_clones_plus_reco_by_eta", "Triplet Clones by #eta;#eta;Clone Rate",
    nEtaBins, etaBinEdges.data());
  TProfile* moduleOverlapClonesRecoByEta = new TProfile(
    "module_overlap_clones_reco_by_eta", "Module Overlap Clones by #eta;#eta;Clone Rate",
    nEtaBins, etaBinEdges.data());
  TProfile* otherClonesRecoByEta = new TProfile(
    "other_clones_reco_by_eta", "Other Clones by #eta;#eta;Clone Rate",
    nEtaBins, etaBinEdges.data());
  // clone rates and absolute numbers of clones by the size/length of MC track
  TH1D* nClonesByMCPSize_forward = new TH1D(
    "n_clones_by_mc_p_size_forward", "Number of Clones by number of MCP hits;N_{hits};Number of Clones",
    nHitBins, hitBinEdges.data());
  TH1D* nClonesByMCPSize_backward = new TH1D(
    "n_clones_by_mc_p_size_backward", "Number of Clones by number of MCP hits;N_{hits};Number of Clones",
    nHitBins, hitBinEdges.data());
  TProfile* cloneRateByMCPSize_forward = new TProfile(
    "clone_rate_by_mc_p_size_forward", "Clone Rate by number of MCP hits;N_{hits};Clone Rate",
    nHitBins, hitBinEdges.data());
  TProfile* cloneRateByMCPSize_backward = new TProfile(
    "clone_rate_by_mc_p_size_backward", "Clone Rate by number of MCP hits;N_{hits};Clone Rate",
    nHitBins, hitBinEdges.data());
  // following four are the same as above but scaled down by number of clones
  // i.e. it does not matter how many clones there are, it's either clone or no clone
  TH1D* nClonesByMCPSize_forward_scaled = new TH1D(
    "n_clones_by_mc_p_size_forward_scaled", "Number of Clones by number of MCP hits;N_{hits};Number of Clones",
    nHitBins, hitBinEdges.data());
  TH1D* nClonesByMCPSize_backward_scaled = new TH1D(
    "n_clones_by_mc_p_size_backward_scaled", "Number of Clones by number of MCP hits;N_{hits};Number of Clones",
    nHitBins, hitBinEdges.data());
  TProfile* cloneRateByMCPSize_forward_scaled = new TProfile(
    "clone_rate_by_mc_p_size_forward_scaled", "Clone Rate by number of MCP hits;N_{hits};Clone Rate",
    nHitBins, hitBinEdges.data());
  TProfile* cloneRateByMCPSize_backward_scaled = new TProfile(
    "clone_rate_by_mc_p_size_backward_scaled", "Clone Rate by number of MCP hits;N_{hits};Clone Rate",
    nHitBins, hitBinEdges.data());
  // also do Histograms for the number of hits per clone type
  TH1D* seedingCloneMCHitDistr = new TH1D(
    "seeding_clone_mc_hit_distribution", "Seeding Clone Hit Distribution;N_{hits};Frequency",
    nHitBins, hitBinEdges.data());
  TH1D* seedingClone2MCHitDistr = new TH1D(
    "seeding_clone2_mc_hit_distribution", "Seeding Clone Hit Distribution;N_{hits};Frequency",
    nHitBins, hitBinEdges.data());
  TH1D* seedingClone3MCHitDistr = new TH1D(
    "seeding_clone3_mc_hit_distribution", "Seeding Clone Hit Distribution;N_{hits};Frequency",
    nHitBins, hitBinEdges.data());
  TH1D* seedingClone4MCHitDistr = new TH1D(
    "seeding_clone4_mc_hit_distribution", "Seeding Clone Hit Distribution;N_{hits};Frequency",
    nHitBins, hitBinEdges.data());
  TH1D* splitTrackCloneMCHitDistr = new TH1D(
    "split_track_clone_mc_hit_distribution", "Split Track Clone Hit Distribution;N_{hits};Frequency",
    nHitBins, hitBinEdges.data());
  TH1D* splitTrackClone_1MissedMCHitDistr = new TH1D(
    "split_track_clone_1missed_mc_hit_distribution", "Split Track Clone Hit Distribution;N_{hits};Frequency",
    nHitBins, hitBinEdges.data());
  TH1D* splitTrackClone_2MissedMCHitDistr = new TH1D(
    "split_track_clone_2missed_mc_hit_distribution", "Split Track Clone Hit Distribution;N_{hits};Frequency",
    nHitBins, hitBinEdges.data());
  TH1D* tripletCloneMCHitDistr = new TH1D(
    "triplet_clone_mc_hit_distribution", "Triplet Clone Hit Distribution;N_{hits};Frequency",
    nHitBins, hitBinEdges.data());
  TH1D* tripletClonePlusMCHitDistr = new TH1D(
    "triplet_clone_plus_mc_hit_distribution", "Triplet Clone Hit Distribution;N_{hits};Frequency",
    nHitBins, hitBinEdges.data());
  TH1D* moduleOverlapCloneMCHitDistr = new TH1D(
    "module_overlap_clone_mc_hit_distribution", "Module Overlap Clone Hit Distribution;N_{hits};Frequency",
    nHitBins, hitBinEdges.data());
  TH1D* otherCloneMCHitDistr = new TH1D(
    "other_clone_mc_hit_distribution", "Other Clone Hit Distribution;N_{hits};Frequency",
    nHitBins, hitBinEdges.data());
  TH1D* seedingCloneRecoHitDistr = new TH1D(
    "seeding_clone_reco_hit_distribution", "Seeding Clone Hit Distribution;N_{hits};Frequency",
    nHitBins, hitBinEdges.data());
  TH1D* seedingClone2RecoHitDistr = new TH1D(
    "seeding_clone2_reco_hit_distribution", "Seeding Clone Hit Distribution;N_{hits};Frequency",
    nHitBins, hitBinEdges.data());
  TH1D* seedingClone3RecoHitDistr = new TH1D(
    "seeding_clone3_reco_hit_distribution", "Seeding Clone Hit Distribution;N_{hits};Frequency",
    nHitBins, hitBinEdges.data());
  TH1D* seedingClone4RecoHitDistr = new TH1D(
    "seeding_clone4_reco_hit_distribution", "Seeding Clone Hit Distribution;N_{hits};Frequency",
    nHitBins, hitBinEdges.data());
  TH1D* splitTrackCloneRecoHitDistr = new TH1D(
    "split_track_clone_reco_hit_distribution", "Split Track Clone Hit Distribution;N_{hits};Frequency",
    nHitBins, hitBinEdges.data());
  TH1D* splitTrackClone_1MissedRecoHitDistr = new TH1D(
    "split_track_clone_1missed_reco_hit_distribution", "Split Track Clone Hit Distribution;N_{hits};Frequency",
    nHitBins, hitBinEdges.data());
  TH1D* splitTrackClone_2MissedRecoHitDistr = new TH1D(
    "split_track_clone_2missed_reco_hit_distribution", "Split Track Clone Hit Distribution;N_{hits};Frequency",
    nHitBins, hitBinEdges.data());
  TH1D* tripletCloneRecoHitDistr = new TH1D(
    "triplet_clone_reco_hit_distribution", "Triplet Clone Hit Distribution;N_{hits};Frequency",
    nHitBins, hitBinEdges.data());
  TH1D* tripletClonePlusRecoHitDistr = new TH1D(
    "triplet_clone_plus_reco_hit_distribution", "Triplet Clone Hit Distribution;N_{hits};Frequency",
    nHitBins, hitBinEdges.data());
  TH1D* moduleOverlapCloneRecoHitDistr = new TH1D(
    "module_overlap_clone_reco_hit_distribution", "Module Overlap Clone Hit Distribution;N_{hits};Frequency",
    nHitBins, hitBinEdges.data());
  TH1D* otherCloneRecoHitDistr = new TH1D(
    "other_clone_reco_hit_distribution", "Other Clone Hit Distribution;N_{hits};Frequency",
    nHitBins, hitBinEdges.data());

  printf("Going through %u MC particles now.\n", nMCTracks);
  unsigned nSeedingClones = 0, nTripletClones = 0, nSplitAndOverlapClones = 0,
           nSplitClones = 0, nSplit1MissedClones = 0, nSplit2MissedClones = 0,
           nAnySplitClones = 0, nSplitTrack2ndTrackFirst = 0;
  unsigned nTotalForward = 0, nTotalBackward = 0,  // total number of MC tracks
           nTotalRecodForward = 0, nTotalRecodBackward = 0,
           nTotalClonesForward = 0, nTotalClonesBackward = 0;

  unsigned nInteresting = 0;
  for (unsigned i_mct = 0; i_mct < nMCTracks; i_mct++) {
    mcTrackTree->GetEntry(i_mct);
    int mcPID_abs = mcPID < 0 ? -mcPID : mcPID;
    bool isInAcceptance = (mcEta >= -5 && mcEta <= -2) ||
                          (mcEta >= 2 && mcEta <= 5);
    // hits in at least 3 different modules and not an electron, as in LHCb's MCTrackInfo
    unsigned nDistinctModules = count_distinct_modules(mcTrackLHCbIDs);
    bool isReconstructible = isInAcceptance && (nDistinctModules >= 3) && (mcPID_abs != 11);
    if (!isReconstructible) continue;  // skip if not reconstructible
    // reconstructed means the MC track was matched to at least one reco track
    unsigned effRegion = (mcEta < 0 ? 1 : 0);
    std::vector<float> effValues = {mcPT, mcEta, mcDocaZ};
    for (unsigned t = 0; t < kEffTypes.size(); t++) {
      unsigned effIdx = effRegion * kEffTypes.size() + t;
      effReconstructible[effIdx]->Fill(effValues[t]);
      if (nMatches) effReconstructed[effIdx]->Fill(effValues[t]);
    }
    // Fill eff and clone numbers
    if (mcEta < 0) {
      nTotalBackward++;
      if (nMatches) {  // only add to reconstructed or clone if matched to mc
        nTotalRecodBackward++;
        nTotalClonesBackward += (nMatches - 1);
        cloneRateByMCPSize_backward->Fill(nMCVeloHits, nMatches - 1);
        cloneRateByMCPSize_backward_scaled->Fill(nMCVeloHits, nMatches > 1);
        nClonesByMCPSize_backward->Fill(nMCVeloHits, nMatches - 1);
        nClonesByMCPSize_backward_scaled->Fill(nMCVeloHits, nMatches > 1);
      }
    } else {
      nTotalForward++;
      if (nMatches) {
        nTotalRecodForward += (nMatches > 0);
        nTotalClonesForward += (nMatches - 1);
        cloneRateByMCPSize_forward->Fill(nMCVeloHits, nMatches - 1);
        cloneRateByMCPSize_forward_scaled->Fill(nMCVeloHits, nMatches > 1);
        nClonesByMCPSize_forward->Fill(nMCVeloHits, nMatches - 1);
        nClonesByMCPSize_forward_scaled->Fill(nMCVeloHits, nMatches > 1);
      }
    }
    if (!nMatches || !matchedRecoTrackIndices) continue;
    std::unordered_set<unsigned> uniqueIDs;
    unsigned longestTrackSize = 0;
    unsigned kTotalIDs = 0;
    unsigned longestTrackTagged = 1;  // represent bool as int for histogram
    unsigned nTriplets = 0;
    unsigned lastHitModuleNumber = 64;  // out of range
    bool tracksRightOrder = false;
    // To check if MC track to be classified as split
    std::unordered_set<unsigned> mcTrackLHCbIDsSet(mcTrackLHCbIDs->begin(),
                                                   mcTrackLHCbIDs->end());

    // track the lhcbIDs of the first three hits of all matched tracks to check for seeding clones
    std::unordered_set<unsigned> firstHitLHCbIDs;
    std::unordered_set<unsigned> secondHitLHCbIDs;
    std::unordered_set<unsigned> thirdHitLHCbIDs;
    for (unsigned i_rt = 0; i_rt < nMatches; i_rt++) {
      unsigned matchIdx = matchedRecoTrackIndices->at(i_rt);
      unsigned evIdx = runNumberToRunIndex[mcTrackRunNo - 1] * kEventsPerRun + (mcTrackEvNo - 1);
      recoEventTree->GetEntry(evIdx);
      unsigned recoTrackIdx = recoEvOffset + matchIdx;
      recoTrackTree->GetEntry(recoTrackIdx);
      if (!recoTrackLHCbIDs) {
        std::cout << "Nothing found in reco track tree for lhcbid at index " << recoTrackIdx << ".\n";
        continue;  // if nothing there, this is weird...
      }
      unsigned nRecoHits = recoTrackLHCbIDs->size();
      nTriplets += (nRecoHits == 3);
      kTotalIDs += nRecoHits;
      if (nRecoHits > longestTrackSize) {
        // if longest track size not zero, that means it has already been changed
        // hence a later track is bigger, which shouldn't be
        if (longestTrackSize) longestTrackTagged = 0;
        longestTrackSize = nRecoHits;
      }
      for (unsigned i_hit_reco = 0; i_hit_reco < nRecoHits; i_hit_reco++) {
        unsigned reco_lhcbid = recoTrackLHCbIDs->at(i_hit_reco);
        mcTrackLHCbIDsSet.erase(reco_lhcbid);
        uniqueIDs.insert(reco_lhcbid);
        if (i_hit_reco == 0) {
          firstHitLHCbIDs.insert(reco_lhcbid);
        } else if (i_hit_reco == 1) {
          secondHitLHCbIDs.insert(reco_lhcbid);
        } else if (i_hit_reco == 2) {
          thirdHitLHCbIDs.insert(reco_lhcbid);
        }
        unsigned moduleNumber = (reco_lhcbid >> 12) & 0x3F;
        if (i_rt == 0 && i_hit_reco == (nRecoHits - 1)) {
          lastHitModuleNumber = moduleNumber;
        } else if (i_rt == 1 && i_hit_reco == 0) {
          // check order of tracks, go from 63->0, so see if next track starts smaller
          // than previous one ended
          tracksRightOrder = (moduleNumber < lastHitModuleNumber);
        }
      }
    }
    // -------------- CLONES --------------
    // all MC hits reconstructed within all reco matched tracks
    bool isSplitTrack = mcTrackLHCbIDsSet.size() == 0 && kTotalIDs == nMCVeloHits &&
                        nMatches > 1;
    // Either one hit missing between, or one overlap/extra hit
    bool isSplitTrack_1Missed = ( (mcTrackLHCbIDsSet.size() == 1 && kTotalIDs == (nMCVeloHits - 1)) ||
                                  (mcTrackLHCbIDsSet.size() == 0 && kTotalIDs == (nMCVeloHits + 1)) ) &&
                                   nMatches > 1;
    // Two hits missing, two overlap (or two extra), or one missing and one extra
    bool isSplitTrack_2Missed = ( (mcTrackLHCbIDsSet.size() == 2 && kTotalIDs == (nMCVeloHits - 2)) ||
                                  (mcTrackLHCbIDsSet.size() == 0 && kTotalIDs == (nMCVeloHits + 2)) ||
                                  (mcTrackLHCbIDsSet.size() == 1 && kTotalIDs == (nMCVeloHits + 1)) ) &&
                                   nMatches > 1;
    nInteresting += (isSplitTrack_2Missed && nMatches == 2 && nMCVeloHits == 4);
    bool isAnySplitTrack = isSplitTrack || isSplitTrack_1Missed || isSplitTrack_2Missed;
    nSplitTrack2ndTrackFirst += tracksRightOrder && isSplitTrack;
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
    nSplitClones += isSplitTrack;
    nSplit1MissedClones += isSplitTrack_1Missed;
    nSplit2MissedClones += isSplitTrack_2Missed;
    nAnySplitClones += isAnySplitTrack;
    nSplitAndOverlapClones += (isAnySplitTrack && containsEvenModules && containsOddModules);
    // can only be tagged as a module overlap clone if not already assigned to something else
    bool isCloneOfModuleOverlap = (nMatches > 1) &&
                                  !(isAnySplitTrack || isAnySeedingClone || isAnyTripletClone) &&
                                  containsEvenModules && containsOddModules;

    float meanNumberOfRecoHits = (float) kTotalIDs / (float) nMatches;
    // fill the clone types and distribution only if at least one match
    if (nMatches && Utils::Functions::inEtaAcceptance(mcEta)) {
      bool isOtherClone = !(isAnySplitTrack || isAnySeedingClone ||
                            isAnyTripletClone || isCloneOfModuleOverlap)
                          && nMatches > 1;
      // fill clone distributions
      hCloneDistribution->Fill(nMatches - 1);
      if (mcEta > 3.5 && mcEta < 5) hCloneDistributionEta35_5->Fill(nMatches - 1);
      // and those for specific types of clones
      // note that max number of clones is 10 (10+) for the distribution
      unsigned nClones = (nMatches > 10) ? 10 : nMatches - 1;
      if (isSplitTrack) {
        hSplitTrackCloneDistribution->Fill(nClones);
        // fill the number of hits for mc and reco
        if (nMCVeloHits < nHitMax) {
          splitTrackCloneMCHitDistr->Fill(nMCVeloHits);
        } else {
          splitTrackCloneMCHitDistr->Fill(nHitMax);
        }
        if (meanNumberOfRecoHits < nHitMax) {
          splitTrackCloneRecoHitDistr->Fill(meanNumberOfRecoHits);
        } else {
          splitTrackCloneRecoHitDistr->Fill(nHitMax);
        }
      } else if (isSplitTrack_1Missed) {
        hSplitTrackClone_1MissedDistribution->Fill(nClones);
        // fill the number of hits for mc and reco
        if (nMCVeloHits < nHitMax) {
          splitTrackClone_1MissedMCHitDistr->Fill(nMCVeloHits);
        } else {
          splitTrackClone_1MissedMCHitDistr->Fill(nHitMax);
        }
        if (meanNumberOfRecoHits < nHitMax) {
          splitTrackClone_1MissedRecoHitDistr->Fill(meanNumberOfRecoHits);
        } else {
          splitTrackClone_1MissedRecoHitDistr->Fill(nHitMax);
        }
      } else if (isSplitTrack_2Missed) {
        hSplitTrackClone_2MissedDistribution->Fill(nClones);
        // fill the number of hits for mc and reco
        if (nMCVeloHits < nHitMax) {
          splitTrackClone_2MissedMCHitDistr->Fill(nMCVeloHits);
        } else {
          splitTrackClone_2MissedMCHitDistr->Fill(nHitMax);
        }
        if (meanNumberOfRecoHits < nHitMax) {
          splitTrackClone_2MissedRecoHitDistr->Fill(meanNumberOfRecoHits);
        } else {
          splitTrackClone_2MissedRecoHitDistr->Fill(nHitMax);
        }
      } else if (isSeedingClone) {
        nSeedingClones++;
        hSeedingCloneDistribution->Fill(nClones);
        // fill the number of hits for mc and reco
        if (nMCVeloHits < nHitMax) {
          seedingCloneMCHitDistr->Fill(nMCVeloHits);
        } else {
          seedingCloneMCHitDistr->Fill(nHitMax);
        }
        if (meanNumberOfRecoHits < nHitMax) {
          seedingCloneRecoHitDistr->Fill(meanNumberOfRecoHits);
        } else {
          seedingCloneRecoHitDistr->Fill(nHitMax);
        }
      } else if (isSeedingClone2) {
        hSeedingClone2Distribution->Fill(nClones);
        // fill the number of hits for mc and reco
        if (nMCVeloHits < nHitMax) {
          seedingClone2MCHitDistr->Fill(nMCVeloHits);
        } else {
          seedingClone2MCHitDistr->Fill(nHitMax);
        }
        if (meanNumberOfRecoHits < nHitMax) {
          seedingClone2RecoHitDistr->Fill(meanNumberOfRecoHits);
        } else {
          seedingClone2RecoHitDistr->Fill(nHitMax);
        }
      } else if (isSeedingClone3) {
        hSeedingClone3Distribution->Fill(nClones);
        // fill the number of hits for mc and reco
        if (nMCVeloHits < nHitMax) {
          seedingClone3MCHitDistr->Fill(nMCVeloHits);
        } else {
          seedingClone3MCHitDistr->Fill(nHitMax);
        }
        if (meanNumberOfRecoHits < nHitMax) {
          seedingClone3RecoHitDistr->Fill(meanNumberOfRecoHits);
        } else {
          seedingClone3RecoHitDistr->Fill(nHitMax);
        }
      } else if (isSeedingClone4) {
        hSeedingClone4Distribution->Fill(nClones);
        // fill the number of hits for mc and reco
        if (nMCVeloHits < nHitMax) {
          seedingClone4MCHitDistr->Fill(nMCVeloHits);
        } else {
          seedingClone4MCHitDistr->Fill(nHitMax);
        }
        if (meanNumberOfRecoHits < nHitMax) {
          seedingClone4RecoHitDistr->Fill(meanNumberOfRecoHits);
        } else {
          seedingClone4RecoHitDistr->Fill(nHitMax);
        }
      } else if (isCloneOfTriplets) {
        nTripletClones++;
        hTripletCloneDistribution->Fill(nClones);
        // fill the number of hits for mc and reco
        if (nMCVeloHits < nHitMax) {
          tripletCloneMCHitDistr->Fill(nMCVeloHits);
        } else {
          tripletCloneMCHitDistr->Fill(nHitMax);
        }
        if (meanNumberOfRecoHits < nHitMax) {
          tripletCloneRecoHitDistr->Fill(meanNumberOfRecoHits);
        } else {
          tripletCloneRecoHitDistr->Fill(nHitMax);
        }
      } else if (isCloneOfTripletsPlus) {
        hTripletClonePlusDistribution->Fill(nClones);
        // fill the number of hits for mc and reco
        if (nMCVeloHits < nHitMax) {
          tripletClonePlusMCHitDistr->Fill(nMCVeloHits);
        } else {
          tripletClonePlusMCHitDistr->Fill(nHitMax);
        }
        if (meanNumberOfRecoHits < nHitMax) {
          tripletClonePlusRecoHitDistr->Fill(meanNumberOfRecoHits);
        } else {
          tripletClonePlusRecoHitDistr->Fill(nHitMax);
        }
      } else if (isCloneOfModuleOverlap) {
        hOverlapCloneDistribution->Fill(nClones);
        // fill the number of hits for mc and reco
        if (nMCVeloHits < nHitMax) {
          moduleOverlapCloneMCHitDistr->Fill(nMCVeloHits);
        } else {
          moduleOverlapCloneMCHitDistr->Fill(nHitMax);
        }
        if (meanNumberOfRecoHits < nHitMax) {
          moduleOverlapCloneRecoHitDistr->Fill(meanNumberOfRecoHits);
        } else {
          moduleOverlapCloneRecoHitDistr->Fill(nHitMax);
        }
      } else if (isOtherClone) {
        hOtherCloneDistribution->Fill(nClones);
        // fill the number of hits for mc and reco
        if (nMCVeloHits < nHitMax) {
          otherCloneMCHitDistr->Fill(nMCVeloHits);
        } else {
          otherCloneMCHitDistr->Fill(nHitMax);
        }
        if (meanNumberOfRecoHits < nHitMax) {
          otherCloneRecoHitDistr->Fill(meanNumberOfRecoHits);
        } else {
          otherCloneRecoHitDistr->Fill(nHitMax);
        }
      }
      // for MC just fill whether it has clones or not
      splitTrackClonesMCByEta->Fill(mcEta, isSplitTrack);
      splitTrackClones_1MissedMCByEta->Fill(mcEta, isSplitTrack_1Missed);
      splitTrackClones_2MissedMCByEta->Fill(mcEta, isSplitTrack_2Missed);
      seedingClonesMCByEta->Fill(mcEta, isSeedingClone);
      seedingClones2MCByEta->Fill(mcEta, isSeedingClone2);
      seedingClones3MCByEta->Fill(mcEta, isSeedingClone3);
      seedingClones4MCByEta->Fill(mcEta, isSeedingClone4);
      tripletClonesMCByEta->Fill(mcEta, isCloneOfTriplets);
      tripletClonesPlusMCByEta->Fill(mcEta, isCloneOfTripletsPlus);
      moduleOverlapClonesMCByEta->Fill(mcEta, isCloneOfModuleOverlap);
      otherClonesMCByEta->Fill(mcEta, isOtherClone);
      // weight the contribution by the number of clones for reco clone plots
      splitTrackClonesRecoByEta->Fill(mcEta, (nMatches - 1) * isSplitTrack);
      splitTrackClones_1MissedRecoByEta->Fill(mcEta, (nMatches - 1) * isSplitTrack_1Missed);
      splitTrackClones_2MissedRecoByEta->Fill(mcEta, (nMatches - 1) * isSplitTrack_2Missed);
      seedingClonesRecoByEta->Fill(mcEta, (nMatches - 1) * isSeedingClone);
      seedingClones2RecoByEta->Fill(mcEta, (nMatches - 1) * isSeedingClone2);
      seedingClones3RecoByEta->Fill(mcEta, (nMatches - 1) * isSeedingClone3);
      seedingClones4RecoByEta->Fill(mcEta, (nMatches - 1) * isSeedingClone4);
      tripletClonesRecoByEta->Fill(mcEta, (nMatches - 1) * isCloneOfTriplets);
      tripletClonesPlusRecoByEta->Fill(mcEta, (nMatches - 1) * isCloneOfTripletsPlus);
      moduleOverlapClonesRecoByEta->Fill(mcEta, (nMatches - 1) * isCloneOfModuleOverlap);
      otherClonesRecoByEta->Fill(mcEta, (nMatches - 1) * isOtherClone);
      // have an if (!nMatches) guard already earlier
      cloneRate->Fill(mcEta, nMatches - 1);
      cloneRatePT->Fill(mcPT, nMatches - 1);
    }  // if nMatches && inEtaAcceptance(mcEta)

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
    hNMatches->Fill(mcEta, nMatches);

    if (i_mct && (i_mct % 100000 == 0)) {
      printf("Finished with %d%% of MC Tracks.\n", (int) (100. * (float) i_mct / nMCTracks));
    }
  }  // MC Particles loop
  printf("Finished with all MC tracks.\n");
  printf("Number of interesting (right now 2 clones, 4 hits, NNLO split) tracks: %u\n", nInteresting);
  printf("With %u seeding clones, %u are also triplet clones.\n",
         nSeedingClones, nTripletClones);
  printf("With %u out of %u (%u & %u) split (1 & 2 missed) clones having module overlap.\n",
         nSplitAndOverlapClones, nSplitClones, nSplit1MissedClones, nSplit2MissedClones);
  printf("Sum of different split clones: %u, while total (ORed) split clones: %u.\n",
         nSplitClones + nSplit1MissedClones + nSplit2MissedClones, nAnySplitClones);
  printf("Sizes of split clone clone histograms: %f, %f, %f.\n",
         hSplitTrackCloneDistribution->GetEntries(),
         hSplitTrackClone_1MissedDistribution->GetEntries(),
         hSplitTrackClone_2MissedDistribution->GetEntries());
  printf("Out of %u split track clones, %u have the second track first.\n",
          nSplitClones, nSplitTrack2ndTrackFirst);
  // now lastly go through all reco tracks and find the ghost rates
  unsigned nRecoTracks = recoTrackTree->GetEntries();
  printf("Going through %u reco tracks now.\n", nRecoTracks);
  unsigned nTotalGhostsForward = 0, nTotalGhostsBackward = 0,
           nTotalRecoTracksForward = 0, nTotalRecoTracksBackward = 0,
           nTotalRecoClonesForward = 0, nTotalRecoClonesBackward = 0;
  for (unsigned i_recot = 0; i_recot < nRecoTracks; i_recot++) {
    recoTrackTree->GetEntry(i_recot);
    unsigned wasTagged = mcMatchIdxForReco != -1;
    ghostRates->Fill(recoEta, 1. - wasTagged);
    if (recoEta < 0) {  // backward
      nTotalRecoTracksBackward++;
      nTotalGhostsBackward += (wasTagged == 0);
      nTotalRecoClonesBackward += isClone;
    } else {
      nTotalRecoTracksForward++;
      nTotalGhostsForward += (wasTagged == 0);
      nTotalRecoClonesForward += isClone;
    }
  }
  printf("\nTotal summary:\n");
  printf("\tForward: Eff %.2f%%, Ghost rate %.2f%%, Clone rate %.2f%%\n",
         100. * (float) nTotalRecodForward / nTotalForward,
         100. * (float) nTotalGhostsForward / nTotalRecoTracksForward,
         100. * (float) nTotalClonesForward / (nTotalClonesForward + nTotalRecodForward));
  printf("\tBackward: Eff %.2f%%, Ghost rate %.2f%%, Clone rate %.2f%%\n",
         100. * (float) nTotalRecodBackward / nTotalBackward,
         100. * (float) nTotalGhostsBackward / nTotalRecoTracksBackward,
         100. * (float) nTotalClonesBackward / (nTotalClonesBackward + nTotalRecodBackward));
  printf("\tThus overall: Eff %.2f%%, Ghost rate %.2f%%, Clone rate %.2f%%\n",
         100. * (float) (nTotalRecodForward + nTotalRecodBackward) /
                        (nTotalForward + nTotalBackward),
         100. * (float) (nTotalGhostsForward + nTotalGhostsBackward) /
                        (nTotalRecoTracksForward + nTotalRecoTracksBackward),
         100. * (float) (nTotalClonesForward + nTotalClonesBackward) /
                        ((nTotalClonesForward + nTotalRecodForward) +
                         (nTotalClonesBackward + nTotalRecodBackward)));
  // write histograms
  TString outDir = (Utils::Definitions::analysisRoot + "/hists/eff_clone_ghosts").c_str();
  gSystem->mkdir(outDir, true);  // TFile does not create the directory itself
  TFile* outFile = new TFile(outDir + "/mc_hists" + input_suffix, "RECREATE");
  hMCHitDistribution->Write();
  hDuplicateIDRates->Write();
  hUniqueIDRates->Write();
  hLongestMatchedTrackRate->Write();
  hLongestMatchedTrackRateClones->Write();
  hLongestMatchedTrackRate5Clones->Write();
  hCloneDistribution->Write();
  hCloneDistributionEta35_5->Write();
  hSplitTrackCloneDistribution->Write();
  hSplitTrackClone_1MissedDistribution->Write();
  hSplitTrackClone_2MissedDistribution->Write();
  hSeedingCloneDistribution->Write();
  hSeedingClone2Distribution->Write();
  hSeedingClone3Distribution->Write();
  hSeedingClone4Distribution->Write();
  hTripletCloneDistribution->Write();
  hTripletClonePlusDistribution->Write();
  hOverlapCloneDistribution->Write();
  hOtherCloneDistribution->Write();
  hNMatches->Write();
  ghostRates->Write();
  cloneRate->Write();
  cloneRatePT->Write();
  // turn the numerator/denominator pairs into efficiencies, named as the plotting macro expects
  for (unsigned i = 0; i < effReconstructible.size(); i++) {
    if (!TEfficiency::CheckConsistency(*effReconstructed[i], *effReconstructible[i])) {
      std::cerr << "Inconsistent histograms for " << effReconstructible[i]->GetName() << ".\n";
      continue;
    }
    TEfficiency* eff = new TEfficiency(*effReconstructed[i], *effReconstructible[i]);
    std::string name = effReconstructible[i]->GetName();
    eff->SetName(name.substr(0, name.find("_reconstructible")).c_str());
    eff->SetTitle(effReconstructible[i]->GetTitle());
    eff->Write();
  }
  splitTrackClonesMCByEta->Write();
  splitTrackClones_1MissedMCByEta->Write();
  splitTrackClones_2MissedMCByEta->Write();
  seedingClonesMCByEta->Write();
  seedingClones2MCByEta->Write();
  seedingClones3MCByEta->Write();
  seedingClones4MCByEta->Write();
  tripletClonesMCByEta->Write();
  tripletClonesPlusMCByEta->Write();
  moduleOverlapClonesMCByEta->Write();
  otherClonesMCByEta->Write();
  splitTrackClonesRecoByEta->Write();
  splitTrackClones_1MissedRecoByEta->Write();
  splitTrackClones_2MissedRecoByEta->Write();
  seedingClonesRecoByEta->Write();
  seedingClones2RecoByEta->Write();
  seedingClones3RecoByEta->Write();
  seedingClones4RecoByEta->Write();
  tripletClonesRecoByEta->Write();
  tripletClonesPlusRecoByEta->Write();
  moduleOverlapClonesRecoByEta->Write();
  otherClonesRecoByEta->Write();
  splitTrackCloneMCHitDistr->Write();
  splitTrackClone_1MissedMCHitDistr->Write();
  splitTrackClone_2MissedMCHitDistr->Write();
  nClonesByMCPSize_forward->Write();
  nClonesByMCPSize_backward->Write();
  nClonesByMCPSize_forward_scaled->Write();
  nClonesByMCPSize_backward_scaled->Write();
  cloneRateByMCPSize_forward->Write();
  cloneRateByMCPSize_backward->Write();
  cloneRateByMCPSize_forward_scaled->Write();
  cloneRateByMCPSize_backward_scaled->Write();
  seedingCloneMCHitDistr->Write();
  seedingClone2MCHitDistr->Write();
  seedingClone3MCHitDistr->Write();
  seedingClone4MCHitDistr->Write();
  tripletCloneMCHitDistr->Write();
  tripletClonePlusMCHitDistr->Write();
  moduleOverlapCloneMCHitDistr->Write();
  otherCloneMCHitDistr->Write();
  splitTrackCloneRecoHitDistr->Write();
  splitTrackClone_1MissedRecoHitDistr->Write();
  splitTrackClone_2MissedRecoHitDistr->Write();
  seedingCloneRecoHitDistr->Write();
  seedingClone2RecoHitDistr->Write();
  seedingClone3RecoHitDistr->Write();
  seedingClone4RecoHitDistr->Write();
  tripletCloneRecoHitDistr->Write();
  tripletClonePlusRecoHitDistr->Write();
  moduleOverlapCloneRecoHitDistr->Write();
  otherCloneRecoHitDistr->Write();
  // Clean up outfile
  outFile->Close();
  delete outFile;

  // clean up file
  file->Close();
  delete file;
}