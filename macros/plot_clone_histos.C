#include <TFile.h>
#include <TH1F.h>
#include <TEfficiency.h>
#include <TCanvas.h>

#include <cstdlib>
#include <iostream>
#include <string>

#include "utils/definitions.h"
#include "utils/basic_functions.h"


void plot_clone_histos(unsigned nEvents=5000, unsigned max_scatter=80000, 
                       unsigned max_dt=0, TString mc_suffix="") {
  gROOT->SetBatch();  // so stuff isn't autoplotted
  TString suffix;
  if (mc_suffix.IsNull()) {
    suffix = Utils::Functions::get_suffix(nEvents, max_scatter, max_dt);
  } else {
    suffix = Form("_%uev_%s", nEvents, mc_suffix.Data());
  }
  TString mcHistParent = (Utils::Definitions::analysisRoot + "/hists/eff_clone_ghosts/mc_hists").c_str();
  TString histosPath = mcHistParent + suffix + ".root";

  TFile* file = TFile::Open(histosPath);
  std::cout << "Opened file: " << histosPath << std::endl;

  TH1D* hMCHitDistribution = (TH1D*) file->Get("mc_hit_distribution");
  TProfile* hDuplicateIDRates = (TProfile*) file->Get("duplicate_match_id");
  TProfile* hUniqueIDRates = (TProfile*) file->Get("unique_id_match_rate");
  TProfile* hLongestMatchedTrackRate = (TProfile*) file->Get("longest_matched_track_fraction");
  TProfile* hLongestMatchedTrackRateClones = (TProfile*) file->Get("longest_matched_track_fraction_1_clones");
  TProfile* hLongestMatchedTrackRate5Clones = (TProfile*) file->Get("longest_matched_track_fraction_5_clones");
  TH1D* hCloneDistribution = (TH1D*) file->Get("clone_distribution");
  TH1D* hCloneDistributionEta35_5 = (TH1D*) file->Get("clone_distribution_eta35_5");
  TH1D* hSplitTrackCloneDistribution = (TH1D*) file->Get("split_track_clone_distribution");
  TH1D* hSplitTrackClone_1MissedDistribution = (TH1D*) file->Get("split_track_clone_1missed_distribution");
  TH1D* hSplitTrackClone_2MissedDistribution = (TH1D*) file->Get("split_track_clone_2missed_distribution");
  printf("Sizes of split track clone histograms: %f, %f, %f.\n",
         hSplitTrackCloneDistribution->GetEntries(),
         hSplitTrackClone_1MissedDistribution->GetEntries(),
         hSplitTrackClone_2MissedDistribution->GetEntries());
  TH1D* hSeedingCloneDistribution = (TH1D*) file->Get("seeding_clone_distribution");
  TH1D* hSeedingClone2Distribution = (TH1D*) file->Get("seeding_clone2_distribution");
  TH1D* hSeedingClone3Distribution = (TH1D*) file->Get("seeding_clone3_distribution");
  TH1D* hSeedingClone4Distribution = (TH1D*) file->Get("seeding_clone4_distribution");
  TH1D* hTripletCloneDistribution = (TH1D*) file->Get("triplet_clone_distribution");
  TH1D* hTripletClonePlusDistribution = (TH1D*) file->Get("triplet_clone_plus_distribution");
  TH1D* hOverlapCloneDistribution = (TH1D*) file->Get("overlap_clone_distribution");
  TH1D* hOtherCloneDistribution = (TH1D*) file->Get("other_clone_distribution");
  TProfile* hNMatches = (TProfile*) file->Get("number_of_matches_of_MC_track");
  TProfile* p_cloneRate = (TProfile*) file->Get("clone_rate");
  printf("Got all basic histograms from file.\n");

  // Clone types by MC and Reco, wrt eta
  TProfile* splitTrackClonesMCByEta = (TProfile*) file->Get("split_track_clones_mc_by_eta");
  TProfile* splitTrackClones_1MissedMCByEta = (TProfile*) file->Get("split_track_clones_1missed_mc_by_eta");
  TProfile* splitTrackClones_2MissedMCByEta = (TProfile*) file->Get("split_track_clones_2missed_mc_by_eta");
  printf("Sizes of split track clone MC by Eta Distributions: %f, %f, %f.\n",
         splitTrackClonesMCByEta->GetEntries(),
         splitTrackClones_1MissedMCByEta->GetEntries(),
         splitTrackClones_2MissedMCByEta->GetEntries());
  TProfile* seedingClonesMCByEta = (TProfile*) file->Get("seeding_clones_mc_by_eta");
  TProfile* seedingClones2MCByEta = (TProfile*) file->Get("seeding_clones2_mc_by_eta");
  TProfile* seedingClones3MCByEta = (TProfile*) file->Get("seeding_clones3_mc_by_eta");
  TProfile* seedingClones4MCByEta = (TProfile*) file->Get("seeding_clones4_mc_by_eta");
  TProfile* tripletClonesMCByEta = (TProfile*) file->Get("triplet_clones_mc_by_eta");
  TProfile* tripletClonesPlusMCByEta = (TProfile*) file->Get("triplet_clones_plus_mc_by_eta");
  TProfile* moduleOverlapClonesMCByEta = (TProfile*) file->Get("module_overlap_clones_mc_by_eta");
  TProfile* otherClonesMCByEta = (TProfile*) file->Get("other_clones_mc_by_eta");
  TProfile* splitTrackClonesRecoByEta = (TProfile*) file->Get("split_track_clones_reco_by_eta");
  TProfile* splitTrackClones_1MissedRecoByEta = (TProfile*) file->Get("split_track_clones_1missed_reco_by_eta");
  TProfile* splitTrackClones_2MissedRecoByEta = (TProfile*) file->Get("split_track_clones_2missed_reco_by_eta");
  printf("Sizes of split track clone Reco by Eta Distributions: %f, %f, %f.\n",
         splitTrackClonesRecoByEta->GetEntries(),
         splitTrackClones_1MissedRecoByEta->GetEntries(),
         splitTrackClones_2MissedRecoByEta->GetEntries());
  TProfile* seedingClonesRecoByEta = (TProfile*) file->Get("seeding_clones_reco_by_eta");
  TProfile* seedingClones2RecoByEta = (TProfile*) file->Get("seeding_clones2_reco_by_eta");
  TProfile* seedingClones3RecoByEta = (TProfile*) file->Get("seeding_clones3_reco_by_eta");
  TProfile* seedingClones4RecoByEta = (TProfile*) file->Get("seeding_clones4_reco_by_eta");
  TProfile* tripletClonesRecoByEta = (TProfile*) file->Get("triplet_clones_reco_by_eta");
  TProfile* tripletClonesPlusRecoByEta = (TProfile*) file->Get("triplet_clones_plus_reco_by_eta");
  TProfile* moduleOverlapClonesRecoByEta = (TProfile*) file->Get("module_overlap_clones_reco_by_eta");
  TProfile* otherClonesRecoByEta = (TProfile*) file->Get("other_clones_reco_by_eta");
  TH1D* nClonesByMCPSize_forward = (TH1D*) file->Get("n_clones_by_mc_p_size_forward");
  TH1D* nClonesByMCPSize_backward = (TH1D*) file->Get("n_clones_by_mc_p_size_backward");
  TH1D* nClonesByMCPSize_forward_scaled = (TH1D*) file->Get("n_clones_by_mc_p_size_forward_scaled");
  TH1D* nClonesByMCPSize_backward_scaled = (TH1D*) file->Get("n_clones_by_mc_p_size_backward_scaled");
  TProfile* cloneRateByMCPSize_forward = (TProfile*) file->Get("clone_rate_by_mc_p_size_forward");
  TProfile* cloneRateByMCPSize_backward = (TProfile*) file->Get("clone_rate_by_mc_p_size_backward");
  TProfile* cloneRateByMCPSize_forward_scaled = (TProfile*) file->Get("clone_rate_by_mc_p_size_forward_scaled");
  TProfile* cloneRateByMCPSize_backward_scaled = (TProfile*) file->Get("clone_rate_by_mc_p_size_backward_scaled");
  // Convert Profiles to Histograms
  TH1D* h_splitTrackClonesMCByEta = Utils::Functions::profileToHist(
    splitTrackClonesMCByEta, "splitTrackClonesMCByEta");
  TH1D* h_splitTrackClones_1MissedMCByEta = Utils::Functions::profileToHist(
    splitTrackClones_1MissedMCByEta, "splitTrackClones_1MissedMCByEta");
  TH1D* h_splitTrackClones_2MissedMCByEta = Utils::Functions::profileToHist(
    splitTrackClones_2MissedMCByEta, "splitTrackClones_2MissedMCByEta");
  TH1D* h_seedingClonesMCByEta = Utils::Functions::profileToHist(
    seedingClonesMCByEta, "seedingClonesMCByEta");
  TH1D* h_seedingClones2MCByEta = Utils::Functions::profileToHist(
    seedingClones2MCByEta, "seedingClones2MCByEta");
  TH1D* h_seedingClones3MCByEta = Utils::Functions::profileToHist(
    seedingClones3MCByEta, "seedingClones3MCByEta");
  TH1D* h_seedingClones4MCByEta = Utils::Functions::profileToHist(
    seedingClones4MCByEta, "seedingClones4MCByEta");
  TH1D* h_tripletClonesMCByEta = Utils::Functions::profileToHist(
    tripletClonesMCByEta, "tripletClonesMCByEta");
  TH1D* h_tripletClonesPlusMCByEta = Utils::Functions::profileToHist(
    tripletClonesPlusMCByEta, "tripletClonesPlusMCByEta");
  TH1D* h_moduleOverlapClonesMCByEta = Utils::Functions::profileToHist(
    moduleOverlapClonesMCByEta, "moduleOverlapClonesMCByEta");
  TH1D* h_otherClonesMCByEta = Utils::Functions::profileToHist(
    otherClonesMCByEta, "otherClonesMCByEta");
  TH1D* h_splitTrackClonesRecoByEta = Utils::Functions::profileToHist(
    splitTrackClonesRecoByEta, "splitTrackClonesRecoByEta");
  TH1D* h_splitTrackClones_1MissedRecoByEta = Utils::Functions::profileToHist(
    splitTrackClones_1MissedRecoByEta, "splitTrackClones_1MissedRecoByEta");
  TH1D* h_splitTrackClones_2MissedRecoByEta = Utils::Functions::profileToHist(
    splitTrackClones_2MissedRecoByEta, "splitTrackClones_2MissedRecoByEta");
  TH1D* h_seedingClonesRecoByEta = Utils::Functions::profileToHist(
    seedingClonesRecoByEta, "seedingClonesRecoByEta");
  TH1D* h_seedingClones2RecoByEta = Utils::Functions::profileToHist(
    seedingClones2RecoByEta, "seedingClones2RecoByEta");
  TH1D* h_seedingClones3RecoByEta = Utils::Functions::profileToHist(
    seedingClones3RecoByEta, "seedingClones3RecoByEta");
  TH1D* h_seedingClones4RecoByEta = Utils::Functions::profileToHist(
    seedingClones4RecoByEta, "seedingClones4RecoByEta");
  TH1D* h_tripletClonesRecoByEta = Utils::Functions::profileToHist(
    tripletClonesRecoByEta, "tripletClonesRecoByEta");
  TH1D* h_tripletClonesPlusRecoByEta = Utils::Functions::profileToHist(
    tripletClonesPlusRecoByEta, "tripletClonesPlusRecoByEta");
  TH1D* h_moduleOverlapClonesRecoByEta = Utils::Functions::profileToHist(
    moduleOverlapClonesRecoByEta, "moduleOverlapClonesRecoByEta");
  TH1D* h_otherClonesRecoByEta = Utils::Functions::profileToHist(
    otherClonesRecoByEta, "otherClonesRecoByEta");
  TH1D* h_cloneRateByMCPSize_forward = Utils::Functions::profileToHist(
    cloneRateByMCPSize_forward, "cloneRateByMCPSize_forward");
  TH1D* h_cloneRateByMCPSize_backward = Utils::Functions::profileToHist(
    cloneRateByMCPSize_backward, "cloneRateByMCPSize_backward");
  TH1D* h_cloneRateByMCPSize_forward_scaled = Utils::Functions::profileToHist(
    cloneRateByMCPSize_forward_scaled, "cloneRateByMCPSize_forward_scaled");
  TH1D* h_cloneRateByMCPSize_backward_scaled = Utils::Functions::profileToHist(
    cloneRateByMCPSize_backward_scaled, "cloneRateByMCPSize_backward_scaled");
  printf("Got all clone type profiles wrt eta from file and converted to histograms.\n");

  // clone hit distributions
  TH1D* h_splitTrackCloneMCHitDistr = (TH1D*) file->Get("split_track_clone_mc_hit_distribution");
  TH1D* h_splitTrackClone_1MissedMCHitDistr = (TH1D*) file->Get("split_track_clone_1missed_mc_hit_distribution");
  TH1D* h_splitTrackClone_2MissedMCHitDistr = (TH1D*) file->Get("split_track_clone_2missed_mc_hit_distribution");
  printf("Sizes of split track hit clone MC Distributions: %f, %f, %f.\n",
         h_splitTrackCloneMCHitDistr->GetEntries(),
         h_splitTrackClone_1MissedMCHitDistr->GetEntries(),
         h_splitTrackClone_2MissedMCHitDistr->GetEntries());
  TH1D* h_seedingCloneMCHitDistr = (TH1D*) file->Get("seeding_clone_mc_hit_distribution");
  TH1D* h_seedingClone2MCHitDistr = (TH1D*) file->Get("seeding_clone2_mc_hit_distribution");
  TH1D* h_seedingClone3MCHitDistr = (TH1D*) file->Get("seeding_clone3_mc_hit_distribution");
  TH1D* h_seedingClone4MCHitDistr = (TH1D*) file->Get("seeding_clone4_mc_hit_distribution");
  TH1D* h_tripletCloneMCHitDistr = (TH1D*) file->Get("triplet_clone_mc_hit_distribution");
  TH1D* h_tripletClonePlusMCHitDistr = (TH1D*) file->Get("triplet_clone_plus_mc_hit_distribution");
  TH1D* h_moduleOverlapCloneMCHitDistr = (TH1D*) file->Get("module_overlap_clone_mc_hit_distribution");
  TH1D* h_otherCloneMCHitDistr = (TH1D*) file->Get("other_clone_mc_hit_distribution");
  TH1D* h_splitTrackCloneRecoHitDistr = (TH1D*) file->Get("split_track_clone_reco_hit_distribution");
  TH1D* h_splitTrackClone_1MissedRecoHitDistr = (TH1D*) file->Get("split_track_clone_1missed_reco_hit_distribution");
  TH1D* h_splitTrackClone_2MissedRecoHitDistr = (TH1D*) file->Get("split_track_clone_2missed_reco_hit_distribution");
  printf("Sizes of split track hit clone Reco Distributions: %f, %f, %f.\n",
         h_splitTrackCloneRecoHitDistr->GetEntries(),
         h_splitTrackClone_1MissedRecoHitDistr->GetEntries(),
         h_splitTrackClone_2MissedRecoHitDistr->GetEntries());
  TH1D* h_seedingCloneRecoHitDistr = (TH1D*) file->Get("seeding_clone_reco_hit_distribution");
  TH1D* h_seedingClone2RecoHitDistr = (TH1D*) file->Get("seeding_clone2_reco_hit_distribution");
  TH1D* h_seedingClone3RecoHitDistr = (TH1D*) file->Get("seeding_clone3_reco_hit_distribution");
  TH1D* h_seedingClone4RecoHitDistr = (TH1D*) file->Get("seeding_clone4_reco_hit_distribution");
  TH1D* h_tripletCloneRecoHitDistr = (TH1D*) file->Get("triplet_clone_reco_hit_distribution");
  TH1D* h_tripletClonePlusRecoHitDistr = (TH1D*) file->Get("triplet_clone_plus_reco_hit_distribution");
  TH1D* h_moduleOverlapCloneRecoHitDistr = (TH1D*) file->Get("module_overlap_clone_reco_hit_distribution");
  TH1D* h_otherCloneRecoHitDistr = (TH1D*) file->Get("other_clone_reco_hit_distribution");
  printf("Got all clone hit distributions from file.\n");

  // make stack plots from both MC and reco clone types
  THStack* clonesMCSum = new THStack("clonesMCSum", "MC Clone Types;#eta;Clone Rate");
  // Make all seeding clones red initially before splitting them later
  h_splitTrackClonesMCByEta->SetFillColor(kYellow);
  h_splitTrackClones_1MissedMCByEta->SetFillColor(kOrange);
  h_splitTrackClones_2MissedMCByEta->SetFillColor(kOrange - 7);
  TH1D* h_seedingClonesMCByEtaSummed =
    (TH1D*) h_seedingClonesMCByEta->Clone("h_seedingClonesMCByEtaSummed");
  h_seedingClonesMCByEtaSummed->Add(h_seedingClones2MCByEta);
  h_seedingClonesMCByEtaSummed->Add(h_seedingClones3MCByEta);
  h_seedingClonesMCByEtaSummed->Add(h_seedingClones4MCByEta);
  h_seedingClonesMCByEtaSummed->SetFillColor(kRed);
  TH1D* h_tripletClonesMCByEtaSummed =
    (TH1D*) h_tripletClonesMCByEta->Clone("h_tripletClonesMCByEtaSummed");
  h_tripletClonesMCByEtaSummed->Add(h_tripletClonesPlusMCByEta);
  h_tripletClonesMCByEtaSummed->SetFillColor(kBlue);
  h_moduleOverlapClonesMCByEta->SetFillColor(kGreen);
  h_otherClonesMCByEta->SetFillColor(kMagenta);
  clonesMCSum->Add(h_seedingClonesMCByEtaSummed);
  clonesMCSum->Add(h_tripletClonesMCByEtaSummed);
  clonesMCSum->Add(h_moduleOverlapClonesMCByEta);
  clonesMCSum->Add(h_splitTrackClones_2MissedMCByEta);
  clonesMCSum->Add(h_splitTrackClones_1MissedMCByEta);
  clonesMCSum->Add(h_splitTrackClonesMCByEta);
  clonesMCSum->Add(h_otherClonesMCByEta);
  // Add legend for MC clone stack
  TLegend* mcCloneLegend = new TLegend(0.35, 0.5, 0.55, 0.7);
  // remove edges from legend
  mcCloneLegend->SetBorderSize(0);
  mcCloneLegend->AddEntry(h_splitTrackClonesMCByEta, "LO Split Track", "f");
  mcCloneLegend->AddEntry(h_splitTrackClones_1MissedMCByEta, "NLO Split Track", "f");
  mcCloneLegend->AddEntry(h_splitTrackClones_2MissedMCByEta, "NNLO Split Track", "f");
  mcCloneLegend->AddEntry(h_seedingClonesMCByEtaSummed, "Seeding", "f");
  mcCloneLegend->AddEntry(h_tripletClonesMCByEtaSummed, "Triplet", "f");
  mcCloneLegend->AddEntry(h_moduleOverlapClonesMCByEta, "Module overlap", "f");
  mcCloneLegend->AddEntry(h_otherClonesMCByEta, "Other", "f");

  THStack* clonesRecoSum = new THStack("clonesRecoSum", "Reco Clone Types;#eta;Clone Rate");
  h_splitTrackClonesRecoByEta->SetFillColor(kYellow);
  h_splitTrackClones_1MissedRecoByEta->SetFillColor(kOrange);
  h_splitTrackClones_2MissedRecoByEta->SetFillColor(kOrange - 7);
  TH1D* h_seedingClonesRecoByEtaSummed =
    (TH1D*) h_seedingClonesRecoByEta->Clone("h_seedingClonesRecoByEtaSummed");
  h_seedingClonesRecoByEtaSummed->Add(h_seedingClones2RecoByEta);
  h_seedingClonesRecoByEtaSummed->Add(h_seedingClones3RecoByEta);
  h_seedingClonesRecoByEtaSummed->Add(h_seedingClones4RecoByEta);
  h_seedingClonesRecoByEtaSummed->SetFillColor(kRed);
  TH1D* h_tripletClonesRecoByEtaSummed =
    (TH1D*) h_tripletClonesRecoByEta->Clone("h_tripletClonesRecoByEtaSummed");
  h_tripletClonesRecoByEtaSummed->Add(h_tripletClonesPlusRecoByEta);
  h_tripletClonesRecoByEtaSummed->SetFillColor(kBlue);
  h_moduleOverlapClonesRecoByEta->SetFillColor(kGreen);
  h_otherClonesRecoByEta->SetFillColor(kMagenta);
  clonesRecoSum->Add(h_seedingClonesRecoByEtaSummed);
  clonesRecoSum->Add(h_tripletClonesRecoByEtaSummed);
  clonesRecoSum->Add(h_moduleOverlapClonesRecoByEta);
  clonesRecoSum->Add(h_splitTrackClones_2MissedRecoByEta);
  clonesRecoSum->Add(h_splitTrackClones_1MissedRecoByEta);
  clonesRecoSum->Add(h_splitTrackClonesRecoByEta);
  clonesRecoSum->Add(h_otherClonesRecoByEta);
  // Add legend for Reco clone stack
  TLegend* recoCloneLegend = new TLegend(0.35, 0.5, 0.55, 0.7);
  recoCloneLegend->SetBorderSize(0);
  recoCloneLegend->AddEntry(h_splitTrackClonesRecoByEta, "LO Split Track", "f");
  recoCloneLegend->AddEntry(h_splitTrackClones_1MissedRecoByEta, "NLO Split Track", "f");
  recoCloneLegend->AddEntry(h_splitTrackClones_2MissedRecoByEta, "NNLO Split Track", "f");
  recoCloneLegend->AddEntry(h_seedingClonesRecoByEtaSummed, "Seeding", "f");
  recoCloneLegend->AddEntry(h_tripletClonesRecoByEtaSummed, "Triplet", "f");
  recoCloneLegend->AddEntry(h_moduleOverlapClonesRecoByEta, "Module overlap", "f");
  recoCloneLegend->AddEntry(h_otherClonesRecoByEta, "Other", "f");
  printf("Read in all histograms and created Stacks for MC and Reco clone types wrt eta.\n");

  // -------------- Also make stack plots for the clone distributions
  THStack* cloneDistrStack = new THStack("cloneDistrStack", "Clone distributions;N_{clones};Frequency");
  hSplitTrackCloneDistribution->SetFillColor(kYellow);
  hSplitTrackClone_1MissedDistribution->SetFillColor(kOrange);
  hSplitTrackClone_2MissedDistribution->SetFillColor(kOrange - 7);
  TH1D* hSeedingCloneDistributionSummed = 
    (TH1D*) hSeedingCloneDistribution->Clone("hSeedingCloneDistributionSummed");
  hSeedingCloneDistributionSummed->Add(hSeedingClone2Distribution);
  hSeedingCloneDistributionSummed->Add(hSeedingClone3Distribution);
  hSeedingCloneDistributionSummed->Add(hSeedingClone4Distribution);
  hSeedingCloneDistributionSummed->SetFillColor(kRed);
  TH1D* hTripletCloneDistributionSummed =
    (TH1D*) hTripletCloneDistribution->Clone("hTripletCloneDistributionSummed");
  hTripletCloneDistributionSummed->Add(hTripletClonePlusDistribution);
  hTripletCloneDistributionSummed->SetFillColor(kBlue);
  hOverlapCloneDistribution->SetFillColor(kGreen);
  hOtherCloneDistribution->SetFillColor(kMagenta);
  cloneDistrStack->Add(hOtherCloneDistribution);
  cloneDistrStack->Add(hTripletCloneDistributionSummed);
  cloneDistrStack->Add(hSeedingCloneDistributionSummed);
  cloneDistrStack->Add(hSplitTrackClone_2MissedDistribution);
  cloneDistrStack->Add(hSplitTrackClone_1MissedDistribution);
  cloneDistrStack->Add(hSplitTrackCloneDistribution);
  cloneDistrStack->Add(hOverlapCloneDistribution);
  // Add legend for clone distribution stack
  TLegend* cloneDistrLegend = new TLegend(0.6, 0.4, 0.8, 0.6);
  cloneDistrLegend->SetBorderSize(0);
  cloneDistrLegend->AddEntry(hSplitTrackCloneDistribution, "LO Split Track", "f");
  cloneDistrLegend->AddEntry(hSplitTrackClone_1MissedDistribution, "NLO Split Track", "f");
  cloneDistrLegend->AddEntry(hSplitTrackClone_2MissedDistribution, "NNLO Split Track", "f");
  cloneDistrLegend->AddEntry(hSeedingCloneDistributionSummed, "Seeding", "f");
  cloneDistrLegend->AddEntry(hTripletCloneDistributionSummed, "Triplet", "f");
  cloneDistrLegend->AddEntry(hOverlapCloneDistribution, "Module overlap", "f");
  cloneDistrLegend->AddEntry(hOtherCloneDistribution, "Other", "f");
  printf("Made stack plots for clone number distributions.\n");
  // ------------- get ratios of hit distributions of various clone types
  // First sum up the various types of clones
  TH1D* hSeedingClonesHitMCDistrSummed = 
    (TH1D*) h_seedingCloneMCHitDistr->Clone("hSeedingClonesHitMCDistrSummed");
  hSeedingClonesHitMCDistrSummed->Add(h_seedingClone2MCHitDistr);
  hSeedingClonesHitMCDistrSummed->Add(h_seedingClone3MCHitDistr);
  hSeedingClonesHitMCDistrSummed->Add(h_seedingClone4MCHitDistr);
  TH1D* hTripletClonesHitMCDistrSummed =
    (TH1D*) h_tripletCloneMCHitDistr->Clone("hTripletClonesHitMCDistrSummed");
  hTripletClonesHitMCDistrSummed->Add(h_tripletClonePlusMCHitDistr);
  TH1D* hSeedingClonesHitRecoDistrSummed = 
    (TH1D*) h_seedingCloneRecoHitDistr->Clone("hSeedingClonesHitRecoDistrSummed");
  hSeedingClonesHitRecoDistrSummed->Add(h_seedingClone2RecoHitDistr);
  hSeedingClonesHitRecoDistrSummed->Add(h_seedingClone3RecoHitDistr);
  hSeedingClonesHitRecoDistrSummed->Add(h_seedingClone4RecoHitDistr);
  TH1D* hTripletClonesHitRecoDistrSummed =
    (TH1D*) h_tripletCloneRecoHitDistr->Clone("hTripletClonesHitRecoDistrSummed");
  hTripletClonesHitRecoDistrSummed->Add(h_tripletClonePlusRecoHitDistr);
  // Now get ratios
  TH1D* hSplitTrackCloneHitDistrRatio =
    (TH1D*) h_splitTrackCloneMCHitDistr->Clone("hSplitTrackCloneHitDistrRatio");
  hSplitTrackCloneHitDistrRatio->Divide(h_splitTrackCloneRecoHitDistr);
  TH1D* hSplitTrackClone_1MissedHitDistrRatio =
    (TH1D*) h_splitTrackClone_1MissedMCHitDistr->Clone("hSplitTrackClone_1MissedHitDistrRatio");
  hSplitTrackClone_1MissedHitDistrRatio->Divide(h_splitTrackClone_1MissedRecoHitDistr);
  TH1D* hSplitTrackClone_2MissedHitDistrRatio =
    (TH1D*) h_splitTrackClone_2MissedMCHitDistr->Clone("hSplitTrackClone_2MissedHitDistrRatio");
  hSplitTrackClone_2MissedHitDistrRatio->Divide(h_splitTrackClone_2MissedRecoHitDistr);
  TH1D* hTotalSeedingCloneHitDistrRatio =
    (TH1D*) hSeedingClonesHitMCDistrSummed->Clone("hSeedingCloneMCHitDistrRatio");
  hTotalSeedingCloneHitDistrRatio->Divide(hSeedingClonesHitRecoDistrSummed);
  TH1D* hTotalTripletCloneHitDistrRatio =
    (TH1D*) hTripletClonesHitMCDistrSummed->Clone("hTripletCloneMCHitDistrRatio");
  hTotalTripletCloneHitDistrRatio->Divide(hTripletClonesHitRecoDistrSummed);
  TH1D* hModuleOverlapCloneHitDistrRatio =
    (TH1D*) h_moduleOverlapCloneMCHitDistr->Clone("hModuleOverlapCloneMCHitDistrRatio");
  hModuleOverlapCloneHitDistrRatio->Divide(h_moduleOverlapCloneRecoHitDistr);
  TH1D* hOtherCloneHitDistrRatio =
    (TH1D*) h_otherCloneMCHitDistr->Clone("hOtherCloneMCHitDistrRatio");
  hOtherCloneHitDistrRatio->Divide(h_otherCloneRecoHitDistr);
  // Together with this, get the max of these distributions for visualisation
  double maxHitNum = hSeedingClonesHitMCDistrSummed->GetMaximum();
  if (hTripletClonesHitMCDistrSummed->GetMaximum() > maxHitNum) {
    maxHitNum = hTripletClonesHitMCDistrSummed->GetMaximum();
  }
  if (h_moduleOverlapCloneMCHitDistr->GetMaximum() > maxHitNum) {
    maxHitNum = h_moduleOverlapCloneMCHitDistr->GetMaximum();
  }
  if (h_otherCloneMCHitDistr->GetMaximum() > maxHitNum) {
    maxHitNum = h_otherCloneMCHitDistr->GetMaximum();
  }
  if (h_splitTrackCloneMCHitDistr->GetMaximum() > maxHitNum) {
    maxHitNum = h_splitTrackCloneMCHitDistr->GetMaximum();
  }
  if (h_splitTrackClone_1MissedMCHitDistr->GetMaximum() > maxHitNum) {
    maxHitNum = h_splitTrackClone_1MissedMCHitDistr->GetMaximum();
  }
  if (h_splitTrackClone_2MissedMCHitDistr->GetMaximum() > maxHitNum) {
    maxHitNum = h_splitTrackClone_2MissedMCHitDistr->GetMaximum();
  }
  // Now do the same with the ratio histograms
  double maxHitNumRatio = hTotalSeedingCloneHitDistrRatio->GetMaximum(),
         minHitNumRatio = hTotalSeedingCloneHitDistrRatio->GetMinimum();
  if (hTotalTripletCloneHitDistrRatio->GetMaximum() > maxHitNumRatio) {
    maxHitNumRatio = hTotalTripletCloneHitDistrRatio->GetMaximum();
  }
  if (hModuleOverlapCloneHitDistrRatio->GetMaximum() > maxHitNumRatio) {
    maxHitNumRatio = hModuleOverlapCloneHitDistrRatio->GetMaximum();
  }
  if (hOtherCloneHitDistrRatio->GetMaximum() > maxHitNumRatio) {
    maxHitNumRatio = hOtherCloneHitDistrRatio->GetMaximum();
  }
  if (hSplitTrackCloneHitDistrRatio->GetMaximum() > maxHitNumRatio) {
    maxHitNumRatio = hSplitTrackCloneHitDistrRatio->GetMaximum();
  }
  if (hSplitTrackClone_1MissedHitDistrRatio->GetMaximum() > maxHitNumRatio) {
    maxHitNumRatio = hSplitTrackClone_1MissedHitDistrRatio->GetMaximum();
  }
  if (hSplitTrackClone_2MissedHitDistrRatio->GetMaximum() > maxHitNumRatio) {
    maxHitNumRatio = hSplitTrackClone_2MissedHitDistrRatio->GetMaximum();
  }
  if (hTotalSeedingCloneHitDistrRatio->GetMinimum() < minHitNumRatio) {
    minHitNumRatio = hTotalSeedingCloneHitDistrRatio->GetMinimum();
  }
  if (hTotalTripletCloneHitDistrRatio->GetMinimum() < minHitNumRatio) {
    minHitNumRatio = hTotalTripletCloneHitDistrRatio->GetMinimum();
  }
  if (hModuleOverlapCloneHitDistrRatio->GetMinimum() < minHitNumRatio) {
    minHitNumRatio = hModuleOverlapCloneHitDistrRatio->GetMinimum();
  }
  if (hOtherCloneHitDistrRatio->GetMinimum() < minHitNumRatio) {
    minHitNumRatio = hOtherCloneHitDistrRatio->GetMinimum();
  }
  if (hSplitTrackCloneHitDistrRatio->GetMinimum() < minHitNumRatio) {
    minHitNumRatio = hSplitTrackCloneHitDistrRatio->GetMinimum();
  }
  if (hSplitTrackClone_1MissedHitDistrRatio->GetMinimum() < minHitNumRatio) {
    minHitNumRatio = hSplitTrackClone_1MissedHitDistrRatio->GetMinimum();
  }
  if (hSplitTrackClone_2MissedHitDistrRatio->GetMinimum() < minHitNumRatio) {
    minHitNumRatio = hSplitTrackClone_2MissedHitDistrRatio->GetMinimum();
  }
  // Set the ranges of the histograms
  hSeedingClonesHitMCDistrSummed->GetYaxis()->SetRangeUser(0., maxHitNum * 1.1);
  hSeedingClonesHitMCDistrSummed->GetXaxis()->SetRangeUser(3, 15);
  hTotalSeedingCloneHitDistrRatio->GetYaxis()->SetRangeUser(minHitNumRatio/1.1, maxHitNumRatio*1.1);
  hTotalSeedingCloneHitDistrRatio->GetXaxis()->SetRangeUser(3, 15);
  printf("Made ratio plots of clone hit distributions and set ranges.\n");


  // Now plot
  TCanvas* canvas = new TCanvas("c", "MC Clone Information", 1200, 1000);
  canvas->Divide(2, 2);

  gStyle->SetOptStat(0);  // remove the info box
  canvas->cd(1);
  // hDuplicateIDRates->Draw();
  nClonesByMCPSize_forward->SetLineColor(kRed);
  nClonesByMCPSize_backward->SetLineColor(kBlue);
  nClonesByMCPSize_forward_scaled->SetLineColor(kRed);
  nClonesByMCPSize_backward_scaled->SetLineColor(kBlue);
  nClonesByMCPSize_forward->SetMarkerColor(kRed);
  nClonesByMCPSize_backward->SetMarkerColor(kBlue);
  nClonesByMCPSize_forward_scaled->SetMarkerColor(kRed);
  nClonesByMCPSize_backward_scaled->SetMarkerColor(kBlue);
  nClonesByMCPSize_forward->SetMarkerStyle(21);
  nClonesByMCPSize_backward->SetMarkerStyle(21);
  nClonesByMCPSize_forward_scaled->SetMarkerStyle(22);
  nClonesByMCPSize_backward_scaled->SetMarkerStyle(22);
  // scale by the number of events to get clones/event (consistent)
  nClonesByMCPSize_forward->Scale(1./nEvents);
  nClonesByMCPSize_backward->Scale(1./nEvents);
  nClonesByMCPSize_forward_scaled->Scale(1./nEvents);
  nClonesByMCPSize_backward_scaled->Scale(1./nEvents);
  nClonesByMCPSize_forward->SetTitle("N_{clones}/Event by N_{MCP hits}");
  nClonesByMCPSize_forward->Draw("P");
  nClonesByMCPSize_backward->Draw("P SAME");
  nClonesByMCPSize_forward_scaled->Draw("P SAME");
  nClonesByMCPSize_backward_scaled->Draw("P SAME");
  nClonesByMCPSize_forward->GetYaxis()->SetRangeUser(0.001, 10.);
  gPad->SetLogy(1);
  // Add legend
  TLegend* leg_Nclones_wrt_mcp_size = new TLegend(0.2, 0.2, 0.4, 0.4);
  leg_Nclones_wrt_mcp_size->AddEntry(nClonesByMCPSize_forward, "Forward", "P");
  leg_Nclones_wrt_mcp_size->AddEntry(nClonesByMCPSize_backward, "Backward", "P");
  leg_Nclones_wrt_mcp_size->AddEntry(nClonesByMCPSize_forward_scaled, "Forward scaled", "P");
  leg_Nclones_wrt_mcp_size->AddEntry(nClonesByMCPSize_backward_scaled, "Backward scaled", "P");
  leg_Nclones_wrt_mcp_size->Draw();

  canvas->cd(2);
  // hUniqueIDRates->Draw(); <-- this is just the 1- of the above duplicate rate
  // set line colours and plot all histograms
  // hLongestMatchedTrackRate->SetLineColor(kBlack);
  // hLongestMatchedTrackRate->Draw();

  // hLongestMatchedTrackRateClones->SetLineColor(kMagenta);
  // hLongestMatchedTrackRateClones->Draw("SAME");

  // hLongestMatchedTrackRate5Clones->SetLineColor(kBlue);
  // hLongestMatchedTrackRate5Clones->Draw("SAME");
  // // Add legend
  // TLegend* leg_lmtRate = new TLegend(0.35, 0.2, 0.55, 0.4);
  // leg_lmtRate->AddEntry(hLongestMatchedTrackRate, "No clones", "l");
  // leg_lmtRate->AddEntry(hLongestMatchedTrackRateClones, "1 #leq N_{clones} #leq 5", "l");
  // leg_lmtRate->AddEntry(hLongestMatchedTrackRate5Clones, "N_{clones} > 5", "l");
  // leg_lmtRate->Draw();

  // plot the generic hit distribution
  hMCHitDistribution->SetLineColor(kBlack);
  hMCHitDistribution->SetTitle("Hit distribution of all MCs");
  // hMCHitDistribution->GetXaxis()->SetRangeUser(3, 15);
  hMCHitDistribution->Draw();

  // clone rate wrt MCP track length
  canvas->cd(3);
  h_cloneRateByMCPSize_forward->SetLineColor(kRed);
  h_cloneRateByMCPSize_backward->SetLineColor(kBlue);
  h_cloneRateByMCPSize_forward_scaled->SetLineColor(kRed);
  h_cloneRateByMCPSize_backward_scaled->SetLineColor(kBlue);
  h_cloneRateByMCPSize_forward->SetMarkerColor(kRed);
  h_cloneRateByMCPSize_backward->SetMarkerColor(kBlue);
  h_cloneRateByMCPSize_forward_scaled->SetMarkerColor(kRed);
  h_cloneRateByMCPSize_backward_scaled->SetMarkerColor(kBlue);
  h_cloneRateByMCPSize_forward->SetMarkerStyle(21);
  h_cloneRateByMCPSize_backward->SetMarkerStyle(21);
  h_cloneRateByMCPSize_forward_scaled->SetMarkerStyle(22);
  h_cloneRateByMCPSize_backward_scaled->SetMarkerStyle(22);
  h_cloneRateByMCPSize_forward->Draw("P");
  h_cloneRateByMCPSize_backward->Draw("P SAME");
  h_cloneRateByMCPSize_forward_scaled->Draw("P SAME");
  h_cloneRateByMCPSize_backward_scaled->Draw("P SAME");
  h_cloneRateByMCPSize_forward->GetYaxis()->SetRangeUser(0.01, 0.8);
  // set y axis to log scale
  gPad->SetLogy(1);

  TLegend* leg_cloneRate_wrt_mcp_size = new TLegend(0.4, 0.5, 0.6, 0.7);
  leg_cloneRate_wrt_mcp_size->AddEntry(nClonesByMCPSize_forward, "Forward", "P");
  leg_cloneRate_wrt_mcp_size->AddEntry(nClonesByMCPSize_backward, "Backward", "P");
  leg_cloneRate_wrt_mcp_size->AddEntry(nClonesByMCPSize_forward_scaled, "Forward scaled", "P");
  leg_cloneRate_wrt_mcp_size->AddEntry(nClonesByMCPSize_backward_scaled, "Backward scaled", "P");
  leg_cloneRate_wrt_mcp_size->Draw();

  canvas->cd(4);
  p_cloneRate->SetLineColor(kMagenta);
  p_cloneRate->GetYaxis()->SetRangeUser(0., 0.2);
  p_cloneRate->Draw();

  TString outputBase = (Utils::Definitions::analysisRoot + "/output/clones/").c_str();
  TString cloneOutputBase = outputBase + TString("mc_clone_plots");
  canvas->SaveAs(cloneOutputBase + suffix + ".pdf");
  delete canvas;
  // delete leg_lmtRate;
  delete leg_cloneRate_wrt_mcp_size;
  delete leg_Nclones_wrt_mcp_size;

  // make another canvas for the clone distributions
  TCanvas* canvas2 = new TCanvas("c2", "Clone Distributions", 1200, 1000);
  canvas2->Divide(2, 2);
  canvas2->cd(1);
  gPad->SetLogy(1);
  // hCloneDistribution->SetTitle("Clone distributions");
  // hCloneDistribution->SetLineColor(kBlack);
  // hCloneDistribution->Draw();
  // hCloneDistributionEta35_5->SetLineColor(kRed);
  // hCloneDistributionEta35_5->Draw("SAME");
  // // Add legend
  // TLegend* cloneDistrLegend = new TLegend(0.6, 0.6, 0.8, 0.8);
  // cloneDistrLegend->AddEntry(hCloneDistribution, "All clones", "l");
  // cloneDistrLegend->AddEntry(hCloneDistributionEta35_5, "3.5 < #eta < 5", "l");
  // cloneDistrLegend->Draw();
  cloneDistrStack->Draw("HIST");
  cloneDistrLegend->Draw();

  canvas2->cd(2);
  // The y axis here runs to ~77000, so the labels are 5 digits wide. The
  // default left margin of 0.1 cannot fit those labels and the rotated y axis
  // title, which then falls outside the pad and gets clipped.
  gPad->SetLeftMargin(0.15);
  // TPad* basepad = (TPad*) gPad;
  // basepad->cd();  // Activate subpad (essential for not overwriting each other)

  // TPad* padTop = new TPad("padTop", "Main Plot", 0.0, 0.3, 1.0, 1.0);
  // TPad* padRatio = new TPad("padRatio", "Ratio Plot", 0.0, 0.0, 1.0, 0.3);

  // padTop->SetBottomMargin(0.0);
  // padTop->Draw();

  // padRatio->SetTopMargin(0.0);
  // padRatio->SetBottomMargin(0.3);
  // padRatio->Draw();

  // padTop->cd();
  hSeedingClonesHitMCDistrSummed->SetTitle("Hit distribution of MCPs with at least one clone");
  hSeedingClonesHitMCDistrSummed->SetLineColor(kRed);
  hSeedingClonesHitMCDistrSummed->Draw("HIST");
  h_splitTrackCloneMCHitDistr->SetLineColor(kYellow);
  h_splitTrackCloneMCHitDistr->Draw("SAME");
  h_splitTrackClone_1MissedMCHitDistr->SetLineColor(kOrange);
  h_splitTrackClone_1MissedMCHitDistr->Draw("SAME");
  h_splitTrackClone_2MissedMCHitDistr->SetLineColor(kOrange - 7);
  h_splitTrackClone_2MissedMCHitDistr->Draw("SAME");
  hTripletClonesHitMCDistrSummed->SetLineColor(kBlue);
  hTripletClonesHitMCDistrSummed->Draw("SAME");
  h_moduleOverlapCloneMCHitDistr->SetLineColor(kGreen);
  h_moduleOverlapCloneMCHitDistr->Draw("SAME");
  h_otherCloneMCHitDistr->SetLineColor(kMagenta);
  h_otherCloneMCHitDistr->Draw("SAME");
  // Add legend
  TLegend* cloneHitDistrLegend = new TLegend(0.62, 0.6, 0.82, 0.85);
  cloneHitDistrLegend->SetBorderSize(0);
  cloneHitDistrLegend->SetFillStyle(0);  // transparent, so the frame box shows through
  cloneHitDistrLegend->AddEntry(h_splitTrackCloneMCHitDistr, "LO Split Track", "l");
  cloneHitDistrLegend->AddEntry(h_splitTrackClone_1MissedMCHitDistr, "NLO Split Track", "l");
  cloneHitDistrLegend->AddEntry(h_splitTrackClone_2MissedMCHitDistr, "NNLO Split Track", "l");
  cloneHitDistrLegend->AddEntry(hSeedingClonesHitMCDistrSummed, "Seeding", "l");
  cloneHitDistrLegend->AddEntry(hTripletClonesHitMCDistrSummed, "Triplet", "l");
  cloneHitDistrLegend->AddEntry(h_moduleOverlapCloneMCHitDistr, "Module overlap", "l");
  cloneHitDistrLegend->AddEntry(h_otherCloneMCHitDistr, "Other", "l");
  cloneHitDistrLegend->Draw();
  // Move to ratio plots
  // padRatio->cd();
  // hTotalSeedingCloneHitDistrRatio->SetTitle(""); // remove title
  // hTotalSeedingCloneHitDistrRatio->SetLineColor(kRed);
  // hTotalSeedingCloneHitDistrRatio->GetYaxis()->SetTitle("MC/Reco");
  // hTotalSeedingCloneHitDistrRatio->GetXaxis()->SetTitleSize(0.12);  // Title font size
  // hTotalSeedingCloneHitDistrRatio->GetXaxis()->SetTitleOffset(1.0);  // Distance from axis
  // hTotalSeedingCloneHitDistrRatio->GetXaxis()->SetLabelSize(0.10);  // label font size
  // hTotalSeedingCloneHitDistrRatio->GetYaxis()->SetTitleSize(0.10);
  // hTotalSeedingCloneHitDistrRatio->GetYaxis()->SetTitleOffset(0.4);
  // hTotalSeedingCloneHitDistrRatio->GetYaxis()->SetLabelSize(0.08);
  // hTotalSeedingCloneHitDistrRatio->GetYaxis()->SetNdivisions(305);
  // hTotalSeedingCloneHitDistrRatio->Draw("EP");
  // hSplitTrackCloneHitDistrRatio->SetLineColor(kYellow);
  // hSplitTrackCloneHitDistrRatio->Draw("SAME");
  // hSplitTrackClone_1MissedHitDistrRatio->SetLineColor(kOrange);
  // hSplitTrackClone_1MissedHitDistrRatio->Draw("SAME");
  // hSplitTrackClone_2MissedHitDistrRatio->SetLineColor(kOrange - 7);
  // hSplitTrackClone_2MissedHitDistrRatio->Draw("SAME");
  // hTotalTripletCloneHitDistrRatio->SetLineColor(kBlue);
  // hTotalTripletCloneHitDistrRatio->Draw("SAME");
  // hModuleOverlapCloneHitDistrRatio->SetLineColor(kGreen);
  // hModuleOverlapCloneHitDistrRatio->Draw("SAME");
  // hOtherCloneHitDistrRatio->SetLineColor(kMagenta);
  // hOtherCloneHitDistrRatio->Draw("SAME");

  canvas2->cd(3);
  // Same as pad 2: for low-clone-rate datasets the labels run 0.005..0.035,
  // which is too wide for the default left margin, clipping the y axis title.
  gPad->SetLeftMargin(0.15);
  clonesMCSum->Draw("HIST");
  mcCloneLegend->Draw();
  canvas2->cd(4);
  gPad->SetLeftMargin(0.15);
  clonesRecoSum->Draw("HIST");
  recoCloneLegend->Draw();

  
  TString cloneDistrOutputBase = outputBase + TString("clone_distributions");
  canvas2->SaveAs(cloneDistrOutputBase + suffix + ".pdf");

  // make a canvas to check individual clone types more in detail
  TCanvas* canvas3 = new TCanvas("c3", "Clone Types", 1200, 1000);
  canvas3->Divide(2, 2);
  
  canvas3->cd(1);
  // check seeding clones with a stack
  THStack* seedingCloneStack = new THStack("seedingCloneStack", "Seeding Clone Types (MC) by #eta;#eta;Clone Rate");
  h_seedingClonesMCByEta->SetFillColor(kRed);
  h_seedingClones2MCByEta->SetFillColor(kBlue);
  h_seedingClones3MCByEta->SetFillColor(kGreen);
  h_seedingClones4MCByEta->SetFillColor(kMagenta);
  seedingCloneStack->Add(h_seedingClonesMCByEta);
  seedingCloneStack->Add(h_seedingClones2MCByEta);
  seedingCloneStack->Add(h_seedingClones3MCByEta);
  seedingCloneStack->Add(h_seedingClones4MCByEta);
  seedingCloneStack->Draw("HIST");
  // Add legend for seeding clone stack
  TLegend* seedingCloneLegend = new TLegend(0.35, 0.3, 0.55, 0.5);
  seedingCloneLegend->AddEntry(h_seedingClonesMCByEta, "1 -M- 1", "f");
  seedingCloneLegend->AddEntry(h_seedingClones2MCByEta, "#leq2 -M- #leq2", "f");
  seedingCloneLegend->AddEntry(h_seedingClones3MCByEta, "#leq3 -M- #leq3", "f");
  seedingCloneLegend->AddEntry(h_seedingClones4MCByEta, "#leq4 -M- #leq4", "f");
  seedingCloneLegend->Draw();

  canvas3->cd(2);
  // now check reco seeding clones 
  THStack* seedingCloneRecoStack = new THStack("seedingCloneRecoStack", "Seeding Clone Types (Reco) by #eta;#eta;Clone Rate");
  h_seedingClonesRecoByEta->SetFillColor(kRed);
  h_seedingClones2RecoByEta->SetFillColor(kBlue);
  h_seedingClones3RecoByEta->SetFillColor(kGreen);
  h_seedingClones4RecoByEta->SetFillColor(kMagenta);
  seedingCloneRecoStack->Add(h_seedingClonesRecoByEta);
  seedingCloneRecoStack->Add(h_seedingClones2RecoByEta);
  seedingCloneRecoStack->Add(h_seedingClones3RecoByEta);
  seedingCloneRecoStack->Add(h_seedingClones4RecoByEta);
  seedingCloneRecoStack->Draw("HIST");
  // Add legend for seeding clone stack
  TLegend* seedingCloneRecoLegend = new TLegend(0.35, 0.3, 0.55, 0.5);
  seedingCloneRecoLegend->AddEntry(h_seedingClonesRecoByEta, "1 -M- 1", "f");
  seedingCloneRecoLegend->AddEntry(h_seedingClones2RecoByEta, "#leq2 -M- #leq2", "f");
  seedingCloneRecoLegend->AddEntry(h_seedingClones3RecoByEta, "#leq3 -M- #leq3", "f");
  seedingCloneRecoLegend->AddEntry(h_seedingClones4RecoByEta, "#leq4 -M- #leq4", "f");
  seedingCloneRecoLegend->Draw();

  canvas3->cd(3);
  // check triplet clones with a stack
  THStack* tripletCloneStack = new THStack("tripletCloneStack", "Triplet Clone Types (MC) by #eta;#eta;Clone Rate");
  h_tripletClonesMCByEta->SetFillColor(kRed);
  h_tripletClonesPlusMCByEta->SetFillColor(kBlue);
  tripletCloneStack->Add(h_tripletClonesMCByEta);
  tripletCloneStack->Add(h_tripletClonesPlusMCByEta);
  tripletCloneStack->Draw("HIST");
  // Add legend for triplet clone stack
  TLegend* tripletCloneLegend = new TLegend(0.35, 0.3, 0.55, 0.5);
  tripletCloneLegend->AddEntry(h_tripletClonesMCByEta, "All Triplets", "f");
  tripletCloneLegend->AddEntry(h_tripletClonesPlusMCByEta, "One not triplet", "f");
  tripletCloneLegend->Draw();

  canvas3->cd(4);
  // now do the same with triplets for reco
  THStack* tripletCloneRecoStack = new THStack("tripletCloneRecoStack", "Triplet Clone Types (Reco) by #eta;#eta;Clone Rate");
  h_tripletClonesRecoByEta->SetFillColor(kRed);
  h_tripletClonesPlusRecoByEta->SetFillColor(kBlue);
  tripletCloneRecoStack->Add(h_tripletClonesRecoByEta);
  tripletCloneRecoStack->Add(h_tripletClonesPlusRecoByEta);
  tripletCloneRecoStack->Draw("HIST");
  TLegend* tripletCloneRecoLegend = new TLegend(0.35, 0.3, 0.55, 0.5);
  tripletCloneRecoLegend->AddEntry(h_tripletClonesRecoByEta, "All Triplets", "f");
  tripletCloneRecoLegend->AddEntry(h_tripletClonesPlusRecoByEta, "One not triplet", "f");
  tripletCloneRecoLegend->Draw();

  TString cloneTypesOutputBase = outputBase + TString("clone_types_detailed");
  canvas3->SaveAs(cloneTypesOutputBase + suffix + ".pdf");

  // clean up
  file->Close();
  delete file;
  delete canvas2;
  delete mcCloneLegend;
  delete recoCloneLegend;
  delete cloneDistrLegend;
}