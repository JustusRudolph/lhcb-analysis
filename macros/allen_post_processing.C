#include <TFile.h>
#include <TH1F.h>
#include <TEfficiency.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <string>

std::string stackRoot = std::getenv("STACK_ROOT");
std::string analysisRoot = std::getenv("ANALYSIS_ROOT");

void eff_pur_plots(bool withElectrons=true) {
  // Open ROOT file and retrieve histograms
  TFile *file = TFile::Open((stackRoot + "/output/PrCheckerPlots.root").c_str());
  if (!file || file->IsZombie()) {
      std::cerr << "Error: Could not open ROOT file." << std::endl;
      return;
  }

  // set up canvases for forward and backward region, and types of plots
  std::vector<std::string> types = {"Pt", "Eta", "docaz"};
  std::vector<std::string> x_axes = {"pT (MeV)", "eta", "docaz (mm)"};
  std::string base_string = "velo_kalman_validator/VeloTracks_eta";
  TCanvas *canvas = new TCanvas("c", "Efficiencies", 1000, 800);
  canvas->Divide(2, 3);
  int c_idx = 1;


  for (int isForward = 0; isForward < 2; isForward++) {
    for (int typeIndex = 0; typeIndex < types.size(); typeIndex++) {
      std::string type = types[typeIndex];
      std::string base_string_region = (isForward ? base_string + "25" : base_string + "-5-2");
      
      TH1F *nonElectronsReconstructible =
        (TH1F*)file->Get((base_string_region + "_notElectrons_" + type + "_reconstructible").c_str());
      TH1F *electronsReconstructible =
        (TH1F*)file->Get((base_string_region + "_electrons_" + type + "_reconstructible").c_str());
      TH1F *nonElectronsReconstructed =
        (TH1F*)file->Get((base_string_region + "_notElectrons_" + type + "_reconstructed").c_str());
      TH1F *electronsReconstructed =
        (TH1F*)file->Get((base_string_region + "_electrons_" + type + "_reconstructed").c_str());

      // Check if histograms are properly loaded
      if (!nonElectronsReconstructible) {
        std::cerr << "Error: Hist "
                  << (base_string_region + "_notElectrons_" + type + "_reconstructible").c_str()
                  << "not found in root file ../../output/PrCheckerPlots.root" << std::endl;
        file->Close();
        return;
      } else if (!electronsReconstructible) {
        std::cerr << "Error: Hist "
                  << (base_string_region + "_electrons_" + type + "_reconstructible").c_str()
                  << "not found in root file ../../output/PrCheckerPlots.root" << std::endl;
        file->Close();
        return;
      } else if (!nonElectronsReconstructed) {
        std::cerr << "Error: Hist "
                  << (base_string_region + "_notElectrons_" + type + "_reconstructed").c_str()
                  << "not found in root file ../../output/PrCheckerPlots.root" << std::endl;
        file->Close();
        return;
      } else if (!electronsReconstructed) {
        std::cerr << "Error: Hist "
                  << (base_string_region + "_electrons_" + type + "_reconstructed").c_str()
                  << "not found in root file ../../output/PrCheckerPlots.root" << std::endl;
        file->Close();
        return;
      }

      // Add the two numerators and denominators
      TH1F *reconstructedTotal = (TH1F*)nonElectronsReconstructed->Clone("reconstructedTotal");
      if (withElectrons) reconstructedTotal->Add(electronsReconstructed);
      
      TH1F *reconstructibleTotal = (TH1F*)nonElectronsReconstructible->Clone("reconstructibleTotal");
      if (withElectrons) reconstructibleTotal->Add(electronsReconstructible);

      // Create TEfficiency object
      if (TEfficiency::CheckConsistency(*reconstructedTotal, *reconstructibleTotal)) {
        TEfficiency *eff = new TEfficiency(*reconstructedTotal, *reconstructibleTotal);
        eff->SetLineColor(kMagenta);
        eff->SetMarkerColor(kMagenta);
        // Draw the efficiency plot
        canvas->cd(c_idx);
        std::string region = (isForward ? "forward" : "backward");
        std::string title = "Efficiency in " + region + " region wrt " + type;
        eff->SetTitle((title + ";" + x_axes[typeIndex] + ";Efficiency").c_str());
        eff->Draw("AP");

        gPad->Update();
        // for pT and eta, cut at 0.8 to have the same
        // TODO: Change to c_idx < 5 for electrons too when better?
        if (withElectrons) {
          if (c_idx < 3 || c_idx == 4) {  // top two are pT
            eff->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.8, 1.01);
            // graph->GetYaxis()->SetNdivisions(10, "N");  // steps of 0.02
          }
        } else {
          if (c_idx < 5) {  // top two are pT
            eff->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.8, 1.01);
            // graph->GetYaxis()->SetNdivisions(10, "N");  // steps of 0.02
          }
        }
        c_idx += 2;  // move to next row
      } else {
        std::cerr << "Histograms are not consistent for TEfficiency calculation." << std::endl;
      }
    }  // types
    c_idx = 2;  // move to second column
  }  // isForward
  canvas->SaveAs((analysisRoot + "/output/efficiency_plots.pdf").c_str());
  // Clean up
  file->Close();
  delete file;
  delete canvas;
}

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
  float mcEta;
  std::vector<unsigned>* matchedRecoTrackIndices = nullptr;
  std::vector<unsigned>* recoTrackLHCbIDs = nullptr;
  // MC Event Tree
  // MC Track Tree
  mcTrackTree->SetBranchAddress("matchedTracks", &matchedRecoTrackIndices);
  mcTrackTree->SetBranchAddress("nMatches", &nMatches);
  mcTrackTree->SetBranchAddress("evNo", &mcTrackEvNo);
  mcTrackTree->SetBranchAddress("runNo", &mcTrackRunNo);
  mcTrackTree->SetBranchAddress("eta", &mcEta);
  mcTrackTree->SetBranchAddress("nHitsVelo", &nMCVeloHits);
  // Reco Event Tree
  recoEventTree->SetBranchAddress("globalTrackOffset", &recoEvOffset);
  // Reco Track Tree
  recoTrackTree->SetBranchAddress("lhcbid", &recoTrackLHCbIDs);

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
  TProfile* hLongestTrackTagged = new TProfile(
    "longest_match_tag_rate", "Longest match tag rate;#eta;Longest Track Tag Rate",
    nBins, etaBinEdges.data());
  TProfile* hNMatches = new TProfile(
    "number_of_matches_of_MC_track", "Number of Matches of MC Track;#eta;N_{matches}",
    nBins, etaBinEdges.data());

  for (unsigned i_mct = 0; i_mct < nMCTracks; i_mct++) {
    mcTrackTree->GetEntry(i_mct);
    if (!nMatches || !matchedRecoTrackIndices) continue;
    std::unordered_set<unsigned> uniqueIDs;
    unsigned longestTrackSize = 0;
    unsigned kTotalIDs = 0;
    unsigned longestTrackTagged = 1;  // represent bool as int for histogram
    for (unsigned matchIdx : *matchedRecoTrackIndices) {
      unsigned evIdx = (mcTrackRunNo - 1) * kEventsPerRun + (mcTrackEvNo - 1);
      recoEventTree->GetEntry(evIdx);
      unsigned recoTrackIdx = recoEvOffset + matchIdx;
      recoTrackTree->GetEntry(recoTrackIdx);
      if (!recoTrackLHCbIDs) {
        std::cout << "Nothing found in reco track tree for lhcbid at index " << recoTrackIdx << ".\n";
        continue;  // if nothing there, this is weird...
      }
      kTotalIDs += recoTrackLHCbIDs->size();
      if (recoTrackLHCbIDs->size() > longestTrackSize) {
        // if longest track size not zero, that means it has already been changed
        // hence a later track is bigger, which shouldn't be
        if (longestTrackSize) longestTrackTagged = 0;
        longestTrackSize = recoTrackLHCbIDs->size();
      }
      for (unsigned lhcbid : *recoTrackLHCbIDs) {
        uniqueIDs.insert(lhcbid);
      }
    }
    float idDuplRate = 1. - (float) uniqueIDs.size() / (float) kTotalIDs;
    float idMatchRate = (float) uniqueIDs.size() / (float) nMCVeloHits;
    float longestTrackRate = (float) longestTrackSize / (float) nMCVeloHits;

    hDuplicateIDRates->Fill(mcEta, idDuplRate);
    hUniqueIDRates->Fill(mcEta, idMatchRate);

    if (nMatches > 5) {  // at least 5 clones
      hLongestMatchedTrackRate5Clones->Fill(mcEta, longestTrackRate);
    } else if (nMatches > 1) {  // at least one clone
      hLongestMatchedTrackRateClones->Fill(mcEta, longestTrackRate);
    } else {  // will be exactly one match now because already checked not zero
      hLongestMatchedTrackRate->Fill(mcEta, longestTrackRate);
    }

    hLongestTrackTagged->Fill(mcEta, longestTrackTagged);
    hNMatches->Fill(mcEta, nMatches);

    if (i_mct && (i_mct % 20000 == 0))
      printf("Finished with %d %% of MC Tracks.\n", (int) (100. * (float) i_mct / nMCTracks));
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
  // Clean up outfile
  outFile->Close();
  delete outFile;

  // clean up file
  file->Close();
  delete file;
}

void allen_post_processing(bool withElectrons=true) {
  gROOT->SetBatch();
  eff_pur_plots(withElectrons);
  get_clone_rates();
}

int main() {
    allen_post_processing();
    return 0;
}