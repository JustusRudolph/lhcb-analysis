#include <TFile.h>
#include <TH1F.h>
#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>

#include "utils/definitions.h"
#include "utils/basic_functions.h"
#include "utils/allen_functions.h"
#include "utils/Hit.h"

void check_ghosts(unsigned nEvents=5000, unsigned max_scatter=80,
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
  // Trees that will be relevant
  TTree* recoTrackTree = (TTree*) file->Get("RecoTrackData");
  if (!recoTrackTree) {
      std::cerr << "RecoTrackData TTree not found!" << std::endl;
      return;
  }
  TTree* mcTrackTree = (TTree*) file->Get("MCTrackData");
  if (!mcTrackTree) {
      std::cerr << "MCTrackData TTree not found!" << std::endl;
      return;
  }
  TTree* mcEventTree = (TTree*) file->Get("MCEventData");
  if (!mcEventTree) {
      std::cerr << "MCEventData TTree not found!" << std::endl;
      return;
  }
  std::vector<double> moduleToZ = Utils::Functions::get_module_to_z_vector();
  // ----------- BINS ---------------
  // eta bins
  int nEtaBins = 50;
  std::vector<float> etaBinEdges(nEtaBins + 1);
  float etaMax{6.}, etaMin{-5};
  float etaStep = (etaMax - etaMin) / (nEtaBins);
  for (unsigned i = 0; i <= nEtaBins; i++) {
    etaBinEdges[i] = etaMin + etaStep * i;
  }
  // dr bins (deflection)
  int nDeflectionBins = 50;
  std::vector<float> deflectionBinEdges(nDeflectionBins + 1);
  float deflectionMax{0.5}, deflectionMin{0.};
  float deflectionStep = (deflectionMax - deflectionMin) / (nDeflectionBins);
  for (unsigned i = 0; i <= nDeflectionBins; i++) {
    deflectionBinEdges[i] = deflectionMin + deflectionStep * i;
  }
  // dt bins (also deflection but in time)
  int nDTBins = 50;
  std::vector<float> dtBinEdges(nDTBins + 1);
  std::vector<float> dtSqBinEdges(nDTBins + 1);
  std::vector<float> chiSqDTBinEdges(nDTBins + 1);  // ps^2 per degree of freedom (nHits - 2)
  float dtMax{500}, dtMin{-500.};  // in ps for visualisation
  float dtSqMax{dtMax * dtMax}, dtSqMin{0.};  // in ps^2 for visualisation
  float chiSqMax{50}, chiSqMin{0.};
  float dtStep = (dtMax - dtMin) / (nDTBins);
  float dtSqStep = (dtSqMax - dtSqMin) / (nDTBins);
  float chiSqStep = (chiSqMax - chiSqMin) / (nDTBins);
  for (unsigned i = 0; i <= nDTBins; i++) {
    dtBinEdges[i] = dtMin + dtStep * i;
    dtSqBinEdges[i] = dtSqMin + dtSqStep * i;
    chiSqDTBinEdges[i] = chiSqMin + chiSqStep * i;
  }
  // deflection bins scaled by dz
  std::vector<float> deflectionBinEdgesDZScaled(nDeflectionBins + 1);
  std::vector<float> deflectionBinEdgesDZSQScaled(nDeflectionBins + 1);
  for (unsigned i = 0; i <= nDeflectionBins; i++) {
    // roughly 40mm per module pair but scaled more for larger displaced ones
    deflectionBinEdgesDZScaled[i] = deflectionBinEdges[i] / 100.;
    deflectionBinEdgesDZSQScaled[i] = deflectionBinEdges[i] / (100. * 100.);
  }
  // hit bins
  int nHitMax = 15, nHitMin = 3;
  int nHitBins = nHitMax - nHitMin + 1;
  float epsilon = 1e-6;
  std::vector<float> hitBinEdges(nHitBins+1);
  for (unsigned i = 0; i <= nHitBins; i++) {
    hitBinEdges[i] = nHitMin + i - epsilon;  // move low bin edge left to avoid floating point errors
  }
  // fraction of MC sizes bins
  int nFractionBins = 50;
  std::vector<float> fractionBinEdges(nFractionBins + 1);
  float fractionMax{1.}, fractionMin{0.};  // should end at 0.7 since 70% is enough to be counted as matched
  float fractionStep = (fractionMax - fractionMin) / (nFractionBins);
  for (unsigned i = 0; i <= nFractionBins; i++) {
    fractionBinEdges[i] = fractionMin + fractionStep * i;
  }
  // ----------------- BRANCH DATA ----------------
  int mcMatchIdx;  // to check if fake or not
  unsigned nHits, recoTrackRunNo, recoTrackEvNo, nMatches, mcEvOffset;
  float recoEta;
  std::vector<float>* x_reco = nullptr;
  std::vector<float>* y_reco = nullptr;
  std::vector<float>* z_reco = nullptr;
  std::vector<float>* t_reco = nullptr;
  std::vector<int>* lhcbID_reco = nullptr;
  std::vector<unsigned>* hit_mc_indices = nullptr; // MC particle of each hit
  // Set addresses
  recoTrackTree->SetBranchAddress("runNo", &recoTrackRunNo);
  recoTrackTree->SetBranchAddress("evNo", &recoTrackEvNo);
  recoTrackTree->SetBranchAddress("lhcbid", &lhcbID_reco);
  recoTrackTree->SetBranchAddress("mcMatchIdx", &mcMatchIdx);
  recoTrackTree->SetBranchAddress("nHits", &nHits);
  recoTrackTree->SetBranchAddress("eta", &recoEta);
  recoTrackTree->SetBranchAddress("x", &x_reco);
  recoTrackTree->SetBranchAddress("y", &y_reco);
  recoTrackTree->SetBranchAddress("z", &z_reco);
  recoTrackTree->SetBranchAddress("t", &t_reco);
  recoTrackTree->SetBranchAddress("mcIndices", &hit_mc_indices);
  // MC Track Tree
  mcTrackTree->SetBranchAddress("nMatches", &nMatches);
  // MC Event Tree
  mcEventTree->SetBranchAddress("globalTrackOffset", &mcEvOffset);

  // Here we want to go through specifically the split tracks
  TH1D* h_dt_0_matched = new TH1D("h_dt_0_matched", "time scatter of hits of matched tracks (0th order);#Delta t [ps];Counts",
                                nDTBins, dtBinEdges.data());
  TH1D* h_dt_0_fake = new TH1D("h_dt_0_fake", "time scatter of hits of fake tracks (0th order);#Delta t [ps];Counts",
                                nDTBins, dtBinEdges.data());
  TH1D* h_dt_1_matched = new TH1D("h_dt_1_matched", "time scatter of hits of matched tracks (1st order);#Delta t [ps];Counts",
                                nDTBins, dtBinEdges.data());
  TH1D* h_dt_1_fake = new TH1D("h_dt_1_fake", "time scatter of hits of fake tracks (1st order);#Delta t [ps];Counts",
                                nDTBins, dtBinEdges.data());
  TH1D* h_dt_seeds_matched = new TH1D("h_dt_seeds_matched", "time scatter of 0th & 2nd hits in seeds of matched tracks;#Delta t [ps];Counts",
                                nDTBins, dtBinEdges.data());
  TH1D* h_dt_seeds_fake = new TH1D("h_dt_seeds_fake", "time scatter of 0th & 2nd hits in seeds of fake tracks;#Delta t [ps];Counts",
                                nDTBins, dtBinEdges.data());
  TH1D* h_dt_seeds_total_matched = new TH1D("h_dt_seeds_total_matched",
    "time scatter of 0th to 1st & 2nd hits in seeds of matched tracks (#chi^2_{t});#Delta t^2 [ps^2];Counts",
    nDTBins, dtSqBinEdges.data());
  TH1D* h_dt_seeds_total_fake = new TH1D("h_dt_seeds_total_fake",
    "time scatter of 0th to 1st & 2nd hits in seeds of fake tracks (#chi^2_{t});#Delta t^2 [ps^2];Counts",
    nDTBins, dtSqBinEdges.data());
  TH1D* h_dt_seeds_total_scaled_matched = new TH1D("h_dt_seeds_total_scaled_matched",
    "time scatter of 0th to 1st & 2nd hits in seeds of matched tracks (#chi^2_{t});#Delta t^2 [ps^2];Counts",
    nDTBins, dtSqBinEdges.data());
  TH1D* h_dt_seeds_total_scaled_fake = new TH1D("h_dt_seeds_total_scaled_fake",
    "time scatter of 0th to 1st & 2nd hits in seeds of fake tracks (#chi^2_{t});#Delta t^2 [ps^2];Counts",
    nDTBins, dtSqBinEdges.data());
  TH1D* h_chi_sq_matched = new TH1D("h_chi_sq_matched", "chi^2 of matched tracks;#chi^{2};Counts",
                                    nDTBins, chiSqDTBinEdges.data());
  TH1D* h_chi_sq_fake = new TH1D("h_chi_sq_fake", "chi^2 of fake tracks;#chi^{2};Counts",
                                 nDTBins, chiSqDTBinEdges.data());
  TH1D* h_dr_sq_dz_sq_matched = new TH1D("h_dr_sq_dz_sq_matched", "deflection of matched tracks (wrt z^{2});#Delta r^{2}/z^{2} [unitless];Counts",
                                          nDeflectionBins, deflectionBinEdgesDZSQScaled.data());
  TH1D* h_dr_sq_dz_sq_fake = new TH1D("h_dr_sq_dz_dq_fake", "deflection of fake tracks (wrt z^{2});#Delta r^{2}/z^{2} [unitless];Counts",
                                          nDeflectionBins, deflectionBinEdgesDZSQScaled.data());
  TH1D* h_dr_sq_dz_matched = new TH1D("h_dr_sq_dz_matched", "deflection of matched tracks;#Delta r^{2}/z [mm];Counts",
                                      nDeflectionBins, deflectionBinEdgesDZScaled.data());
  TH1D* h_dr_sq_dz_fake = new TH1D("h_dr_sq_dz_fake", "deflection of fake tracks;#Delta r^{2}/z [mm];Counts",
                                   nDeflectionBins, deflectionBinEdgesDZScaled.data());
  TH1D* h_eta_matched = new TH1D("h_eta_matched", "eta of matched tracks;#eta;Counts",
                                nEtaBins, etaBinEdges.data());
  TH1D* h_eta_fake = new TH1D("h_eta_fake", "eta of fake tracks;#eta;Counts",
                                nEtaBins, etaBinEdges.data());
  TH1D* h_size_matched = new TH1D("h_size_matched", "size of matched tracks;size;Counts",
                                  nHitBins, hitBinEdges.data());
  TH1D* h_size_fake = new TH1D("h_size_fake", "size of fake tracks;size;Counts",
                               nHitBins, hitBinEdges.data());
  TH1D* h_number_mcs_matched = new TH1D("h_number_mcs_matched", "Number of MC particles used for matched track;N_{MCs};Counts",
                                        10, 1, 11);
  TH1D* h_number_mcs_fake = new TH1D("h_number_mcs_fake", "Number of MC particles used for fake track;N_{MCs};Counts",
                                      10, 1, 11);
  TH1D* h_number_mc_splits_matched = new TH1D("h_number_mc_splits_matched", "N_{splits}+1 in reco track;N_{MCs};Counts",
                                              10, 1, 11);
  TH1D* h_number_mc_splits_fake = new TH1D("h_number_mc_splits_fake", "N_{splits}+1 in reco track;N_{MCs};Counts",
                                           10, 1, 11);
  TProfile* p_recod_hits_fake_by_size = new TProfile("p_recod_hits_fake_by_size", "Fraction of hits in fake tracks belonging to non-recod MCs;N_{reco}/N_{Hits};N_{Hits}",
                                                     nHitBins, hitBinEdges.data());
  // where are the splits in the tracks? The first two are the sums, then 2, 3, 4, 5 are the numbers of splits + 1
  TH1D* h_split_index_fake = new TH1D("h_split_index_fake", "Distribution of where fake tracks split;index in track;Counts",
                                      nHitMax, 1.-epsilon, nHitMax-epsilon);  // define split to occur when the next is different, i.e. start at 1
  TH1D* h_split_index_matched = new TH1D("h_split_index_matched", "Distribution of where matched tracks split;index in track;Counts",
                                         nHitMax, 1.-epsilon, nHitMax-epsilon);
  TH1D* h_split_index_fake_2 = new TH1D("h_split_index_fake_2", "Distribution of where fake tracks split (2 Splits);index in track;Counts",
                                        nHitMax, 1.-epsilon, nHitMax-epsilon);
  TH1D* h_split_index_matched_2 = new TH1D("h_split_index_matched_2", "Distribution of where matched tracks split (2 Splits);index in track;Counts",
                                         nHitMax, 1.-epsilon, nHitMax-epsilon);
  TH1D* h_split_index_fake_3 = new TH1D("h_split_index_fake_3", "Distribution of where fake tracks split (3 Splits);index in track;Counts",
                                        nHitMax, 1.-epsilon, nHitMax-epsilon);
  TH1D* h_split_index_matched_3 = new TH1D("h_split_index_matched_3", "Distribution of where matched tracks split (3 Splits);index in track;Counts",
                                           nHitMax, 1.-epsilon, nHitMax-epsilon);
  TH1D* h_split_index_fake_4 = new TH1D("h_split_index_fake_4", "Distribution of where fake tracks split (4 Splits);index in track;Counts",
                                        nHitMax, 1.-epsilon, nHitMax-epsilon);
  TH1D* h_split_index_fake_5 = new TH1D("h_split_index_fake_5", "Distribution of where fake tracks split (5 Splits);index in track;Counts",
                                        nHitMax, 1.-epsilon, nHitMax-epsilon);
  TH1D* h_dt_splits_matched = new TH1D("h_dt_splits_matched", "time scatter of split hits of matched tracks (0th order);#Delta t [ps];Counts",
                                       nDTBins, dtBinEdges.data());
  TH1D* h_dt_splits_fake = new TH1D("h_dt_splits_fake", "time scatter of split hits of fake tracks (0th order);#Delta t [ps];Counts",
                                    nDTBins, dtBinEdges.data());
  TH1D* h_fraction_mc_sizes_2 = new TH1D("h_fraction_mc_sizes_2", "Distribution of MC hit usage (normalised), 2 MCs;n_{hits}^{MC} / n_{hits}^{track};Counts",
                                          nFractionBins, fractionBinEdges.data());
  TH1D* h_fraction_mc_sizes_3 = new TH1D("h_fraction_mc_sizes_3", "Distribution of MC hit usage (normalised), 3 MCs;n_{hits}^{MC} / n_{hits}^{track};Counts",
                                          nFractionBins, fractionBinEdges.data());
  TH1D* h_fraction_mc_sizes_4 = new TH1D("h_fraction_mc_sizes_4", "Distribution of MC hit usage (normalised), 4 MCs;n_{hits}^{MC} / n_{hits}^{track};Counts",
                                          nFractionBins, fractionBinEdges.data());

  // Go through reco tracks
  unsigned nFakes = 0;
  unsigned nRecoTracks = recoTrackTree->GetEntries();
  printf("Going through %u reconstructed tracks now.\n", nRecoTracks);
  for (unsigned i_reco = 0; i_reco < nRecoTracks; i_reco++) {
    // print progress
    if (i_reco && (i_reco % 1000000 == 0))
      printf("Finished with %d%% of tracks.\n", (int) (100. * (float) i_reco / nRecoTracks));
    recoTrackTree->GetEntry(i_reco);
    // fill sizes & eta immediately since on track basis
    if (mcMatchIdx == -1) {
      nFakes++;
      h_eta_fake->Fill(recoEta);
      if (nHits > nHitMax)
        h_size_fake->Fill(nHitMax);
      else
        h_size_fake->Fill(nHits);
    } else {
      h_eta_matched->Fill(recoEta);
      if (nHits > nHitMax)
        h_size_matched->Fill(nHitMax);
      else
        h_size_matched->Fill(nHits);
    }

    Hit::BaseHit h0, h1;  // changed throughout the track
    float time_scatter_0;  // 0th order extrapolation
    float time_scatter_1;  // 1st order extrapolation
    float dr_sq_scatter_dz_sq;
    float dr_sq_scatter_dz;

    // create unordered map of unique mc indices in track to frequency
    std::unordered_map<unsigned, unsigned> unique_mc_indices;
    unsigned nMatchedHits = 0;
    // check for number of splits
    unsigned current_mc_index{ hit_mc_indices->at(0) }, nMCSplits{0};
    std::vector<unsigned> track_split_indices;
    for (unsigned i_hit = 0; i_hit < nHits; i_hit++) {
      if (i_hit == 0) {
        h0 = Hit::BaseHit(lhcbID_reco->at(i_hit), x_reco->at(i_hit),
                          y_reco->at(i_hit), z_reco->at(i_hit), t_reco->at(i_hit));
      } else if (i_hit != (nHits - 1)) {  // not at the end or start, can still extrapolate to next layer
        h1 = Hit::BaseHit(lhcbID_reco->at(i_hit), x_reco->at(i_hit),
                          y_reco->at(i_hit), z_reco->at(i_hit), t_reco->at(i_hit));
        // estimate the absolute time difference to next hit
        float next_module_z = Utils::Functions::get_next_module_z(h0, h1, moduleToZ);
        float abs_dt_estimate_0 = Utils::Functions::estimate_abs_dt(h0, h1, next_module_z, 0);
        float abs_dt_estimate_1 = Utils::Functions::estimate_abs_dt(h0, h1, next_module_z, 1);
        // check measured time difference from midpoint of h0 and h1
        float abs_dt = std::abs( t_reco->at(i_hit + 1) - (h0.t + h1.t) / 2 );
        time_scatter_0 = abs_dt - abs_dt_estimate_0;
        time_scatter_1 = abs_dt - abs_dt_estimate_1;

        // fill the split histograms
        if (hit_mc_indices->at(i_hit) != hit_mc_indices->at(i_hit + 1)) {
          if (mcMatchIdx == -1) {
            h_dt_splits_fake->Fill(1000 * time_scatter_0); // in picoseconds
          } else {
            h_dt_splits_matched->Fill(1000 * time_scatter_0);
          }
        }
        // and check if mc index was reconstructed
        unsigned evIdx = (recoTrackRunNo - 1) * kEventsPerRun + (recoTrackEvNo - 1);
        mcEventTree->GetEntry(evIdx);
        unsigned mcTrackIdx = mcEvOffset + hit_mc_indices->at(i_hit);
        mcTrackTree->GetEntry(mcTrackIdx);
        nMatchedHits += (nMatches > 0);

        // estimate the deflection (dr^2) to next hit
        std::tuple<float, float, float> pos_estimate =
          Utils::Functions::estimateNextPhi(h0, h1, next_module_z);
        float dx = std::get<0>(pos_estimate) - x_reco->at(i_hit + 1);
        float dy = std::get<1>(pos_estimate) - y_reco->at(i_hit + 1);
        float dr_sq_scatter = dx * dx + dy * dy;
        float dz_abs = std::abs(next_module_z - h0.z);
        dr_sq_scatter_dz_sq = dr_sq_scatter / (dz_abs * dz_abs);
        dr_sq_scatter_dz = dr_sq_scatter / dz_abs;
        h0 = Hit::BaseHit(h1); // update h0 to be h1 for next iteration
      }
      if (i_hit == 2) {  // check the seed
        // first get difference in r between the two layers (not just z)
        float dx = x_reco->at(i_hit) - x_reco->at(i_hit - 2);
        float dy = y_reco->at(i_hit) - y_reco->at(i_hit - 2);
        float dz = z_reco->at(i_hit) - z_reco->at(i_hit - 2);
        float dr_sq = dx * dx + dy * dy + dz * dz;
        float dt_sq = dr_sq / (Utils::Definitions::c_f_float_mmns * Utils::Definitions::c_f_float_mmns);
        float dt_estimate_seed = std::sqrt(dt_sq);
        float dt_observed_seed = std::abs(t_reco->at(i_hit) - t_reco->at(i_hit - 2));
        float dt_scatter_2_ps = 1000 * (dt_observed_seed - dt_estimate_seed);
        // fill the seed histograms
        if (mcMatchIdx == -1) {
          h_dt_seeds_fake->Fill(dt_scatter_2_ps); // in picoseconds
        } else {
          h_dt_seeds_matched->Fill(dt_scatter_2_ps);
        }
        // now do the same between the 0th and 1st hits
        dx = x_reco->at(i_hit) - x_reco->at(i_hit - 1);
        dy = y_reco->at(i_hit) - y_reco->at(i_hit - 1);
        dz = z_reco->at(i_hit) - z_reco->at(i_hit - 1);
        dr_sq = dx * dx + dy * dy + dz * dz;
        dt_sq = dr_sq / (Utils::Definitions::c_f_float_mmns * Utils::Definitions::c_f_float_mmns);
        dt_estimate_seed = std::sqrt(dt_sq);
        dt_observed_seed = std::abs(t_reco->at(i_hit) - t_reco->at(i_hit - 1));
        float dt_scatter_1_ps = 1000 * (dt_observed_seed - dt_estimate_seed);

        float dt_scatter_sq = dt_scatter_1_ps * dt_scatter_1_ps + dt_scatter_2_ps * dt_scatter_2_ps;
        float scale_factor = 2.0;  // change this how you want (to be tuned)
        float dt_scatter_sq_scaled = dt_scatter_sq;
        if (dt_scatter_1_ps * dt_scatter_2_ps < 0) {
          // if the two time scatterings are in opposite directions, this will typically make
          // the chi_sq better, even though the 1st and 2nd are very far apart (indicating split)
          dt_scatter_sq_scaled *= scale_factor;
        }
        if (mcMatchIdx == -1) {
          h_dt_seeds_total_fake->Fill(dt_scatter_sq);
          h_dt_seeds_total_scaled_fake->Fill(dt_scatter_sq_scaled);
        } else {
          h_dt_seeds_total_matched->Fill(dt_scatter_sq);
          h_dt_seeds_total_scaled_matched->Fill(dt_scatter_sq_scaled);
        }

      }
      // fill the unique mc indices map
      unsigned mcIndex = hit_mc_indices->at(i_hit);
      unique_mc_indices[mcIndex]++;  // [] inserts the key already & initialises to zero
      if (mcIndex != current_mc_index) { // for hit 0 this will guaranteed be the same
        nMCSplits++;
        track_split_indices.push_back(i_hit);
        // fill total split histogram to avoid having to loop over vector later
        if (mcMatchIdx == -1) {
          h_split_index_fake->Fill(i_hit);
        } else {
          h_split_index_matched->Fill(i_hit);
        }
      }
      current_mc_index = mcIndex;

      if (mcMatchIdx == -1) {
        // fake track
        h_dr_sq_dz_fake->Fill(dr_sq_scatter_dz);
        h_dr_sq_dz_sq_fake->Fill(dr_sq_scatter_dz_sq);
        h_dt_0_fake->Fill(1000 * time_scatter_0); // in picoseconds
        h_dt_1_fake->Fill(1000 * time_scatter_1);

      } else {
        // matched track
        h_dr_sq_dz_matched->Fill(dr_sq_scatter_dz);
        h_dr_sq_dz_sq_matched->Fill(dr_sq_scatter_dz_sq);
        h_dt_0_matched->Fill(1000 * time_scatter_0);
        h_dt_1_matched->Fill(1000 * time_scatter_1);
      }
    } // hits in track
    // fill the p_recod_hits_fake_by_size histogram
    if (nHits > nHitMax)
      p_recod_hits_fake_by_size->Fill(nHitMax, (float) nMatchedHits / nHits);
    else
      p_recod_hits_fake_by_size->Fill(nHits, (float) nMatchedHits / nHits);
    // get chi^2 fit for t/z
    TGraphErrors* track_t_z_graph = new TGraphErrors(nHits);
    for (int i = 0; i < nHits; ++i) {
        track_t_z_graph->SetPoint(i, z_reco->at(i), 1000*t_reco->at(i));  // ps/mm
        track_t_z_graph->SetPointError(i, 0.0, 50);  // constant 50ps resolution
    }
    TFitResultPtr fit = track_t_z_graph->Fit("pol1", "SQ");  // linear fit
    double chi2 = fit->Chi2();
    unsigned ndf = fit->Ndf();
    assert(ndf == nHits - 2);  // should just be two parameters, i.e. -2
    if (mcMatchIdx == -1) {
      h_chi_sq_fake->Fill(chi2 / ndf);
    } else {
      h_chi_sq_matched->Fill(chi2 / ndf);
    }
    delete track_t_z_graph;
    
    // I'm aware this is super inefficient with all the -1 checks, but it's nicer
    // to read this instead of all sets of histograms at once (compartmentalise!!)

    // fill number of unique mc indices & fractions, as well as splits
    if (mcMatchIdx == -1) {
      h_number_mc_splits_fake->Fill(nMCSplits + 1);  // to be able to directly compare to number of unique mcs
      h_number_mcs_fake->Fill(unique_mc_indices.size());
      if (unique_mc_indices.size() == 2) {
        for (auto& [mcIdx, nHitsFromMC] : unique_mc_indices) {
          h_fraction_mc_sizes_2->Fill((float) nHitsFromMC / nHits);
        }
      } else if (unique_mc_indices.size() == 3) {
        for (auto& [mcIdx, nHitsFromMC] : unique_mc_indices) {
          h_fraction_mc_sizes_3->Fill((float) nHitsFromMC / nHits);
        }
      } else if (unique_mc_indices.size() == 4) {
        for (auto& [mcIdx, nHitsFromMC] : unique_mc_indices) {
          h_fraction_mc_sizes_4->Fill((float) nHitsFromMC / nHits);
        }
      }
    } else {
      h_number_mcs_matched->Fill(unique_mc_indices.size());
      h_number_mc_splits_matched->Fill(nMCSplits + 1);  // to be able to directly compare to number of unique mcs
    }

    // Now do the indices with 2, 3, 4, 5 MCs
    // first make a vector to access the right histogram
    std::vector<TH1D*> split_histograms_fake = {
      h_split_index_fake_2, h_split_index_fake_3, h_split_index_fake_4, h_split_index_fake_5
    };
    std::vector<TH1D*> split_histograms_matched = { h_split_index_matched_2, h_split_index_matched_3 };
    unsigned hist_vec_index = nMCSplits - 1;  // 1 split is 2 MCs, i.e. index 0
    // make sure that we have at least one split (this should be given for fakes but whatever still added it)
    if (mcMatchIdx == -1  && hist_vec_index < 4 && hist_vec_index >= 0) {
      for (unsigned split_index : track_split_indices) {
        // printf("Filling hist %u with split index %u for fake track with %u splits.\n",
        //        hist_vec_index, split_index, nMCSplits);
        split_histograms_fake[hist_vec_index]->Fill(split_index);
      }
    } else if (mcMatchIdx != -1 && hist_vec_index < 2 && hist_vec_index >= 0) {
      for (unsigned split_index : track_split_indices) {
        // printf("Filling hist %u with split index %u for matched track with %u splits.\n",
        //        hist_vec_index, split_index, nMCSplits);
        split_histograms_matched[hist_vec_index]->Fill(split_index);
      }
    }
    // printf("\tAfter filling:\n");
    for (unsigned int splIdx = 0; splIdx < split_histograms_fake.size(); splIdx++) {
      TH1D* hist = split_histograms_fake[splIdx];
      // printf("\t\tHist %u has %f entries\n", splIdx, hist->GetEntries());
    }
  } // tracks
  printf("Out of %u reconstructed tracks, %u are fakes.\n", nRecoTracks, nFakes);
  
  // Write histograms and clean up
  TString outPrefix =
    (Utils::Definitions::analysisRoot + "/hists/ghosts/basic_hists").c_str();
  TFile* outFile = new TFile(outPrefix + input_suffix, "RECREATE");
  h_dt_0_matched->Write();
  h_dt_0_fake->Write();
  h_dt_1_matched->Write();
  h_dt_1_fake->Write();
  h_dt_seeds_matched->Write();
  h_dt_seeds_fake->Write();
  h_dt_seeds_total_matched->Write();
  h_dt_seeds_total_fake->Write();
  h_dt_seeds_total_scaled_matched->Write();
  h_dt_seeds_total_scaled_fake->Write();
  h_chi_sq_matched->Write();
  h_chi_sq_fake->Write();
  h_dr_sq_dz_matched->Write();
  h_dr_sq_dz_fake->Write();
  h_dr_sq_dz_sq_matched->Write();
  h_dr_sq_dz_sq_fake->Write();
  h_eta_matched->Write();
  h_eta_fake->Write();
  h_size_matched->Write();
  h_size_fake->Write();
  h_split_index_matched->Write();
  h_split_index_fake->Write();
  h_split_index_matched_2->Write();
  h_split_index_fake_2->Write();
  h_split_index_matched_3->Write();
  h_split_index_fake_3->Write();
  h_split_index_fake_4->Write();
  h_split_index_fake_5->Write();
  h_number_mcs_fake->Write();
  h_number_mcs_matched->Write();
  h_number_mc_splits_matched->Write();
  h_number_mc_splits_fake->Write();
  p_recod_hits_fake_by_size->Write();
  h_dt_splits_matched->Write();
  h_dt_splits_fake->Write();
  h_fraction_mc_sizes_2->Write();
  h_fraction_mc_sizes_3->Write();
  h_fraction_mc_sizes_4->Write();

  outFile->Close();
  delete outFile;
  file->Close();
  delete file;
}