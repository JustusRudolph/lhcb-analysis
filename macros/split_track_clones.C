#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <THStack.h>
#include <TLegend.h>
#include <TMath.h>
#include <iostream>
#include <unordered_set>
#include <vector>

#include "utils/MCDataAnalysis.h"
#include "utils/definitions.h"
#include "utils/basic_functions.h"
#include "utils/Hit.h"

/*
 * Split track clones: an MC track that has been reconstructed as several reco
 * tracks, at leading order without losing a single hit (LO), or with one/two
 * hits missing or extra (NLO/NNLO).
 */
class SplitTrackClones : public Utils::MCDataAnalysis {
 public:
  using Utils::MCDataAnalysis::MCDataAnalysis;

  void book() override {
    auto phiBins = binEdges(nPhiBins, -0.1, 0.1);  // short range for granularity
    auto ptBins = binEdges(nPTBins, 0., 5000.);  // only up to TeV
    auto deflectionBins = binEdges(nDeflectionBins, 0., 0.5);
    // deflection bins scaled by dz, roughly 40mm per module pair but scaled
    // more for larger displaced ones
    std::vector<float> deflectionBinsDZScaled(nDeflectionBins + 1);
    std::vector<float> deflectionBinsDZSQScaled(nDeflectionBins + 1);
    for (unsigned i = 0; i <= nDeflectionBins; i++) {
      deflectionBinsDZScaled[i] = deflectionBins[i] / 100.;
      deflectionBinsDZSQScaled[i] = deflectionBins[i] / (100. * 100.);
    }
    // one set of histograms per order, "" is LO, _1Missed NLO, _2Missed NNLO
    for (unsigned order = 0; order < kOrders; order++) {
      TString tag = orderTag(order);
      bindHistogram<TH1D>("delta_phi" + tag, "#Delta#phi Distribution;#eta;#Delta#phi",
                 nPhiBins, phiBins.data());
      bindHistogram<TH1D>("deflection" + tag, "Deflection Distribution;#eta;Deflection",
                 nDeflectionBins, deflectionBins.data());
      bindHistogram<TH1D>("deflection_per_z" + tag,
                 "Deflection per z Distribution;#eta;Deflection per z",
                 nDeflectionBins, deflectionBinsDZScaled.data());
      bindHistogram<TH1D>("deflection_per_z_sq" + tag,
                 "Deflection per z^2 Distribution;#eta;Deflection per z^2",
                 nDeflectionBins, deflectionBinsDZSQScaled.data());
      bindHistogram<TH1D>("pT" + tag, "pT Distribution;#eta;pT", nPTBins, ptBins.data());
      bindHistogram<TH2D>("delta_phi_pT" + tag,
                 ("pT vs phi (" + orderName(order) + " Split Track);#phi;pT").Data(),
                 nPhiBins, phiBins.data(), nPTBins, ptBins.data());
    }
    bindHistogram<TH1D>("pT_reference", "pT Distribution;#eta;pT", nPTBins, ptBins.data());
  }

  void fill() override {
    if (!openInput()) return;
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
    // MC Track Tree
    mcTracks()->SetBranchAddress("runNo", &mcTrackRunNo);
    mcTracks()->SetBranchAddress("evNo", &mcTrackEvNo);
    mcTracks()->SetBranchAddress("nHitsVelo", &nMCVeloHits);
    mcTracks()->SetBranchAddress("nMatches", &nMatches);
    mcTracks()->SetBranchAddress("pt", &mcPT);
    mcTracks()->SetBranchAddress("matchedTracks", &matchedRecoTrackIndices);
    mcTracks()->SetBranchAddress("lhcbid", &mcTrackLHCbIDs);
    mcTracks()->SetBranchAddress("x", &x_mc);
    mcTracks()->SetBranchAddress("y", &y_mc);
    // Reco Trees
    recoEvents()->SetBranchAddress("globalTrackOffset", &recoEvOffset);
    recoTracks()->SetBranchAddress("lhcbid", &recoTrackLHCbIDs);
    recoTracks()->SetBranchAddress("mcMatchIdx", &mcMatchIdxForReco);
    recoTracks()->SetBranchAddress("eta", &recoEta);
    recoTracks()->SetBranchAddress("x", &x_reco);
    recoTracks()->SetBranchAddress("y", &y_reco);

    // Go through MC particles
    unsigned nMCTracks = mcTracks()->GetEntries();
    printf("Going through %u MC particles now.\n", nMCTracks);
    for (unsigned i_mct = 0; i_mct < nMCTracks; i_mct++) {
      printProgress(i_mct, nMCTracks, "MC Tracks");

      mcTracks()->GetEntry(i_mct);
      getHistFromFile<TH1D>("pT_reference")->Fill(mcPT);  // fill always, not just with clones
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
        recoEvents()->GetEntry(recoEventEntry(mcTrackRunNo, mcTrackEvNo));
        unsigned recoTrackIdx = recoEvOffset + matchIdx;
        recoTracks()->GetEntry(recoTrackIdx);
        if (!recoTrackLHCbIDs) {
          std::cout << "Nothing found in reco track tree for lhcbid at index "
                    << recoTrackIdx << ".\n";
          continue;  // if nothing there, this is weird...
        }
        unsigned nRecoHits = recoTrackLHCbIDs->size();
        kTotalIDs += nRecoHits;
        for (unsigned i_hit_reco = 0; i_hit_reco < nRecoHits; i_hit_reco++) {
          unsigned reco_lhcbid = recoTrackLHCbIDs->at(i_hit_reco);
          mcTrackLHCbIDsSet.erase(reco_lhcbid);
          unsigned moduleNumber = (reco_lhcbid >> 12) & 0x3F;
          float z = moduleZ(moduleNumber);
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
              // track moving forwards, trace out detector backwards (this is what we should always do)
              // jump by steps of 2 since dealing with module pairs
              dz_to_next_module = moduleZ(moduleNumber - 2) - t0_penultimateHit.z;
            } else {
              // detector being traced forwards
              dz_to_next_module = moduleZ(moduleNumber + 2) - t0_penultimateHit.z;
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
      unsigned nUnmatched = mcTrackLHCbIDsSet.size();
      // all MC hits reconstructed within all reco matched tracks
      bool isSplitTrack = nUnmatched == 0 && kTotalIDs == nMCVeloHits;
      // Either one hit missing between, one overlap, or one extra hit (covered by case 2)
      bool isSplitTrack_1Missed = (nUnmatched == 1 && kTotalIDs == (nMCVeloHits - 1)) ||
                                  (nUnmatched == 0 && kTotalIDs == (nMCVeloHits + 1));
      // Two hits missing, two overlap (or two extra), or one missing and one extra
      bool isSplitTrack_2Missed = (nUnmatched == 2 && kTotalIDs == (nMCVeloHits - 2)) ||
                                  (nUnmatched == 0 && kTotalIDs == (nMCVeloHits + 2)) ||
                                  (nUnmatched == 1 && kTotalIDs == (nMCVeloHits + 1));

      // if not split track, continue (nMatches > 1 checked above already)
      if (!(isSplitTrack || isSplitTrack_1Missed || isSplitTrack_2Missed)) continue;
      unsigned order = isSplitTrack ? 0 : (isSplitTrack_1Missed ? 1 : 2);
      TString tag = orderTag(order);

      float dr_sq = dx * dx + dy * dy;
      float dr_sq_per_z = dr_sq / abs(dz_to_next_module);
      float dr_sq_per_z_sq = dr_sq / (dz_to_next_module * dz_to_next_module);
      // Fill histograms
      getHistFromFile<TH1D>("delta_phi" + tag)->Fill(dPhi);
      getHistFromFile<TH1D>("deflection" + tag)->Fill(dr_sq);
      getHistFromFile<TH1D>("deflection_per_z" + tag)->Fill(dr_sq_per_z);
      getHistFromFile<TH1D>("deflection_per_z_sq" + tag)->Fill(dr_sq_per_z_sq);
      getHistFromFile<TH1D>("pT" + tag)->Fill(mcPT);
      getHistFromFile<TH2D>("delta_phi_pT" + tag)->Fill(dPhi, mcPT);
    }  // MC Particles loop
  }

  void plot() override {
    // Make stacks for delta phi, deflection and pT distributions
    THStack* deltaPhiStack = stack(
      "delta_phi", "#Delta#Phi Distributions;#Delta#Phi (rad);Counts");
    THStack* deflectionStack = stack(
      "deflection", "Deflection in xy (squared);#Deltar^{2} (mm^{2});Counts");
    THStack* deflectionPerZStack = stack(
      "deflection_per_z", "Deflection by unit z;#Deltar^{2} / z (mm);Counts");
    THStack* deflectionPerZSqStack = stack(
      "deflection_per_z_sq", "Deflection by unit z^{2};#Deltar^{2} / z^{2} (unitless);Counts");
    THStack* pTStack = stack("pT", "pT Distributions");

    // Make one legend for all (will be the same)
    TLegend* splitTrackLegend = new TLegend(0.6, 0.6, 0.8, 0.8);
    for (unsigned order = 0; order < kOrders; order++)
      splitTrackLegend->AddEntry(getHistFromFile<TH1D>("delta_phi" + orderTag(order)),
                                 (orderName(order) + " Split Track").Data(), "f");

    // Create a canvas to draw the histograms
    TCanvas* canvas_1d = new TCanvas("canvas_1d", "Split Track Clones", 800, 800);
    canvas_1d->Divide(2, 2);
    canvas_1d->cd(1);  // (0,0)
    deltaPhiStack->Draw("HIST");
    splitTrackLegend->Draw();
    canvas_1d->cd(2);  // (0,1)
    deflectionStack->Draw("HIST");
    splitTrackLegend->Draw();
    canvas_1d->cd(3);  // (1,0)
    deflectionPerZStack->Draw("HIST");
    splitTrackLegend->Draw();
    canvas_1d->cd(4);  // (1,1)
    deflectionPerZSqStack->Draw("HIST");
    splitTrackLegend->Draw();
    saveCanvas(canvas_1d, fPlotBase + "_1D.pdf");
    delete canvas_1d;

    // scale the reference down to the largest split track distribution
    TH1D* h_pT_reference = getHistFromFile<TH1D>("pT_reference");
    double maxValpTSplit = 0.;
    for (unsigned order = 0; order < kOrders; order++)
      maxValpTSplit = std::max(maxValpTSplit, getHistFromFile<TH1D>("pT" + orderTag(order))->GetMaximum());
    h_pT_reference->Scale(maxValpTSplit / h_pT_reference->GetMaximum());
    h_pT_reference->SetLineColor(kRed);
    // This will need a special legend
    TLegend* pTLegend = new TLegend(0.6, 0.6, 0.8, 0.8);
    for (unsigned order = 0; order < kOrders; order++)
      pTLegend->AddEntry(getHistFromFile<TH1D>("pT" + orderTag(order)),
                         (orderName(order) + " Split Track").Data(), "f");
    pTLegend->AddEntry(h_pT_reference, "pT Reference", "l");

    TCanvas* canvas_2d = new TCanvas("canvas_2d", "Split Track Clones 2D", 800, 800);
    canvas_2d->Divide(2, 2);
    for (unsigned order = 0; order < kOrders; order++) {
      canvas_2d->cd(order + 1);
      getHistFromFile<TH2D>("delta_phi_pT" + orderTag(order))->Draw("COLZ");
    }
    canvas_2d->cd(4);  // (1,1)
    pTStack->Draw("HIST");
    h_pT_reference->Draw("SAME");
    pTLegend->Draw();
    saveCanvas(canvas_2d, fPlotBase + "_2D.pdf");

    // clean up, the histograms themselves belong to the base class
    delete canvas_2d;
    delete deltaPhiStack;
    delete deflectionStack;
    delete deflectionPerZStack;
    delete deflectionPerZSqStack;
    delete pTStack;
    delete splitTrackLegend;
    delete pTLegend;
  }

 public:
  // Plots are saved as <base>_1D.pdf and <base>_2D.pdf.
  void setPlotBase(TString base) { fPlotBase = base; }

 private:
  TString fPlotBase;
  static const unsigned kOrders = 3;
  static const unsigned nPhiBins = 50, nPTBins = 50, nDeflectionBins = 50;

  static TString orderTag(unsigned order) {
    const char* tags[kOrders] = {"", "_1Missed", "_2Missed"};
    return tags[order];
  }
  static TString orderName(unsigned order) {
    const char* names[kOrders] = {"LO", "NLO", "NNLO"};
    return names[order];
  }

  // Stack the three orders of a distribution on top of each other.
  THStack* stack(TString name, TString title) {
    const int colours[kOrders] = {kYellow, kOrange, kOrange - 7};
    THStack* stack = new THStack((name + "Stack").Data(), title.Data());
    for (unsigned order = 0; order < kOrders; order++) {
      TH1D* hist = getHistFromFile<TH1D>(name + orderTag(order));
      hist->SetFillColor(colours[order]);
      stack->Add(hist);
    }
    return stack;
  }
};

/*
 * mode: kAuto reuses the histogram file if it is there, kRefill always runs
 * over the data again, kPlotOnly only redraws, kFillOnly skips the plots.
 */
void split_track_clones(unsigned nEvents=5000, unsigned max_scatter=80000, unsigned max_dt=0,
                        Utils::AnalysisBase::Mode mode=Utils::AnalysisBase::kAuto) {
  TString suffix = Utils::Functions::get_suffix(nEvents, max_scatter, max_dt);
  TString stackRoot = Utils::Definitions::stackRoot.c_str();
  TString analysisRoot = Utils::Definitions::analysisRoot.c_str();

  // TODO: writing to clones_test/ for now so the existing files are left
  // alone, change back to clones/ once this is confirmed to work
  SplitTrackClones analysis(analysisRoot + "hists/clones_test/split_track_hists" + suffix + ".root",
                            stackRoot + "output/MCData_Checking" + suffix + ".root");
  analysis.setPlotBase(analysisRoot + "output/clones_test/split_track_clones" + suffix);
  analysis.run(mode);
}
