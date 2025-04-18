#include <TFile.h>
#include <TH1F.h>

#include "utils/definitions.h"
#include "utils/basic_functions.h"
#include "utils/Hit.h"

void plot_split_track_clones() {
  gROOT->SetBatch();  // so stuff isn't autoplotted
  TFile* file = TFile::Open(
    (Utils::Definitions::analysisRoot + "hists/clones/split_track_hists.root").c_str());
  if (!file || file->IsZombie()) {
    std::cerr << "Error: Could not open ROOT file." << std::endl;
    return;
  }

  TH1D* h_deltaPhi_splitTrackClone = (TH1D*) file->Get("delta_phi");
  TH1D* h_deltaPhi_splitTrackClone_1Missed = (TH1D*) file->Get("delta_phi_1Missed");
  TH1D* h_deltaPhi_splitTrackClone_2Missed = (TH1D*) file->Get("delta_phi_2Missed");
  TH1D* h_pT_splitTrackClone = (TH1D*) file->Get("pT");
  TH1D* h_pT_splitTrackClone_1Missed = (TH1D*) file->Get("pT_1Missed");
  TH1D* h_pT_splitTrackClone_2Missed = (TH1D*) file->Get("pT_2Missed");
  TH2D* h_phi_pT_splitTrackClone = (TH2D*) file->Get("delta_phi_pT");
  TH2D* h_phi_pT_splitTrackClone_1Missed = (TH2D*) file->Get("delta_phi_pT_1Missed");
  TH2D* h_phi_pT_splitTrackClone_2Missed = (TH2D*) file->Get("delta_phi_pT_2Missed");
  
  // Make stacks for delta phi and pT distributions
  THStack* deltaPhiStack = new THStack("deltaPhiStack", "#Delta#Phi Distributions");
  h_deltaPhi_splitTrackClone->SetFillColor(kYellow);
  h_deltaPhi_splitTrackClone_1Missed->SetFillColor(kOrange);
  h_deltaPhi_splitTrackClone_2Missed->SetFillColor(kOrange - 7);
  deltaPhiStack->Add(h_deltaPhi_splitTrackClone);
  deltaPhiStack->Add(h_deltaPhi_splitTrackClone_1Missed);
  deltaPhiStack->Add(h_deltaPhi_splitTrackClone_2Missed);
  THStack* pTStack = new THStack("pTStack", "pT Distributions");
  h_pT_splitTrackClone->SetFillColor(kYellow);
  h_pT_splitTrackClone_1Missed->SetFillColor(kOrange);
  h_pT_splitTrackClone_2Missed->SetFillColor(kOrange - 7);
  pTStack->Add(h_pT_splitTrackClone);
  pTStack->Add(h_pT_splitTrackClone_1Missed);
  pTStack->Add(h_pT_splitTrackClone_2Missed);
  // Make a legend for both (will be the same)
  TLegend* splitTrackLegend = new TLegend(0.6, 0.6, 0.8, 0.8);
  splitTrackLegend->AddEntry(h_deltaPhi_splitTrackClone, "LO Split Track", "f");
  splitTrackLegend->AddEntry(h_deltaPhi_splitTrackClone_1Missed, "NLO Split Track", "f");
  splitTrackLegend->AddEntry(h_deltaPhi_splitTrackClone_2Missed, "NNLO Split Track", "f");
  
  // Create a canvas to draw the histograms
  TCanvas* canvas = new TCanvas("canvas", "Split Track Clones", 800, 1200);
  canvas->Divide(2, 3);
  canvas->cd(1);  // (0,0)
  deltaPhiStack->Draw("HIST");
  splitTrackLegend->Draw();
  canvas->cd(2);  // (0,1)
  pTStack->Draw("HIST");
  splitTrackLegend->Draw();
  canvas->cd(3);  // (1,0)
  h_phi_pT_splitTrackClone->Draw("COLZ");
  canvas->cd(4);  // (1,1)
  h_phi_pT_splitTrackClone_1Missed->Draw("COLZ");
  canvas->cd(5);  // (2,1)
  h_phi_pT_splitTrackClone_2Missed->Draw("COLZ");

  // Save the canvas as a PDF
  canvas->SaveAs(
    (Utils::Definitions::analysisRoot + "/output/split_track_clones.pdf").c_str());
  // Clean up
  delete canvas;
  delete deltaPhiStack;
  delete pTStack;
  delete splitTrackLegend;
  delete file;
}