#include <TFile.h>
#include <TH1F.h>
#include <TString.h>

#include "utils/definitions.h"
#include "utils/basic_functions.h"
#include "utils/Hit.h"

void plot_split_track_clones(unsigned nEvents=5000, unsigned max_scatter=80000, unsigned max_dt=0) {
  gROOT->SetBatch();  // so stuff isn't autoplotted
  gStyle->SetOptStat(0); // remove the info box
  gStyle->SetAxisMaxDigits(3);
  TString suffix = Utils::Functions::get_suffix(nEvents, max_scatter, max_dt);
  TString input_suffix = suffix + ".root";
  TFile* file = TFile::Open(
    TString((Utils::Definitions::analysisRoot + "hists/clones/split_track_hists").c_str()) +
    input_suffix);
  if (!file || file->IsZombie()) {
    std::cerr << "Error: Could not open ROOT file." << std::endl;
    return;
  }

  TH1D* h_deltaPhi_splitTrackClone = (TH1D*) file->Get("delta_phi");
  TH1D* h_deltaPhi_splitTrackClone_1Missed = (TH1D*) file->Get("delta_phi_1Missed");
  TH1D* h_deltaPhi_splitTrackClone_2Missed = (TH1D*) file->Get("delta_phi_2Missed");
  TH1D* h_deltaPhi_reference = (TH1D*) file->Get("delta_phi_reference");
  TH1D* h_deflection_splitTrackClone = (TH1D*) file->Get("deflection");
  TH1D* h_deflection_splitTrackClone_1Missed = (TH1D*) file->Get("deflection_1Missed");
  TH1D* h_deflection_splitTrackClone_2Missed = (TH1D*) file->Get("deflection_2Missed");
  TH1D* h_deflection_reference = (TH1D*) file->Get("deflection_reference");
  TH1D* h_deflection_per_z_splitTrackClone = (TH1D*) file->Get("deflection_per_z");
  TH1D* h_deflection_per_z_splitTrackClone_1Missed = (TH1D*) file->Get("deflection_per_z_1Missed");
  TH1D* h_deflection_per_z_splitTrackClone_2Missed = (TH1D*) file->Get("deflection_per_z_2Missed");
  TH1D* h_deflection_per_z_reference = (TH1D*) file->Get("deflection_per_z_reference");
  TH1D* h_deflection_per_z_sq_splitTrackClone = (TH1D*) file->Get("deflection_per_z_sq");
  TH1D* h_deflection_per_z_sq_splitTrackClone_1Missed = (TH1D*) file->Get("deflection_per_z_sq_1Missed");
  TH1D* h_deflection_per_z_sq_splitTrackClone_2Missed = (TH1D*) file->Get("deflection_per_z_sq_2Missed");
  TH1D* h_deflection_per_z_sq_reference = (TH1D*) file->Get("deflection_per_z_sq_reference");
  TH1D* h_pT_reference = (TH1D*) file->Get("pT_reference");
  TH1D* h_pT_splitTrackClone = (TH1D*) file->Get("pT");
  TH1D* h_pT_splitTrackClone_1Missed = (TH1D*) file->Get("pT_1Missed");
  TH1D* h_pT_splitTrackClone_2Missed = (TH1D*) file->Get("pT_2Missed");
  TH2D* h_phi_pT_splitTrackClone = (TH2D*) file->Get("delta_phi_pT");
  TH2D* h_phi_pT_splitTrackClone_1Missed = (TH2D*) file->Get("delta_phi_pT_1Missed");
  TH2D* h_phi_pT_splitTrackClone_2Missed = (TH2D*) file->Get("delta_phi_pT_2Missed");
  
  // Make stacks for delta phi, deflection and pT distributions
  THStack* deltaPhiStack = new THStack(
    "deltaPhiStack", "#Delta#Phi Distributions;#Delta#Phi (rad);Counts");
  h_deltaPhi_splitTrackClone->SetFillColor(kYellow);
  h_deltaPhi_splitTrackClone_1Missed->SetFillColor(kOrange);
  h_deltaPhi_splitTrackClone_2Missed->SetFillColor(kOrange - 7);
  deltaPhiStack->Add(h_deltaPhi_splitTrackClone);
  deltaPhiStack->Add(h_deltaPhi_splitTrackClone_1Missed);
  deltaPhiStack->Add(h_deltaPhi_splitTrackClone_2Missed);

  THStack* deflectionStack = new THStack(
    "deflectionStack", "Deflection in xy (squared);#Deltar^{2} (mm^{2});Counts");
  h_deflection_splitTrackClone->SetFillColor(kYellow);
  h_deflection_splitTrackClone_1Missed->SetFillColor(kOrange);
  h_deflection_splitTrackClone_2Missed->SetFillColor(kOrange - 7);
  deflectionStack->Add(h_deflection_splitTrackClone);
  deflectionStack->Add(h_deflection_splitTrackClone_1Missed);
  deflectionStack->Add(h_deflection_splitTrackClone_2Missed);

  THStack* deflectionPerZStack = new THStack(
    "deflectionPerZStack", "Deflection by unit z;#Deltar^{2} / z (mm);Counts");
  h_deflection_per_z_splitTrackClone->SetFillColor(kYellow);
  h_deflection_per_z_splitTrackClone_1Missed->SetFillColor(kOrange);
  h_deflection_per_z_splitTrackClone_2Missed->SetFillColor(kOrange - 7);
  deflectionPerZStack->Add(h_deflection_per_z_splitTrackClone);
  deflectionPerZStack->Add(h_deflection_per_z_splitTrackClone_1Missed);
  deflectionPerZStack->Add(h_deflection_per_z_splitTrackClone_2Missed);

  THStack* deflectionPerZSqStack = new THStack(
    "deflectionPerZSqStack", "Deflection by unit z^{2};#Deltar^{2} / z^{2} (unitless);Counts");
  h_deflection_per_z_sq_splitTrackClone->SetFillColor(kYellow);
  h_deflection_per_z_sq_splitTrackClone_1Missed->SetFillColor(kOrange);
  h_deflection_per_z_sq_splitTrackClone_2Missed->SetFillColor(kOrange - 7);
  deflectionPerZSqStack->Add(h_deflection_per_z_sq_splitTrackClone);
  deflectionPerZSqStack->Add(h_deflection_per_z_sq_splitTrackClone_1Missed);
  deflectionPerZSqStack->Add(h_deflection_per_z_sq_splitTrackClone_2Missed);

  // Reference distributions are drawn as red points on top of the stacks, scaled to the
  // summed height of the stack so that the shapes can be compared directly.
  auto scaleReferenceToStack = [](TH1D* reference, THStack* stack) {
    // last entry of the stack is the cumulative sum of all its histograms
    TH1* stackSum = (TH1*) stack->GetStack()->Last();
    double maxReference = reference->GetMaximum();
    if (maxReference > 0) reference->Scale(stackSum->GetMaximum() / maxReference);
    reference->SetLineColor(kRed);
    reference->SetMarkerColor(kRed);
    reference->SetMarkerStyle(20);
    reference->SetMarkerSize(0.5);
  };
  scaleReferenceToStack(h_deltaPhi_reference, deltaPhiStack);
  scaleReferenceToStack(h_deflection_reference, deflectionStack);
  scaleReferenceToStack(h_deflection_per_z_reference, deflectionPerZStack);
  scaleReferenceToStack(h_deflection_per_z_sq_reference, deflectionPerZSqStack);

  THStack* pTStack = new THStack("pTStack", "pT Distributions");
  h_pT_splitTrackClone->SetFillColor(kYellow);
  h_pT_splitTrackClone_1Missed->SetFillColor(kOrange);
  h_pT_splitTrackClone_2Missed->SetFillColor(kOrange - 7);
  pTStack->Add(h_pT_splitTrackClone);
  pTStack->Add(h_pT_splitTrackClone_1Missed);
  pTStack->Add(h_pT_splitTrackClone_2Missed);

  unsigned maxValpTReference = h_pT_reference->GetMaximum();
  unsigned maxValpTSplitTrack = h_pT_splitTrackClone->GetMaximum();
  unsigned maxValpTSplitTrack_1Missed = h_pT_splitTrackClone_1Missed->GetMaximum();
  unsigned maxValpTSplitTrack_2Missed = h_pT_splitTrackClone_2Missed->GetMaximum();
  unsigned maxValpTSplit = std::max(maxValpTSplitTrack, std::max(maxValpTSplitTrack_1Missed, maxValpTSplitTrack_2Missed));
  float reference_ratio =  float(maxValpTSplit) / float(maxValpTReference);
  h_pT_reference->Scale(reference_ratio);
  h_pT_reference->SetLineColor(kRed);
  // This will need a special legend
  TLegend* pTLegend = new TLegend(0.6, 0.6, 0.8, 0.8);
  pTLegend->AddEntry(h_pT_splitTrackClone, "LO Split Track", "f");
  pTLegend->AddEntry(h_pT_splitTrackClone_1Missed, "NLO Split Track", "f");
  pTLegend->AddEntry(h_pT_splitTrackClone_2Missed, "NNLO Split Track", "f");
  pTLegend->AddEntry(h_pT_reference, "pT Reference", "l");

  // Make one legend for all (will be the same)
  TLegend* splitTrackLegend = new TLegend(0.6, 0.6, 0.8, 0.8);
  splitTrackLegend->AddEntry(h_deltaPhi_splitTrackClone, "LO Split Track", "f");
  splitTrackLegend->AddEntry(h_deltaPhi_splitTrackClone_1Missed, "NLO Split Track", "f");
  splitTrackLegend->AddEntry(h_deltaPhi_splitTrackClone_2Missed, "NNLO Split Track", "f");
  splitTrackLegend->AddEntry(h_deltaPhi_reference, "Reference (scaled)", "p");

  // Create a canvas to draw the histograms
  TCanvas* canvas_1d = new TCanvas("canvas_1d", "Split Track Clones", 800, 800);
  canvas_1d->Divide(2, 2);
  canvas_1d->cd(1);  // (0,0)
  deltaPhiStack->Draw("HIST");
  h_deltaPhi_reference->Draw("P SAME");
  splitTrackLegend->Draw();
  canvas_1d->cd(2);  // (0,1)
  deflectionStack->Draw("HIST");
  h_deflection_reference->Draw("P SAME");
  splitTrackLegend->Draw();
  canvas_1d->cd(3);  // (1,0)
  deflectionPerZStack->Draw("HIST");
  // deflectionPerZStack->GetXaxis()->SetMaxDigits(1);
  h_deflection_per_z_reference->Draw("P SAME");
  splitTrackLegend->Draw();
  canvas_1d->cd(4);  // (1,1)
  deflectionPerZSqStack->Draw("HIST");
  h_deflection_per_z_sq_reference->Draw("P SAME");
  splitTrackLegend->Draw();


  // Save the canvas_1d as a PDF
  TString outputBase =
    (Utils::Definitions::analysisRoot + "/output/clones/split_track_clones").c_str();
  canvas_1d->SaveAs(outputBase + "_1D" + suffix + ".pdf");
  // Clean up canvas
  delete canvas_1d;

  TCanvas* canvas_2d = new TCanvas("canvas_2d", "Split Track Clones 2D", 800, 800);
  canvas_2d->Divide(2, 2);
  canvas_2d->cd(1);  // (0,0)
  h_phi_pT_splitTrackClone->Draw("COLZ");
  canvas_2d->cd(2);  // (0,1)
  h_phi_pT_splitTrackClone_1Missed->Draw("COLZ");
  canvas_2d->cd(3);  // (1,0)
  h_phi_pT_splitTrackClone_2Missed->Draw("COLZ");
  canvas_2d->cd(4);  // (1,1)
  pTStack->Draw("HIST");
  h_pT_reference->Draw("SAME");
  pTLegend->Draw();

  // Save the canvas_2d as a PDF
  canvas_2d->SaveAs(outputBase + "_2D" + suffix + ".pdf");

  // clean up everything
  delete canvas_2d;
  delete deltaPhiStack;
  delete pTStack;
  delete deflectionStack;
  delete splitTrackLegend;
  delete file;
}