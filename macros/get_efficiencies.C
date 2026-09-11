#include <TEfficiency.h>
#include <TFile.h>
#include <TH1F.h>
#include <TProfile.h>
#include <TString.h>
#include <TSystem.h>

#include <iostream>
#include <string>
#include <vector>

#include "utils/definitions.h"
#include "utils/basic_functions.h"

// variables the efficiency is calculated against, with their axis titles
const std::vector<std::string> kEfficiencyTypes = {"Pt", "Eta", "docaz"};
const std::vector<std::string> kEfficiencyAxes = {"p_{T} (MeV)", "#eta", "DOCA_{z} (mm)"};

TH1F* get_prchecker_hist(TFile* file, const std::string& name) {
  TH1F* hist = (TH1F*) file->Get(name.c_str());
  if (!hist) {
    std::cerr << "Error: Hist " << name << " not found in " << file->GetName() << std::endl;
  }
  return hist;
}

/*
 * Build the efficiency of one variable in one region out of the PrChecker histograms.
 * Returns nullptr if any input is missing or the numerator/denominator do not match.
 */
TEfficiency* build_efficiency(TFile* file, unsigned typeIndex, bool isForward,
                              bool withElectrons) {
  const std::string& type = kEfficiencyTypes[typeIndex];
  std::string region = (isForward ? "forward" : "backward");
  std::string base_string_region =
    std::string("velo_kalman_validator/VeloTracks_eta") + (isForward ? "25" : "-5-2");

  TH1F* nonElectronsReconstructible =
    get_prchecker_hist(file, base_string_region + "_notElectrons_" + type + "_reconstructible");
  TH1F* electronsReconstructible =
    get_prchecker_hist(file, base_string_region + "_electrons_" + type + "_reconstructible");
  TH1F* nonElectronsReconstructed =
    get_prchecker_hist(file, base_string_region + "_notElectrons_" + type + "_reconstructed");
  TH1F* electronsReconstructed =
    get_prchecker_hist(file, base_string_region + "_electrons_" + type + "_reconstructed");
  if (!nonElectronsReconstructible || !electronsReconstructible ||
      !nonElectronsReconstructed || !electronsReconstructed) return nullptr;

  // Add the two numerators and denominators
  TH1F* reconstructedTotal =
    (TH1F*) nonElectronsReconstructed->Clone(("reconstructedTotal_" + region + type).c_str());
  reconstructedTotal->SetDirectory(nullptr);  // only an intermediate, don't write it
  if (withElectrons) reconstructedTotal->Add(electronsReconstructed);

  TH1F* reconstructibleTotal =
    (TH1F*) nonElectronsReconstructible->Clone(("reconstructibleTotal_" + region + type).c_str());
  reconstructibleTotal->SetDirectory(nullptr);
  if (withElectrons) reconstructibleTotal->Add(electronsReconstructible);

  if (!TEfficiency::CheckConsistency(*reconstructedTotal, *reconstructibleTotal)) {
    std::cerr << "Histograms are not consistent for TEfficiency calculation (" << region
              << ", " << type << ")." << std::endl;
    delete reconstructedTotal;
    delete reconstructibleTotal;
    return nullptr;
  }
  TEfficiency* eff = new TEfficiency(*reconstructedTotal, *reconstructibleTotal);
  // name is what the comparison macro looks for, keep in sync with compare_efficiencies.C
  eff->SetName(("efficiency_" + region + "_" + type).c_str());
  std::string title = "Efficiency vs " + kEfficiencyAxes[typeIndex];
  eff->SetTitle((title + ";" + kEfficiencyAxes[typeIndex] + ";Efficiency").c_str());

  delete reconstructedTotal;
  delete reconstructibleTotal;
  return eff;
}

/*
 * Write the efficiencies (and ghost rates) of one dataset into hists/efficiency, so that
 * compare_efficiencies.C can overlay several of them. Give outName to choose the file name
 * the comparison macro will be pointed at, it defaults to the usual data suffix.
 */
void get_efficiencies(unsigned nEvents=5000, unsigned max_scatter=80000,
                      unsigned max_dt=0, bool withElectrons=false, TString outName="") {
  gROOT->SetBatch();  // so stuff isn't autoplotted
  TString suffix = Utils::Functions::get_suffix(nEvents, max_scatter, max_dt);
  TString input_suffix = suffix + ".root";
  TFile* file = TFile::Open(
    TString(Utils::Definitions::stackRoot + "/output/PrCheckerPlots") + input_suffix);
  if (!file || file->IsZombie()) {
    std::cerr << "Error: Could not open ROOT file." << std::endl;
    return;
  }
  // ghost and clone rates are calculated separately by get_clone_rates.C, which writes them
  // into the mc_hists file. They are optional here.
  std::vector<std::string> rateNames = {"ghost_rates", "clone_rate", "clone_rate_pt"};
  std::vector<TProfile*> rates(rateNames.size(), nullptr);
  TString mcHistPath =
    TString(Utils::Definitions::analysisRoot + "/hists/clones/mc_hists") + input_suffix;
  TFile* mc_histo_file = TFile::Open(mcHistPath);
  if (!mc_histo_file || mc_histo_file->IsZombie()) {
    std::cerr << "Warning: " << mcHistPath << " not found, skipping ghost and clone rates."
              << "\nRun get_clone_rates.C to generate them first." << std::endl;
  } else {
    for (unsigned i = 0; i < rateNames.size(); i++) {
      rates[i] = (TProfile*) mc_histo_file->Get(rateNames[i].c_str());
      if (!rates[i]) std::cerr << "Warning: no " << rateNames[i] << " in " << mcHistPath << ".\n";
    }
  }

  if (outName.IsNull()) outName = TString("efficiencies") + suffix;
  TString outDir = (Utils::Definitions::analysisRoot + "/hists/efficiency").c_str();
  gSystem->mkdir(outDir, true);  // TFile does not create the directory itself
  TString outPath = outDir + "/" + outName + ".root";
  TFile* outFile = new TFile(outPath, "RECREATE");
  outFile->cd();

  // both forward and backward region into the same file
  for (bool isForward : {true, false}) {
    for (unsigned typeIndex = 0; typeIndex < kEfficiencyTypes.size(); typeIndex++) {
      TEfficiency* eff = build_efficiency(file, typeIndex, isForward, withElectrons);
      if (!eff) continue;
      eff->Write();
      delete eff;
    }
  }
  // the eta rates cover both regions, they only differ in the drawn x range
  for (unsigned i = 0; i < rateNames.size(); i++) {
    if (rates[i]) rates[i]->Write(rateNames[i].c_str());
  }
  std::cout << "Wrote efficiencies to " << outPath << std::endl;

  outFile->Close();
  delete outFile;
  if (mc_histo_file) {
    mc_histo_file->Close();
    delete mc_histo_file;
  }
  file->Close();
  delete file;
}
