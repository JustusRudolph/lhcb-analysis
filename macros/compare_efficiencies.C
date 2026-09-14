#include <TFile.h>
#include <TString.h>

#include <iostream>
#include <vector>

#include "utils/definitions.h"
#include "utils/basic_functions.h"
#include "utils/efficiency_plots.h"

/*
 * Compare the efficiencies, ghost and clone rates of two or more datasets written by
 * get_eff_clone_ghosts.C. histNames are the part after mc_hists_ of the files inside
 * hists/eff_clone_ghosts, labels are what ends up in the legend and default to those names.
 * One PDF is written per region. For a single dataset use plot_efficiencies.C instead.
 *
 * e.g. root -l 'compare_efficiencies.C({"5000ev_20nm_0ps", "5000ev_80000nm_0ps"},
 *                                      {"20nm", "80000nm"})'
 */
void compare_efficiencies(std::vector<TString> histNames, std::vector<TString> labels={},
                          TString outName="") {
  gROOT->SetBatch();  // so stuff isn't autoplotted
  gStyle->SetOptStat(0);  // remove the info box
  if (histNames.size() < 2) {
    std::cerr << "Error: give at least two hist names to compare." << std::endl;
    return;
  }
  std::vector<TFile*> files{};
  std::vector<TString> usedLabels{};
  for (unsigned i = 0; i < histNames.size(); i++) {
    TFile* file = Utils::EfficiencyPlots::openHistFile(histNames[i]);
    if (!file) continue;
    files.push_back(file);
    // fall back to the file name if no (or not enough) labels were given
    usedLabels.push_back(i < labels.size() ? labels[i] : histNames[i]);
  }
  if (files.size() < 2) {
    std::cerr << "Error: less than two of the given hists could be opened." << std::endl;
    return;
  }
  if (outName.IsNull()) {
    outName = histNames[0];
    for (unsigned i = 1; i < histNames.size(); i++) outName += "_vs_" + histNames[i];
  }
  TString outputBase =
    TString((Utils::Definitions::analysisRoot + "/output/efficiency/").c_str()) + outName;

  Utils::EfficiencyPlots::drawAndSave(files, usedLabels, outputBase);

  for (TFile* file : files) {
    file->Close();
    delete file;
  }
}
