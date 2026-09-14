#include <TFile.h>
#include <TString.h>
#include <TSystem.h>

#include <iostream>
#include <vector>

#include "utils/definitions.h"
#include "utils/basic_functions.h"
#include "utils/efficiency_plots.h"

/*
 * Draw the efficiencies, ghost and clone rates of a single dataset written by
 * get_eff_clone_ghosts.C. The dataset is picked the same way as there, label is what ends up
 * in the legend and defaults to the data suffix. One PDF is written per region.
 *
 * e.g. root -l 'plot_efficiencies.C(5000, 20, 0)'
 */
// max_dt is in picoseconds and scatter in micrometers
void plot_efficiencies(unsigned nEvents=5000, unsigned max_scatter=80000, unsigned max_dt=0,
                       TString mc_file_suffix="", TString label="", TString outName="") {
  gROOT->SetBatch();  // so stuff isn't autoplotted
  gStyle->SetOptStat(0);  // remove the info box
  TString suffix;
  if (mc_file_suffix.IsNull()) {
    suffix = Utils::Functions::get_suffix(nEvents, max_scatter, max_dt);
  } else {
    suffix = Form("_%uev_%s", nEvents, mc_file_suffix.Data());
  }
  // the file name already has the underscore after mc_hists, so drop the leading one
  TString histName = suffix;
  histName.Remove(0, 1);

  TFile* file = Utils::EfficiencyPlots::openHistFile(histName);
  if (!file) return;
  std::vector<TString> labelToPass{};  // default is empty
  if (!label.IsNull()) labelToPass.push_back(label);
  if (outName.IsNull()) outName = histName;
  TString outputBase =
    TString((Utils::Definitions::analysisRoot + "/output/efficiency/").c_str()) + outName;

  Utils::EfficiencyPlots::drawAndSave({file}, labelToPass, outputBase);

  file->Close();
  delete file;
}
