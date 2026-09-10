// Everything that is specific to the LHCb MC checking files, so that the
// AnalysisBase itself stays usable for other experiments. Derive from this
// instead of from AnalysisBase to get the four trees of an MCData_Checking
// file, the event indexing and the module z positions.
//
// The input file is opened the first time one of the tree accessors is used
// and closed again in the destructor, so a run that only reads the histogram
// file back never touches the data.
#pragma once

#include <TBranch.h>
#include <TFile.h>
#include <TH1.h>
#include <TLeaf.h>
#include <TString.h>
#include <TTree.h>

#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "AnalysisBase.h"
#include "definitions.h"

namespace Utils {

class MCDataAnalysis : public AnalysisBase {
 public:
  MCDataAnalysis(TString histFilePath, TString inputFilePath = "")
      : AnalysisBase(histFilePath), fInputFilePath(inputFilePath) {}

  ~MCDataAnalysis() override { closeInput(); }

  // -------------------------------------------------------------- input trees
  /*
   * Open the MC checking file and grab the four trees. Called on demand, but
   * worth calling explicitly at the top of fill() to bail out early if the
   * file isn't there.
   */
  bool openInput() {
    if (fInFile) return true;
    if (fInputFilePath.IsNull()) {
      std::cerr << "Error: No input file set." << std::endl;
      return false;
    }
    fInFile = TFile::Open(fInputFilePath);
    if (!fInFile || fInFile->IsZombie()) {
      std::cerr << "Error: Could not open ROOT file " << fInputFilePath << "." << std::endl;
      return false;
    }
    fMCTrackTree = getTree("MCTrackData");
    fMCEventTree = getTree("MCEventData");
    fRecoTrackTree = getTree("RecoTrackData");
    fRecoEventTree = getTree("RecoEventData");
    if (!fMCTrackTree || !fMCEventTree || !fRecoTrackTree || !fRecoEventTree) {
      closeInput();
      return false;
    }
    // built before anything sets its own branch addresses
    fMCEventIndex = buildEventIndex(fMCEventTree);
    fRecoEventIndex = buildEventIndex(fRecoEventTree);
    return true;
  }

  void closeInput() {
    if (!fInFile) return;
    fInFile->Close();
    delete fInFile;
    fInFile = nullptr;
    fMCTrackTree = fMCEventTree = fRecoTrackTree = fRecoEventTree = nullptr;
  }

  TFile* inputFile() { openInput(); return fInFile; }
  TTree* mcTracks() { openInput(); return fMCTrackTree; }
  TTree* mcEvents() { openInput(); return fMCEventTree; }
  TTree* recoTracks() { openInput(); return fRecoTrackTree; }
  TTree* recoEvents() { openInput(); return fRecoEventTree; }

  // Entry of an event in the reco/MC event tree, the runs are not in order.
  Long64_t recoEventEntry(unsigned runNo, unsigned evNo) {
    openInput();
    return eventEntry(fRecoEventIndex, runNo, evNo);
  }
  Long64_t mcEventEntry(unsigned runNo, unsigned evNo) {
    openInput();
    return eventEntry(fMCEventIndex, runNo, evNo);
  }

  // z position of a module, 0 if the module doesn't exist (as GetBinContent).
  double moduleZ(int module) {
    if (fModuleToZ.empty()) loadModuleToZ();
    if (module < 0 || module >= (int) fModuleToZ.size()) return 0.;
    return fModuleToZ[module];
  }
  const std::vector<double>& moduleToZ() {
    if (fModuleToZ.empty()) loadModuleToZ();
    return fModuleToZ;
  }

  TString inputFilePath() const { return fInputFilePath; }
  void setInputFilePath(TString path) { closeInput(); fInputFilePath = path; }

  // Only used to guess an event entry if the event tree has no run/event
  // branches to index on, see eventEntry().
  void setEventsPerRun(unsigned eventsPerRun) { fEventsPerRun = eventsPerRun; }

 protected:
  TTree* getTree(const char* name) const {
    TTree* tree = (TTree*) fInFile->Get(name);
    if (!tree) std::cerr << name << " TTree not found!" << std::endl;
    return tree;
  }

 private:
  typedef std::map<std::pair<unsigned, unsigned>, Long64_t> EventIndex;

  /*
   * (runNo, evNo) -> entry, since neither the runs nor the events within them
   * are guaranteed to be in order. Reads only the two branches it needs, so
   * branch addresses set by the macro are left alone.
   */
  static EventIndex buildEventIndex(TTree* tree) {
    EventIndex index;
    TBranch* runBranch = tree->GetBranch("runNo");
    TBranch* evBranch = tree->GetBranch("evNo");
    if (!runBranch || !evBranch) return index;  // fall back to the ordering
    TLeaf* runLeaf = runBranch->GetLeaf("runNo");
    TLeaf* evLeaf = evBranch->GetLeaf("evNo");
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++) {
      runBranch->GetEntry(i);
      evBranch->GetEntry(i);
      index[{(unsigned) runLeaf->GetValue(), (unsigned) evLeaf->GetValue()}] = i;
    }
    return index;
  }

  Long64_t eventEntry(const EventIndex& index, unsigned runNo, unsigned evNo) const {
    auto entry = index.find({runNo, evNo});
    if (entry != index.end()) return entry->second;
    // no run/event branches to index on, assume everything is in order
    return (Long64_t) (runNo - 1) * fEventsPerRun + (evNo - 1);
  }

  void loadModuleToZ() {
    TFile* moduleFile = TFile::Open(
      (Utils::Definitions::analysisRoot + "hists/module_mc_info.root").c_str());
    if (!moduleFile || moduleFile->IsZombie()) {
      std::cerr << "Error: Could not open module info ROOT file." << std::endl;
      return;
    }
    TH1D* moduleToZ = (TH1D*) moduleFile->Get("module_to_z");
    if (moduleToZ) {
      fModuleToZ.resize(moduleToZ->GetNbinsX());
      for (int i = 0; i < moduleToZ->GetNbinsX(); i++)
        fModuleToZ[i] = moduleToZ->GetBinContent(i + 1);  // +1 because ROOT 1 indexed
    }
    moduleFile->Close();
    delete moduleFile;
  }

  TString fInputFilePath;
  unsigned fEventsPerRun{200};

  TFile* fInFile{nullptr};
  TTree* fMCTrackTree{nullptr};
  TTree* fMCEventTree{nullptr};
  TTree* fRecoTrackTree{nullptr};
  TTree* fRecoEventTree{nullptr};
  EventIndex fMCEventIndex, fRecoEventIndex;
  std::vector<double> fModuleToZ;
};

}  // namespace Utils
