#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <string>


void MCDataCheck(unsigned run = 5, unsigned event = 133, unsigned moduleNo = 60,
                 float phiMin = 2.76, float phiMax = 2.764, unsigned kEventsPerRun=200,
                 std::string mcFilePath="MCData_Checking.root") {
  std::string stackRoot = std::getenv("STACK_ROOT");
  std::string analysisRoot = std::getenv("ANALYSIS_ROOT");
  
  TFile *file = TFile::Open((stackRoot + mcFilePath).c_str());
  if (!file || file->IsZombie()) {
      std::cerr << "Error: Could not open ROOT file." << std::endl;
      return;
  }
  // get mc track tree
  TTree* mcTrackTree = (TTree*) file->Get("MCTrackData");
  if (!mcTrackTree) {
      std::cerr << "TTree not found!" << std::endl;
      return;
  }
  // for offsets also get event tree
  TTree* mcEventTree = (TTree*) file->Get("MCEventData");
  if (!mcEventTree) {
      std::cerr << "TTree not found!" << std::endl;
      return;
  }

  // for now, just check pid of the concerned cluster
  unsigned runNo, evNo, nMatches, mcTrackOffset;
  int pid;
  std::vector<unsigned>* lhcbid = nullptr;
  std::vector<float>* x = nullptr;
  std::vector<float>* y = nullptr;
  std::vector<float>* t = nullptr;
  float mcEta;
  
  mcTrackTree->SetBranchAddress("runNo", &runNo);
  mcTrackTree->SetBranchAddress("evNo", &evNo);
  mcTrackTree->SetBranchAddress("nMatches", &nMatches);
  mcTrackTree->SetBranchAddress("lhcbid", &lhcbid);
  mcTrackTree->SetBranchAddress("pid", &pid);
  mcTrackTree->SetBranchAddress("x", &x);
  mcTrackTree->SetBranchAddress("y", &y);
  mcTrackTree->SetBranchAddress("t", &t);
  mcTrackTree->SetBranchAddress("eta", &mcEta);

  mcEventTree->SetBranchAddress("globalTrackOffset", &mcTrackOffset);

  // get start and end points of tracks to check for faster access
  unsigned trackEndIdx = mcTrackTree->GetEntries();
  unsigned evIdx = (run - 1) * kEventsPerRun + (event - 1);
  mcEventTree->GetEntry(evIdx);
  unsigned startIdx = mcTrackOffset;
  if (evIdx != mcEventTree->GetEntries()) {
    // if not the last event, get the next one for the right offset
    mcEventTree->GetEntry(evIdx + 1);
    trackEndIdx = mcTrackOffset;
  }

  for (unsigned i_track = startIdx; i_track < trackEndIdx; i_track++) {
    mcTrackTree->GetEntry(i_track);
    // skip if not the requested run and event (shouldn't matter anyway because of offset)
    if (runNo != run || evNo != event) continue;
    for (unsigned i_hit = 0; i_hit < lhcbid->size(); i_hit++) {
      unsigned modID = (lhcbid->at(i_hit) >> 12) & 0x3F;
      // skip if not the requested module
      if (modID != moduleNo) continue;
      float phi = TMath::ATan2(y->at(i_hit), x->at(i_hit));
      // skip if phi out of bounds
      if (phi < phiMin || phi > phiMax) continue;
      printf("MCP %u at eta=%.3f with hit id %u at (x=%.3f, y=%.3f, phi=%.3f, t=%.3f) in Module %u with PID: %d\n",
             i_track, mcEta, lhcbid->at(i_hit), x->at(i_hit), y->at(i_hit), phi, t->at(i_hit), modID, pid);
      printf("\tThis MCP with %u matches has hits in the following %lu modules: {",
             nMatches, lhcbid->size());
      for (unsigned i_hit2 = 0; i_hit2 < lhcbid->size(); i_hit2++) {
        printf("%u", (lhcbid->at(i_hit2) >> 12) & 0x3F);
        if (i_hit2 != lhcbid->size() - 1) printf(", ");
      }
      printf("}\n\n");
    }
  }

  return;
}