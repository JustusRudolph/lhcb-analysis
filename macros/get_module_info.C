#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TMath.h>

#include <cstdlib>
#include <iostream>
#include <string>

// tMin and tMax bound the timestamp histogram, the out of range counts are printed so that
// it is clear how much of the sample falls outside them
void get_module_info(int runNoCheck=-1, int evNoCheck=-1, int moduleNoCheck=-1,
                     float tMin=-5., float tMax=15., int nTimeBins=1000) {
  gROOT->SetBatch();  // don't show hists while running
  std::string stackRoot = std::getenv("STACK_ROOT");
  std::string analysisRoot = std::getenv("ANALYSIS_ROOT");

  TFile *file = TFile::Open((stackRoot + "/fake_clusters.root").c_str());
  if (!file || file->IsZombie()) {
      std::cerr << "Error: Could not open ROOT file." << std::endl;
      return;
  }
  // Get the bins for module edges in case of plotting
  int nBins = 64;
  std::vector<float> moduleBinEdges(nBins+1);
  float step = 1;
  for (unsigned i = 0; i <= nBins; i++) {
    moduleBinEdges[i] = -0.5 + step * i;
  }

  TH1D* moduleDensity = new TH1D("moduleDensity", "Occupancy by module;module ID;Counts",
                                 nBins, moduleBinEdges.data());
  TProfile* moduleToZ = new TProfile("module_to_z", "Z positions of modules",
                                     nBins, moduleBinEdges.data());
  TH1D* phiDensity = new TH1D("phiDensity", "Occupancy by #phi;#phi;Counts", 100, -3.15, 3.15);
  TH2D* xyDensity = new TH2D("xyDensity", "Occupancy by x-y;x;y", 20, -50, 50, 20, -50, 50);
  TProfile* moduleTimes = new TProfile("module_times",
                                       "Mean hit timestamp by module;module ID;t (ns)",
                                       nBins, moduleBinEdges.data());
  TH1D* hitTimes = new TH1D("hit_times", "Hit timestamps;t (ns);Hits",
                            nTimeBins, tMin, tMax);
  
  // now get the data from tree
  unsigned lhcbid, evNo, runNo;
  float x, y, z, t;
  TTree* fakeClusTree = (TTree*) file->Get("FakeClusters");
  if (!fakeClusTree) {
      std::cerr << "TTree not found!" << std::endl;
      return;
  }
  fakeClusTree->SetBranchAddress("runNo", &runNo);
  fakeClusTree->SetBranchAddress("evNo", &evNo);
  fakeClusTree->SetBranchAddress("id", &lhcbid);
  fakeClusTree->SetBranchAddress("x", &x);
  fakeClusTree->SetBranchAddress("y", &y);
  fakeClusTree->SetBranchAddress("z", &z);
  fakeClusTree->SetBranchAddress("t", &t);

  float tMinSeen = 1e9f, tMaxSeen = -1e9f;
  unsigned nEntries = fakeClusTree->GetEntries();
  for (unsigned i = 0; i < nEntries; i++) {
    fakeClusTree->GetEntry(i);
    // only take from specific run and event if specified
    if (runNoCheck >=0 && runNo != runNoCheck) continue;
    if (evNoCheck >=0 && evNo != evNoCheck) continue;
    unsigned tvMod = (lhcbid >> 12) & 0x3F;
    moduleToZ->Fill(tvMod, z);
    moduleDensity->Fill(tvMod);
    // times are not restricted to the requested module, they describe the whole detector
    moduleTimes->Fill(tvMod, t);
    hitTimes->Fill(t);
    if (t < tMinSeen) tMinSeen = t;
    if (t > tMaxSeen) tMaxSeen = t;

    // phi only filled for specific module that's requested
    if (moduleNoCheck >= 0 && tvMod != moduleNoCheck) continue;
    float phi = TMath::ATan2(y, x);
    phiDensity->Fill(phi);
    xyDensity->Fill(x, y);

    // print status
    if (i && (i % 1000000 == 0))
      printf("Finished with %d %% of fake clusters.\n", (int) (100. * (float) i / nEntries));
  }

  printf("Hit timestamps: seen %.3f .. %.3f ns, histogram covers %.3f .. %.3f ns\n",
         tMinSeen, tMaxSeen, tMin, tMax);
  printf("  %.0f below range, %.0f above range, out of %.0f hits\n",
         hitTimes->GetBinContent(0), hitTimes->GetBinContent(nTimeBins + 1),
         hitTimes->GetEntries());

  // Finally write the profile to a file
  TFile* outFile = new TFile((analysisRoot + "/hists/module_mc_info.root").c_str(), "RECREATE");
  moduleToZ->Write();
  moduleDensity->Write();
  phiDensity->Write();
  xyDensity->Write();
  moduleTimes->Write();
  hitTimes->Write();

  // clean up
  outFile->Close();
  delete outFile;

  file->Close();
  delete file;

  return;
}