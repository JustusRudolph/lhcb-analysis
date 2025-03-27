#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>

#include <cstdlib>
#include <iostream>
#include <string>

void get_module_z_positions() {
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

  TProfile* moduleToZ = new TProfile("module_to_z", "Z positions of modules",
                                     nBins, moduleBinEdges.data());
  
  // now get the data from tree
  unsigned lhcbid;
  float z;
  TTree* fakeClusTree = (TTree*) file->Get("FakeClusters");
  if (!fakeClusTree) {
      std::cerr << "TTree not found!" << std::endl;
      return;
  }
  fakeClusTree->SetBranchAddress("id", &lhcbid);
  fakeClusTree->SetBranchAddress("z", &z);

  unsigned nEntries = fakeClusTree->GetEntries();
  for (unsigned i = 0; i < nEntries; i++) {
    fakeClusTree->GetEntry(i);
    unsigned tvMod = (lhcbid >> 12) & 0x3F;
    moduleToZ->Fill(tvMod, z);

    // print status
    if (i && (i % 1000000 == 0))
      printf("Finished with %d %% of fake clusters.\n", (int) (100. * (float) i / nEntries));
  }

  // Finally write the profile to a file
  TFile* outFile = new TFile((analysisRoot + "/hists/module_z_positions.root").c_str(), "RECREATE");
  moduleToZ->Write();

  // clean up
  outFile->Close();
  delete outFile;

  file->Close();
  delete file;

  return;
}