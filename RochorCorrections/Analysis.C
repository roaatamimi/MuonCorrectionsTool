#include "Analysis.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include <vector>

void Analysis::main()
{
  //Creating a TTree from a ROOT-file
  TFile *f = new TFile("ObjectInfoNtuple.root");
  TTree *t1 = (TTree*)f->Get("mtree");

  //Variables to hold values read from the tree
  Int_t numbermuon = 0;
  vector<float>* muon_pt = NULL;
  vector<float>* muon_eta = NULL;
  vector<float>* muon_phi = NULL;
  vector<float>* muon_ch = NULL;

  //Set addresses to make the tree populate the variables when reading an entry
  t1->SetBranchAddress("numbermuon", &numbermuon);
  t1->SetBranchAddress("muon_pt", &muon_pt);
  t1->SetBranchAddress("muon_eta", &muon_eta);
  t1->SetBranchAddress("muon_phi", &muon_phi);
  t1->SetBranchAddress("muon_ch", &muon_ch);

  //Variable for muon mass
  const Float_t muon_m = 0.1056583745; //GeV

  // Number of entries
  Int_t nentries = (Int_t)t1->GetEntries();
  std::cout << "Number of entries: " <<nentries << std::endl << std::flush;
  if (nentries == 0) return;

  rochcor2012 *rmcor = new rochcor2012(); // make the pointer of rochcor class

  //for-loop of the event
  for (Int_t k=0; k<nentries; ++k){
    t1->GetEntry(k);

    //Run the correction if the event contains muons
    if (numbermuon > 0){
      TLorentzVector mu; //TLorentzVector of the reconstructed muon

      //Set TLorentzVector of muon object
      mu.SetPtEtaPhiM(muon_pt->at(0), muon_eta->at(0), muon_phi->at(0), muon_m);

      float qter = 1.0; // added it by Higgs group’s request to propagate the uncertainty

      //If you run MC, apply the muon momentum correction, “momcor_mc( )” function (only for MC)
      //ntrk (number of track layer) is one of input and it can slightly improved the extra smearing
      rmcor->momcor_mc(mu, muon_ch->at(0), 0, qter); //ntrk is the third parameter

      //If you run data, apply the muon momentum correction, "momcor_data()" function (only for Data)
      // No run dependence for 2012 data, so default of “runopt=0”
      rmcor->momcor_data(mu, muon_ch->at(0), 0, qter); //runopt is the third parameter
    }
  }
}
