#include "Analysis.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include <vector>
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "Math/Vector4Dfwd.h"
#include "Math/Vector4D.h"
#include "TStreamerInfo.h"
#include "ROOT/RDF/RInterface.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TStyle.h"

using namespace ROOT::VecOps;

// Compute the invariant mass of two muon four-vectors
float computeInvariantMass(RVec<float>& pt, RVec<float>& eta, RVec<float>& phi, RVec<float>& mass) {
  ROOT::Math::PtEtaPhiMVector m1(pt[0], eta[0], phi[0], mass[0]);
  ROOT::Math::PtEtaPhiMVector m2(pt[1], eta[1], phi[1], mass[1]);
  return (m1 + m2).mass();
}

// Compute the corrected invariant mass of two muon four-vectors
float correctedInvariantMass(RVec<TLorentzVector> cor_muons) {
  ROOT::Math::PtEtaPhiMVector m1(cor_muons[0].Pt(), cor_muons[0].Eta(), cor_muons[0].Phi(), cor_muons[0].M());
  ROOT::Math::PtEtaPhiMVector m2(cor_muons[1].Pt(), cor_muons[1].Eta(), cor_muons[1].Phi(), cor_muons[1].M());
  return (m1 + m2).mass();
}

// Create TLorentzVectors
RVec<TLorentzVector> createVector(RVec<float>& pt, RVec<float>& eta, RVec<float>& phi, RVec<float>& mass) {
  TLorentzVector mu1;
  TLorentzVector mu2;
  mu1.SetPtEtaPhiM(pt[0], eta[0], phi[0], mass[0]);
  mu2.SetPtEtaPhiM(pt[1], eta[1], phi[1], mass[1]);
  RVec<TLorentzVector> vectors {mu1, mu2};

  return vectors;
}

// Add corrections to MC muons
RVec<TLorentzVector> correcMCMuon(RVec<TLorentzVector> muons, RVec<int>& charge) {
  rochcor2012 rmcor; // make the pointer of rochcor class
  float ntrk = 0; //ntrk (number of track layer) is one of input and it can slightly improved the extra smearing
  float qter = 1.0; // added it by Higgs group’s request to propagate the uncertainty

  rmcor.momcor_mc(muons[0], charge[0], ntrk, qter);
  rmcor.momcor_mc(muons[1], charge[1], ntrk, qter);
  RVec<TLorentzVector> vectors {muons[0], muons[1]};

  return vectors;
}

// Add corrections to data muons
RVec<TLorentzVector> correctDataMuon(RVec<TLorentzVector> muons, RVec<int>& charge) {
  rochcor2012 rmcor; // make the notpointer of rochcor class
  float runopt = 0; //No run dependence for 2012 data, so default of “runopt=0”
  float qter = 1.0; // added it by Higgs group’s request to propagate the uncertainty

  rmcor.momcor_mc(muons[0], charge[0], runopt, qter);
  rmcor.momcor_mc(muons[1], charge[1], runopt, qter);
  RVec<TLorentzVector> vectors {muons[0], muons[1]};

  return vectors;
}

void Analysis::main()
{
  // Create dataframe from NanoAOD files
  //ROOT::RDataFrame df("Events", "root://eospublic.cern.ch//eos/opendata/cms/derived-data/AOD2NanoAODOutreachTool/Run2012BC_DoubleMuParked_Muons.root");
  ROOT::RDataFrame df("Events", "Run2012BC_DoubleMuParked_Muons.root"); // when saved locally

  // Select events with exactly two muons
  auto df_2mu = df.Filter("nMuon == 2", "Events with exactly two muons");

  // Select events with two muons of opposite charge
  auto df_os = df_2mu.Filter("Muon_charge[0] != Muon_charge[1]", "Muons with opposite charge");

  // Compute invariant mass of the dimuon system
  auto df_mass = df_os.Define("Dimuon_mass", computeInvariantMass, {"Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass"});

  // Add TLorentzVectors
  auto df_tlv = df_mass.Define("TLVectors", createVector, {"Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass"});

  // Run the correctios and add corrected muons as a new column

  // If you run MC, apply the muon momentum correction "correctMCMuon" function (only for MC)
  //auto df_tlv = df_cor.Define("CorrectedMuons", correctMCMuon, {"TLVectors","Muon_charge"});

  // If you run data, apply the muon momentum correction "correctDataMuon" function (only for Data)
  auto df_cor = df_tlv.Define("CorrectedMuons", correctDataMuon, {"TLVectors","Muon_charge"});

  // Compute invariant mass of the corrected dimuon system
  auto df_mass_cor = df_cor.Define("Dimuon_mass_cor", correctedInvariantMass, {"CorrectedMuons"});

  // Check columns
  auto colNames = df_mass_cor.GetColumnNames();
  for (auto &&colName : colNames) std::cout << colName << std::endl;

  // Save the new dataframe to an output-file
  // Does not include TLVectors or corrected TLVectors
  std::cout << "Saving dataframe" << std::endl;
  df_mass_cor.Snapshot("Events", "Run2012BC_DoubleMuParked_Muons_RochCor.root", {"Dimuon_mass", "Dimuon_mass_cor", "nMuon", "Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass", "Muon_charge"});

}
