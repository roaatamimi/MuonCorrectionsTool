#include "Analysis.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "Math/Vector4D.h"

using namespace ROOT::VecOps;

// Compute the invariant mass of two muon four-vectors
float computeInvariantMass(RVec<float>& pt, RVec<float>& eta, RVec<float>& phi, RVec<float>& mass) {
  ROOT::Math::PtEtaPhiMVector m1(pt[0], eta[0], phi[0], mass[0]);
  ROOT::Math::PtEtaPhiMVector m2(pt[1], eta[1], phi[1], mass[1]);
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
RVec<TLorentzVector> correctMCMuon(RVec<TLorentzVector> muons, RVec<int>& charge) {
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

// Change datatype from vector to float
float vecToFloat(RVec<float>& value) {
  return value[0];
}

// Extract corrected variables from TLorentzVectors
// Pt
RVec<float> correctedPt(RVec<TLorentzVector> muons){
  RVec<float> values {static_cast<float>(muons[0].Pt()), static_cast<float>(muons[1].Pt())};
  return values;
}
// Eta
RVec<float> correctedEta(RVec<TLorentzVector> muons){
  RVec<float> values {static_cast<float>(muons[0].Eta()), static_cast<float>(muons[1].Eta())};
  return values;
}
// Phi
RVec<float> correctedPhi(RVec<TLorentzVector> muons){
  RVec<float> values {static_cast<float>(muons[0].Phi()), static_cast<float>(muons[1].Phi())};
  return values;
}
// Mass
RVec<float> correctedMass(RVec<TLorentzVector> muons){
  RVec<float> values {static_cast<float>(muons[0].M()), static_cast<float>(muons[1].M())};
  return values;
}


void Analysis::main()
{
  // Create dataframe from NanoAOD files
  ROOT::RDataFrame df("Events", "root://eospublic.cern.ch//eos/opendata/cms/derived-data/AOD2NanoAODOutreachTool/Run2012BC_DoubleMuParked_Muons.root");
  //ROOT::RDataFrame df("Events", "Run2012BC_DoubleMuParked_Muons.root"); // use when saved locally

  // Select events with exactly two muons
  auto df_2mu = df.Filter("nMuon == 2", "Events with exactly two muons");

  // Select events with two muons of opposite charge
  auto df_os = df_2mu.Filter("Muon_charge[0] != Muon_charge[1]", "Muons with opposite charge");

  // Compute invariant mass of the dimuon system
  auto df_mass = df_os.Define("Dimuon_mass", computeInvariantMass, {"Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass"});

  // Eta of positive charge muons
  auto df_eta_pos_vec = df_mass.Define("mask", "Muon_charge > 0")
                           .Define("Muon_eta_pos_vec", "Muon_eta[mask]");

  // Eta of negative charge muons
  auto df_eta_neg_vec = df_eta_pos_vec.Define("mask2", "Muon_charge < 0")
                              .Define("Muon_eta_neg_vec", "Muon_eta[mask2]");

  // Eta values as floats instead of vectors
  auto df_eta_pos = df_eta_neg_vec.Define("Muon_eta_pos", vecToFloat, {"Muon_eta_pos_vec"});
  auto df_eta_neg = df_eta_pos.Define("Muon_eta_neg", vecToFloat, {"Muon_eta_neg_vec"});

  // Add TLorentzVectors
  auto df_tlv = df_eta_neg.Define("TLVectors", createVector, {"Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass"});

  // Run the correctios and add corrected muons as a new column
  // CHOOSE THE CORRECT FUNCTION FOR DATA OR MC AND COMMENT THE OTHER

  // If you run MC, apply the muon momentum correction "correctMCMuon" function (only for MC)
  auto df_cor = df_tlv.Define("CorrectedMuons", correctMCMuon, {"TLVectors","Muon_charge"});

  // If you run data, apply the muon momentum correction "correctDataMuon" function (only for Data)
  auto df_cor = df_tlv.Define("CorrectedMuons", correctDataMuon, {"TLVectors","Muon_charge"});

  // Add columns for corrected values
  auto df_cor_pt = df_cor.Define("Muon_pt_cor", correctedPt, {"CorrectedMuons"});
  auto df_cor_eta = df_cor_pt.Define("Muon_eta_cor", correctedEta, {"CorrectedMuons"});
  auto df_cor_phi = df_cor_eta.Define("Muon_phi_cor", correctedPhi, {"CorrectedMuons"});
  auto df_cor_mass = df_cor_phi.Define("Muon_mass_cor", correctedMass, {"CorrectedMuons"});

  // Compute invariant mass of the corrected dimuon system
  auto df_mass_cor = df_cor_mass.Define("Dimuon_mass_cor", computeInvariantMass, {"Muon_pt_cor", "Muon_eta_cor", "Muon_phi_cor", "Muon_mass_cor"});

  // Check columns
  auto colNames = df_mass_cor.GetColumnNames();
  for (auto &&colName : colNames) std::cout << colName << std::endl;

  // Save the new dataframe to a root-file
  std::cout << "Saving dataframe" << std::endl;
  df_mass_cor.Snapshot("Events", "Run2012BC_DoubleMuParked_Muons_RochCor.root", {"Dimuon_mass", "Dimuon_mass_cor", "nMuon", "Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass", "Muon_charge", "Muon_eta_pos", "Muon_eta_neg", "Muon_pt_cor", "Muon_eta_cor", "Muon_phi_cor", "Muon_mass_cor"});

}
