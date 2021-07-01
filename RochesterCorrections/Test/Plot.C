#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooFitResult.h"
#include "TH1F.h"

using namespace RooFit;

double* doFit(string condition, double* init_conditions, string invariantMass)
{
  // Read data to a TTree
  TFile *file0 = TFile::Open("./RochesterCorrections/Test/Run2012BC_DoubleMuParked_Muons_RochCor.root");
  TTree *DataTree = (TTree*)file0->Get(("Events"));

  // Create canvas for plot
  TCanvas* c = new TCanvas;

  // Invariant mass
  RooRealVar Dimuon_mass(invariantMass.c_str(), invariantMass.c_str(), 88, 94);

  // Frame for plotting
  RooPlot *frame = Dimuon_mass.frame(RooFit::Title("Invariant Mass"));

  // Limits for eta
  RooRealVar quantity("Muon_eta_pos", "Muon_eta_pos", -3, 3);

  // Select the bin associated to the condition
  RooFormulaVar* redeuce = new RooFormulaVar("EtaBin", condition.c_str(), RooArgList(quantity));
  RooDataSet *Data_ALL = new RooDataSet("DATA", "DATA", DataTree, RooArgSet(Dimuon_mass, quantity), *redeuce);

  // Binned RooDataHist for performing the fit
  RooDataHist* dh_ALL = Data_ALL->binnedClone();

  // Gaussian variables
  RooRealVar sigma("sigma","sigma",init_conditions[1], 0, 4);
  RooRealVar mean1("mean1","mean1",init_conditions[0], 90, 92);

  // Fit function
  RooGaussian *gaussian1 = new RooGaussian("signal1", "signal1", Dimuon_mass, mean1, sigma);

  // Fit result
  RooFitResult* fitres = new RooFitResult;
  fitres = gaussian1->fitTo(*dh_ALL, RooFit::Save());

  // Output
  double* output = new double[2];

  output[0] = mean1.getVal();
  output[1] = mean1.getError();

 // Plotting the fit
 frame->SetXTitle("#mu^{+}#mu^{-} invariant mass [GeV/c^{2}]");
 Data_ALL->plotOn(frame);
 gaussian1->plotOn(frame);
 gaussian1->paramOn(frame);

 c->cd();
 frame->Draw("");

// Save plots to Fit Result
 if (invariantMass == "Dimuon_mass_cor") {
     c->SaveAs(("Fit Result/" + condition + "Cor.pdf").c_str());
 } else {
     c->SaveAs(("Fit Result/" + condition + ".pdf").c_str());
 }

  // DELETING ALLOCATED MEMORY
  delete file0;
  delete Data_ALL;
  delete dh_ALL;
  delete c;

  return output;

}

// Bin conditions
string* get_conditions(int bin_n, double* bins, string quantity)
{
    string* conditions = new string[bin_n];
    for (int i = 0; i < bin_n; i++)
    {
        conditions[i] = quantity + ">" + to_string(bins[i]) + " && " + quantity + "<" + to_string(bins[i+1]);
    }
    return conditions;
}

// Make final histogram
TH1F* make_hist(string name, double** values, int qnt, int bin_n, Double_t* binning)
{
    TH1F* hist = new TH1F(name.c_str(), name.c_str(), bin_n, binning);
    hist->SetMaximum(91);
    hist->SetMinimum(90.6);

    for (int i = 0; i < bin_n; i++)
    {
        hist->SetBinContent(i, values[i][qnt]);
        hist->SetBinError(i, values[i][qnt+1]);
    }

    return hist;
}

int main()
{
  // Eta bins
  double bins[] = {-2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0};
  int bin_n = 20;

  // Initial parameters
  double *init_conditions = new double[2];
  init_conditions[0] = 91; // Z peak
  init_conditions[1] = 2; //sigma

  // Bin conditions
  string* conditions = get_conditions(bin_n, bins, "Muon_eta_pos");

  // Arrays for fit results
  double ** yields_n_errs = new double*[bin_n];
  double ** yields_n_errs_cor = new double*[bin_n];

  // Loop over bins
  for (int i = 0; i < bin_n; i++)
  {
      yields_n_errs[i] = doFit(conditions[i], init_conditions, "Dimuon_mass");
      yields_n_errs_cor[i] = doFit(conditions[i], init_conditions, "Dimuon_mass_cor");
  }

  // Make histograms
  TH1F *yield_ALL  = make_hist("ALL" , yields_n_errs, 0, bin_n, bins);
  TH1F *yield_ALL_cor  = make_hist("ALL" , yields_n_errs_cor, 0, bin_n, bins);

  // Canvas for plotting results
  TCanvas *c1 = new TCanvas("c1","c1");

  // Plot results
  yield_ALL->GetXaxis()->SetTitle("#eta of #mu^{+}");
  yield_ALL->GetYaxis()->SetTitle("Mean of M(#mu^{+}#mu^{-})");
  yield_ALL->Draw();
  yield_ALL_cor->Draw("SAME PLC");

  // Save plot
  c1->SaveAs("./RochesterCorrections/Test/dimuon_mass_eta.pdf");
}
