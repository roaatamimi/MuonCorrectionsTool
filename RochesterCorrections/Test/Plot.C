#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TChain.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooFormulaVar.h"
#include "RooVoigtian.h"
#include "RooFitResult.h"
#include "TH1F.h"
#include "TList.h"
#include "TFileMerger.h"

using namespace RooFit;

double* doFit(string condition, string invariantMass, string mu_charge, bool isData)
{  
  // Create a TTree
  const char *pathToFile = "";
  
  if(isData) {
    pathToFile = "./RochesterCorrections/Test/Run2012BC_DoubleMuParked_Muons_Cor.root";
  } else {
    pathToFile = "./RochesterCorrections/Test/ZZTo2e2mu_Cor.root";
  }
  
  TFile *f1 = TFile::Open(pathToFile);
  TTree *DataTree = (TTree*)f1->Get("Events");

  // Create canvas for plot
  TCanvas* c = new TCanvas;

  // Invariant mass
  RooRealVar Dimuon_mass(invariantMass.c_str(), invariantMass.c_str(), 88, 94);

  // Frame for plotting
  RooPlot *frame = Dimuon_mass.frame(RooFit::Title("Invariant Mass"));

  // Eta
  RooRealVar quantity(mu_charge.c_str(), mu_charge.c_str(), -3, 3);

  // Select the bin associated to the condition
  RooFormulaVar* binSelection = new RooFormulaVar("EtaBin", condition.c_str(), RooArgList(quantity));
  RooDataSet *data = new RooDataSet("data", "data", DataTree, RooArgSet(Dimuon_mass, quantity), *binSelection);

  // Binned RooDataHist for performing the fit
  RooDataHist* dataHist = data->binnedClone();

  //Voigtian variables
  RooRealVar sigma("sigma","sigma", 2, 0, 4);
  RooRealVar mean("mean","mean", 91, 90, 92);
  RooRealVar gamma("gamma", "gamma", 2.5, 0, 10);

  // Fit function
  RooVoigtian *voigtian = new RooVoigtian("voigtianFit", "voigtianFit", Dimuon_mass, mean, gamma, sigma);

  // Fit result
  RooFitResult* fitres = new RooFitResult;
  fitres = voigtian->fitTo(*dataHist, RooFit::Save());

  // Output
  double* output = new double[2];

  output[0] = mean.getVal();
  output[1] = mean.getError();

  // Plotting the fit
  frame->SetXTitle("#mu^{+}#mu^{-} invariant mass [GeV/c^{2}]");
  data->plotOn(frame);
  voigtian->plotOn(frame);
  voigtian->paramOn(frame);

  c->cd();
  frame->Draw("");

  // Save plots to Fit Result
  if (mu_charge == "Muon_eta_pos") { // positive
    if (isData) { // data
      if (invariantMass == "Dimuon_mass_cor") { // cor
        c->SaveAs(("./RochesterCorrections/Test/Fit Result/" + condition + "_PDCor.pdf").c_str());
      } else { // not cor
        c->SaveAs(("./RochesterCorrections/Test/Fit Result/" + condition + "_PD.pdf").c_str());
      }
    } else { // MC
      if (invariantMass == "Dimuon_mass_cor") { // cor
        c->SaveAs(("./RochesterCorrections/Test/Fit Result/" + condition + "_PMCCor.pdf").c_str());
      } else { // not cor
        c->SaveAs(("./RochesterCorrections/Test/Fit Result/" + condition + "_PMC.pdf").c_str());
      }
    }
  } else { // negative
    if (isData) { // data
      if (invariantMass == "Dimuon_mass_cor") { // cor
        c->SaveAs(("./RochesterCorrections/Test/Fit Result/" + condition + "_NDCor.pdf").c_str());
      } else { // not cor
        c->SaveAs(("./RochesterCorrections/Test/Fit Result/" + condition + "_ND.pdf").c_str());
      }
    } else { // MC
      if (invariantMass == "Dimuon_mass_cor") { // cor
        c->SaveAs(("./RochesterCorrections/Test/Fit Result/" + condition + "_NMCCor.pdf").c_str());
      } else { // not cor
        c->SaveAs(("./RochesterCorrections/Test/Fit Result/" + condition + "_NMC.pdf").c_str());
      }
    }

  }

  // DELETING ALLOCATED MEMORY
  delete f1;
  delete data;
  delete dataHist;
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
    hist->SetMaximum(91.2);
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

  // Bin conditions for positive muon
  string pos = "Muon_eta_pos";
  string* conditions_pos = get_conditions(bin_n, bins, pos);

  // Bin conditions for negative muon
  string neg = "Muon_eta_neg";
  string* conditions_neg = get_conditions(bin_n, bins, neg);

  // Arrays for fit results
  double ** res_pos_data = new double*[bin_n];
  double ** res_pos_data_cor = new double*[bin_n];
  double ** res_pos_MC = new double*[bin_n];
  double ** res_pos_MC_cor = new double*[bin_n];

  double ** res_neg_data = new double*[bin_n];
  double ** res_neg_data_cor = new double*[bin_n];
  double ** res_neg_MC = new double*[bin_n];
  double ** res_neg_MC_cor = new double*[bin_n];

  // Loop over bins
  for (int i = 0; i < bin_n; i++)
  {
    res_pos_data[i] = doFit(conditions_pos[i], "Dimuon_mass", pos, true);
    res_pos_data_cor[i] = doFit(conditions_pos[i], "Dimuon_mass_cor", pos, true);

    res_pos_MC[i] = doFit(conditions_pos[i], "Dimuon_mass", pos, false);
    res_pos_MC_cor[i] = doFit(conditions_pos[i], "Dimuon_mass_cor", pos, false);

    res_neg_data[i] = doFit(conditions_neg[i], "Dimuon_mass", neg, true);
    res_neg_data_cor[i] = doFit(conditions_neg[i], "Dimuon_mass_cor", neg, true);

    res_neg_MC[i] = doFit(conditions_neg[i], "Dimuon_mass", neg, false);
    res_neg_MC_cor[i] = doFit(conditions_neg[i], "Dimuon_mass_cor", neg, false);
  }

  // Make histograms
  TH1F *hist1 = make_hist("hist1" , res_pos_data, 0, bin_n, bins);
  TH1F *hist2 = make_hist("hist2" , res_pos_data_cor, 0, bin_n, bins);
  TH1F *hist3 = make_hist("hist3" , res_pos_MC, 0, bin_n, bins);
  TH1F *hist4 = make_hist("hist4" , res_pos_MC_cor, 0, bin_n, bins);
  TH1F *hist5 = make_hist("hist5" , res_neg_data, 0, bin_n, bins);
  TH1F *hist6 = make_hist("hist6" , res_neg_data_cor, 0, bin_n, bins);
  TH1F *hist7 = make_hist("hist7" , res_neg_MC, 0, bin_n, bins);
  TH1F *hist8 = make_hist("hist8" , res_neg_MC_cor, 0, bin_n, bins);

  // Canvas for plotting results
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->Divide(2,2);

  // Plot non-corrected positive
  c1->cd(1);
  hist1->SetTitle("Before corrections");
  hist1->GetXaxis()->SetTitle("#eta of #mu^{+}");
  hist1->GetXaxis()->SetTitleSize(0.04);
  hist1->GetYaxis()->SetTitle("Mean of M(#mu^{+}#mu^{-})");
  hist1->GetYaxis()->SetTitleSize(0.04);
  hist1->SetStats(0);
  hist1->Draw("PLC");
  hist3->Draw("SAME PLC");

  // Plot corrected positive
  c1->cd(2);
  hist2->SetTitle("After corrections");
  hist2->GetXaxis()->SetTitle("#eta of #mu^{+}");
  hist2->GetXaxis()->SetTitleSize(0.04);
  hist2->GetYaxis()->SetTitle("Mean of M(#mu^{+}#mu^{-})");
  hist2->GetYaxis()->SetTitleSize(0.04);
  hist2->SetStats(0);
  hist2->Draw("PLC");
  hist4->Draw("SAME PLC");

  // Plot non-corrected negative
  c1->cd(3);
  hist5->SetTitle("Before corrections");
  hist5->GetXaxis()->SetTitle("#eta of #mu^{-}");
  hist5->GetXaxis()->SetTitleSize(0.04);
  hist5->GetYaxis()->SetTitle("Mean of M(#mu^{+}#mu^{-})");
  hist5->GetYaxis()->SetTitleSize(0.04);
  hist5->SetStats(0);
  hist5->Draw("PLC");
  hist7->Draw("SAME PLC");

  // Plot corrected negative
  c1->cd(4);
  hist6->SetTitle("After corrections");
  hist6->GetXaxis()->SetTitle("#eta of #mu^{-}");
  hist6->GetXaxis()->SetTitleSize(0.04);
  hist6->GetYaxis()->SetTitle("Mean of M(#mu^{+}#mu^{-})");
  hist6->GetYaxis()->SetTitleSize(0.04);
  hist6->SetStats(0);
  hist6->Draw("PLC");
  hist8->Draw("SAME PLC");

  // Save plot
  c1->SaveAs("./RochesterCorrections/Test/dimuon_mass_eta_all.pdf");
}
