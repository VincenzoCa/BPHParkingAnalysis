//https://github.com/ICBPHCMS/BPHParkingAnalysis/blob/master/NtupleProducer/macro/fitBmass_fromHistos.C

//plot and fit BToKEE and BToKMuMu

//example to run
//root fitBmass_fromHistos.C'(1, "/my_path/KEE_histos_file.root")'
//root fitBmass_fromHistos.C'(1, "/my_path/KMuMu_histos_file.root")'


#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>

#include "TROOT.h"
#include "TSystem.h"
#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
//#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TString.h"
#include "TCut.h"
#include "TMath.h"
#include "TApplication.h"
#include "TError.h"
#include "TGraphErrors.h"
#include "TPad.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TChain.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"
#include "RooFitResult.h"
#include "TText.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"

using namespace RooFit;


void fitBmass_fromHistos(int isEleFinalState, std::string inFile){

    gROOT->Reset();
    gROOT->Macro("setStyle.C");

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    std::cout << " isEleFinalState = " << isEleFinalState << std::endl;


    TFile* inF = TFile::Open(inFile.c_str());


    TH1D* h_Bmass;
    TH1D* h_Bmass_PF;
    TH1D* h_Bmass_mix;
    TH1D* h_Bmass_Low;
  
    h_Bmass = (TH1D*)inF->Get("KEE_mass_B3r")->Clone("h_KEE_mass_B3r");    
    h_Bmass_PF = (TH1D*)inF->Get("KEE_mass_PF_B3r")->Clone("h_KEE_mass_PF_B3r");
    h_Bmass_mix = (TH1D*)inF->Get("KEE_mass_mix_B3r")->Clone("h_KEE_mass_mix_B3r");
    h_Bmass_Low = (TH1D*)inF->Get("KEE_mass_Low_B3r")->Clone("h_KEE_mass_Low_B3r");


    //now fitting
    float nEv_postFit = 0.;
    float nEvError_postFit = 0.;
    float nBkg_postFit = 0.;
    float nBkgError_postFit = 0.;
    float nBkgInt_postFit = 0.;
    float nBkgIntError_postFit = 0.;
    float chi2 = 0.;

    RooWorkspace w("w");    
    w.factory("x[4.7, 6.]");//4.5, 6.

    w.factory("nbackground[1000, 0, 10000]");//10000, 0, 100000   
    //w.factory("nbackgroundR[10, 0, 100]");   
    w.factory("nsignal[100, 0.0, 100000]");//100, 0.0, 10000
    
    w.factory("Gaussian::smodel(x,mu[5.3,4.5,6],sigma[0.05,0,0.15])");
    //w.factory("Gaussian::smodel(x,mu[5.3,4.5,6],sigma[0.02,0,0.03])");
    //w.factory("RooCBShape::smodel(x,m[5.3,4.5,6],s[0.1,0.,1.],a[1.2,0.,3.],n[1,0.1,6.])");
    //w.factory("RooCBShape::CBall(x[0,15], mean[11000,13000], sigma[5000,200000], alpha[0,10000],n[0,100000])");
    RooAbsPdf * smodel = w.pdf("smodel");
    
    //Background: exponential + gaussian
    //
    //w.factory("Exponential::bmodel1(x,tau[-2,-3,0])");
    //RooAbsPdf * bmodel1 = w.pdf("bmodel1");
    //w.factory("Gaussian::bmodel2(x,mub[5.3,4.5,6],sigmab[0.15,0.05,2.])");
    //RooAbsPdf * bmodel2 = w.pdf("bmodel2");
    //w.factory("SUM::modelb(nbackground * bmodel1, -nbackgroundR * bmodel1, nbackgroundR * bmodel2)");
    //RooAbsPdf * modelb = w.pdf("modelb");

    //Background: exponential
    //
    w.factory("Exponential::modelb(x,tau[-2,-3,0])");
    RooAbsPdf * modelb = w.pdf("modelb");

    w.factory("SUM::model(nbackground * modelb, nsignal * smodel)");
    RooAbsPdf * model = w.pdf("model");
    
    RooDataHist hBMass("hBMass", "hBMass", *w.var("x"), Import(*(h_Bmass)));
    RooDataHist hBMass_PF("hBMass_PF", "hBMass_PF", *w.var("x"), Import(*(h_Bmass_PF)));
    RooDataHist hBMass_mix("hBMass_mix", "hBMass_mix", *w.var("x"), Import(*(h_Bmass_mix)));
    RooDataHist hBMass_Low("hBMass_Low", "hBMass_Low", *w.var("x"), Import(*(h_Bmass_Low)));
    w.Print();

    RooFitResult * r = model->fitTo(hBMass, Minimizer("Minuit2"),Save(true));
    std::cout << " fit status = " << r->status() << std::endl;

    RooPlot * plot = w.var("x")->frame();
    if(isEleFinalState){
      plot->SetXTitle("K(JPsi)ee mass (GeV)");
      //else plot->SetXTitle("Kee mass (GeV)");
    }
    else{
      plot->SetXTitle("K(JPsi)#mu#mu mass (GeV)");
      //else plot->SetXTitle("K#mu#mu mass (GeV)");
    }
    plot->SetTitle("");
    plot->SetAxisRange(4.,6);
    hBMass.plotOn(plot);
    model->plotOn(plot);
    model->plotOn(plot, Components("modelb"),LineStyle(kDashed));
    model->plotOn(plot, Components("smodel"),LineColor(kOrange));
    hBMass_PF.plotOn(plot,LineColor(kRed),MarkerColor(kRed)) ;
    hBMass_mix.plotOn(plot,LineColor(kGreen+1),MarkerColor(kGreen+1)) ;
    hBMass_Low.plotOn(plot,LineColor(kMagenta),MarkerColor(kMagenta)) ;
    chi2 = plot->chiSquare();

    RooRealVar* parS = (RooRealVar*) r->floatParsFinal().find("nsignal");
    RooRealVar* parB = (RooRealVar*) r->floatParsFinal().find("nbackground");
    nEv_postFit = parS->getValV();
    nEvError_postFit = parS->getError();
    nBkg_postFit = parB->getValV();
    nBkgError_postFit = parB->getError();

    std::cout << " **** JPsi selection signal events = \t " << parS->getValV() << " error = " << parS->getError() 
	      << " bkg events = " << parB->getValV() << " error = " << parB->getError() << std::endl;

    RooRealVar* parMean = (RooRealVar*) r->floatParsFinal().find("mu");
    RooRealVar* parSigma = (RooRealVar*) r->floatParsFinal().find("sigma");

    float meanVal = parMean->getValV();
    float sigmaVal = parSigma->getValV();

    std::cout << "\n  parMean = " << parMean->getValV() << " parSigma = " << parSigma->getValV() << std::endl;
    
    w.var("x")->setRange("signalRange", meanVal - 3.*sigmaVal, meanVal + 3.*sigmaVal);
    //RooAbsReal* bkgIntegral = w.pdf("bmodel")->createIntegral(*(w.var("x")), NormSet(*(w.var("x"))), Range("signalRange")) ;
    RooAbsReal* bkgIntegral = w.pdf("modelb")->createIntegral(*(w.var("x")), NormSet(*(w.var("x"))), Range("signalRange")) ;
    // w.var("x")->setRange("signalRange", 4.5, 6.);
    // RooAbsReal* bkgIntegral = w.pdf("bmodel")->createIntegral(*(w.var("x")), NormSet(*(w.var("x"))), Range("signalRange")) ;
    std::cout << "\n  bkgIntegral = " << bkgIntegral->getVal() << std::endl;
    nBkgInt_postFit = bkgIntegral->getVal() * nBkg_postFit;
    nBkgIntError_postFit = nBkgError_postFit * bkgIntegral->getVal();

    
    TCanvas * cc = new TCanvas();
    cc->SetLogy(0);
    plot->Draw();

    TLatex tL;
    tL.SetNDC();
    tL.SetTextSize(0.05);
    tL.SetTextFont(42);
    tL.DrawLatex(0.65,0.9, Form("S %.1f +/- %1.f",nEv_postFit, nEvError_postFit));
    TLatex tL2;
    tL2.SetNDC();
    tL2.SetTextSize(0.05);
    tL2.SetTextFont(42);
    tL2.DrawLatex(0.65,0.85, Form("B %.1f +/- %1.f",nBkgInt_postFit, nBkgIntError_postFit));


    if(isEleFinalState) cc->Print(Form("plots/KEE_%s.png",h_Bmass->GetName()), "png");
    else cc->Print(Form("plots/KMuMu_%s.png",h_Bmass->GetName()), "png");
  
    
    std::cout << " ***** summary ***** "<< std::endl;

    std::cout   << "\n \n  category: " << h_Bmass->GetName()
                << "\n \t signal = " << nEv_postFit << "+/-" << nEvError_postFit
                << "\n \t bkg = " << nBkg_postFit << "+/-" << nBkgError_postFit 
                << "\n \t bkg in +/- 3sigma = " << nBkgInt_postFit << " +/- " << nBkgIntError_postFit
                << " chi2 = " << chi2 << std::endl;

} 
