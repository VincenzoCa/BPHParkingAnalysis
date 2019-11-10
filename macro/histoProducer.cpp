/*
g++ -Wall -o histoProducer `root-config --cflags --glibs` -lRooFitCore histoProducer.cpp
./histoProducer --ntupleList (list.txt) --JOBid (1,2..) --outputFolder ("outfolder") --testFile ("path")
*/

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>

#include "TChain.h"
#include "TCut.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

int main(int argc, char **argv){
    
    if(argc < 2) {
        std::cout << " Missing arguments " << std::endl;
        return -1;
    }

    std::string ntupleList = "-1";
    std::string JOBid = "-1";
    std::string outputFolder = "-1";
    std::string testFile = "-1";
 
    for (int i = 1; i < argc; ++i) {
        if(std::string(argv[i]) == "--ntupleList") {
            if (i + 1 < argc) {
                ntupleList = argv[i+1];
                break;
            } 
            else {
                std::cerr << " --ntupleList option requires one argument " << std::endl;
                return 1;
            }
        }
    }  
    for (int i = 1; i < argc; ++i) {
        if(std::string(argv[i]) == "--JOBid") {
            if (i + 1 < argc) {
            JOBid = argv[i+1];
            break;
        } 
        else {
            std::cerr << " --JOBid option requires one argument " << std::endl;
            return 1;
        }
        }
    }  
    for (int i = 1; i < argc; ++i) {
        if(std::string(argv[i]) == "--outputFolder") {
            if (i + 1 < argc) {
            outputFolder = argv[i+1];
            break;
            } 
            else {
                std::cerr << " --outputFolder option requires one argument " << std::endl;
                return 1;
            }
        }
    }
    for (int i = 1; i < argc; ++i) {
        if(std::string(argv[i]) == "--testFile") {
            if (i + 1 < argc) {
                testFile = argv[i+1];
                break;
            } 
            else {
                std::cerr << " --testFile option requires one argument " << std::endl;
                return 1;
            }
        }
    }

  
    if(ntupleList != "-1" && (JOBid == "-1" || outputFolder == "-1")){
        std::cout << " configuration ERROR => splitting file based but missing JOBid and output folder " << std::endl;
        return -1;
    }

    if(ntupleList == "-1" && testFile == "-1"){
        std::cout << " configuration ERROR => need a file list or a test file " << std::endl;
        return -1;
    }

    std::cout << " ntupleList = " << ntupleList << " JOBid = " << JOBid << " outputFolder = " << outputFolder << std::endl;


    TChain* ch = new TChain("skimTree");
  
    if(ntupleList != "-1"){
        std::string rootFileName;
        std::ifstream inFileLong;
        inFileLong.open(ntupleList.c_str(), std::ios::in);
        while(!inFileLong.eof()){
            if(inFileLong >> rootFileName){
                ch -> Add(rootFileName.c_str());
                std::cout << " adding " << rootFileName << std::endl;
            }
        }
    }
    else{
        ch -> Add(testFile.c_str());
    }    
    
  
    //KEE
    TH1F *KEE_mass_B3f = new TH1F("KEE_mass_B3f","", 75, 4.5, 6.);
    KEE_mass_B3f -> Sumw2();
    KEE_mass_B3f -> SetLineColor(kRed);
    KEE_mass_B3f -> SetLineWidth(1);    
    //
    TH1F *KEE_fitMass_B3f = new TH1F("KEE_fitMass_B3f","", 75, 4.5, 6.);
    KEE_fitMass_B3f -> Sumw2();  
    KEE_fitMass_B3f -> SetLineColor(kRed);
    KEE_fitMass_B3f -> SetLineWidth(1);       
    //
    TH1F *KEE_mass_B3r = new TH1F("KEE_mass_B3r","", 75, 4.5, 6.);
    KEE_mass_B3r -> Sumw2();
    KEE_mass_B3r -> SetLineColor(kRed);
    KEE_mass_B3r -> SetLineWidth(1);       
    //
    TH1F *KEE_fitMass_B3r = new TH1F("KEE_fitMass_B3r","", 75, 4.5, 6.);
    KEE_fitMass_B3r -> Sumw2(); 
    KEE_fitMass_B3r -> SetLineColor(kRed);
    KEE_fitMass_B3r -> SetLineWidth(1);       
    //
    TH1F *KEE_mass_PF_B3f = new TH1F("KEE_mass_PF_B3f","", 75, 4.5, 6.);
    KEE_mass_PF_B3f -> Sumw2();
    KEE_mass_PF_B3f -> SetLineColor(kRed);
    KEE_mass_PF_B3f -> SetLineWidth(1);       
    //
    TH1F *KEE_fitMass_PF_B3f = new TH1F("KEE_fitMass_PF_B3f","", 75, 4.5, 6.);
    KEE_fitMass_PF_B3f -> Sumw2();  
    KEE_fitMass_PF_B3f -> SetLineColor(kRed);
    KEE_fitMass_PF_B3f -> SetLineWidth(1);       
    //
    TH1F *KEE_mass_PF_B3r = new TH1F("KEE_mass_PF_B3r","", 75, 4.5, 6.);
    KEE_mass_PF_B3r -> Sumw2();
    KEE_mass_PF_B3r -> SetLineColor(kRed);
    KEE_mass_PF_B3r -> SetLineWidth(1);       
    //
    TH1F *KEE_fitMass_PF_B3r = new TH1F("KEE_fitMass_PF_B3r","", 75, 4.5, 6.);
    KEE_fitMass_PF_B3r -> Sumw2();
    KEE_fitMass_PF_B3r -> SetLineColor(kRed);
    KEE_fitMass_PF_B3r -> SetLineWidth(1);       
    //
    TH1F *KEE_mass_mix_B3f = new TH1F("KEE_mass_mix_B3f","", 75, 4.5, 6.);
    KEE_mass_mix_B3f -> Sumw2();
    KEE_mass_mix_B3f -> SetLineColor(kRed);
    KEE_mass_mix_B3f -> SetLineWidth(1);       
    //
    TH1F *KEE_fitMass_mix_B3f = new TH1F("KEE_fitMass_mix_B3f","", 75, 4.5, 6.);
    KEE_fitMass_mix_B3f -> Sumw2();  
    KEE_fitMass_mix_B3f -> SetLineColor(kRed);
    KEE_fitMass_mix_B3f -> SetLineWidth(1);       
    //
    TH1F *KEE_mass_mix_B3r = new TH1F("KEE_mass_mix_B3r","", 75, 4.5, 6.);
    KEE_mass_mix_B3r -> Sumw2();
    KEE_mass_mix_B3r -> SetLineColor(kRed);
    KEE_mass_mix_B3r -> SetLineWidth(1);       
    //
    TH1F *KEE_fitMass_mix_B3r = new TH1F("KEE_fitMass_mix_B3r","", 75, 4.5, 6.);
    KEE_fitMass_mix_B3r -> Sumw2();
    KEE_fitMass_mix_B3r -> SetLineColor(kRed);
    KEE_fitMass_mix_B3r -> SetLineWidth(1);       
    //
    TH1F *KEE_mass_Low_B3f = new TH1F("KEE_mass_Low_B3f","", 75, 4.5, 6.);
    KEE_mass_Low_B3f -> Sumw2();
    KEE_mass_Low_B3f -> SetLineColor(kRed);
    KEE_mass_Low_B3f -> SetLineWidth(1);   
    //
    TH1F *KEE_fitMass_Low_B3f = new TH1F("KEE_fitMass_Low_B3f","", 75, 4.5, 6.);
    KEE_fitMass_Low_B3f -> Sumw2();
    KEE_fitMass_Low_B3f -> SetLineColor(kRed);
    KEE_fitMass_Low_B3f -> SetLineWidth(1);       
    //
    TH1F *KEE_mass_Low_B3r = new TH1F("KEE_mass_Low_B3r","", 75, 4.5, 6.);
    KEE_mass_Low_B3r -> Sumw2();
    KEE_mass_Low_B3r -> SetLineColor(kRed);
    KEE_mass_Low_B3r -> SetLineWidth(1);       
    //
    TH1F *KEE_fitMass_Low_B3r = new TH1F("KEE_fitMass_Low_B3r","", 75, 4.5, 6.);
    KEE_fitMass_Low_B3r -> Sumw2();
    KEE_fitMass_Low_B3r -> SetLineColor(kRed);
    KEE_fitMass_Low_B3r -> SetLineWidth(1);
    
    //KMuMu
    TH1F *KMM_mass_PF_B3f = new TH1F("KMM_mass_PF_B3f","", 75, 4.5, 6.);
    KMM_mass_PF_B3f -> Sumw2();
    KMM_mass_PF_B3f -> SetLineColor(kRed);
    KMM_mass_PF_B3f -> SetLineWidth(1);       
    //
    TH1F *KMM_fitMass_PF_B3f = new TH1F("KMM_fitMass_PF_B3f","", 75, 4.5, 6.);
    KMM_fitMass_PF_B3f -> Sumw2();  
    KMM_fitMass_PF_B3f -> SetLineColor(kRed);
    KMM_fitMass_PF_B3f -> SetLineWidth(1);       
    //
    TH1F *KMM_mass_PF_B3r = new TH1F("KMM_mass_PF_B3r","", 75, 4.5, 6.);
    KMM_mass_PF_B3r -> Sumw2();
    KMM_mass_PF_B3r -> SetLineColor(kRed);
    KMM_mass_PF_B3r -> SetLineWidth(1);       
    //
    TH1F *KMM_fitMass_PF_B3r = new TH1F("KMM_fitMass_PF_B3r","", 75, 4.5, 6.);
    KMM_fitMass_PF_B3r -> Sumw2();
    KMM_fitMass_PF_B3r -> SetLineColor(kRed);
    KMM_fitMass_PF_B3r -> SetLineWidth(1);    

    
    //KEE
    TCut e_PFov = "e1_isPFov_sk == 0 && e2_isPFov_sk == 0";
    TCut e_PF = "e1_isPF_sk == 1 && e2_isPF_sk == 1";
    TCut e_mix = "((e1_isPF_sk == 1 && e2_isLow_sk == 1) || (e1_isLow_sk == 1 && e2_isPF_sk == 1))"; 
    TCut e_Low = "e1_isLow_sk == 1 && e2_isLow_sk == 1";  
    TCut e_mllraw = "KEE_mll_raw_sk > 2.9 && KEE_mll_raw_sk < 3.3";
    TCut e_mllff = "KEE_mll_fullfit_sk > 2.9 && KEE_mll_fullfit_sk < 3.3";
    //KMuMu
    TCut mu_PF = "mu1_isPF_sk == 1 && mu2_isPF_sk == 1";
    TCut mu_mllraw = "KMM_mll_raw_sk > 2.9 && KMM_mll_raw_sk < 3.3";
    TCut mu_mllff = "KMM_mll_fullfit_sk > 2.9 && KMM_mll_fullfit_sk < 3.3";    

    
    //KEE    
    ch -> Draw("KEE_mass_sk >> KEE_mass_B3f", e_PFov && e_mllff,"goff");
    ch -> Draw("KEE_fit_mass_sk >> KEE_fitMass_B3f", e_PFov && e_mllff,"goff");
    ch -> Draw("KEE_mass_sk >> KEE_mass_B3r", e_PFov && e_mllraw,"goff");
    ch -> Draw("KEE_fit_mass_sk >> KEE_fitMass_B3r", e_PFov && e_mllraw,"goff");
    //    
    ch -> Draw("KEE_mass_sk >> KEE_mass_PF_B3f", e_PF && e_mllff,"goff");
    ch -> Draw("KEE_fit_mass_sk >> KEE_fitMass_PF_B3f", e_PF && e_mllff,"goff");
    ch -> Draw("KEE_mass_sk >> KEE_mass_PF_B3r", e_PF && e_mllraw,"goff");
    ch -> Draw("KEE_fit_mass_sk >> KEE_fitMass_PF_B3r", e_PF && e_mllraw,"goff");
    //    
    ch -> Draw("KEE_mass_sk >> KEE_mass_mix_B3f", e_PFov && e_mix && e_mllff,"goff");
    ch -> Draw("KEE_fit_mass_sk >> KEE_fitMass_mix_B3f", e_PFov && e_mix && e_mllff,"goff");
    ch -> Draw("KEE_mass_sk >> KEE_mass_mix_B3r", e_PFov && e_mix && e_mllraw,"goff");
    ch -> Draw("KEE_fit_mass_sk >> KEE_fitMass_mix_B3r", e_PFov && e_mix && e_mllraw,"goff");  
    //
    ch -> Draw("KEE_mass_sk >> KEE_mass_Low_B3f", e_PFov && e_Low && e_mllff,"goff");
    ch -> Draw("KEE_fit_mass_sk >> KEE_fitMass_Low_B3f", e_PFov && e_Low && e_mllff,"goff");
    ch -> Draw("KEE_mass_sk >> KEE_mass_Low_B3r", e_PFov && e_Low && e_mllraw,"goff");
    ch -> Draw("KEE_fit_mass_sk >> KEE_fitMass_Low_B3r", e_PFov && e_Low && e_mllraw,"goff");  
    
    //KMuMu
    ch -> Draw("KMM_mass_sk >> KMM_mass_PF_B3f", mu_PF && mu_mllff,"goff");
    ch -> Draw("KMM_fit_mass_sk >> KMM_fitMass_PF_B3f", mu_PF && mu_mllff,"goff");
    ch -> Draw("KMM_mass_sk >> KMM_mass_PF_B3r", mu_PF && mu_mllraw,"goff");
    ch -> Draw("KMM_fit_mass_sk >> KMM_fitMass_PF_B3r", mu_PF && mu_mllraw,"goff");    

    
    std::string outName = "histo_output_from_tree";
    if(JOBid != "-1") outName = outputFolder;
    else outName += ".root";
    TFile outHistos(outName.c_str(), "recreate");
    
    outHistos.cd();
  
    //KEE
    KEE_mass_B3f -> Write();
    KEE_fitMass_B3f -> Write();  
    KEE_mass_B3r -> Write();
    KEE_fitMass_B3r -> Write();
    //
    KEE_mass_PF_B3f -> Write();
    KEE_fitMass_PF_B3f -> Write();  
    KEE_mass_PF_B3r -> Write();
    KEE_fitMass_PF_B3r -> Write();
    //
    KEE_mass_mix_B3f -> Write();
    KEE_fitMass_mix_B3f -> Write();  
    KEE_mass_mix_B3r -> Write();
    KEE_fitMass_mix_B3r -> Write();  
    //
    KEE_mass_Low_B3f -> Write();
    KEE_fitMass_Low_B3f -> Write();  
    KEE_mass_Low_B3r -> Write();
    KEE_fitMass_Low_B3r -> Write();
    
    //KMuMu
    KMM_mass_PF_B3f -> Write();
    KMM_fitMass_PF_B3f -> Write();  
    KMM_mass_PF_B3r -> Write();
    KMM_fitMass_PF_B3r -> Write();    
    
    outHistos.Close();
} 
