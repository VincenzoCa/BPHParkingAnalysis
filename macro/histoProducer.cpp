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


    gROOT->Reset();
    gROOT->Macro("./setStyle.C");
    gSystem->Load("libRooFit") ;
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    TChain* ch = new TChain("skimTree");
  
    if(ntupleList != "-1"){
        std::string rootFileName;
        std::ifstream inFileLong;
        inFileLong.open(ntupleList.c_str(), std::ios::in);
        while(!inFileLong.eof()){
            if(inFileLong >> rootFileName){
                ch->Add(rootFileName.c_str());
                std::cout << " adding " << rootFileName << std::endl;
            }
        }
    }
    else{
        ch->Add(testFile.c_str());
    }    
    
    
    TCut PFov = "l1_isPFov_cut == 0 && l2_isPFov_cut == 0";
    TCut PF = "l1_isPF_cut == 1 && l2_isPF_cut == 1";
    TCut mix = "((l1_isPF_cut == 1 && l2_isLow_cut == 1) || (l1_isLow_cut == 1 && l2_isPF_cut == 1))"; 
    TCut Low = "l1_isLow_cut == 1 && l2_isLow_cut == 1";  
    TCut mllraw = "mll_raw_cut > 2.9 && mll_raw_cut < 3.3";
    TCut mllff = "mll_fullfit_cut > 2.9 && mll_fullfit_cut < 3.3";
  
  
    TH1F *Bmass_B3f = new TH1F("Bmass_B3f","", 75, 4.5, 6.);
    Bmass_B3f->Sumw2();
    Bmass_B3f->SetLineColor(kRed);
    Bmass_B3f->SetLineWidth(2);    
    
    TH1F *BfitMass_B3f = new TH1F("BfitMass_B3f","", 75, 4.5, 6.);
    BfitMass_B3f->Sumw2();  
    BfitMass_B3f->SetLineColor(kRed);
    BfitMass_B3f->SetLineWidth(2);       
    
    TH1F *Bmass_B3r = new TH1F("Bmass_B3r","", 75, 4.5, 6.);
    Bmass_B3r->Sumw2();
    Bmass_B3r->SetLineColor(kRed);
    Bmass_B3r->SetLineWidth(2);       
    
    TH1F *BfitMass_B3r = new TH1F("BfitMass_B3r","", 75, 4.5, 6.);
    BfitMass_B3r->Sumw2(); 
    BfitMass_B3r->SetLineColor(kRed);
    BfitMass_B3r->SetLineWidth(2);       
    
    TH1F *Bmass_PF_B3f = new TH1F("Bmass_PF_B3f","", 75, 4.5, 6.);
    Bmass_PF_B3f->Sumw2();
    Bmass_PF_B3f->SetLineColor(kRed);
    Bmass_PF_B3f->SetLineWidth(2);       
    
    TH1F *BfitMass_PF_B3f = new TH1F("BfitMass_PF_B3f","", 75, 4.5, 6.);
    BfitMass_PF_B3f->Sumw2();  
    BfitMass_PF_B3f->SetLineColor(kRed);
    BfitMass_PF_B3f->SetLineWidth(2);       
    
    TH1F *Bmass_PF_B3r = new TH1F("Bmass_PF_B3r","", 75, 4.5, 6.);
    Bmass_PF_B3r->Sumw2();
    Bmass_PF_B3r->SetLineColor(kRed);
    Bmass_PF_B3r->SetLineWidth(2);       
    
    TH1F *BfitMass_PF_B3r = new TH1F("BfitMass_PF_B3r","", 75, 4.5, 6.);
    BfitMass_PF_B3r->Sumw2();
    BfitMass_PF_B3r->SetLineColor(kRed);
    BfitMass_PF_B3r->SetLineWidth(2);       
    
    TH1F *Bmass_mix_B3f = new TH1F("Bmass_mix_B3f","", 75, 4.5, 6.);
    Bmass_mix_B3f->Sumw2();
    Bmass_mix_B3f->SetLineColor(kRed);
    Bmass_mix_B3f->SetLineWidth(2);       
    
    TH1F *BfitMass_mix_B3f = new TH1F("BfitMass_mix_B3f","", 75, 4.5, 6.);
    BfitMass_mix_B3f->Sumw2();  
    BfitMass_mix_B3f->SetLineColor(kRed);
    BfitMass_mix_B3f->SetLineWidth(2);       
    
    TH1F *Bmass_mix_B3r = new TH1F("Bmass_mix_B3r","", 75, 4.5, 6.);
    Bmass_mix_B3r->Sumw2();
    Bmass_mix_B3r->SetLineColor(kRed);
    Bmass_mix_B3r->SetLineWidth(2);       
    
    TH1F *BfitMass_mix_B3r = new TH1F("BfitMass_mix_B3r","", 75, 4.5, 6.);
    BfitMass_mix_B3r->Sumw2();
    BfitMass_mix_B3r->SetLineColor(kRed);
    BfitMass_mix_B3r->SetLineWidth(2);       
    
    TH1F *Bmass_Low_B3f = new TH1F("Bmass_Low_B3f","", 75, 4.5, 6.);
    Bmass_Low_B3f->Sumw2();
    Bmass_Low_B3f->SetLineColor(kRed);
    Bmass_Low_B3f->SetLineWidth(2);   
    
    TH1F *BfitMass_Low_B3f = new TH1F("BfitMass_Low_B3f","", 75, 4.5, 6.);
    BfitMass_Low_B3f->Sumw2();
    BfitMass_Low_B3f->SetLineColor(kRed);
    BfitMass_Low_B3f->SetLineWidth(2);       
    
    TH1F *Bmass_Low_B3r = new TH1F("Bmass_Low_B3r","", 75, 4.5, 6.);
    Bmass_Low_B3r->Sumw2();
    Bmass_Low_B3r->SetLineColor(kRed);
    Bmass_Low_B3r->SetLineWidth(2);       
    
    TH1F *BfitMass_Low_B3r = new TH1F("BfitMass_Low_B3r","", 75, 4.5, 6.);
    BfitMass_Low_B3r->Sumw2();
    BfitMass_Low_B3r->SetLineColor(kRed);
    BfitMass_Low_B3r->SetLineWidth(2);       
   

    ch->Draw("B_mass_cut>>Bmass_B3f", PFov && mllff,"goff");
    ch->Draw("B_fit_mass_cut>>BfitMass_B3f", PFov && mllff,"goff");
    ch->Draw("B_mass_cut>>Bmass_B3r", PFov && mllraw,"goff");
    ch->Draw("B_fit_mass_cut>>BfitMass_B3r", PFov && mllraw,"goff");
    //    
    ch->Draw("B_mass_cut>>Bmass_PF_B3f", PF && mllff,"goff");
    ch->Draw("B_fit_mass_cut>>BfitMass_PF_B3f", PF && mllff,"goff");
    ch->Draw("B_mass_cut>>Bmass_PF_B3r", PF && mllraw,"goff");
    ch->Draw("B_fit_mass_cut>>BfitMass_PF_B3r", PF && mllraw,"goff");
    //    
    ch->Draw("B_mass_cut>>Bmass_mix_B3f", PFov && mix && mllff,"goff");
    ch->Draw("B_fit_mass_cut>>BfitMass_mix_B3f", PFov && mix && mllff,"goff");
    ch->Draw("B_mass_cut>>Bmass_mix_B3r", PFov && mix && mllraw,"goff");
    ch->Draw("B_fit_mass_cut>>BfitMass_mix_B3r", PFov && mix && mllraw,"goff");  
    //
    ch->Draw("B_mass_cut>>Bmass_Low_B3f", PFov && Low && mllff,"goff");
    ch->Draw("B_fit_mass_cut>>BfitMass_Low_B3f", PFov && Low && mllff,"goff");
    ch->Draw("B_mass_cut>>Bmass_Low_B3r", PFov && Low && mllraw,"goff");
    ch->Draw("B_fit_mass_cut>>BfitMass_Low_B3r", PFov && Low && mllraw,"goff");  
    

    std::string outName = "histo_output_from_tree";
    if(JOBid != "-1") outName = outputFolder;
    else outName += ".root";
    TFile outHistos(outName.c_str(), "recreate");
    
    outHistos.cd();
  
    Bmass_B3f->Write();
    BfitMass_B3f->Write();  
    Bmass_B3r->Write();
    BfitMass_B3r->Write();
    //
    Bmass_PF_B3f->Write();
    BfitMass_PF_B3f->Write();  
    Bmass_PF_B3r->Write();
    BfitMass_PF_B3r->Write();
    //
    Bmass_mix_B3f->Write();
    BfitMass_mix_B3f->Write();  
    Bmass_mix_B3r->Write();
    BfitMass_mix_B3r->Write();  
    //
    Bmass_Low_B3f->Write();
    BfitMass_Low_B3f->Write();  
    Bmass_Low_B3r->Write();
    BfitMass_Low_B3r->Write();  
    
    outHistos.Close();
} 
