/*
source settings.sh
g++ -Wall -o analyzer_RDF `root-config --cflags --glibs ` analyzer_RDF.cpp
./analyzer_RDF --inList input_list.txt --JOBid (1,2..) --outFile /path/to/out_file.root --skTree (0, 1) --testFile /path/to/input_file.root
*/

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>

#include "TChain.h"
#include "TFile.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"

#include <time.h> 
 
#include <ROOT/RDataFrame.hxx>
#include "ROOT/RVec.hxx"

using namespace ROOT::VecOps;
using RV_f = RVec<float>;
using RV_ui = RVec<unsigned int>; 
using RV_b = RVec<bool>;


//----------------- Function declarations -----------------

// Returns indices of triplets passing B and K selections
RV_ui Cuts(unsigned int nTr, RV_f& BpT, RV_f& cosAlpha, RV_f& svprob, RV_f& LxySig, RV_f& K_DCASig, RV_f& KpT);
// Returns indices of good KEE triplets
RV_ui EleCuts(RV_ui& inp_idx, RV_f& l1pT, RV_f& l2pT, RV_f& l1mvaId, RV_f& l2mvaId, RV_b& l1Veto, RV_b& l2Veto);
// Returns indices of good KMuMu triplets
RV_ui MuCuts(RV_ui& inp_idx, RV_f& l1pT, RV_f& l2pT, RV_ui& nTrMu);
// Return indices of triplets with  b1 < mll < b2
RV_ui Bin(RV_ui& inp_idx, double b1, double b2, RV_f& mll);
// Return indices of triplets with both leptons not being PFoverlap
RV_ui noPFover(RV_ui& inp_idx, RV_ui& l1isPFover, RV_ui& l2isPFover);
// Return indices of triplets with both leptons being PF (or LowPt)
RV_ui bothX(RV_ui& inp_idx, RV_ui& l1isX, RV_ui& l2isX);
// Return indices of mixed triplets: (l1_is_PF, l2_is_lowPt) or (l1_is_lowPt, l2_is_PF)
RV_ui mix(RV_ui& inp_idx, RV_ui& l1isPF, RV_ui& l2isPF, RV_ui& l1isLow, RV_ui& l2isLow);



int main(int argc, char **argv){

    time_t start = time(NULL);
    printf(" %s\n", ctime(&start)); 
    
    if(argc < 2) {
        std::cout << " Missing arguments " << std::endl;
        return -1;
    }

    std::string inList = "-1";
    std::string JOBid = "-1";
    std::string outFile = "-1";
    int isSkTree = 0;
    std::string testFile = "-1";
 
    for (int i = 1; i < argc; ++i) {
        if(std::string(argv[i]) == "--inList") {
            if (i + 1 < argc) {
                inList = argv[i+1];
                break;
            } 
            else {
                std::cerr << " --inList option requires one argument " << std::endl;
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
        if(std::string(argv[i]) == "--outFile") {
            if (i + 1 < argc) {
            outFile = argv[i+1];
            break;
            } 
            else {
                std::cerr << " --outFile option requires one argument " << std::endl;
                return 1;
            }
        }
    }
    for (int i = 1; i < argc; ++i) {
      if(std::string(argv[i]) == "--skTree") {
	if (i + 1 < argc) {
	  isSkTree = atoi(argv[i+1]);
	  break;
	}
	else {
	  std::cerr << " --skTree option requires one argument " << std::endl;
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

  
    if(inList != "-1" && (JOBid == "-1" || outFile == "-1")){
        std::cout << " configuration ERROR => splitting file based but missing JOBid and output file path " << std::endl;
        return -1;
    }

    if(inList == "-1" && testFile == "-1"){
        std::cout << " configuration ERROR => need a file list or a test file " << std::endl;
        return -1;
    }

    std::cout << " inList = " << inList << " JOBid = " << JOBid << " outFile = " << outFile << " skTree = " << isSkTree << "\n" << std::endl;


    TChain ch("Events");
  
    if(inList != "-1"){
        std::string rootFileName;
        std::ifstream inFileLong;
        inFileLong.open(inList.c_str(), std::ios::in);
        while(!inFileLong.eof()){
            if(inFileLong >> rootFileName){
	      ch.Add(rootFileName.c_str());
                std::cout << " adding " << rootFileName << std::endl;
            }
        }
    }
    else{
      ch.Add(testFile.c_str());
    }
    

    // Enable multi-threading
    ROOT::EnableImplicitMT();
    
    ROOT::RDataFrame df(ch);

    auto entries = df.Count();
    std::cout << "\n " << *entries << " entries" << std::endl;


    std::string outName = "output_RDF_histo";
    if(JOBid != "-1") outName = outFile;
    else outName += ".root";


    // KEE & KMuMu: removing cross-reference
    auto branch_def = df.Define("K_DCASig", "Take(ProbeTracks_DCASig, BToKEE_kIdx)")
                        .Define("K_pT", "Take(ProbeTracks_pt, BToKEE_kIdx)")
                        .Define("e1_pT", "Take(Electron_pt, BToKEE_l1Idx)")
                        .Define("e2_pT", "Take(Electron_pt, BToKEE_l2Idx)")
                        .Define("e1_Veto", "Take(Electron_convVeto, BToKEE_l1Idx)")
                        .Define("e2_Veto", "Take(Electron_convVeto, BToKEE_l2Idx)")
                        .Define("e1_isPFov", "(ROOT::VecOps::RVec<unsigned int>) Take(Electron_isPFoverlap, BToKEE_l1Idx)")
                        .Define("e2_isPFov", "(ROOT::VecOps::RVec<unsigned int>) Take(Electron_isPFoverlap, BToKEE_l2Idx)")
                        .Define("e1_isPF", "(ROOT::VecOps::RVec<unsigned int>) Take(Electron_isPF, BToKEE_l1Idx)")
                        .Define("e2_isPF", "(ROOT::VecOps::RVec<unsigned int>) Take(Electron_isPF, BToKEE_l2Idx)")
                        .Define("e1_isLow", "(ROOT::VecOps::RVec<unsigned int>) Take(Electron_isLowPt, BToKEE_l1Idx)")
                        .Define("e2_isLow", "(ROOT::VecOps::RVec<unsigned int>) Take(Electron_isLowPt, BToKEE_l2Idx)")
                        .Define("e1_mvaId", "Take(Electron_mvaId, BToKEE_l1Idx)")
                        .Define("e2_mvaId", "Take(Electron_mvaId, BToKEE_l2Idx)")
                        .Define("LxySig", "BToKEE_l_xy/BToKEE_l_xy_unc")
                        .Define("W", "ROOT::VecOps::RVec<float>(nBToKEE, 0.999)")
                        .Define("K_DCASig_mu", "Take(ProbeTracks_DCASig, BToKMuMu_kIdx)")
                        .Define("K_pT_mu", "Take(ProbeTracks_pt, BToKMuMu_kIdx)")
                        .Define("mu1_pT", "Take(Muon_pt, BToKMuMu_l1Idx)")
                        .Define("mu2_pT", "Take(Muon_pt, BToKMuMu_l2Idx)")
                        .Define("mu1_isPF", "(ROOT::VecOps::RVec<unsigned int>) Take(Muon_isPFcand, BToKMuMu_l1Idx)")
                        .Define("mu2_isPF", "(ROOT::VecOps::RVec<unsigned int>) Take(Muon_isPFcand, BToKMuMu_l2Idx)")
                        .Define("mu1_isTr", "(ROOT::VecOps::RVec<int>) Take(Muon_isTriggering, BToKMuMu_l1Idx)")
                        .Define("mu2_isTr", "(ROOT::VecOps::RVec<int>) Take(Muon_isTriggering, BToKMuMu_l2Idx)")
                        .Define("nAddTrMu", "(ROOT::VecOps::RVec<unsigned int>) (nTriggerMuon - mu1_isTr - mu2_isTr)")
                        .Define("LxySig_mu", "BToKMuMu_l_xy/BToKMuMu_l_xy_unc")
                        .Define("W_mu", "ROOT::VecOps::RVec<float>(nBToKMuMu, 0.999)")
                        .Define("b1", "2.9")
                        .Define("b2", "3.3");


    // Finding indices of triplets passing cuts for KEE and KMuMu
    auto ind_cuts = branch_def.Define("Idx_KEE_tmp", Cuts, { "nBToKEE", "BToKEE_pt", "BToKEE_cos2D", "BToKEE_svprob", "LxySig", "K_DCASig", "K_pT" } )
                              .Define("Idx_KEE", EleCuts, { "Idx_KEE_tmp", "e1_pT", "e2_pT", "e1_mvaId", "e2_mvaId", "e1_Veto", "e2_Veto" } )
                              .Define("Idx_KMM_tmp", Cuts, { "nBToKMuMu", "BToKMuMu_pt", "BToKMuMu_cos2D", "BToKMuMu_svprob", "LxySig_mu", "K_DCASig_mu", "K_pT_mu" } )
                              .Define("Idx_KMM", MuCuts, { "Idx_KMM_tmp", "mu1_pT", "mu2_pT", "nAddTrMu" } );
    

    if(!isSkTree){
      
      // Finding indices of triplets for different configurations (All, PF, mix, Low)
      auto config = ind_cuts.Define("Idx_noPFov", noPFover, { "Idx_KEE", "e1_isPFov", "e2_isPFov" } )
	                    .Define("Idx_PF", bothX, { "Idx_KEE", "e1_isPF", "e2_isPF" } )
	                    .Define("Idx_mix", mix, { "Idx_noPFov", "e1_isPF", "e2_isPF", "e1_isLow", "e2_isLow" } )
	                    .Define("Idx_Low", bothX, { "Idx_noPFov", "e1_isLow", "e2_isLow" } )
	                    .Define("Idx_PF_mu", bothX, { "Idx_KMM", "mu1_isPF", "mu2_isPF" } );
      

      
      // Bin ranges
      // B0 = [0, 1], B1 = [1, 2.5], B2 = [2.5, 2.9], B3 = [2.9, 3.3], B4 = [3.3, 3.58], B5 = [3.58, 100]
      // Finding indices of triplets with 2.9 < mll_fullfit(_raw) < 3.1
      auto config_bin = config.Define("noPFov_B3f", Bin, { "Idx_noPFov", "b1", "b2", "BToKEE_mll_fullfit" } )
	                      .Define("noPFov_B3r", Bin, { "Idx_noPFov", "b1", "b2", "BToKEE_mll_raw" } )
	                      .Define("PF_B3f", Bin, { "Idx_PF", "b1", "b2", "BToKEE_mll_fullfit" } )
	                      .Define("PF_B3r", Bin, { "Idx_PF", "b1", "b2", "BToKEE_mll_raw" } )
	                      .Define("mix_B3f", Bin, { "Idx_mix", "b1", "b2", "BToKEE_mll_fullfit" } )
	                      .Define("mix_B3r", Bin, { "Idx_mix", "b1", "b2", "BToKEE_mll_raw" } ) 
	                      .Define("Low_B3f", Bin, { "Idx_Low", "b1", "b2", "BToKEE_mll_fullfit" } )
	                      .Define("Low_B3r", Bin, { "Idx_Low", "b1", "b2", "BToKEE_mll_raw" } )
          	              .Define("PF_B3f_mu", Bin, { "Idx_PF_mu", "b1", "b2", "BToKMuMu_mll_fullfit" } )
	                      .Define("PF_B3r_mu", Bin, { "Idx_PF_mu", "b1", "b2", "BToKMuMu_mll_raw" } );       


      // Defining quantities for histograms
      auto config_bin_q = config_bin.Define("KEE_mass_B3f","Take(BToKEE_mass, noPFov_B3f)")
	                            .Define("KEE_fit_mass_B3f","Take(BToKEE_fit_mass, noPFov_B3f)")
	                            .Define("W_B3f","Take(W, noPFov_B3f)")
	                            .Define("KEE_mass_B3r","Take(BToKEE_mass, noPFov_B3r)")
	                            .Define("KEE_fit_mass_B3r","Take(BToKEE_fit_mass, noPFov_B3r)")
	                            .Define("W_B3r","Take(W, noPFov_B3r)")
	                            .Define("KEE_mass_PF_B3f","Take(BToKEE_mass, PF_B3f)")
	                            .Define("KEE_fit_mass_PF_B3f","Take(BToKEE_fit_mass, PF_B3f)")
	                            .Define("W_PF_B3f","Take(W, PF_B3f)")
	                            .Define("KEE_mass_PF_B3r","Take(BToKEE_mass, PF_B3r)")
	                            .Define("KEE_fit_mass_PF_B3r","Take(BToKEE_fit_mass, PF_B3r)")
	                            .Define("W_PF_B3r","Take(W, PF_B3r)")
	                            .Define("KEE_mass_mix_B3f","Take(BToKEE_mass, mix_B3f)")
	                            .Define("KEE_fit_mass_mix_B3f","Take(BToKEE_fit_mass, mix_B3f)")
	                            .Define("W_mix_B3f","Take(W, mix_B3f)")
	                            .Define("KEE_mass_mix_B3r","Take(BToKEE_mass, mix_B3r)")
	                            .Define("KEE_fit_mass_mix_B3r","Take(BToKEE_fit_mass, mix_B3r)")
	                            .Define("W_mix_B3r","Take(W, mix_B3r)")
	                            .Define("KEE_mass_Low_B3f","Take(BToKEE_mass, Low_B3f)")
	                            .Define("KEE_fit_mass_Low_B3f","Take(BToKEE_fit_mass, Low_B3f)")
	                            .Define("W_Low_B3f","Take(W, Low_B3f)")
	                            .Define("KEE_mass_Low_B3r","Take(BToKEE_mass, Low_B3r)")
	                            .Define("KEE_fit_mass_Low_B3r","Take(BToKEE_fit_mass, Low_B3r)")
	                            .Define("W_Low_B3r","Take(W, Low_B3r)")
	                            .Define("KMM_mass_PF_B3f","Take(BToKMuMu_mass, PF_B3f_mu)")
	                            .Define("KMM_fit_mass_PF_B3f","Take(BToKMuMu_fit_mass, PF_B3f_mu)")
	                            .Define("W_PF_B3f_mu","Take(W_mu, PF_B3f_mu)")
	                            .Define("KMM_mass_PF_B3r","Take(BToKMuMu_mass, PF_B3r_mu)")
	                            .Define("KMM_fit_mass_PF_B3r","Take(BToKMuMu_fit_mass, PF_B3r_mu)")
	                            .Define("W_PF_B3r_mu","Take(W_mu, PF_B3r_mu)");

      
      // Save histograms
      // KEE histograms
      auto h_KEE_mass_B3f = config_bin_q.Histo1D( {"KEE_mass_B3f", "", 75, 4.5, 6.0}, "KEE_mass_B3f", "W_B3f");
      auto h_KEE_fitMass_B3f = config_bin_q.Histo1D( {"KEE_fitMass_B3f", "", 75, 4.5, 6.0}, "KEE_fit_mass_B3f", "W_B3f");
      auto h_KEE_mass_B3r = config_bin_q.Histo1D( {"KEE_mass_B3r", "", 75, 4.5, 6.0}, "KEE_mass_B3r", "W_B3r");
      auto h_KEE_fitMass_B3r = config_bin_q.Histo1D( {"KEE_fitMass_B3r", "", 75, 4.5, 6.0}, "KEE_fit_mass_B3r", "W_B3r");
      auto h_KEE_mass_PF_B3f = config_bin_q.Histo1D( {"KEE_mass_PF_B3f", "", 75, 4.5, 6.0}, "KEE_mass_PF_B3f", "W_PF_B3f");
      auto h_KEE_fitMass_PF_B3f = config_bin_q.Histo1D( {"KEE_fitMass_PF_B3f", "", 75, 4.5, 6.0}, "KEE_fit_mass_PF_B3f", "W_PF_B3f");
      auto h_KEE_mass_PF_B3r = config_bin_q.Histo1D( {"KEE_mass_PF_B3r", "", 75, 4.5, 6.0}, "KEE_mass_PF_B3r", "W_PF_B3r");
      auto h_KEE_fitMass_PF_B3r = config_bin_q.Histo1D( {"KEE_fitMass_PF_B3r", "", 75, 4.5, 6.0}, "KEE_fit_mass_PF_B3r", "W_PF_B3r");
      auto h_KEE_mass_mix_B3f = config_bin_q.Histo1D( {"KEE_mass_mix_B3f", "", 75, 4.5, 6.0}, "KEE_mass_mix_B3f", "W_mix_B3f");
      auto h_KEE_fitMass_mix_B3f = config_bin_q.Histo1D( {"KEE_fitMass_mix_B3f", "", 75, 4.5, 6.0}, "KEE_fit_mass_mix_B3f", "W_mix_B3f");
      auto h_KEE_mass_mix_B3r = config_bin_q.Histo1D( {"KEE_mass_mix_B3r", "", 75, 4.5, 6.0}, "KEE_mass_mix_B3r", "W_mix_B3r");
      auto h_KEE_fitMass_mix_B3r = config_bin_q.Histo1D( {"KEE_fitMass_mix_B3r", "", 75, 4.5, 6.0}, "KEE_fit_mass_mix_B3r", "W_mix_B3r");
      auto h_KEE_mass_Low_B3f = config_bin_q.Histo1D( {"KEE_mass_Low_B3f", "", 75, 4.5, 6.0}, "KEE_mass_Low_B3f", "W_Low_B3f");
      auto h_KEE_fitMass_Low_B3f = config_bin_q.Histo1D( {"KEE_fitMass_Low_B3f", "", 75, 4.5, 6.0}, "KEE_fit_mass_Low_B3f", "W_Low_B3f");
      auto h_KEE_mass_Low_B3r = config_bin_q.Histo1D( {"KEE_mass_Low_B3r", "", 75, 4.5, 6.0}, "KEE_mass_Low_B3r", "W_Low_B3r");
      auto h_KEE_fitMass_Low_B3r = config_bin_q.Histo1D( {"KEE_fitMass_Low_B3r", "", 75, 4.5, 6.0}, "KEE_fit_mass_Low_B3r", "W_Low_B3r");
      // KMuMu histograms
      auto h_KMM_mass_PF_B3f = config_bin_q.Histo1D( {"KMM_mass_PF_B3f", "", 75, 4.5, 6.0}, "KMM_mass_PF_B3f", "W_PF_B3f_mu");
      auto h_KMM_fitMass_PF_B3f = config_bin_q.Histo1D( {"KMM_fitMass_PF_B3f", "", 75, 4.5, 6.0}, "KMM_fit_mass_PF_B3f", "W_PF_B3f_mu");
      auto h_KMM_mass_PF_B3r = config_bin_q.Histo1D( {"KMM_mass_PF_B3r", "", 75, 4.5, 6.0}, "KMM_mass_PF_B3r", "W_PF_B3r_mu");
      auto h_KMM_fitMass_PF_B3r = config_bin_q.Histo1D( {"KMM_fitMass_PF_B3r", "", 75, 4.5, 6.0}, "KMM_fit_mass_PF_B3r", "W_PF_B3r_mu");


      TFile outHistos(outName.c_str(), "recreate");

      outHistos.cd();

      h_KEE_mass_B3f->Write();
      h_KEE_fitMass_B3f->Write();
      h_KEE_mass_B3r->Write();
      h_KEE_fitMass_B3r->Write();
      h_KEE_mass_PF_B3f->Write();
      h_KEE_fitMass_PF_B3f->Write();
      h_KEE_mass_PF_B3r->Write();
      h_KEE_fitMass_PF_B3r->Write();
      h_KEE_mass_mix_B3f->Write();
      h_KEE_fitMass_mix_B3f->Write();
      h_KEE_mass_mix_B3r->Write();
      h_KEE_fitMass_mix_B3r->Write();
      h_KEE_mass_Low_B3f->Write();
      h_KEE_fitMass_Low_B3f->Write();
      h_KEE_mass_Low_B3r->Write();
      h_KEE_fitMass_Low_B3r->Write();
      h_KMM_mass_PF_B3f->Write();
      h_KMM_fitMass_PF_B3f->Write();
      h_KMM_mass_PF_B3r->Write();
      h_KMM_fitMass_PF_B3r->Write();
      
      outHistos.Close();
      
    }
    else{

      // Skimmed quantities
      auto tree_q = ind_cuts.Define("e1_isPF_sk","Take(e1_isPF, Idx_KEE)")
	                    .Define("e2_isPF_sk","Take(e2_isPF, Idx_KEE)")
       	                    .Define("e1_isPFov_sk","Take(e1_isPFov, Idx_KEE)")
	                    .Define("e2_isPFov_sk","Take(e2_isPFov, Idx_KEE)")
	                    .Define("e1_isLow_sk","Take(e1_isLow, Idx_KEE)")
	                    .Define("e2_isLow_sk","Take(e2_isLow, Idx_KEE)")
	                    .Define("KEE_mll_fullfit_sk","Take(BToKEE_mll_fullfit, Idx_KEE)")
	                    .Define("KEE_mll_raw_sk","Take(BToKEE_mll_raw, Idx_KEE)")
	                    .Define("KEE_mass_sk","Take(BToKEE_mass, Idx_KEE)")
  	                    .Define("KEE_fit_mass_sk","Take(BToKEE_fit_mass, Idx_KEE)")
	                    .Define("mu1_isPF_sk","Take(mu1_isPF, Idx_KMM)")
	                    .Define("mu2_isPF_sk","Take(mu2_isPF, Idx_KMM)")
	                    .Define("KMM_mll_fullfit_sk","Take(BToKMuMu_mll_fullfit, Idx_KMM)")
	                    .Define("KMM_mll_raw_sk","Take(BToKMuMu_mll_raw, Idx_KMM)")
	                    .Define("KMM_mass_sk","Take(BToKMuMu_mass, Idx_KMM)")
	                    .Define("KMM_fit_mass_sk","Take(BToKMuMu_fit_mass, Idx_KMM)");

      tree_q.Snapshot("skimTree", outName, "\\b([^ ]*)(_sk)|(Idx_)([^ ]*)"); 
    }

    time_t finish = time(NULL);
    printf("\n %s\n", ctime(&finish));

}  



//----------------- Function definitions -----------------

// Returns indices of triplets passing B and K selections                                                     
RV_ui Cuts(unsigned int nTr,
           RV_f& BpT,
           RV_f& cosAlpha,
           RV_f& svprob,
           RV_f& LxySig,
           RV_f& K_DCASig,
           RV_f& KpT){

  RV_ui out_idx;
  for (auto i=0; i<nTr; ++i){
    if( BpT[i] < 3. || cosAlpha[i] < 0.999 || svprob[i] < 0.1 ||
        LxySig[i] < 6.0 || K_DCASig[i] < 2. || KpT[i] < 3.) continue;
    out_idx.push_back(i);
  }
  return out_idx;
} 


// Returns indices of good KEE triplets                                                 
RV_ui EleCuts(RV_ui& inp_idx,
              RV_f& l1pT,
              RV_f& l2pT,
              RV_f& l1mvaId,
              RV_f& l2mvaId,
              RV_b& l1Veto,
              RV_b& l2Veto){

  RV_ui out_idx;
  for (auto i : inp_idx){
    if( l1mvaId[i] > 3.94 && l2mvaId[i] > 3.94 &&
                l1Veto[i] && l2Veto[i] &&
        l1pT[i] > 1.5 && l2pT[i] > 0.5) out_idx.push_back(i);
  }

  return out_idx;
}


// Returns indices of good KMuMu triplets
RV_ui MuCuts(RV_ui& inp_idx, 
	     RV_f& l1pT,
	     RV_f& l2pT, 
	     RV_ui& nTrMu){

  RV_ui out_idx;
  for (auto i : inp_idx){
    if( l1pT[i] > 1.5 && l2pT[i] > 0.5 && nTrMu[i] > 0) out_idx.push_back(i);
  }
  return out_idx;
}


// Return indices of triplets with  b1 < mll < b2
RV_ui Bin(RV_ui& inp_idx,
          double b1,double b2,
          RV_f& mll){

  RV_ui out_idx;
  for (auto i : inp_idx){
    if( mll[i] >= b1 && mll[i] < b2 )out_idx.push_back(i);
  }
  return out_idx;
}


// Return indices of triplets with both leptons not being PFoverlap                                                    
RV_ui noPFover(RV_ui& inp_idx,
               RV_ui& l1isPFover,
               RV_ui& l2isPFover){

  RV_ui out_idx;
  for (auto i : inp_idx){
    if( l1isPFover[i] == 0 && l2isPFover[i] == 0 ) out_idx.push_back(i);
  }
  return out_idx;
}


// Return indices of triplets with both leptons being PF (or LowPt)                                                    
RV_ui bothX(RV_ui& inp_idx,
            RV_ui& l1isX,
            RV_ui& l2isX){

  RV_ui out_idx;
  for (auto i : inp_idx){
    if( l1isX[i] == 1 && l2isX[i] == 1 ) out_idx.push_back(i);
  }
  return out_idx;
}


// Return indices of mixed triplets: (l1_is_PF, l2_is_lowPt) or (l1_is_lowPt, l2_is_PF)                    
RV_ui mix(RV_ui& inp_idx,
          RV_ui& l1isPF,
          RV_ui& l2isPF,
          RV_ui& l1isLow,
          RV_ui& l2isLow){

  RV_ui out_idx;
  for (auto i : inp_idx){
    if( (l1isPF[i] == 1 && l2isLow[i] ==1) || (l1isLow[i] ==1 && l2isPF[i] == 1) ) out_idx.push_back(i);
  }
  return out_idx;
}
