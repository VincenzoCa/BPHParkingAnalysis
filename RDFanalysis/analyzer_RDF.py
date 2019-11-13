# source settings.sh
# python analyzer_RDF.py --inList input_list.txt --JOBid (1,2,...) --outFile /path/output.root
#
# By default, the output will be stored as histograms. Please add --tree if you want a tree as output
# The output consists of both KEE and KMuMu quantities

import ROOT
import module as Mod
# Enable multi-threading
ROOT.ROOT.EnableImplicitMT()
import argparse
import numpy as np

Mod.print_time('Start')

parser = argparse.ArgumentParser()
parser.add_argument("--inList", "-l", help="set list")
parser.add_argument("--JOBid", "-j", help="set job ID")
parser.add_argument("--outFile", "-f", help="set output folder")
parser.add_argument('--tree', action='store_true')
args = parser.parse_args()

listName = args.inList
outFileName = args.outFile
isTree = args.tree
print('\n --list {} --JOBid {} --outFile {} --tree {}'.format(listName, args.JOBid, outFileName, str(isTree)))
    
# Create dataframe from NanoAOD files 
files = ROOT.std.vector("string")()
mylist = [line.rstrip('\n') for line in open(listName)]
print('\n Input files:\n'+'\n'.join(map(str, mylist))+'\n')
for f in mylist : files.push_back(f)
df = ROOT.ROOT.RDataFrame("Events", files)

entries = df.Count()
print('\n {} entries '.format(entries.GetValue()))


ROOT.gInterpreter.Declare(Mod.Cuts_code)
ROOT.gInterpreter.Declare(Mod.EleCuts_code)
ROOT.gInterpreter.Declare(Mod.MuCuts_code)
ROOT.gInterpreter.Declare(Mod.Bin_code)
ROOT.gInterpreter.Declare(Mod.noPFover_code)
ROOT.gInterpreter.Declare(Mod.bothX_code)
ROOT.gInterpreter.Declare(Mod.mix_code)


# KEE & KMuMu: removing cross-reference
branch_def = df.Define("K_DCASig", "Take(ProbeTracks_DCASig, BToKEE_kIdx)") \
               .Define("K_pT", "Take(ProbeTracks_pt, BToKEE_kIdx)") \
               .Define("e1_pT", "Take(Electron_pt, BToKEE_l1Idx)") \
               .Define("e2_pT", "Take(Electron_pt, BToKEE_l2Idx)") \
               .Define("e1_Veto", "Take(Electron_convVeto, BToKEE_l1Idx)") \
               .Define("e2_Veto", "Take(Electron_convVeto, BToKEE_l2Idx)") \
               .Define("e1_isPFov", "(RVec<unsigned int>) Take(Electron_isPFoverlap, BToKEE_l1Idx)") \
               .Define("e2_isPFov", "(RVec<unsigned int>) Take(Electron_isPFoverlap, BToKEE_l2Idx)") \
               .Define("e1_isPF", "(RVec<unsigned int>) Take(Electron_isPF, BToKEE_l1Idx)") \
               .Define("e2_isPF", "(RVec<unsigned int>) Take(Electron_isPF, BToKEE_l2Idx)") \
               .Define("e1_isLow", "(RVec<unsigned int>) Take(Electron_isLowPt, BToKEE_l1Idx)") \
               .Define("e2_isLow", "(RVec<unsigned int>) Take(Electron_isLowPt, BToKEE_l2Idx)") \
               .Define("e1_mvaId", "Take(Electron_mvaId, BToKEE_l1Idx)") \
               .Define("e2_mvaId", "Take(Electron_mvaId, BToKEE_l2Idx)") \
               .Define("LxySig", "BToKEE_l_xy/BToKEE_l_xy_unc") \
               .Define("W", "RVec<float>(nBToKEE, 0.999)") \
               .Define("K_DCASig_mu", "Take(ProbeTracks_DCASig, BToKMuMu_kIdx)") \
               .Define("K_pT_mu", "Take(ProbeTracks_pt, BToKMuMu_kIdx)") \
               .Define("mu1_pT", "Take(Muon_pt, BToKMuMu_l1Idx)") \
               .Define("mu2_pT", "Take(Muon_pt, BToKMuMu_l2Idx)") \
               .Define("mu1_isPF", "(RVec<unsigned int>) Take(Muon_isPFcand, BToKMuMu_l1Idx)") \
               .Define("mu2_isPF", "(RVec<unsigned int>) Take(Muon_isPFcand, BToKMuMu_l2Idx)") \
               .Define("mu1_isTr", "(RVec<int>) Take(Muon_isTriggering, BToKMuMu_l1Idx)") \
               .Define("mu2_isTr", "(RVec<int>) Take(Muon_isTriggering, BToKMuMu_l2Idx)") \
               .Define("nAddTrMu", "(RVec<unsigned int>) (nTriggerMuon - mu1_isTr - mu2_isTr)") \
               .Define("LxySig_mu", "BToKMuMu_l_xy/BToKMuMu_l_xy_unc") \
               .Define("W_mu", "RVec<float>(nBToKMuMu, 0.999)")


# Finding indices of triplets passing cuts for KEE and KMuMu
ind_cuts = branch_def.Define("tmp_Idx_KEE", "Cuts( nBToKEE, BToKEE_pt, BToKEE_cos2D, BToKEE_svprob, LxySig, K_DCASig, K_pT )") \
                     .Define("Idx_KEE", "EleCuts( tmp_Idx_KEE, e1_pT, e2_pT, e1_mvaId, e2_mvaId, e1_Veto, e2_Veto )") \
                     .Define("tmp_Idx_KMM", "Cuts( nBToKMuMu, BToKMuMu_pt, BToKMuMu_cos2D, BToKMuMu_svprob, LxySig_mu, K_DCASig_mu, K_pT_mu )") \
                     .Define("Idx_KMM", "MuCuts( tmp_Idx_KMM, mu1_pT, mu2_pT, nAddTrMu )")


if  np.logical_not(isTree):


    # Finding indices of triplets for different configurations (All, PF, mix, Low)
    config = ind_cuts.Define("Idx_noPFov", "noPFover( Idx_KEE, e1_isPFov, e2_isPFov )" ) \
                     .Define("Idx_PF", "bothX( Idx_KEE, e1_isPF, e2_isPF )" ) \
                     .Define("Idx_mix", "mix( Idx_noPFov, e1_isPF, e2_isPF, e1_isLow, e2_isLow )" ) \
                     .Define("Idx_Low", "bothX( Idx_noPFov, e1_isLow, e2_isLow )" ) \
                     .Define("Idx_PF_mu", "bothX( Idx_KMM, mu1_isPF, mu2_isPF )" )


    # Bin ranges
    # B0 = [0, 1], B1 = [1, 2.5], B2 = [2.5, 2.9], B3 = [2.9, 3.3], B4 = [3.3, 3.58], B5 = [3.58, 100]
    # Finding indices of triplets with 2.9 < mll_fullfit(_raw) < 3.1
    config_bin = config.Define("noPFov_B3f", "Bin( Idx_noPFov, 2.9, 3.3, BToKEE_mll_fullfit )" ) \
                       .Define("noPFov_B3r", "Bin( Idx_noPFov, 2.9, 3.3, BToKEE_mll_raw )" ) \
                       .Define("PF_B3f", "Bin( Idx_PF, 2.9, 3.3, BToKEE_mll_fullfit )" ) \
                       .Define("PF_B3r", "Bin( Idx_PF, 2.9, 3.3, BToKEE_mll_raw )" ) \
                       .Define("mix_B3f", "Bin( Idx_mix, 2.9, 3.3, BToKEE_mll_fullfit )" ) \
                       .Define("mix_B3r", "Bin( Idx_mix, 2.9, 3.3, BToKEE_mll_raw )" ) \
                       .Define("Low_B3f", "Bin( Idx_Low, 2.9, 3.3, BToKEE_mll_fullfit )" ) \
                       .Define("Low_B3r", "Bin( Idx_Low, 2.9, 3.3, BToKEE_mll_raw )" ) \
                       .Define("PF_B3f_mu", "Bin( Idx_PF_mu, 2.9, 3.3, BToKMuMu_mll_fullfit )" ) \
                       .Define("PF_B3r_mu", "Bin( Idx_PF_mu, 2.9, 3.3, BToKMuMu_mll_raw )" ) 


    # Defining quantities for histograms
    config_bin_q = config_bin.Define("KEE_mass_B3f","Take(BToKEE_mass, noPFov_B3f)") \
                             .Define("KEE_fit_mass_B3f","Take(BToKEE_fit_mass, noPFov_B3f)") \
                             .Define("W_B3f","Take(W, noPFov_B3f)") \
                             .Define("KEE_mass_B3r","Take(BToKEE_mass, noPFov_B3r)") \
                             .Define("KEE_fit_mass_B3r","Take(BToKEE_fit_mass, noPFov_B3r)") \
                             .Define("W_B3r","Take(W, noPFov_B3r)") \
                             .Define("KEE_mass_PF_B3f","Take(BToKEE_mass, PF_B3f)") \
                             .Define("KEE_fit_mass_PF_B3f","Take(BToKEE_fit_mass, PF_B3f)") \
                             .Define("W_PF_B3f","Take(W, PF_B3f)") \
                             .Define("KEE_mass_PF_B3r","Take(BToKEE_mass, PF_B3r)") \
                             .Define("KEE_fit_mass_PF_B3r","Take(BToKEE_fit_mass, PF_B3r)") \
                             .Define("W_PF_B3r","Take(W, PF_B3r)") \
                             .Define("KEE_mass_mix_B3f","Take(BToKEE_mass, mix_B3f)") \
                             .Define("KEE_fit_mass_mix_B3f","Take(BToKEE_fit_mass, mix_B3f)") \
                             .Define("W_mix_B3f","Take(W, mix_B3f)") \
                             .Define("KEE_mass_mix_B3r","Take(BToKEE_mass, mix_B3r)") \
                             .Define("KEE_fit_mass_mix_B3r","Take(BToKEE_fit_mass, mix_B3r)") \
                             .Define("W_mix_B3r","Take(W, mix_B3r)") \
                             .Define("KEE_mass_Low_B3f","Take(BToKEE_mass, Low_B3f)") \
                             .Define("KEE_fit_mass_Low_B3f","Take(BToKEE_fit_mass, Low_B3f)") \
                             .Define("W_Low_B3f","Take(W, Low_B3f)") \
                             .Define("KEE_mass_Low_B3r","Take(BToKEE_mass, Low_B3r)") \
                             .Define("KEE_fit_mass_Low_B3r","Take(BToKEE_fit_mass, Low_B3r)") \
                             .Define("W_Low_B3r","Take(W, Low_B3r)") \
                             .Define("KMM_mass_PF_B3f","Take(BToKMuMu_mass, PF_B3f_mu)") \
                             .Define("KMM_fit_mass_PF_B3f","Take(BToKMuMu_fit_mass, PF_B3f_mu)") \
                             .Define("W_PF_B3f_mu","Take(W_mu, PF_B3f_mu)") \
                             .Define("KMM_mass_PF_B3r","Take(BToKMuMu_mass, PF_B3r_mu)") \
                             .Define("KMM_fit_mass_PF_B3r","Take(BToKMuMu_fit_mass, PF_B3r_mu)") \
                             .Define("W_PF_B3r_mu","Take(W_mu, PF_B3r_mu)")


    # Save histograms
    Mod.print_time('Before histo')
    # KEE histograms
    h_KEE_mass_B3f = config_bin_q.Histo1D( ("KEE_mass_B3f", "", 75, 4.5, 6.0), "KEE_mass_B3f", "W_B3f")
    h_KEE_fitMass_B3f = config_bin_q.Histo1D( ("KEE_fitMass_B3f", "", 75, 4.5, 6.0), "KEE_fit_mass_B3f", "W_B3f")
    h_KEE_mass_B3r = config_bin_q.Histo1D( ("KEE_mass_B3r", "", 75, 4.5, 6.0), "KEE_mass_B3r", "W_B3r")
    h_KEE_fitMass_B3r = config_bin_q.Histo1D( ("KEE_fitMass_B3r", "", 75, 4.5, 6.0), "KEE_fit_mass_B3r", "W_B3r")
    h_KEE_mass_PF_B3f = config_bin_q.Histo1D( ("KEE_mass_PF_B3f", "", 75, 4.5, 6.0), "KEE_mass_PF_B3f", "W_PF_B3f")
    h_KEE_fitMass_PF_B3f = config_bin_q.Histo1D( ("KEE_fitMass_PF_B3f", "", 75, 4.5, 6.0), "KEE_fit_mass_PF_B3f", "W_PF_B3f")
    h_KEE_mass_PF_B3r = config_bin_q.Histo1D( ("KEE_mass_PF_B3r", "", 75, 4.5, 6.0), "KEE_mass_PF_B3r", "W_PF_B3r")
    h_KEE_fitMass_PF_B3r = config_bin_q.Histo1D( ("KEE_fitMass_PF_B3r", "", 75, 4.5, 6.0), "KEE_fit_mass_PF_B3r", "W_PF_B3r")
    h_KEE_mass_mix_B3f = config_bin_q.Histo1D( ("KEE_mass_mix_B3f", "", 75, 4.5, 6.0), "KEE_mass_mix_B3f", "W_mix_B3f")
    h_KEE_fitMass_mix_B3f = config_bin_q.Histo1D( ("KEE_fitMass_mix_B3f", "", 75, 4.5, 6.0), "KEE_fit_mass_mix_B3f", "W_mix_B3f")
    h_KEE_mass_mix_B3r = config_bin_q.Histo1D( ("KEE_mass_mix_B3r", "", 75, 4.5, 6.0), "KEE_mass_mix_B3r", "W_mix_B3r")
    h_KEE_fitMass_mix_B3r = config_bin_q.Histo1D( ("KEE_fitMass_mix_B3r", "", 75, 4.5, 6.0), "KEE_fit_mass_mix_B3r", "W_mix_B3r")
    h_KEE_mass_Low_B3f = config_bin_q.Histo1D( ("KEE_mass_Low_B3f", "", 75, 4.5, 6.0), "KEE_mass_Low_B3f", "W_Low_B3f")
    h_KEE_fitMass_Low_B3f = config_bin_q.Histo1D( ("KEE_fitMass_Low_B3f", "", 75, 4.5, 6.0), "KEE_fit_mass_Low_B3f", "W_Low_B3f")
    h_KEE_mass_Low_B3r = config_bin_q.Histo1D( ("KEE_mass_Low_B3r", "", 75, 4.5, 6.0), "KEE_mass_Low_B3r", "W_Low_B3r")
    h_KEE_fitMass_Low_B3r = config_bin_q.Histo1D( ("KEE_fitMass_Low_B3r", "", 75, 4.5, 6.0), "KEE_fit_mass_Low_B3r", "W_Low_B3r")
    # KMuMu histograms
    h_KMM_mass_PF_B3f = config_bin_q.Histo1D( ("KMM_mass_PF_B3f", "", 75, 4.5, 6.0), "KMM_mass_PF_B3f", "W_PF_B3f_mu")
    h_KMM_fitMass_PF_B3f = config_bin_q.Histo1D( ("KMM_fitMass_PF_B3f", "", 75, 4.5, 6.0), "KMM_fit_mass_PF_B3f", "W_PF_B3f_mu")
    h_KMM_mass_PF_B3r = config_bin_q.Histo1D( ("KMM_mass_PF_B3r", "", 75, 4.5, 6.0), "KMM_mass_PF_B3r", "W_PF_B3r_mu")
    h_KMM_fitMass_PF_B3r = config_bin_q.Histo1D( ("KMM_fitMass_PF_B3r", "", 75, 4.5, 6.0), "KMM_fit_mass_PF_B3r", "W_PF_B3r_mu")

    outHistFile = ROOT.TFile.Open(outFileName,"RECREATE")
    outHistFile.cd()

    h_KEE_mass_B3f.Write()
    h_KEE_fitMass_B3f.Write()
    h_KEE_mass_B3r.Write()
    h_KEE_fitMass_B3r.Write()
    h_KEE_mass_PF_B3f.Write()
    h_KEE_fitMass_PF_B3f.Write()
    h_KEE_mass_PF_B3r.Write()
    h_KEE_fitMass_PF_B3r.Write()
    h_KEE_mass_mix_B3f.Write()
    h_KEE_fitMass_mix_B3f.Write()
    h_KEE_mass_mix_B3r.Write()
    h_KEE_fitMass_mix_B3r.Write()
    h_KEE_mass_Low_B3f.Write()
    h_KEE_fitMass_Low_B3f.Write()
    h_KEE_mass_Low_B3r.Write()
    h_KEE_fitMass_Low_B3r.Write()
    h_KMM_mass_PF_B3f.Write()
    h_KMM_fitMass_PF_B3f.Write()
    h_KMM_mass_PF_B3r.Write()
    h_KMM_fitMass_PF_B3r.Write()

    outHistFile.Close()


else:


    # Skimmed quantities
    tree_q = ind_cuts.Define("e1_isPF_sk","Take(e1_isPF, Idx_KEE)") \
                     .Define("e2_isPF_sk","Take(e2_isPF, Idx_KEE)") \
                     .Define("e1_isPFov_sk","Take(e1_isPFov, Idx_KEE)") \
                     .Define("e2_isPFov_sk","Take(e2_isPFov, Idx_KEE)") \
                     .Define("e1_isLow_sk","Take(e1_isLow, Idx_KEE)") \
                     .Define("e2_isLow_sk","Take(e2_isLow, Idx_KEE)") \
                     .Define("KEE_mll_fullfit_sk","Take(BToKEE_mll_fullfit, Idx_KEE)") \
                     .Define("KEE_mll_raw_sk","Take(BToKEE_mll_raw, Idx_KEE)") \
                     .Define("KEE_mass_sk","Take(BToKEE_mass, Idx_KEE)") \
                     .Define("KEE_fit_mass_sk","Take(BToKEE_fit_mass, Idx_KEE)") \
                     .Define("mu1_isPF_sk","Take(mu1_isPF, Idx_KMM)") \
                     .Define("mu2_isPF_sk","Take(mu2_isPF, Idx_KMM)") \
                     .Define("KMM_mll_fullfit_sk","Take(BToKMuMu_mll_fullfit, Idx_KMM)") \
                     .Define("KMM_mll_raw_sk","Take(BToKMuMu_mll_raw, Idx_KMM)") \
                     .Define("KMM_mass_sk","Take(BToKMuMu_mass, Idx_KMM)") \
                     .Define("KMM_fit_mass_sk","Take(BToKMuMu_fit_mass, Idx_KMM)")
    
    Mod.print_time('Before_snapshot')
    
    # Save skimmed branches in skimTree
    tree_q.Snapshot("skimTree", outFileName, "\\b([^ ]*)(_sk)|^(Idx_)([^ ]*)")


Mod.print_time('End')
