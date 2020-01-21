# source settings.sh
# python analyzer_RDF.py --outFile /path/output.root
#
#    # --inList input_list.txt | for input list
#    # --testFile testFile.root | to test a single input file
#    # --JOBid (1, 2, 3, ...) | for batch jobs
#
# MC: 
#
#    # Add --isMC
#    # Add --isEE | for electron channel. Muon channel by default
#    # Add --isResonant | for resonant channel
#
# By default, the output will be stored as histograms. Please add --tree if you want a tree as output
# The output consists of both KEE and KMuMu quantities for data

import ROOT
import module as Mod
# Enable multi-threading
ROOT.ROOT.EnableImplicitMT()
import argparse
import numpy as np

Mod.print_time('Start')

parser = argparse.ArgumentParser()
parser.add_argument("--inList", "-l", help="set list")
parser.add_argument("--testFile", "-tf", help="set test file")
parser.add_argument("--JOBid", "-j", help="set job ID")
parser.add_argument("--outFile", "-f", help="set output folder")
parser.add_argument('--isMC', action='store_true')
parser.add_argument('--isEE', action='store_true')
parser.add_argument('--isResonant', action='store_true')
parser.add_argument('--tree', action='store_true')
args = parser.parse_args()

listName = args.inList
tFile = args.testFile
outFileName = args.outFile
isMC = args.isMC
isEE = args.isEE
isResonant = args.isResonant
isTree = args.tree
print('\n --inList {} --testFile {} --JOBid {} --outFile {} --isEE {} --isMC {} --isResonant {} --tree {}'.format(listName, tFile, args.JOBid, outFileName, str(isEE), str(isMC), str(isResonant), str(isTree)))
    
# Create dataframe from NanoAOD files 
files = ROOT.std.vector("string")() if listName else ROOT.std.vector("string")(1)
if listName:
    #files = ROOT.std.vector("string")()
    mylist = [line.rstrip('\n') for line in open(listName)]
    print('\n Input files:\n'+'\n'.join(map(str, mylist))+'\n')
    for f in mylist : files.push_back(f)
else:
    files.push_back(tFile)
df = ROOT.ROOT.RDataFrame("Events", files)

entries = df.Count()
print('\n {} entries '.format(entries.GetValue()))


ROOT.gInterpreter.Declare(Mod.Cuts_code)
ROOT.gInterpreter.Declare(Mod.EleCuts_code)
ROOT.gInterpreter.Declare(Mod.MuCuts_code)
ROOT.gInterpreter.Declare(Mod.Bin_code)
ROOT.gInterpreter.Declare(Mod.Bin_MC_code)
ROOT.gInterpreter.Declare(Mod.noPFover_code)
ROOT.gInterpreter.Declare(Mod.bothX_code)
ROOT.gInterpreter.Declare(Mod.mix_code)
ROOT.gInterpreter.Declare(Mod.flagGenMatchExt_code)
ROOT.gInterpreter.Declare(Mod.computedR_code)
ROOT.gInterpreter.Declare(Mod.bestRank_code)

isRes_def =  'true' if isResonant else 'false'

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
               .Define("W_mu", "RVec<float>(nBToKMuMu, 0.999)") \
               .Define("isResonant", isRes_def) \
               .Define("ele_JPsi_1", "2.969") \
               .Define("ele_JPsi_2", "3.227") \
               .Define("mu_JPsi_1", "3.0189") \
               .Define("mu_JPsi_2", "3.1743")                  

# https://amartell.web.cern.ch/amartell/Analysis/BParking/NANOstudies_2019/Bmass_MC_JPsi_ee_Binclusive_unbinned/ee_All.png
# Ele JPsi bin: 2.969, 3.227

# https://amartell.web.cern.ch/amartell/Analysis/BParking/NANOstudies_2019/Bmass_MC_JPsi_mm_Binclusive_unbinned/mm_All.png
# Mu JPsi bin: 3.0189, 3.1743

# Finding indices of triplets passing cuts for KEE and KMuMu
ind_cuts = branch_def.Define("tmp_Idx_KEE", "Cuts( nBToKEE, BToKEE_pt, BToKEE_cos2D, BToKEE_svprob, LxySig, K_DCASig, K_pT )") \
                     .Define("Idx_KEE", "EleCuts( tmp_Idx_KEE, e1_pT, e2_pT, e1_mvaId, e2_mvaId, e1_Veto, e2_Veto )") \
                     .Define("tmp_Idx_KMM", "Cuts( nBToKMuMu, BToKMuMu_pt, BToKMuMu_cos2D, BToKMuMu_svprob, LxySig_mu, K_DCASig_mu, K_pT_mu )") \
                     .Define("Idx_KMM", "MuCuts( tmp_Idx_KMM, mu1_pT, mu2_pT, nAddTrMu )")
    
#config = ind_cuts

if np.logical_not(isMC):

    # Finding indices of triplets for different configurations (All, PF, mix, Low)
    config = ind_cuts.Define("Idx_noPFov", "noPFover( Idx_KEE, e1_isPFov, e2_isPFov )" ) \
                     .Define("Idx_PF", "bothX( Idx_KEE, e1_isPF, e2_isPF )" ) \
                     .Define("Idx_mix", "mix( Idx_noPFov, e1_isPF, e2_isPF, e1_isLow, e2_isLow )" ) \
                     .Define("Idx_Low", "bothX( Idx_noPFov, e1_isLow, e2_isLow )" ) \
                     .Define("Idx_PF_mu", "bothX( Idx_KMM, mu1_isPF, mu2_isPF )" )
else:
            
    # 443 = JPsi    521 = B+  321 = K
    GenPart_l1_idx = "Take(Electron_genPartIdx, BToKEE_l1Idx)" if isEE else "Take(Muon_genPartIdx, BToKMuMu_l1Idx)"
    GenPart_l2_idx = "Take(Electron_genPartIdx, BToKEE_l2Idx)" if isEE else "Take(Muon_genPartIdx, BToKMuMu_l2Idx)"
    GenPart_k_idx = "Take(ProbeTracks_genPartIdx, BToKEE_kIdx)" if isEE else "Take(ProbeTracks_genPartIdx, BToKMuMu_kIdx)"
        
    l1_eta = "Take(Electron_eta, BToKEE_l1Idx)" if isEE else "Take(Muon_eta, BToKMuMu_l1Idx)"
    l2_eta = "Take(Electron_eta, BToKEE_l2Idx)" if isEE else "Take(Muon_eta, BToKMuMu_l2Idx)"
    k_eta = "Take(ProbeTracks_eta, BToKEE_kIdx)" if isEE else "Take(ProbeTracks_eta, BToKMuMu_kIdx)"
    l1_phi = "Take(Electron_phi, BToKEE_l1Idx)" if isEE else "Take(Muon_phi, BToKMuMu_l1Idx)"
    l2_phi = "Take(Electron_phi, BToKEE_l2Idx)" if isEE else "Take(Muon_phi, BToKMuMu_l2Idx)"
    k_phi = "Take(ProbeTracks_phi, BToKEE_kIdx)" if isEE else "Take(ProbeTracks_phi, BToKMuMu_kIdx)"     
            
    n_Btriplet = "nBToKEE" if isEE else "nBToKMuMu"
    Idx_tmp = "noPFover( Idx_KEE, e1_isPFov, e2_isPFov )" if isEE else "bothX( Idx_KMM, mu1_isPF, mu2_isPF )"
    Idx_Bin = "Bin_MC( Idx_KEE, ele_JPsi_1, ele_JPsi_2, BToKEE_mll_fullfit )" if isEE else "Bin( Idx_KMM, 2.9, 3.3, BToKMuMu_mll_fullfit )"
    Idx_noPFov = "Idx_best_dR" if isEE else "RVec<unsigned int>(n_Btriplet, 0)"
    Idx_PF = "bothX( Idx_noPFov, e1_isPF, e2_isPF )" if isEE else "RVec<unsigned int>(n_Btriplet, 0)"
    Idx_mix = "mix( Idx_noPFov, e1_isPF, e2_isPF, e1_isLow, e2_isLow )" if isEE else "RVec<unsigned int>(n_Btriplet, 0)"
    Idx_Low = "bothX( Idx_noPFov, e1_isLow, e2_isLow )" if isEE else "RVec<unsigned int>(n_Btriplet, 0)"
    Idx_PF_mu = "RVec<unsigned int>(n_Btriplet, 0)" if isEE else "Idx_best_dR"
        
    config  = ind_cuts.Define("GenPart_l1_idx", GenPart_l1_idx) \
                      .Define("GenPart_l2_idx", GenPart_l2_idx) \
                      .Define("GenPart_k_idx", GenPart_k_idx) \
                      .Define("GenPart_l1_pdgId", "Take(GenPart_pdgId, GenPart_l1_idx)") \
                      .Define("GenPart_l2_pdgId", "Take(GenPart_pdgId, GenPart_l2_idx)") \
                      .Define("GenPart_k_pdgId", "Take(GenPart_pdgId, GenPart_k_idx)") \
                      .Define("GenMothPart_l1_idx", "Take(GenPart_genPartIdxMother, GenPart_l1_idx)") \
                      .Define("GenMothPart_l2_idx", "Take(GenPart_genPartIdxMother, GenPart_l2_idx)") \
                      .Define("GenMothPart_k_idx", "Take(GenPart_genPartIdxMother, GenPart_k_idx)") \
                      .Define("GenMothPart_l1_pdgId", "Take(GenPart_pdgId, GenMothPart_l1_idx)") \
                      .Define("GenMothPart_l2_pdgId", "Take(GenPart_pdgId, GenMothPart_l2_idx)") \
                      .Define("GenMothPart_k_pdgId", "Take(GenPart_pdgId, GenMothPart_k_idx)") \
                      .Define("GenGMothPart_l1_idx", "Take(GenPart_genPartIdxMother, GenMothPart_l1_idx)") \
                      .Define("GenGMothPart_l2_idx", "Take(GenPart_genPartIdxMother, GenMothPart_l2_idx)") \
                      .Define("GenGMothPart_k_idx", "Take(GenPart_genPartIdxMother, GenMothPart_k_idx)") \
                      .Define("GenGMothPart_l1_pdgId", "Take(GenPart_pdgId, GenGMothPart_l1_idx)") \
                      .Define("GenGMothPart_l2_pdgId", "Take(GenPart_pdgId, GenGMothPart_l2_idx)") \
                      .Define("GenGMothPart_k_pdgId", "Take(GenPart_pdgId, GenGMothPart_k_idx)") \
                      .Define("GenPart_l1_eta", "Take(GenPart_eta, GenPart_l1_idx)") \
                      .Define("GenPart_l2_eta", "Take(GenPart_eta, GenPart_l2_idx)") \
                      .Define("GenPart_k_eta", "Take(GenPart_eta, GenPart_k_idx)") \
                      .Define("GenPart_l1_phi", "Take(GenPart_phi, GenPart_l1_idx)") \
                      .Define("GenPart_l2_phi", "Take(GenPart_phi, GenPart_l2_idx)") \
                      .Define("GenPart_k_phi", "Take(GenPart_phi, GenPart_k_idx)") \
                      .Define("l1_eta", l1_eta) \
                      .Define("l2_eta", l2_eta) \
                      .Define("k_eta", k_eta) \
                      .Define("l1_phi", l1_phi) \
                      .Define("l2_phi", l2_phi) \
                      .Define("k_phi", k_phi) \
                      .Define("isGenMatched", "flagGenMatchExt(isResonant, GenPart_l1_pdgId, GenPart_l2_pdgId, GenPart_k_pdgId, GenMothPart_l1_pdgId, GenMothPart_l2_pdgId, GenMothPart_k_pdgId, GenGMothPart_l1_pdgId, GenGMothPart_l2_pdgId, GenGMothPart_k_pdgId)" ) \
                      .Define("dRwithGen", "computedR(isGenMatched, GenPart_l1_eta, GenPart_l2_eta, GenPart_k_eta, l1_eta, l2_eta, k_eta, GenPart_l1_phi, GenPart_l2_phi, GenPart_k_phi, l1_phi, l2_phi, k_phi)" ) \
                      .Define("n_Btriplet", n_Btriplet) \
                      .Define("Idx_tmp", Idx_tmp) \
                      .Define("Idx_Bin", Idx_Bin) \
                      .Define("Idx_best_dR", "bestRank( n_Btriplet, Idx_tmp, Idx_Bin, isGenMatched, dRwithGen )" ) \
                      .Define("Idx_noPFov", Idx_noPFov) \
                      .Define("Idx_PF", Idx_PF) \
                      .Define("Idx_mix", Idx_mix) \
                      .Define("Idx_Low", Idx_Low) \
                      .Define("Idx_PF_mu", Idx_PF_mu)

# Bin ranges
# B0 = [0, 1], B1 = [1, 2.5], B2 = [2.5, 2.9], B3 = [JPsi_1, JPsi_2], B4 = [3.3, 3.58], B5 = [3.58, 100]
# Finding indices of triplets with JPsi_1 < mll_fullfit < JPsi_2
config_bin = config.Define("noPFov_B3f", "Bin( Idx_noPFov, ele_JPsi_1, ele_JPsi_2, BToKEE_mll_fullfit )" ) \
                   .Define("PF_B3f", "Bin( Idx_PF, ele_JPsi_1, ele_JPsi_2, BToKEE_mll_fullfit )" ) \
                   .Define("mix_B3f", "Bin( Idx_mix, ele_JPsi_1, ele_JPsi_2, BToKEE_mll_fullfit )" ) \
                   .Define("Low_B3f", "Bin( Idx_Low, ele_JPsi_1, ele_JPsi_2, BToKEE_mll_fullfit )" ) \
                   .Define("PF_B3f_mu", "Bin( Idx_PF_mu, mu_JPsi_1, mu_JPsi_2, BToKMuMu_mll_fullfit )" )

# Defining final quantities
if np.logical_not(isMC):

    config_bin_q = config_bin.Define("KEE_fit_mass_B3f_T","Take(BToKEE_fit_mass, noPFov_B3f)") \
                             .Define("W_B3f","Take(W, noPFov_B3f)") \
                             .Define("KEE_fit_mass_PF_B3f_T","Take(BToKEE_fit_mass, PF_B3f)") \
                             .Define("W_PF_B3f","Take(W, PF_B3f)") \
                             .Define("KEE_fit_mass_mix_B3f_T","Take(BToKEE_fit_mass, mix_B3f)") \
                             .Define("W_mix_B3f","Take(W, mix_B3f)") \
                             .Define("KEE_fit_mass_Low_B3f_T","Take(BToKEE_fit_mass, Low_B3f)") \
                             .Define("W_Low_B3f","Take(W, Low_B3f)") \
                             .Define("KMM_fit_mass_PF_B3f_T","Take(BToKMuMu_fit_mass, PF_B3f_mu)") \
                             .Define("W_PF_B3f_mu","Take(W_mu, PF_B3f_mu)") 
        
else:

    config_bin_q = config_bin.Define("KEE_fit_mass_B3f_T","Take(BToKEE_fit_mass, noPFov_B3f)") \
                             .Define("e1_pT_B3f_T","Take(e1_pT, noPFov_B3f)") \
                             .Define("e2_pT_B3f_T","Take(e2_pT, noPFov_B3f)") \
                             .Define("W_B3f","Take(W, noPFov_B3f)") \
                             .Define("KEE_fit_mass_PF_B3f_T","Take(BToKEE_fit_mass, PF_B3f)") \
                             .Define("e1_pT_PF_B3f_T","Take(e1_pT, PF_B3f)") \
                             .Define("e2_pT_PF_B3f_T","Take(e2_pT, PF_B3f)") \
                             .Define("W_PF_B3f","Take(W, PF_B3f)") \
                             .Define("KEE_fit_mass_mix_B3f_T","Take(BToKEE_fit_mass, mix_B3f)") \
                             .Define("e1_pT_mix_B3f_T","Take(e1_pT, mix_B3f)") \
                             .Define("e2_pT_mix_B3f_T","Take(e2_pT, mix_B3f)") \
                             .Define("W_mix_B3f","Take(W, mix_B3f)") \
                             .Define("KEE_fit_mass_Low_B3f_T","Take(BToKEE_fit_mass, Low_B3f)") \
                             .Define("e1_pT_Low_B3f_T","Take(e1_pT, Low_B3f)") \
                             .Define("e2_pT_Low_B3f_T","Take(e2_pT, Low_B3f)") \
                             .Define("W_Low_B3f","Take(W, Low_B3f)") \
                             .Define("KMM_fit_mass_PF_B3f_T","Take(BToKMuMu_fit_mass, PF_B3f_mu)") \
                             .Define("mu1_pT_PF_B3f_T","Take(mu1_pT, PF_B3f_mu)") \
                             .Define("mu2_pT_PF_B3f_T","Take(mu2_pT, PF_B3f_mu)") \
                             .Define("W_PF_B3f_mu","Take(W_mu, PF_B3f_mu)")
                             
if  np.logical_not(isTree):

    # Save histograms
    Mod.print_time('Before histo')
    # KEE histograms
    h_KEE_fitMass_B3f = config_bin_q.Histo1D( ("KEE_fitMass_B3f", "", 75, 4.5, 6.0), "KEE_fit_mass_B3f_T", "W_B3f")
    h_KEE_fitMass_PF_B3f = config_bin_q.Histo1D( ("KEE_fitMass_PF_B3f", "", 75, 4.5, 6.0), "KEE_fit_mass_PF_B3f_T", "W_PF_B3f")
    h_KEE_fitMass_mix_B3f = config_bin_q.Histo1D( ("KEE_fitMass_mix_B3f", "", 75, 4.5, 6.0), "KEE_fit_mass_mix_B3f_T", "W_mix_B3f")
    h_KEE_fitMass_Low_B3f = config_bin_q.Histo1D( ("KEE_fitMass_Low_B3f", "", 75, 4.5, 6.0), "KEE_fit_mass_Low_B3f_T", "W_Low_B3f")
    # KMuMu histograms
    h_KMM_fitMass_PF_B3f = config_bin_q.Histo1D( ("KMM_fitMass_PF_B3f", "", 75, 4.5, 6.0), "KMM_fit_mass_PF_B3f_T", "W_PF_B3f_mu")
    if isMC:
        h_e1_pT_B3f = config_bin_q.Histo1D( ("e1_pT_B3f", "", 200, 0., 20.), "e1_pT_B3f_T", "W_B3f")
        h_e2_pT_B3f = config_bin_q.Histo1D( ("e2_pT_B3f", "", 200, 0., 20.), "e2_pT_B3f_T", "W_B3f")
        h_e1_pT_PF_B3f = config_bin_q.Histo1D( ("e1_pT_PF_B3f", "", 200, 0., 20.), "e1_pT_PF_B3f_T", "W_PF_B3f")
        h_e2_pT_PF_B3f = config_bin_q.Histo1D( ("e2_pT_PF_B3f", "", 200, 0., 20.), "e2_pT_PF_B3f_T", "W_PF_B3f")
        h_e1_pT_mix_B3f = config_bin_q.Histo1D( ("e1_pT_mix_B3f", "", 200, 0., 20.), "e1_pT_mix_B3f_T", "W_mix_B3f")
        h_e2_pT_mix_B3f = config_bin_q.Histo1D( ("e2_pT_mix_B3f", "", 200, 0., 20.), "e2_pT_mix_B3f_T", "W_mix_B3f")        
        h_e1_pT_Low_B3f = config_bin_q.Histo1D( ("e1_pT_Low_B3f", "", 200, 0., 20.), "e1_pT_Low_B3f_T", "W_Low_B3f")
        h_e2_pT_Low_B3f = config_bin_q.Histo1D( ("e2_pT_Low_B3f", "", 200, 0., 20.), "e2_pT_Low_B3f_T", "W_Low_B3f")
        h_mu1_pT_PF_B3f = config_bin_q.Histo1D( ("mu1_pT_PF_B3f", "", 200, 0., 20.), "mu1_pT_PF_B3f_T", "W_PF_B3f_mu")
        h_mu2_pT_PF_B3f = config_bin_q.Histo1D( ("mu2_pT_PF_B3f", "", 200, 0., 20.), "mu2_pT_PF_B3f_T", "W_PF_B3f_mu")

    outHistFile = ROOT.TFile.Open(outFileName,"RECREATE")
    outHistFile.cd()

    h_KEE_fitMass_B3f.Write()
    h_KEE_fitMass_PF_B3f.Write()
    h_KEE_fitMass_mix_B3f.Write()
    h_KEE_fitMass_Low_B3f.Write()
    h_KMM_fitMass_PF_B3f.Write()
    if isMC:            
        h_e1_pT_B3f.Write()
        h_e2_pT_B3f.Write()
        h_e1_pT_PF_B3f.Write()
        h_e2_pT_PF_B3f.Write()            
        h_e1_pT_mix_B3f.Write()
        h_e2_pT_mix_B3f.Write()        
        h_e1_pT_Low_B3f.Write()
        h_e2_pT_Low_B3f.Write()
        h_mu1_pT_PF_B3f.Write()
        h_mu2_pT_PF_B3f.Write()
    
    outHistFile.Close()


else:
    
    Mod.print_time('Before_snapshot')
    
    # Save skimmed branches in skimTree
    config_bin_q.Snapshot("newTree", outFileName, "\\b([^ ]*)(_T)")


Mod.print_time('End') 
