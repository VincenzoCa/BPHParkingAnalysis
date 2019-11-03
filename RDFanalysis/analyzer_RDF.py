#source settings.sh
#python analyzer_RDF.py --isEleCh (0,1) --inList input_list.txt --JOBid (1,2,...) --outFile /path/output.root

import ROOT
import module as Mod
# Enable multi-threading
ROOT.ROOT.EnableImplicitMT()
import argparse

Mod.print_time('Start')

parser = argparse.ArgumentParser()
parser.add_argument("--isEleCh", "-c", help="set channel")
parser.add_argument("--inList", "-l", help="set list")
parser.add_argument("--JOBid", "-j", help="set job ID")
parser.add_argument("--outFile", "-f", help="set output folder")
args = parser.parse_args()

isEle = args.isEleCh
listName = args.inList
outFileName = args.outFile

print('\n --isEleCh ' + isEle + ' --list ' + listName + ' --JOBid ' + args.JOBid + ' --outFile ' + outFileName)
    
# Create dataframe from NanoAOD files 
files = ROOT.std.vector("string")()
mylist = [line.rstrip('\n') for line in open(listName)]
print('\n' + ' Input files:\n'+'\n'.join(map(str, mylist))+'\n')
for f in mylist : files.push_back(f)
df = ROOT.ROOT.RDataFrame("Events", files)

entries = df.Count()
print("\n %s entries " %entries.GetValue())


ROOT.gInterpreter.Declare(Mod.Cuts_code)
ROOT.gInterpreter.Declare(Mod.Bin_code)
ROOT.gInterpreter.Declare(Mod.noPFover_code)
ROOT.gInterpreter.Declare(Mod.bothPF_code)
ROOT.gInterpreter.Declare(Mod.bothLow_code)
ROOT.gInterpreter.Declare(Mod.mix_code)
ROOT.gInterpreter.Declare(Mod.weights_code)


#Electron or Muon channel
isKEE_ch = 'true' #if isEle else 'false'
nTriplet_ch = 'nBToKEE' #if isEle else 'nBToKMuMu'
B_pT_ch = 'BToKEE_pt' #if isEle else 'BToKMuMu_pt'
cos2D_ch = 'BToKEE_cos2D' #if isEle else 'BToKMuMu_cos2D'
B_svprob_ch = 'BToKEE_svprob' #if isEle else 'BToKMuMu_svprob'
B_Lxy_ch = 'BToKEE_l_xy' #if isEle else 'BToKMuMu_l_xy'
B_Lxy_unc_ch = 'BToKEE_l_xy_unc' #if isEle else 'BToKMuMu_l_xy_unc'
K_DCASig_ch = 'Take(ProbeTracks_DCASig, BToKEE_kIdx)' #if isEle else 'Take(ProbeTracks_DCASig, BToKMuMu_kIdx)'
K_pT_ch = 'Take(ProbeTracks_pt, BToKEE_kIdx)' #if isEle else 'Take(ProbeTracks_pt, BToKMuMu_kIdx)'
l1_pT_ch = 'Take(Electron_pt, BToKEE_l1Idx)' #if isEle else 'Take(Muon_pt, BToKMuMu_l1Idx)'
l2_pT_ch = 'Take(Electron_pt, BToKEE_l2Idx)' #if isEle else 'Take(Muon_pt, BToKMuMu_l2Idx)'
l1_Veto_ch = 'Take(Electron_convVeto, BToKEE_l1Idx)'
l2_Veto_ch = 'Take(Electron_convVeto, BToKEE_l2Idx)'
l1_isPFov_ch = 'Take(Electron_isPFoverlap, BToKEE_l1Idx)'
l2_isPFov_ch = 'Take(Electron_isPFoverlap, BToKEE_l2Idx)'
l1_isPF_ch = 'Take(Electron_isPF, BToKEE_l1Idx)'
l2_isPF_ch = 'Take(Electron_isPF, BToKEE_l2Idx)'
l1_isLow_ch = 'Take(Electron_isLowPt, BToKEE_l1Idx)'
l2_isLow_ch = 'Take(Electron_isLowPt, BToKEE_l2Idx)'
l1_mvaId_ch = 'Take(Electron_mvaId, BToKEE_l1Idx)' 
l2_mvaId_ch = 'Take(Electron_mvaId, BToKEE_l2Idx)' 
mll_fit_ch = 'BToKEE_mll_fullfit' #if isEle else 'BToKMuMu_mll_fullfit'
mll_raw_ch = 'BToKEE_mll_raw' #if isEle else 'BToKMuMu_mll_raw'
B_mass_ch = 'BToKEE_mass' #if isEle else 'BToKMuMu_mass'
B_fit_mass_ch = 'BToKEE_fit_mass' #if isEle else 'BToKMuMu_fit_mass'


#Define useful quantities
branch_def = df.Define("isKEE",isKEE_ch) \
               .Define("nTriplet",nTriplet_ch) \
               .Define("B_pT", B_pT_ch) \
               .Define("cos2D",cos2D_ch) \
               .Define("B_svprob", B_svprob_ch) \
               .Define("B_Lxy", B_Lxy_ch) \
               .Define("B_Lxy_unc", B_Lxy_unc_ch) \
               .Define("K_DCASig", K_DCASig_ch) \
               .Define("K_pT", K_pT_ch) \
               .Define("l1_pT",l1_pT_ch) \
               .Define("l2_pT",l2_pT_ch) \
               .Define("l1_Veto",l1_Veto_ch) \
               .Define("l2_Veto",l2_Veto_ch) \
               .Define("l1_isPFov",l1_isPFov_ch) \
               .Define("l2_isPFov",l2_isPFov_ch) \
               .Define("l1_isPF",l1_isPF_ch) \
               .Define("l2_isPF",l2_isPF_ch) \
               .Define("l1_isLow",l1_isLow_ch) \
               .Define("l2_isLow",l2_isLow_ch) \
               .Define("l1_mvaId",l1_mvaId_ch) \
               .Define("l2_mvaId",l2_mvaId_ch) \
               .Define("mll_fit",mll_fit_ch) \
               .Define("mll_raw",mll_raw_ch) \
               .Define("B_mass",B_mass_ch) \
               .Define("B_fit_mass",B_fit_mass_ch) \
               .Define("W","weights(nTriplet)")


#Finding indices of triplets passing cuts
ind_cuts = branch_def.Define("Idx", "Cuts( nTriplet, B_pT, cos2D, B_svprob, B_Lxy, B_Lxy_unc, K_DCASig, K_pT, l1_pT, l2_pT, l1_mvaId, l2_mvaId, l1_Veto, l2_Veto )")


#Finding indices of triplets for different configurations
config = ind_cuts.Define("Idx_noPFov", "noPFover( Idx, l1_isPFov, l2_isPFov )" ) \
                 .Define("Idx_PF", "bothPF( Idx, l1_isPF, l2_isPF )" ) \
                 .Define("Idx_mix", "mix( Idx_noPFov, l1_isPF, l2_isPF, l1_isLow, l2_isLow )" ) \
                 .Define("Idx_Low", "bothLow( Idx_noPFov, l1_isLow, l2_isLow )" )


#Bin ranges
#B0 = [0, 1], B1 = [1, 2.5], B2 = [2.5, 2.9], B3 = [2.9, 3.3], B4 = [3.3, 3.58], B5 = [3.58, 100]
#Finding indices of triplets with 2.9 < mll_fullfit(_raw) < 3.1
config_bin = config.Define("noPFov_B3f", "Bin( Idx_noPFov, 2.9, 3.3, mll_fit )" ) \
                   .Define("noPFov_B3r", "Bin( Idx_noPFov, 2.9, 3.3, mll_raw )" ) \
                   .Define("PF_B3f", "Bin( Idx_PF, 2.9, 3.3, mll_fit )" ) \
                   .Define("PF_B3r", "Bin( Idx_PF, 2.9, 3.3, mll_raw )" ) \
                   .Define("mix_B3f", "Bin( Idx_mix, 2.9, 3.3, mll_fit )" ) \
                   .Define("mix_B3r", "Bin( Idx_mix, 2.9, 3.3, mll_raw )" ) \
                   .Define("Low_B3f", "Bin( Idx_Low, 2.9, 3.3, mll_fit )" ) \
                   .Define("Low_B3r", "Bin( Idx_Low, 2.9, 3.3, mll_raw )" )


#Defining quantities
config_bin_q = config_bin.Define("B_mass_B3f","Take(B_mass, noPFov_B3f)") \
                         .Define("B_fit_mass_B3f","Take(B_fit_mass, noPFov_B3f)") \
                         .Define("W_B3f","Take(W, noPFov_B3f)") \
                         .Define("B_mass_B3r","Take(B_mass, noPFov_B3r)") \
                         .Define("B_fit_mass_B3r","Take(B_fit_mass, noPFov_B3r)") \
                         .Define("W_B3r","Take(W, noPFov_B3r)") \
                         .Define("B_mass_PF_B3f","Take(B_mass, PF_B3f)") \
                         .Define("B_fit_mass_PF_B3f","Take(B_fit_mass, PF_B3f)") \
                         .Define("W_PF_B3f","Take(W, PF_B3f)") \
                         .Define("B_mass_PF_B3r","Take(B_mass, PF_B3r)") \
                         .Define("B_fit_mass_PF_B3r","Take(B_fit_mass, PF_B3r)") \
                         .Define("W_PF_B3r","Take(W, PF_B3r)") \
                         .Define("B_mass_mix_B3f","Take(B_mass, mix_B3f)") \
                         .Define("B_fit_mass_mix_B3f","Take(B_fit_mass, mix_B3f)") \
                         .Define("W_mix_B3f","Take(W, mix_B3f)") \
                         .Define("B_mass_mix_B3r","Take(B_mass, mix_B3r)") \
                         .Define("B_fit_mass_mix_B3r","Take(B_fit_mass, mix_B3r)") \
                         .Define("W_mix_B3r","Take(W, mix_B3r)") \
                         .Define("B_mass_Low_B3f","Take(B_mass, Low_B3f)") \
                         .Define("B_fit_mass_Low_B3f","Take(B_fit_mass, Low_B3f)") \
                         .Define("W_Low_B3f","Take(W, Low_B3f)") \
                         .Define("B_mass_Low_B3r","Take(B_mass, Low_B3r)") \
                         .Define("B_fit_mass_Low_B3r","Take(B_fit_mass, Low_B3r)") \
                         .Define("W_Low_B3r","Take(W, Low_B3r)") 


#'''
# SAVE HISTOGRAMS
Mod.print_time('Before histo')

h_Bmass_B3f = config_bin_q.Histo1D( ("Bmass_B3f", "Bmass_B3f", 75, 4.5, 6.0), "B_mass_B3f", "W_B3f")
h_BfitMass_B3f = config_bin_q.Histo1D( ("BfitMass_B3f", "BfitMass_B3f", 75, 4.5, 6.0), "B_fit_mass_B3f", "W_B3f")
h_Bmass_B3r = config_bin_q.Histo1D( ("Bmass_B3r", "Bmass_B3r", 75, 4.5, 6.0), "B_mass_B3r", "W_B3r")
h_BfitMass_B3r = config_bin_q.Histo1D( ("BfitMass_B3r", "BfitMass_B3r", 75, 4.5, 6.0), "B_fit_mass_B3r", "W_B3r")
h_Bmass_PF_B3f = config_bin_q.Histo1D( ("Bmass_PF_B3f", "Bmass_PF_B3f", 75, 4.5, 6.0), "B_mass_PF_B3f", "W_PF_B3f")
h_BfitMass_PF_B3f = config_bin_q.Histo1D( ("BfitMass_PF_B3f", "BfitMass_PF_B3f", 75, 4.5, 6.0), "B_fit_mass_PF_B3f", "W_PF_B3f")
h_Bmass_PF_B3r = config_bin_q.Histo1D( ("Bmass_PF_B3r", "Bmass_PF_B3r", 75, 4.5, 6.0), "B_mass_PF_B3r", "W_PF_B3r")
h_BfitMass_PF_B3r = config_bin_q.Histo1D( ("BfitMass_PF_B3r", "BfitMass_PF_B3r", 75, 4.5, 6.0), "B_fit_mass_PF_B3r", "W_PF_B3r")
h_Bmass_mix_B3f = config_bin_q.Histo1D( ("Bmass_mix_B3f", "Bmass_mix_B3f", 75, 4.5, 6.0), "B_mass_mix_B3f", "W_mix_B3f")
h_BfitMass_mix_B3f = config_bin_q.Histo1D( ("BfitMass_mix_B3f", "BfitMass_mix_B3f", 75, 4.5, 6.0), "B_fit_mass_mix_B3f", "W_mix_B3f")
h_Bmass_mix_B3r = config_bin_q.Histo1D( ("Bmass_mix_B3r", "Bmass_mix_B3r", 75, 4.5, 6.0), "B_mass_mix_B3r", "W_mix_B3r")
h_BfitMass_mix_B3r = config_bin_q.Histo1D( ("BfitMass_mix_B3r", "BfitMass_mix_B3r", 75, 4.5, 6.0), "B_fit_mass_mix_B3r", "W_mix_B3r")
h_Bmass_Low_B3f = config_bin_q.Histo1D( ("Bmass_Low_B3f", "Bmass_Low_B3f", 75, 4.5, 6.0), "B_mass_Low_B3f", "W_Low_B3f")
h_BfitMass_Low_B3f = config_bin_q.Histo1D( ("BfitMass_Low_B3f", "BfitMass_Low_B3f", 75, 4.5, 6.0), "B_fit_mass_Low_B3f", "W_Low_B3f")
h_Bmass_Low_B3r = config_bin_q.Histo1D( ("Bmass_Low_B3r", "Bmass_Low_B3r", 75, 4.5, 6.0), "B_mass_Low_B3r", "W_Low_B3r")
h_BfitMass_Low_B3r = config_bin_q.Histo1D( ("BfitMass_Low_B3r", "BfitMass_Low_B3r", 75, 4.5, 6.0), "B_fit_mass_Low_B3r", "W_Low_B3r")

outHistFile = ROOT.TFile.Open(outFileName,"RECREATE")
outHistFile.cd()

h_Bmass_B3f.Write()
h_BfitMass_B3f.Write()
h_Bmass_B3r.Write()
h_BfitMass_B3r.Write()
h_Bmass_PF_B3f.Write()
h_BfitMass_PF_B3f.Write()
h_Bmass_PF_B3r.Write()
h_BfitMass_PF_B3r.Write()
h_Bmass_mix_B3f.Write()
h_BfitMass_mix_B3f.Write()
h_Bmass_mix_B3r.Write()
h_BfitMass_mix_B3r.Write()
h_Bmass_Low_B3f.Write()
h_BfitMass_Low_B3f.Write()
h_Bmass_Low_B3r.Write()
h_BfitMass_Low_B3r.Write()

outHistFile.Close()


'''
# SAVE THE NEW BRANCHES IN newtree
Mod.print_time('Before_snapshot')

brList = ROOT.vector('string')()
for brName in ["B_mass_B3f","B_fit_mass_B3f", "B_mass_B3r", "B_fit_mass_B3r", "B_mass_PF_B3f", "B_fit_mass_PF_B3f", "B_mass_PF_B3r", "B_fit_mass_PF_B3r", "B_mass_noPFovLow_B3f", "B_fit_mass_noPFovLow_B3f", "B_mass_noPFovLow_B3r", "B_fit_mass_noPFovLow_B3r"]:
    brList.push_back(brName)
config_bin_q.Snapshot("newtree", outFileName, brList)
'''

Mod.print_time('End')
