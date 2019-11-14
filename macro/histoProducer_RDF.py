# A Python script to produce histograms from skimmed nanoAOD files
# The output consists of both KEE and KMuMu histograms
#
# source settings.sh
# python histoProducer_RDF.py --inList skimNano_list.txt --JOBid (1,2,...) --outFile /path/output.root

import ROOT
import Hmodule as Hmod
import sys
sys.path.insert(1, '../RDFanalysis')
import module as Mod
# Enable multi-threading
ROOT.ROOT.EnableImplicitMT()
import argparse

Mod.print_time('Start')

parser = argparse.ArgumentParser()
parser.add_argument("--inList", "-l", help="set list")
parser.add_argument("--JOBid", "-j", help="set job ID")
parser.add_argument("--outFile", "-f", help="set output folder")
args = parser.parse_args()

listName = args.inList
outFileName = args.outFile
print('\n --list {} --JOBid {} --outFile {}'.format(listName, args.JOBid, outFileName))
    
# Create dataframe from NanoAOD files 
files = ROOT.std.vector("string")()
mylist = [line.rstrip('\n') for line in open(listName)]
print('\n Input files:\n'+'\n'.join(map(str, mylist))+'\n')
for f in mylist : files.push_back(f)
df = ROOT.ROOT.RDataFrame("skimTree", files)

entries = df.Count()
print('\n {} entries '.format(entries.GetValue()))


ROOT.gInterpreter.Declare(Hmod.noPFover_code)
ROOT.gInterpreter.Declare(Hmod.bothPF_code)
ROOT.gInterpreter.Declare(Mod.bothX_code)
ROOT.gInterpreter.Declare(Mod.mix_code)
ROOT.gInterpreter.Declare(Mod.Bin_code)


# Finding indices of triplets for different configurations (All, PF, mix, Low)
config = df.Define("Idx_noPFov", "noPFover( Idx_KEE, e1_isPFov_sk, e2_isPFov_sk )" ) \
           .Define("Idx_PF", "bothPF( Idx_KEE, e1_isPF_sk, e2_isPF_sk )" ) \
           .Define("Idx_mix", "mix( Idx_noPFov, e1_isPF_sk, e2_isPF_sk, e1_isLow_sk, e2_isLow_sk )" ) \
           .Define("Idx_Low", "bothX( Idx_noPFov, e1_isLow_sk, e2_isLow_sk )" ) \
           .Define("Idx_PF_mu", "bothPF( Idx_KMM, mu1_isPF_sk, mu2_isPF_sk )" ) \
           .Define("W", "RVec<float>((int)Idx_KEE.size(), 0.999)" ) \
           .Define("W_mu", "RVec<float>((int)Idx_KMM.size(), 0.999)" )


# Bin ranges
# B0 = [0, 1], B1 = [1, 2.5], B2 = [2.5, 2.9], B3 = [2.9, 3.3], B4 = [3.3, 3.58], B5 = [3.58, 100]
# Finding indices of triplets with 2.9 < mll_fullfit(_raw) < 3.1
config_bin = config.Define("noPFov_B3f", "Bin( Idx_noPFov, 2.9, 3.3, KEE_mll_fullfit_sk )" ) \
                   .Define("noPFov_B3r", "Bin( Idx_noPFov, 2.9, 3.3, KEE_mll_raw_sk )" ) \
                   .Define("PF_B3f", "Bin( Idx_PF, 2.9, 3.3, KEE_mll_fullfit_sk )" ) \
                   .Define("PF_B3r", "Bin( Idx_PF, 2.9, 3.3, KEE_mll_raw_sk )" ) \
                   .Define("mix_B3f", "Bin( Idx_mix, 2.9, 3.3, KEE_mll_fullfit_sk )" ) \
                   .Define("mix_B3r", "Bin( Idx_mix, 2.9, 3.3, KEE_mll_raw_sk )" ) \
                   .Define("Low_B3f", "Bin( Idx_Low, 2.9, 3.3, KEE_mll_fullfit_sk )" ) \
                   .Define("Low_B3r", "Bin( Idx_Low, 2.9, 3.3, KEE_mll_raw_sk )" ) \
                   .Define("PF_B3f_mu", "Bin( Idx_PF_mu, 2.9, 3.3, KMM_mll_fullfit_sk )" ) \
                   .Define("PF_B3r_mu", "Bin( Idx_PF_mu, 2.9, 3.3, KMM_mll_raw_sk )" ) 


# Defining quantities for histograms
config_bin_q = config_bin.Define("KEE_mass_B3f","Take(KEE_mass_sk, noPFov_B3f)") \
                         .Define("KEE_fit_mass_B3f","Take(KEE_fit_mass_sk, noPFov_B3f)") \
                         .Define("W_B3f","Take(W, noPFov_B3f)") \
                         .Define("KEE_mass_B3r","Take(KEE_mass_sk, noPFov_B3r)") \
                         .Define("KEE_fit_mass_B3r","Take(KEE_fit_mass_sk, noPFov_B3r)") \
                         .Define("W_B3r","Take(W, noPFov_B3r)") \
                         .Define("KEE_mass_PF_B3f","Take(KEE_mass_sk, PF_B3f)") \
                         .Define("KEE_fit_mass_PF_B3f","Take(KEE_fit_mass_sk, PF_B3f)") \
                         .Define("W_PF_B3f","Take(W, PF_B3f)") \
                         .Define("KEE_mass_PF_B3r","Take(KEE_mass_sk, PF_B3r)") \
                         .Define("KEE_fit_mass_PF_B3r","Take(KEE_fit_mass_sk, PF_B3r)") \
                         .Define("W_PF_B3r","Take(W, PF_B3r)") \
                         .Define("KEE_mass_mix_B3f","Take(KEE_mass_sk, mix_B3f)") \
                         .Define("KEE_fit_mass_mix_B3f","Take(KEE_fit_mass_sk, mix_B3f)") \
                         .Define("W_mix_B3f","Take(W, mix_B3f)") \
                         .Define("KEE_mass_mix_B3r","Take(KEE_mass_sk, mix_B3r)") \
                         .Define("KEE_fit_mass_mix_B3r","Take(KEE_fit_mass_sk, mix_B3r)") \
                         .Define("W_mix_B3r","Take(W, mix_B3r)") \
                         .Define("KEE_mass_Low_B3f","Take(KEE_mass_sk, Low_B3f)") \
                         .Define("KEE_fit_mass_Low_B3f","Take(KEE_fit_mass_sk, Low_B3f)") \
                         .Define("W_Low_B3f","Take(W, Low_B3f)") \
                         .Define("KEE_mass_Low_B3r","Take(KEE_mass_sk, Low_B3r)") \
                         .Define("KEE_fit_mass_Low_B3r","Take(KEE_fit_mass_sk, Low_B3r)") \
                         .Define("W_Low_B3r","Take(W, Low_B3r)") \
                         .Define("KMM_mass_PF_B3f","Take(KMM_mass_sk, PF_B3f_mu)") \
                         .Define("KMM_fit_mass_PF_B3f","Take(KMM_fit_mass_sk, PF_B3f_mu)") \
                         .Define("W_PF_B3f_mu","Take(W_mu, PF_B3f_mu)") \
                         .Define("KMM_mass_PF_B3r","Take(KMM_mass_sk, PF_B3r_mu)") \
                         .Define("KMM_fit_mass_PF_B3r","Take(KMM_fit_mass_sk, PF_B3r_mu)") \
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


Mod.print_time('End') 
