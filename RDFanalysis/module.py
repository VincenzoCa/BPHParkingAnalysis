import time
import datetime


def print_time(where): 
    ts_start = time.time()
    print("\n {} - Time: {}".format(where,datetime.datetime.fromtimestamp(ts_start).strftime('%Y-%m-%d %H:%M:%S'))) 


# Returns indices of triplets passing B and K selections
Cuts_code = '''
using namespace ROOT::VecOps;
RVec<unsigned int> Cuts(unsigned int nB,
                        RVec<float>& BpT, 
                        RVec<float>& cosAlpha, 
                        RVec<float>& svprob,
                        RVec<float>& LxySig,
                        RVec<float>& K_DCASig,
                        RVec<float>& KpT){
    RVec<unsigned int> out_idx(nB, 0);
    for (auto i=0; i<nB; ++i){
            if( BpT[i] < 3. || cosAlpha[i] < 0.999 || svprob[i] < 0.1 ||
                LxySig[i] < 6.0 || K_DCASig[i] < 2. || KpT[i] < 3.) continue;
            out_idx[i] = 1;            
    }
    return out_idx;
}
'''


# Returns indices of good KEE triplets
EleCuts_code = '''
RVec<unsigned int> EleCuts(RVec<unsigned int>& inp_idx,
                           RVec<float>& l1pT,
                           RVec<float>& l2pT,
                           RVec<float>& l1mvaId,
                           RVec<float>& l2mvaId,
                           RVec<bool>& l1Veto,
                           RVec<bool>& l2Veto){
    RVec<unsigned int> out_idx(inp_idx.size(), 0);
    auto idx = Nonzero(inp_idx);
    for (auto i : idx){
    //old ID: 3.94
            if( l1mvaId[i] > 4.24 && l2mvaId[i] > 4.24 &&
                l1Veto[i] && l2Veto[i] &&
                l1pT[i] > 1.5 && l2pT[i] > 0.5) out_idx[i] = 1;    
    }
    return out_idx;
} 
'''


# Returns indices of good KMuMu triplets
MuCuts_code = '''
RVec<unsigned int> MuCuts(RVec<unsigned int>& inp_idx,
                          RVec<float>& l1pT,
                          RVec<float>& l2pT,
                          RVec<unsigned int>& nTrMu){
    
    RVec<unsigned int> out_idx(inp_idx.size(), 0);
    auto idx = Nonzero(inp_idx);    
    for (auto i : idx){    
            if( l1pT[i] > 1.5 && l2pT[i] > 0.5 && nTrMu[i] > 0) out_idx[i] = 1;
    }
    return out_idx;
}
'''


# Return indices of triplets with  b1 < mll < b2
Bin_code = '''
RVec<unsigned int> Bin(RVec<unsigned int>& inp_idx,
                       float b1,float b2,
                       RVec<float>& mll){    
    RVec<unsigned int> out_idx;
    auto idx = Nonzero(inp_idx);
    for (auto i : idx){
            if( mll[i] >= b1 && mll[i] < b2 )out_idx.push_back(i);
    }
    return out_idx;
}
'''


Bin_MC_code = '''
RVec<unsigned int> Bin_MC(RVec<unsigned int>& inp_idx,
                          float b1,float b2,
                          RVec<float>& mll){
    RVec<unsigned int> out_idx(inp_idx.size(), 0);
    auto idx = Nonzero(inp_idx);
    for (auto i : idx){
            if( mll[i] >= b1 && mll[i] < b2 )out_idx[i] = 1;
    }
    return out_idx; 
} 
'''


# Return indices of triplets with both leptons not being PFoverlap
noPFover_code = '''
RVec<unsigned int> noPFover(RVec<unsigned int>& inp_idx,
                            RVec<unsigned int>& l1isPFover,
                            RVec<unsigned int>& l2isPFover){
    
    RVec<unsigned int> out_idx(inp_idx.size(), 0);
    auto idx = Nonzero(inp_idx);
    for (auto i : idx){    
            if( l1isPFover[i] == 0 && l2isPFover[i] == 0 ) out_idx[i] = 1;
    }
    return out_idx;
}
'''


# Return indices of triplets with both leptons being PF (or LowPt)
bothX_code = '''
RVec<unsigned int> bothX(RVec<unsigned int>& inp_idx,
                         RVec<unsigned int>& l1isX,
                         RVec<unsigned int>& l2isX){
    RVec<unsigned int> out_idx(inp_idx.size(), 0);
    auto idx = Nonzero(inp_idx);    
    for (auto i : idx){
            if( l1isX[i] == 1 && l2isX[i] == 1 ) out_idx[i] = 1;
    }
    return out_idx;
}
'''


# Return indices of mixed triplets: (l1_is_PF, l2_is_lowPt) or (l1_is_lowPt, l2_is_PF)
mix_code = '''
RVec<unsigned int> mix(RVec<unsigned int>& inp_idx,
                       RVec<unsigned int>& l1isPF,
                       RVec<unsigned int>& l2isPF,
                       RVec<unsigned int>& l1isLow,
                       RVec<unsigned int>& l2isLow){
    RVec<unsigned int> out_idx(inp_idx.size(), 0);
    auto idx = Nonzero(inp_idx);   
    for (auto i : idx){
            if( (l1isPF[i] == 1 && l2isLow[i] ==1) || (l1isLow[i] ==1 && l2isPF[i] == 1) ) out_idx[i] = 1;
    }
    return out_idx;
}
'''


# Return indices of gen-matched triplets (MC)
flagGenMatchExt_code = '''
RVec<unsigned int> flagGenMatchExt(bool& isKstar,
                                   bool& isResonant,//int
                                   RVec<int>& l1_pdgId, RVec<int>& l2_pdgId,
                                   RVec<int>& k_pdgId, RVec<int>& Ml1_pdgId,
                                   RVec<int>& Ml2_pdgId, RVec<int>& Mk_pdgId,
                                   RVec<int>& GMl1_pdgId, RVec<int>& GMl2_pdgId,
                                   RVec<int>& GMk_pdgId){
  auto totN = l1_pdgId.size();
  RVec<unsigned int> matched(totN, 0);
  for(unsigned int ij=0; ij<totN; ++ij){
                       
    if(l1_pdgId[ij] == -1 || l2_pdgId[ij] == -1 || k_pdgId[ij] == -1) continue;
    if(l1_pdgId[ij] != -1. * l2_pdgId[ij] || std::abs(k_pdgId[ij]) != 321) continue;
    
    if(!isKstar){
    //B+ = 521, J/Psi(1S) = 443
        if(isResonant){
            if( Ml1_pdgId[ij] == Ml2_pdgId[ij] && Ml1_pdgId[ij] == 443 &&
                GMl2_pdgId[ij] == GMl1_pdgId[ij] && GMl1_pdgId[ij] == Mk_pdgId[ij] &&
                std::abs(GMl1_pdgId[ij]) == 521)
            matched[ij] = 1;
        }
        else if (!isResonant){
            if(Ml1_pdgId[ij] == Ml2_pdgId[ij] && Ml2_pdgId[ij] == Mk_pdgId[ij] &&
               std::abs(Ml1_pdgId[ij]) == 521)
            matched[ij] = 1;
        }
        
    }
    else if(isKstar){
    //B0 = 511, K*(892)0 = 313, J/Psi(1S) = 443
        if(isResonant){
            if( Ml1_pdgId[ij] == Ml2_pdgId[ij] && Ml1_pdgId[ij] == 443 &&
                GMl2_pdgId[ij] == GMl1_pdgId[ij] && GMl1_pdgId[ij] == GMk_pdgId[ij] &&
                std::abs(GMk_pdgId[ij]) == 511 &&
                std::abs(Mk_pdgId[ij]) == 313)
            matched[ij] = 1;
        }
        else if (!isResonant){
            if(Ml1_pdgId[ij] == Ml2_pdgId[ij] && Ml2_pdgId[ij] == GMk_pdgId[ij] &&
               std::abs(GMk_pdgId[ij]) == 511 && 
               std::abs(Mk_pdgId[ij]) == 313)
            matched[ij] = 1;
        }    
    
    }
    
  }
  return matched;
}
'''


# Compute dR for gen-matched triplets (-1 if not gen-matched)
computedR_code = '''
RVec<float> computedR(RVec<unsigned int>& matchedIdxs,
                      RVec<float>& gen_eta1, RVec<float>& gen_eta2, RVec<float>& gen_eta3,
                      RVec<float>& reco_eta1, RVec<float>& reco_eta2, RVec<float>& reco_eta3,
                      RVec<float>& gen_phi1, RVec<float>& gen_phi2, RVec<float>& gen_phi3,
                      RVec<float>& reco_phi1, RVec<float>& reco_phi2, RVec<float>& reco_phi3){
  auto idx = Nonzero(matchedIdxs);
  auto genEta1 = Take(gen_eta1, idx);
  auto genEta2 = Take(gen_eta2, idx);
  auto genEta3 = Take(gen_eta3, idx);
  auto genPhi1 = Take(gen_phi1, idx);
  auto genPhi2 = Take(gen_phi2, idx);
  auto genPhi3 = Take(gen_phi3, idx);
  auto recoEta1 = Take(reco_eta1, idx);
  auto recoEta2 = Take(reco_eta2, idx);
  auto recoEta3 = Take(reco_eta3, idx);
  auto recoPhi1 = Take(reco_phi1, idx);
  auto recoPhi2 = Take(reco_phi2, idx);
  auto recoPhi3 = Take(reco_phi3, idx);
  RVec<float> dR_1 = DeltaR(genEta1, recoEta1, genPhi1, recoPhi1);
  RVec<float> dR_2 = DeltaR(genEta2, recoEta2, genPhi2, recoPhi2);
  RVec<float> dR_3 = DeltaR(genEta3, recoEta3, genPhi3, recoPhi3);
  RVec<float> sum(gen_eta1.size(), -1);
  int selected = 0;
  for(unsigned int ij=0; ij<matchedIdxs.size(); ++ij){
    if(matchedIdxs[ij] == 1){
      sum[ij] = dR_1[selected] + dR_2[selected] + dR_3[selected];
      ++selected;
    }
    else sum[ij] = -1;
  }
  return sum;
}
'''


# Return for each event the index of the triplet with the best dR (MC)
bestRank_code = '''
RVec<unsigned int> bestRank(unsigned int nB,
                            RVec<unsigned int>& goodB,
                            RVec<unsigned int>& Idx_Bin,
                            RVec<unsigned int>& isGenMatch,
                            RVec<float>& dR){
  RVec<unsigned int> rank(nB, 0);
  auto sortIndices = Argsort(dR);
  for(unsigned int ij=0; ij<nB; ++ij){
    if(goodB[sortIndices[ij]] && Idx_Bin[sortIndices[ij]] && isGenMatch[sortIndices[ij]]){
        rank[sortIndices[ij]] = 1;
        break;
    }
  }
  return rank;
}
''' 
