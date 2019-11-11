import time
import datetime


def print_time(where): 
    ts_start = time.time()
    print("\n {} - Time: {}".format(where,datetime.datetime.fromtimestamp(ts_start).strftime('%Y-%m-%d %H:%M:%S'))) 


# Returns indices of triplets passing B and K selections
Cuts_code = '''
using namespace ROOT::VecOps;
RVec<unsigned int> Cuts(unsigned int nTr,
                        RVec<float>& BpT, 
                        RVec<float>& cosAlpha, 
                        RVec<float>& svprob,
                        RVec<float>& LxySig,
                        RVec<float>& K_DCASig,
                        RVec<float>& KpT){

    RVec<unsigned int> out_idx;

    for (auto i=0; i<nTr; ++i){
            if( BpT[i] < 3. || cosAlpha[i] < 0.999 || svprob[i] < 0.1 ||
                LxySig[i] < 6.0 || K_DCASig[i] < 2. || KpT[i] < 3.) continue;
            out_idx.push_back(i);
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

    RVec<unsigned int> out_idx;
    
    for (auto i : inp_idx){
            if( l1mvaId[i] > 3.94 && l2mvaId[i] > 3.94 &&
                l1Veto[i] && l2Veto[i] &&
                l1pT[i] > 1.5 && l2pT[i] > 0.5) out_idx.push_back(i);            
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

    RVec<unsigned int> out_idx;

    for (auto i : inp_idx){
            if( l1pT[i] > 1.5 && l2pT[i] > 0.5 && nTrMu[i] > 0) out_idx.push_back(i);
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

    for (auto i : inp_idx){
            if( mll[i] >= b1 && mll[i] < b2 )out_idx.push_back(i);
    }
    return out_idx;
}
'''


# Return indices of triplets with both leptons not being PFoverlap
noPFover_code = '''
RVec<unsigned int> noPFover(RVec<unsigned int>& inp_idx,
                            RVec<unsigned int>& l1isPFover,
                            RVec<unsigned int>& l2isPFover){

    RVec<unsigned int> out_idx;
    for (auto i : inp_idx){
            if( l1isPFover[i] == 0 && l2isPFover[i] == 0 ) out_idx.push_back(i);
    }
    return out_idx;
}
'''


# Return indices of triplets with both leptons being PF (or LowPt)
bothX_code = '''
RVec<unsigned int> bothX(RVec<unsigned int>& inp_idx,
                         RVec<unsigned int>& l1isX,
                         RVec<unsigned int>& l2isX){

    RVec<unsigned int> out_idx;                                            
    for (auto i : inp_idx){
            if( l1isX[i] == 1 && l2isX[i] == 1 ) out_idx.push_back(i);
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

    RVec<unsigned int> out_idx;
    for (auto i : inp_idx){
            if( (l1isPF[i] == 1 && l2isLow[i] ==1) || (l1isLow[i] ==1 && l2isPF[i] == 1) ) out_idx.push_back(i);
    }
    return out_idx;
}
'''
