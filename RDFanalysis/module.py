import time
import datetime


def print_time(where): 
    ts_start = time.time()
    print("\n {} - Time: {}".format(where,datetime.datetime.fromtimestamp(ts_start).strftime('%Y-%m-%d %H:%M:%S'))) 


#returns indices of triplets passing selections
Cuts_code = '''
using namespace ROOT::VecOps;
RVec<unsigned int> Cuts(unsigned int nTr,
                        RVec<float>& BpT, 
                        RVec<float>& cosAlpha, 
                        RVec<float>& svprob,
                        RVec<float>& Lxy,
                        RVec<float>& Lxy_unc,
                        RVec<float>& K_DCASig,
                        RVec<float>& KpT,
                        RVec<float>& l1pT,
                        RVec<float>& l2pT,
                        RVec<float>& l1mvaId,
                        RVec<float>& l2mvaId,
                        RVec<bool>& l1Veto,
                        RVec<bool>& l2Veto){

    RVec<unsigned int> out_idx;

    for (auto i=0; i<nTr; ++i){
            if( BpT[i] < 3. || cosAlpha[i] < 0.999 || svprob[i] < 0.1 ||
                Lxy[i]/(Lxy_unc[i]) < 6.0 || K_DCASig[i] < 2. || KpT[i] < 3.) continue;
            if( l1mvaId[i] > 3.94 && l2mvaId[i] > 3.94 && 
                l1Veto[i] && l2Veto[i] && 
                l1pT[i] > 1.5 && l2pT[i] > 0.5) out_idx.push_back(i);
    }

    return out_idx;
}
'''


#return indices of triplets with  b1 < mll < b2
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


#return indices of triplets with both leptons not being PFoverlap
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


#return indices of triplets with both leptons being PF
bothPF_code = '''
RVec<unsigned int> bothPF(RVec<unsigned int>& inp_idx,
                          RVec<unsigned int>& l1isPF,
                          RVec<unsigned int>& l2isPF){

    RVec<unsigned int> out_idx;                                            
    for (auto i : inp_idx){
            if( l1isPF[i] == 1 && l2isPF[i] == 1 ) out_idx.push_back(i);
    }
    return out_idx;
}
'''


#return indices of triplets with both leptons being lowPt
bothLow_code = '''
RVec<unsigned int> bothLow(RVec<unsigned int>& inp_idx,
                           RVec<unsigned int>& l1isLow,                                         
                           RVec<unsigned int>& l2isLow){

    RVec<unsigned int> out_idx;
    for (auto i : inp_idx){
            if( l1isLow[i] == 1 && l2isLow[i] == 1 ) out_idx.push_back(i);
    }
    return out_idx;
}
'''


#return indices of mixed triplets: (l1_is_PF, l2_is_lowPt) or (l1_is_lowPt, l2_is_PF)
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


#return weights for histograms
weights_code = '''
RVec<float> weights(unsigned int nTr){
    
    RVec<float> w;
    for (auto i=0; i<nTr; ++i){
            w.push_back(0.999);
    }
    return w;
}
'''
