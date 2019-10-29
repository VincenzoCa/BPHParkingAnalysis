import time
import datetime


def print_time(where): 
    ts_start = time.time()
    print("\n {} - Time: {}".format(where,datetime.datetime.fromtimestamp(ts_start).strftime('%Y-%m-%d %H:%M:%S'))) 


#returns indices of triplets passing selections
Cuts_code = '''
ROOT::VecOps::RVec<unsigned int> Cuts(bool isE, 
                                      unsigned int nTr,
                                      ROOT::VecOps::RVec<float>& BpT, 
                                      ROOT::VecOps::RVec<float>& cosAlpha, 
                                      ROOT::VecOps::RVec<float>& svprob,
                                      ROOT::VecOps::RVec<float>& Lxy,
                                      ROOT::VecOps::RVec<float>& Lxy_unc,
                                      ROOT::VecOps::RVec<float>& KpT,
                                      ROOT::VecOps::RVec<float>& l1pT,
                                      ROOT::VecOps::RVec<float>& l2pT,
                                      ROOT::VecOps::RVec<float>& l1mvaId,
                                      ROOT::VecOps::RVec<float>& l2mvaId,
                                      ROOT::VecOps::RVec<bool>& l1Veto,
                                      ROOT::VecOps::RVec<bool>& l2Veto){

    ROOT::VecOps::RVec<unsigned int> out_idx;

    float pT1, pT2;
    pT1 = isE ? 1.5 : 2.;
    pT2 = isE ? 0.5 : 2;

    for (auto i=0; i<nTr; ++i){
            if( BpT[i] < 3. || cosAlpha[i] < 0.999 || svprob[i] < 0.1 ||
                Lxy[i]/sqrt(Lxy_unc[i]) < 6.0 || KpT[i] < 3. || 
                l1pT[i] < pT1 || l2pT[i] < pT2) continue;
            if( isE && l1mvaId[i] > 3.96 && l2mvaId[i] > 3.96 && 
                l1Veto[i] && l2Veto[i] ) out_idx.push_back(i);
    }

    return out_idx;
}
'''


#return indices of triplets with  b1 < mll < b2
Bin_code = '''
ROOT::VecOps::RVec<unsigned int> Bin(ROOT::VecOps::RVec<unsigned int>& inp_idx,
                                     float b1,float b2,
                                     ROOT::VecOps::RVec<float>& mll){    

    ROOT::VecOps::RVec<unsigned int> out_idx;

    for (auto i : inp_idx){
            if( mll[i] >= b1 && mll[i] < b2 )out_idx.push_back(i);
    }
    return out_idx;
}
'''


#return indices of triplets with both leptons not being PFoverlap
noPFover_code = '''
ROOT::VecOps::RVec<unsigned int> noPFover(ROOT::VecOps::RVec<unsigned int>& inp_idx,
                                          ROOT::VecOps::RVec<bool>& l1isPFover,
                                          ROOT::VecOps::RVec<bool>& l2isPFover){

    ROOT::VecOps::RVec<unsigned int> out_idx;
    for (auto i : inp_idx){
            if( !l1isPFover[i] && !l1isPFover[i] ) out_idx.push_back(i);
    }
    return out_idx;
}
'''


#return indices of triplets with both leptons being PF
bothPF_code = '''
ROOT::VecOps::RVec<unsigned int> bothPF(ROOT::VecOps::RVec<unsigned int>& inp_idx,
                                        ROOT::VecOps::RVec<bool>& l1isPF,
                                        ROOT::VecOps::RVec<bool>& l2isPF){

    ROOT::VecOps::RVec<unsigned int> out_idx;                                            
    for (auto i : inp_idx){
            if( l1isPF[i] && l1isPF[i] ) out_idx.push_back(i);
    }
    return out_idx;
}
'''


#return indices of triplets with both leptons being lowPt
bothLow_code = '''
ROOT::VecOps::RVec<unsigned int> bothLow(ROOT::VecOps::RVec<unsigned int>& inp_idx,
                                         ROOT::VecOps::RVec<bool>& l1isLow,                                         
                                         ROOT::VecOps::RVec<bool>& l2isLow){

    ROOT::VecOps::RVec<unsigned int> out_idx;
    for (auto i : inp_idx){
            if( l1isLow[i] && l1isLow[i] ) out_idx.push_back(i);
    }
    return out_idx;
}
'''
