import time
import datetime


def print_time(where): 
    ts_start = time.time()
    print("\n {} - Time: {}".format(where,datetime.datetime.fromtimestamp(ts_start).strftime('%Y-%m-%d %H:%M:%S'))) 


# Return indices of triplets with both leptons not being PFoverlap
noPFover_code = '''
using namespace ROOT::VecOps;
RVec<unsigned int> noPFover(RVec<unsigned int>& inp_idx,
                            RVec<unsigned int>& l1isPFover,
                            RVec<unsigned int>& l2isPFover){

    unsigned int n = inp_idx.size();
    RVec<unsigned int> out_idx;

    for (auto i=0; i<n; ++i){
            if( l1isPFover[i] == 0 && l2isPFover[i] == 0 ) out_idx.push_back(i);
    }
    return out_idx;
}
'''


# Return indices of triplets with both leptons being PF
bothPF_code = '''
RVec<unsigned int> bothPF(RVec<unsigned int>& inp_idx,
                          RVec<unsigned int>& l1isPF,
                          RVec<unsigned int>& l2isPF){

    unsigned int n = inp_idx.size();
    RVec<unsigned int> out_idx;                                            

    for (auto i=0; i<n; ++i){
            if( l1isPF[i] == 1 && l2isPF[i] == 1 ) out_idx.push_back(i);
    }
    return out_idx;
}
'''


# Return indices of triplets with both leptons being LowPt
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
