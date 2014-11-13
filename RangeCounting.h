#ifndef OP
#define OP
#include "./order_preserving.h"
#endif

#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<math.h>
#include<ctime>

#define UINT32 unsigned int
#define DEFINED_L 6 // DEFINED_L = ceil(sqrt(32)),
                    // each x-coodinate of point is stored in (DEFINES_L+1) bit(s)
#define DEFINED_B 4 // DEFINED_B = 32 / (DEFINED_L + 1)
#define PACK_INDEX(x) ((x) / DEFINED_B)

using namespace std;

int rec_count_inversions(int, int, vector<UINT32>&);
void normalizeQueryPoints(vector< pair<int,bool> >&, vector< pair<int,bool> >&);
void dividingQueryPoints(vector<UINT32>&, int, vector<UINT32>&, int*, vector<UINT32>&, int*);
void packingQueryPoints(vector< pair<int, bool> >&, vector<UINT32>&);
int count_inversions(vector< pair<int, bool> >&);

int rec_count_inversions(int p_size, int l, vector<UINT32>& P){
    int h=1, p0_size, p1_size;
    vector<UINT32> P0, P1;
    dividingQueryPoints(P, p_size, P0, &p0_size, P1, &p1_size);
    printf("p0_size:%d, p1_size:%d\n", p0_size, p1_size);
    exit(1);
    if( l <= DEFINED_L ){
        vector<UINT32> P0, P1;
        dividingQueryPoints(P, p_size, P0, &p0_size, P1, &p1_size);
        rec_count_inversions(p0_size, l-1, P0);
        rec_count_inversions(p1_size, l-1, P1);
    }else{
        h = DEFINED_L;
    }
    return 0;
};

void normalizeQueryPoints(vector< pair<int,bool> >& P, vector< pair<int,bool> >& ret){//{{{
    // normalizetion: make each value unique and carry out bucket sort
    vector<int> charCount(UB_ALPHABET_SIZE, 0),
                insetingIndices(UB_ALPHABET_SIZE);
    //  // count the occurence for each character
    vector< pair<int, bool> >::iterator it_s=P.begin(), end_it_s=P.end();
    for(; it_s!=end_it_s; it_s++){
        charCount[it_s->first]++;
    }
    //  // inseting indices(, that is to be updated)
    for(int i=0, rank=0, end_i=charCount.size(); i<end_i; i++){
        insetingIndices[i] = rank;
        rank += charCount[i];
    }
    //  // bucket sort
    int p_size = P.size();
    vector<int> sortedP(p_size),
                sortingArray(p_size);
    for(int i=0; i<p_size; i++){
        sortedP[ insetingIndices[ P[i].first ] ] = P[i].first;
        sortingArray[ insetingIndices[ P[i].first ] ] = i;
        insetingIndices[ P[i].first ]++;
    }
    //  // create inverseSortingArray and normalizedP
    for(int i=0, sa_i=sortingArray[i]; i<p_size; sa_i=sortingArray[++i]){
        ret[ sa_i ] = pair<int, bool>(i, P[ sa_i ].second);
    }
    // cout << "--- (debug) -------------------------------" << endl;
    // for(int i=0; i<p_size; i++){ cout << P[i].first << " ";}
    // cout << endl;
    // for(int i=0; i<p_size; i++){ cout << sortedP[i] << " ";}
    // cout << endl;
    // for(int i=0; i<p_size; i++){ cout << sortingArray[i] << " ";}
    // cout << endl;
    // for(int i=0; i<p_size; i++){ printf("%d(%d) ", ret[i].first, ret[i].second); }
    // cout << endl;
    // cout << "------------------------------------------" << endl;
};//}}}

void dividingQueryPoints(vector<UINT32>& P, int p_size, vector<UINT32>& P0, int* p0_size, vector<UINT32>& P1, int* p1_size){//{{{
    P0.resize(1, 0);
    P1.resize(1, 0);
    *p0_size = *p1_size = 0;
    const static UINT32 mask = (1 << (DEFINED_L+1)) - 1;
    for(int i=0, end_i=P.size(); i<end_i; i++){
        // printf("P[%d]=%u\n", i, P[i]);
        int end_j = (i != end_i-1) ? DEFINED_B : (p_size % DEFINED_B);
        for(int j=end_j-1; j>=0; j--){
            if(P[i] & ( 1 << ( j * (DEFINED_L+1) ) )){
                P1[ PACK_INDEX(*p1_size) ] <<= (DEFINED_L+1);
                P1[ PACK_INDEX(*p1_size) ] |= (( P[i] >> (j * (DEFINED_L+1)) ) & mask);
                if( !((++*p1_size) % DEFINED_B) ){ P1.push_back(0); }
            }else{
                P0[ PACK_INDEX(*p0_size) ] <<= (DEFINED_L+1);
                P0[ PACK_INDEX(*p0_size) ] |= (( P[i] >> (j * (DEFINED_L+1)) ) & mask);
                if( !((++*p0_size) % DEFINED_B) ){ P0.push_back(0); }
            }
        }
        // printf("P0 (%d-points are stored): ", *p0_size);
        // for(int i=0; i<P0.size(); i++){ printf("%u" ,P0[i]); }
        // cout << endl;
        // printf("P1 (%d-points are stored): ", *p1_size);
        // for(int i=0; i<P1.size(); i++){ printf("%u" ,P1[i]); }
        // cout << endl;
    }
};//}}}

void packingQueryPoints(vector< pair<int, bool> >& P, vector<UINT32>& ret){//{{{
    if( ret.size() < PACK_INDEX( P.size() ) ){
        fprintf(stderr, "need enough region for packing\n");
        exit(1);
    }
    for(int i=0, end_i=P.size(); i<end_i; i++){
        UINT32 p_i = (UINT32)P[i].first;
        p_i = (P[i].second) ? ((p_i << 1) + 1) : (p_i << 1);
        ret[ PACK_INDEX(i) ] <<= (DEFINED_L + 1);
        ret[ PACK_INDEX(i) ] |= p_i;
    }
};//}}}

int count_inversions(vector< pair<int, bool> >& P){
    vector< pair<int, bool> > normalizedP(P.size());
    normalizeQueryPoints(P, normalizedP);
    cout << "normalized P:\t";
    for(int j=0; j<P.size(); j++){
        printf("%d(%d) ", normalizedP[j].first, normalizedP[j].second);
    }
    cout << endl;
    int p_size = P.size(),
        pack_size = (p_size / DEFINED_B) + 1;
    vector<UINT32> packedP(pack_size, 0);
    packingQueryPoints(normalizedP, packedP);

    int l = (int)ceil( log2( P.size() ) );
    return rec_count_inversions(p_size, l, packedP);
};

class RangeCounting{//{{{
private:
    static int const w = 32; // w: 32 bits
    static int const L = 6; // L: ceil( sqrt(w) )
    static int const A = 5; // A: (1 + epsilon) * ceil( log(w) )
    static int const H = 5; // H: epsilon * ceil( log(w) )
    int n, l;
public:
    ~RangeCounting(){};
    RangeCounting(vector<int>&);
    int query(int, int);
};

RangeCounting::RangeCounting(vector<int>& S): n( (int)S.size() ){
};

int RangeCounting::query(int x, int y){
    if( l == 0 ){ // initialization
        l = (int)ceil(log(n));
        vector<int> bucket(UB_ALPHABET_SIZE);
    }

    if( (l <= H) and (n <= (1 << A)) ){ // case: 0'
    }else if( (l <= H) and (n > (1 << A)) ){ // case: 0''
    }else if( (H < l) and (l <= L) ){ // case: 1
    }else if( l > L ){ // case: 2
    }else{
        fprintf(stderr, "not matched for any case\n");
        exit(1);
    }

    return 0;
}//}}}

/* vim:set foldmethod=marker commentstring=//%s : */
