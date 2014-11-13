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
#define DEFINED_w 32
#define DEFINED_L 6 // DEFINED_L = ceil(sqrt(32)),
                    // each x-coodinate of point is stored in (DEFINES_L+1) bit(s)
#define DEFINED_B 4 // DEFINED_B = 32 / (DEFINED_L + 1)
#define NUM_POINTS_IN_WORD(l) ((DEFINED_w)/(l+1)) // 'l' meand # of bits represeinting each p.y
                                                  // and the rest one bit is its color
#define NUM_WORDS(n,l) (((n)/NUM_POINTS_IN_WORD(l))+1)
#define WORD_INDEX(i,l) ((i)/(NUM_POINTS_IN_WORD(l)))

using namespace std;

int rec_count_inversions(int, int, vector<UINT32>&);
int count_inversions_tilde(int, vector<UINT32>&);
void normalizeQueryPoints(vector< pair<int,bool> >&, vector< pair<int,bool> >&);
void dividingQueryPoints(int, vector<UINT32>&, int, vector<UINT32>&, int*, vector<UINT32>&, int*);
void packingQueryPoints(int, vector< pair<int, bool> >&, vector<UINT32>&);
int count_inversions(vector< pair<int, bool> >&);
void showWord(int, vector<UINT32>&);

void showWords(int l, vector<UINT32>& z){//{{{
    for(int i=0; i<z.size(); i++){
        printf("z[%d]: ", i);
        for(int j=DEFINED_w-1; j>=0; j--){
            cout << (bool)( z[i] & (1<<j) ) << " ";
            if( !(j%l) ){ cout << " "; }
        }
        cout << endl;
    }
};//}}}

int rec_count_inversions(int l, vector<UINT32>& P, int p_size){//{{{
    int h=1, p0_size, p1_size;
    if( !(l>0) ){
        return 0;
    }else if( l <= DEFINED_L ){
        vector<UINT32> P0, P1;
        dividingQueryPoints(l, P, p_size, P0, &p0_size, P1, &p1_size);
        printf("l:%d, p0_size:%d, p1_size:%d\n", l, p0_size, p1_size);
        return rec_count_inversions(l-1, P0, p0_size)
                + rec_count_inversions(l-1, P1, p1_size)
                + count_inversions_tilde(l, P);
    }else{
        h = DEFINED_L;
    }
    return 0;
};//}}}

int count_inversions_tilde(int l, vector<UINT32>& P_tilde){//{{{
    if( !(l>0) ){return 0;}
    int ret = 0;
    for(int i=0, end_i=P_tilde.size(); i<end_i; i++){
    }
    return ret;
};//}}}

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

void dividingQueryPoints(int l, vector<UINT32>& P, int p_size, vector<UINT32>& P0, int* p0_size, vector<UINT32>& P1, int* p1_size){//{{{
    P0.resize(1, 0);
    P1.resize(1, 0);
    (*p0_size) = (*p1_size) = 0;
    UINT32 mask = (1 << (l+1)) - 1;
    for(int i=0, end_i=P.size(); i<end_i; i++){
        int end_j = (i != end_i-1) ? NUM_POINTS_IN_WORD(l) : (p_size % NUM_POINTS_IN_WORD(l));
        for(int j=end_j-1, j_shift=j*(l+1); j>=0; j_shift=(--j)*(l+1) ){
            // printf("l:%d, mask:%u, P[%d]=%u, j_shift:%d, masked:%u, border:%u\n",
            //         l, mask, i, P[i], j_shift, (( P[i] >> j_shift ) & mask), (1 << l));
            UINT clipped = (( P[i] >> j_shift ) & mask);
            if( clipped & (1 << l) ){ // divide [0,2^(l-1)), [2^(l-1),2^l)
                // P1[ WORD_INDEX(*p1_size,l) ] <<= (l+1);
                P1[ WORD_INDEX(*p1_size,l) ] <<= l;
                P1[ WORD_INDEX(*p1_size,l) ] |= (clipped & (mask >> 1)); // append clipped point with re-scaling
                if( !((++*p1_size) % NUM_POINTS_IN_WORD(l)) ){ P1.push_back(0); }
            }else{
                // P0[ WORD_INDEX(*p0_size,l) ] <<= (l+1);
                P0[ WORD_INDEX(*p0_size,l) ] <<= l;
                P0[ WORD_INDEX(*p0_size,l) ] |= (clipped & (mask >> 1));
                if( !((++*p0_size) % NUM_POINTS_IN_WORD(l)) ){ P0.push_back(0); }
            }
        }
    }
    // printf("P0 (%d-points are stored):\n", *p0_size);
    // showWords(l, P0);
    // printf("P1 (%d-points are stored):\n", *p1_size);
    // showWords(l, P1);
};//}}}

void packingQueryPoints(int l, vector< pair<int, bool> >& P, vector<UINT32>& ret){//{{{
    if( ret.size() < NUM_WORDS(P.size(),l) ){
        fprintf(stderr, "need enough region for packing\n");
        exit(1);
    }
    for(int i=0, end_i=P.size(); i<end_i; i++){
        UINT32 p_i = (UINT32)P[i].first;
        p_i = (P[i].second) ? ((p_i << 1) + 1) : (p_i << 1);
        ret[ WORD_INDEX(i,l) ] <<= (l + 1);
        ret[ WORD_INDEX(i,l) ] |= p_i;
    }
    cout << endl;
};//}}}

int count_inversions(vector< pair<int, bool> >& P){//{{{
    int p_size = P.size(),
        l = (int)ceil( log2(p_size) );
    printf("P.size():%d, l:%d, # of points in a word:%d, # of words:%d\n", p_size, l, NUM_POINTS_IN_WORD(l), NUM_WORDS(p_size,l));
    vector< pair<int, bool> > normalizedP(p_size);
    vector<UINT32> packedP(NUM_WORDS(p_size,l), 0);
    // normalize given points in P
    normalizeQueryPoints(P, normalizedP);
    cout << "normalized P:\t";
    for(int j=0; j<p_size; j++){
        printf("%d(%d) ", normalizedP[j].first, normalizedP[j].second);
    }
    cout << endl;
    // pack some points into a word
    packingQueryPoints(l, normalizedP, packedP);
    printf("points in P are packed:\n");
    showWords(l+1, packedP);

    return rec_count_inversions(l, packedP, p_size);
};//}}}

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
