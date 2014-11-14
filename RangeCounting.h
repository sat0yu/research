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

#define UINT64 unsigned long long
#define UINT32 unsigned int
#define UINT16 unsigned short
#define DEFINED_w 16
#define NUM_POINTS_IN_WORD(l) ((DEFINED_w)/(l+1)) // 'l' meand # of bits represeinting each p.y
                                                  // and the rest one bit is its color
#define NUM_WORDS(n,l) ((int)(ceil((n)/(double)NUM_POINTS_IN_WORD(l))))
#define WORD_INDEX(i,l) ((i)/(NUM_POINTS_IN_WORD(l)))

using namespace std;

int rec_count_inversions(vector< vector< pair<UINT16, UINT16> > >&, int, int, vector<UINT16>&);
int count_inversions_tilde(int, vector<UINT16>&);
void normalizeQueryPoints(vector< pair<int,bool> >&, vector< pair<int,bool> >&);
void dividingQueryPoints(vector< vector< pair<UINT16, UINT16> > >&, vector<int>&, int, vector<UINT16>&, int, vector<UINT16>&, int*, vector<UINT16>&, int*);
void packingQueryPoints(int, vector< pair<int, bool> >&, vector<UINT16>&);
int count_inversions(vector< pair<int, bool> >&);
void construstWordOperationTable(vector< vector< pair<UINT16, UINT16> > >&);
void showWords(int, vector<UINT16>&);
void showOneWord(int, UINT16);

void showWords(int l, vector<UINT16>& z){//{{{
    for(int i=0; i<z.size(); i++){
        printf("z[%d]: ", i);
        for(int j=DEFINED_w-1; j>=0; j--){
            cout << (bool)( z[i] & (1<<j) ) << " ";
            if( !(j%l) ){ cout << " "; }
        }
        cout << endl;
    }
};//}}}

void showOneWord(int l, UINT16 z){//{{{
    for(int j=DEFINED_w-1; j>=0; j--){
        cout << (bool)( z & (1<<j) ) << " ";
        if( !(j%l) ){ cout << " "; }
    }
    cout << endl;
};//}}}

int rec_count_inversions(vector< vector< pair<UINT16, UINT16> > >& table, vector<int>& countingTable, int l, vector<UINT16>& P, int p_size){//{{{
    const static int L = (int)ceil(sqrt(DEFINED_w));
    int h=1, p0_size, p1_size;
    if( !(l>0) ){
        return 0;
    }else if( l <= L ){
        vector<UINT16> P0, P1;
        dividingQueryPoints(table, countingTable, l, P, p_size, P0, &p0_size, P1, &p1_size);
        printf("l:%d, p0_size:%d, p1_size:%d\n", l, p0_size, p1_size);
        return rec_count_inversions(table, countingTable, l-1, P0, p0_size)
                + rec_count_inversions(table, countingTable, l-1, P1, p1_size)
                + count_inversions_tilde(l, P);
    }else{
        h = L;
    }
    return 0;
};//}}}

int count_inversions_tilde(int l, vector<UINT16>& P_tilde){//{{{
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

void dividingQueryPoints(vector< vector< pair<UINT16, UINT16> > >& table, vector<int>& countingTable, int l, vector<UINT16>& P, int p_size, vector<UINT16>& P0, int* p0_size, vector<UINT16>& P1, int* p1_size){//{{{
    UINT16 mask = (1 << (l+1)) - 1, count_mask = 0;
    for(int i=0, end_i=NUM_POINTS_IN_WORD(l); i<end_i; i++){ // create the mask to be utilized in counting
        count_mask <<= (l+1);
        count_mask |= (1 << l);
    }
    (*p0_size) = (*p1_size) = 0;
    for(int i=0, end_i=P.size(); i<end_i; i++){ // count |P0| and |P1|
        (*p1_size) += countingTable[ P[i] & count_mask ];
        // UINT16 MSBs = P[i] & count_mask;
        // for(int j=0; j<NUM_POINTS_IN_WORD(l); j++){
        //     (*p1_size) += ( (MSBs >> (j * (l+1))) & mask ) ? 1 : 0;
        // }
    }
    (*p0_size) = p_size - (*p1_size);
    P0.resize(NUM_WORDS((*p0_size), l-1), 0); // reserve enough memory regions
    P1.resize(NUM_WORDS((*p1_size), l-1), 0);
    printf("p0_size:%d, p1_size:%d\n", *p0_size, *p1_size);
    printf("NUM_POINTS_IN_WORD(l-1):%d, P0.size():%d, P1.size():%d\n",
            NUM_POINTS_IN_WORD(l-1), NUM_WORDS((*p0_size), l-1), NUM_WORDS((*p1_size), l-1));

    int rest_p0, rest_p1, points_per_word;
    rest_p0 = rest_p1 = points_per_word = NUM_POINTS_IN_WORD(l-1);
    for(int i=0, end_i=P.size(), p0_idx=0, p1_idx=0; i<end_i; i++){
        int num_p0=0, num_p1=0;
        num_p1 = countingTable[ P[i] & count_mask ];
        if( (p_size % NUM_POINTS_IN_WORD(l)) and (i == end_i-1) ){ // the last elem. of P[i]
            num_p0 = (p_size % NUM_POINTS_IN_WORD(l)) - num_p1;
        }else{
            num_p0 = NUM_POINTS_IN_WORD(l) - num_p1;
        }
        // cout << "---------------------------" << endl;
        // showOneWord(l+1, P[i]);
        // showOneWord(l, table[l][P[i]].first);
        // showOneWord(l, table[l][P[i]].second);
        // cout << "---------------------------" << endl;

        // insert divided P[i], for each whose point, p.y is in [0, 2^(l-1))
        if(num_p0 <= rest_p0){
            // printf("case0: num_p0:%d, rest_p0:%d, points_per_word:%d\n", num_p0, rest_p0, points_per_word);
            P0[p0_idx] <<= (num_p0 * l);
            P0[p0_idx] |= table[l][P[i]].first;
            rest_p0 -= num_p0; // update the rest space
            if( rest_p0 == 0 ){ // increment the index of P0 as needed
                p0_idx++;
                rest_p0 = points_per_word;
            }
        }else{
            // printf("case1: num_p0:%d, rest_p0:%d, points_per_word:%d\n", num_p0, rest_p0, points_per_word);
            UINT16 inserting_mask = ((1 << (rest_p0 * l)) - 1), // a mask for the first half of divided P[i]
                   shifted = (table[l][P[i]].first >> ((num_p0 - rest_p0) * l));
            P0[p0_idx] <<= (rest_p0 * l); // shift left to reserve the space
            P0[p0_idx] |= (shifted & rest_p0);
            int remain = num_p0 - rest_p0;

            p0_idx++; // increment the index of P0
            rest_p0 = points_per_word;

            inserting_mask = ((1 << (remain * l)) - 1); // insert the last half of divided P[i]
            P0[p0_idx] |= (table[l][P[i]].first & inserting_mask);
            rest_p0 = points_per_word - remain;
        }

        // insert divided P[i], for each whose point, p.y is in [2^(l-1), 2^l)
        if(num_p1 <= rest_p1){
            // printf("case0: num_p1:%d, rest_p1:%d, points_per_word:%d\n", num_p1, rest_p1, points_per_word);
            P1[p1_idx] <<= (num_p1 * l);
            P1[p1_idx] |= table[l][P[i]].second;
            rest_p1 -= num_p1; // update the rest space
            if( rest_p1 == 0 ){ // increment the index of P1 as needed
                p1_idx++;
                rest_p1 = points_per_word;
            }
        }else{
            // printf("case1: num_p1:%d, rest_p1:%d, points_per_word:%d\n", num_p1, rest_p1, points_per_word);
            UINT16 inserting_mask = ((1 << (rest_p1 * l)) - 1), // a mask for the second half of divided P[i]
                   shifted = (table[l][P[i]].second >> ((num_p1 - rest_p1) * l));
            P1[p1_idx] <<= (rest_p1 * l); // shift left to reserve the space
            P1[p1_idx] |= (shifted & rest_p1);
            int remain = num_p1 - rest_p1;

            p1_idx++; // increment the index of P1
            rest_p1 = points_per_word;

            inserting_mask = ((1 << (remain * l)) - 1); // insert the last half of divided P[i]
            P1[p1_idx] |= (table[l][P[i]].second & inserting_mask);
            rest_p1 = points_per_word - remain;
        }

    }
    printf("P0[] (%d-points are stored):\n", *p0_size);
    showWords(l, P0);
    printf("P1[] (%d-points are stored):\n", *p1_size);
    showWords(l, P1);
    exit(1);
};//}}}

void packingQueryPoints(int l, vector< pair<int, bool> >& P, vector<UINT16>& ret){//{{{
    if( ret.size() < NUM_WORDS(P.size(),l) ){
        fprintf(stderr, "need enough region for packing\n");
        exit(1);
    }
    for(int i=0, end_i=P.size(); i<end_i; i++){
        UINT16 p_i = (UINT16)P[i].first;
        p_i = (P[i].second) ? ((p_i << 1) + 1) : (p_i << 1);
        ret[ WORD_INDEX(i,l) ] <<= (l + 1);
        ret[ WORD_INDEX(i,l) ] |= p_i;
    }
    cout << endl;
};//}}}

void construstWordOperationTable(vector< vector< pair<UINT16, UINT16> > >& table){//{{{
    const static int L = (int)ceil(sqrt(DEFINED_w)) + 1;
    for(int l = L; l > 0; l--){
        UINT16 mask = (1 << l) - 1;
        for(int w = 0, end_w =(1<<(DEFINED_w)); w < end_w; w++){
            // cout << "---------------------------" << endl;
            // showOneWord(l, w);
            // // 'l' meand # of bits represeinting each p.y and the rest one bit is its color
            UINT16 P0 = 0, P1 = 0;
            for(int i=NUM_POINTS_IN_WORD(l), i_shift=i*l; i>=0; i_shift=(--i)*l ){
                // printf("l:%d, mask:%u, w=%u, j_shift:%d, masked:%u, border:%u\n",
                //         l, mask, w, i_shift, (( w >> i_shift ) & mask), (1 << l));
                UINT clipped = (( w >> i_shift ) & mask);
                if( clipped & (1 << (l-1)) ){ // divide [0,2^(l-1)), [2^(l-1),2^l)
                    P1 <<= (l-1);
                    P1 |= (clipped & (mask >> 1)); // append clipped point with re-scaling
                }else{
                    P0 <<= (l-1);
                    P0 |= (clipped & (mask >> 1));
                }
            }
            // showOneWord(l-1, P0);
            // showOneWord(l-1, P1);
            // cout << "---------------------------" << endl;
            table[l-1][w].first = P0;
            table[l-1][w].second = P1;
        }
    }
};//}}}

int count16bit(UINT16 v){//{{{
    unsigned short count = (v & 0x5555) + ((v >> 1) & 0x5555);
    count = (count & 0x3333) + ((count >> 2) & 0x3333);
    count = (count & 0x0f0f) + ((count >> 4) & 0x0f0f);
    return (count & 0x00ff) + ((count >> 8) & 0x00ff);
}//}}}

int count_inversions(vector< pair<int, bool> >& P){//{{{
    const static int L = (int)ceil(sqrt(DEFINED_w)) + 1;
    vector< vector< pair<UINT16, UINT16> > > // the table for O(1) word operation
        table(L, vector< pair<UINT16, UINT16> >((1<<DEFINED_w), pair<UINT16, UINT16>(0,0)));
    construstWordOperationTable(table);
    vector<int> countingTable((1 << DEFINED_w), 0); // the table for O(1) word operation
    for(int i=0, end_i=countingTable.size(); i<end_i; i++){
        countingTable[i] = count16bit( (UINT16)i );
    }

    int p_size = P.size(),
        l = (int)ceil( log2(p_size) );
    printf("P.size():%d, l:%d, # of points in a word:%d, # of words:%d\n",
            p_size, l, NUM_POINTS_IN_WORD(l), NUM_WORDS(p_size,l));
    vector< pair<int, bool> > normalizedP(p_size);
    vector<UINT16> packedP(NUM_WORDS(p_size,l), 0);
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

    return rec_count_inversions(table, countingTable, l, packedP, p_size);
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
