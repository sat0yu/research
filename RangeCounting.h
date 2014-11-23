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

#define WORD_SIZE 16 // w: 16 bits ( w = epsilon * log(n), then 2^16 = 65536 )
#define L 4 // L: ceil( sqrt(WORD_SIZE) )
#define A 6 // A: (1 + epsilon) * ceil( log(WORD_SIZE) ), given epsilon=1/2, A is 6
#define H 2 // H: epsilon * ceil( log(WORD_SIZE) ), given epsilon=1/2, H is 2

#define UINT64 unsigned long long
#define UINT32 unsigned int
#define UINT16 unsigned short
#define BITS_OF(x) (8 * (int)sizeof(x))
#define MIN(x,y) (((x) < (y)) ? (x) : (y))

using namespace std;

struct Word{
private:
public:
    UINT16 bits, head_mask;
    size_t n, l, g, rest;
    Word(int l):
        bits(0), head_mask((1<<l)-1), n(0), l(l), g(WORD_SIZE/l), rest(WORD_SIZE/l){};
    Word(const Word &rhs):
        bits(rhs.bits), head_mask(rhs.head_mask), n(rhs.n), l(rhs.l), g(rhs.g), rest(rhs.rest){};
    void append(int, UINT16);
    void push_back(int);
    void showBits() const;
    int lessThan(int) const;
    int lessThanAt(int) const;
    int lessThanAt(int, int) const;
    int operator[](int) const;
};
void Word::append(int m, UINT16 w){//{{{
    if( m > rest ){
        fprintf(stderr, "there does not remain enough bits in a word\n");
        exit(1);
    }
    UINT16 masked = w & ( ( 1 << (m * l) ) - 1 );
    bits |= (masked << (l * n));
    n += m;
    rest -= m;
};//}}}
void Word::push_back(int i){//{{{
    if( !(rest > 0) ){
        fprintf(stderr, "there does not remain enough bits in a word\n");
        exit(1);
    }
    UINT16 masked = i & head_mask; // slice out useless bits
    bits |= (masked << (l * n++));
    rest--;
};//}}}
void Word::showBits() const{//{{{
    for(int j=WORD_SIZE-1; j>=0; j--){
        cout << (bool)(bits & (1<<j)) << " ";
        if( !(j%l) ){ cout << " "; }
    }
    cout << endl;
};//}}}
int Word::lessThan(int i) const{//{{{
    // // return the number of points p_i in the word,
    // // whose p_i.y is less than i.

    // ==================================================
    // // this function works in O(WORD_SIZE/l)-time
    // ==================================================
    int ret = 0;
    for(int j=0; j<g; j++){
        if( ((bits >> (l*j)) & head_mask) < (UINT16)i ){ ret++; }
    }
    return ret;
};//}}}
int Word::lessThanAt(int pos) const{//{{{
    // // return the number of points p_i in the word,
    // // whose p_i.y is less than an integer at position pos.

    // ==================================================
    // // this function works in O(WORD_SIZE/l)-time
    // ==================================================
    if( !( (pos >= 0) and (pos < n) ) ){
        fprintf(stderr, "invalid query; given pos=%d\n", pos);
        exit(1);
    }
    int ret = 0;
    UINT16 q = ( bits >> (pos * l) ) & head_mask;
    for(int i=0; i<pos; i++){
        if( ( ( bits >> (l * i) ) & head_mask ) < q ){ ret++; }
    }
    return ret;
};//}}}
int Word::lessThanAt(int x, int y) const{//{{{
    // // return the number of points p_i in the word,
    // // whose p_i.y is less than y and is at before x.

    // ==================================================
    // // this function works in O(WORD_SIZE/l)-time
    // ==================================================
    if( !( (x >= 0) and (y >= 0) ) ){
        fprintf(stderr, "invalid query; given (x,y)=(%d,%d)\n", x, y);
        exit(1);
    }
    int ret = 0;
    for(int i=0; i<x; i++){
        if( ( ( bits >> (l * i) ) & head_mask ) < (UINT16)y ){ ret++; }
    }
    return ret;
};//}}}
int Word::operator[](int i) const{//{{{
    if( (i < 0) or (i >= n) ){
        fprintf(stderr, "out of range reference\n");
        exit(1);
    }
    return (int)( ( bits & ( head_mask << (i * l) ) ) >> (i * l) );
}//}}}


class PackedIntegers{
private:
    UINT16 head_mask;
    vector< vector<int> > C;
    void _append(int, UINT16);
    void constructCountingMask();
public:
    int n, l, g;
    UINT16 counting_mask;
    ~PackedIntegers(){};
    PackedIntegers(int l, vector<int>&);
    PackedIntegers(int l);
    PackedIntegers(const PackedIntegers& rhs);
    PackedIntegers();
    vector<Word> words;
    void append(int, UINT16);
    void push_integer(int);
    int get_integer(int);
    int rc_query(int);
    int rc_query(int, int);
    void showWords() const;
    void showOneWord(UINT16) const;
    void constructDataStructure();
    void sequence(vector<int>&);
};
PackedIntegers::PackedIntegers(int l, vector<int>& P)://{{{
n(0), l(l), words(1,Word(l)), head_mask((1 << l) - 1), counting_mask(0), g(WORD_SIZE/l){
    for(int i=0, end_i=P.size(); i<end_i; i++){ // add each element in given P
        if( !(P[i] < (1 << l)) ){
            fprintf(stderr, "too large value in P\n");
            exit(1);
        }
        push_integer(P[i]);
    }
    constructCountingMask();
    // showWords();
};//}}}
PackedIntegers::PackedIntegers(int l)://{{{
n(0), l(l), words(1,Word(l)), head_mask((1 << l) - 1), counting_mask(0), g(WORD_SIZE/l){
    constructCountingMask();
};//}}}
PackedIntegers::PackedIntegers(const PackedIntegers& rhs)://{{{
n(rhs.n), l(rhs.l), words(rhs.words), counting_mask(rhs.counting_mask),
head_mask((1 << rhs.l) - 1),g(rhs.g), C(rhs.C){};//}}}
PackedIntegers::PackedIntegers()://{{{
n(0), l(BITS_OF(UINT16)), words(1,Word(l)), counting_mask(1 << (BITS_OF(UINT16) - 1)),
head_mask((1 << BITS_OF(UINT16)) - 1), g(1){};//}}}
void PackedIntegers::constructCountingMask(){//{{{
    for(int i = 0; i < g; i++){
        // (1 << (l-1)) represents the MSB of the head integer in a word
        counting_mask |= ((1 << (l-1)) << (i*l));
    }
    // cout << "counting_mask:";
    // showOneWord(counting_mask);
};//}}}
void PackedIntegers::showWords() const{//{{{
    cout << "--------------------------------" << endl;
    for(int i=0; i<words.size(); i++){
        printf("z[%d]: ", i);
        words[i].showBits();
    }
    cout << "--------------------------------" << endl;
};//}}}
void PackedIntegers::showOneWord(UINT16 z) const{//{{{
    cout << "--------------------------------" << endl;
    for(int j=WORD_SIZE-1; j>=0; j--){
        cout << (bool)( z & (1<<j) ) << " ";
        if( !(j%l) ){ cout << " "; }
    }
    cout << endl;
    cout << "--------------------------------" << endl;
};//}}}
void PackedIntegers::push_integer(int i){//{{{
    words.back().push_back(i);
    n++;
    if(words.back().rest == 0){
        words.push_back(Word(l));
    }
};//}}}
int PackedIntegers::get_integer(int i){//{{{
    if( i >= n ){
        fprintf(stderr, "out of range reference; given i=%d\n", i);
        exit(1);
    }
    // printf("words[%d][%d]:%d\n", i/g, i%g, words[i/g][i%g]);
    return (int)( words[i/g][i%g] );
};//}}}
int PackedIntegers::rc_query(int x){//{{{
    int i = x / g, j = get_integer(x);
    return C.at(i).at(j) + words.at(i).lessThanAt(x % g);
};//}}}
int PackedIntegers::rc_query(int x, int y){//{{{
    int i = x / g;
    return C.at(i).at(y) + words.at(i).lessThanAt(x % g, y);
};//}}}
void PackedIntegers::_append(int m, UINT16 w){//{{{
    words.back().append(m, w);
    n += m;
    if(words.back().rest == 0){ // reserve a space for the next word
        words.push_back(Word(l));
    }
};//}}}
void PackedIntegers::append(int m, UINT16 w){//{{{
    // showOneWord(w);
    if( m <= words.back().rest ){
        _append(m, w);
    }else{
        int temp = words.back().rest; // append dividing into two
        _append(temp, w);
        _append(m - temp, ( w >> (l * temp) ));
    }
};//}}}
void PackedIntegers::constructDataStructure(){//{{{
    C = vector<vector<int> >(
            // taking care of the case given a greater x (,y) than n (,P[n-1] respectively),
            // reserve spaces in surplus for C~
            words.size() + 1, vector<int>((1 << H) + 1,0)
        );
    for(int i=1, end_i=C.size(); i<end_i; i++){ // skip C_0i, i in [0,2~H)
        for(int j=0, end_j=C[0].size(); j<end_j; j++){
            C[i][j] = C[i-1][j] + words[i-1].lessThan(j);
        }
    }
};//}}}
void PackedIntegers::sequence(vector<int>& res){//{{{
    res.resize(n, 0);
    for(int i=0, end_i=words.size(); i<end_i; i++){
        int idx = i * g;
        for(int j=0, end_j=words[i].n; j<end_j; j++){
            res[idx+j] = (int)words[i][j];
        }
    };
};//}}}


class RangeCounting{
private:
    static vector< vector< pair<UINT16, UINT16> > > *splitingTable;
    int n, l, log2n, pow2A, case00_h;
    vector<int> P;
    vector< vector<int> > case00_pTilde;
    vector< RangeCounting > case1_sublists;
    vector< RangeCounting > case2_sublists;
    RangeCounting *case1_pTilde;
    RangeCounting *case2_pTilde;
    PackedIntegers case0_ds;
    vector< PackedIntegers > case00_sublists;
    void constructInCase0();
    void constructInCase00();
    void constructInCase1();
    void constructInCase2();
    void constructSplitingTable();
    void dividePackedIntegersIntoPow2(int, PackedIntegers&, vector<PackedIntegers>&);
    int count16bit(UINT16) const;
public:
    ~RangeCounting(){
        if( case1_pTilde != NULL ){ delete case1_pTilde; }
        if( case2_pTilde != NULL ){ delete case2_pTilde; }
    };
    RangeCounting(vector<int>&);
    RangeCounting(vector<int>&, int);
    RangeCounting(const RangeCounting&);
    int get_integer(int);
    int queryCase0(int);
    int queryCase0(int, int);
    int queryCase00(int);
    int queryCase00(int, int);
    int queryCase1(int);
    int queryCase1(int, int);
    int queryCase2(int);
    int query(int);
    int query(int, int);
};
vector< vector< pair<UINT16, UINT16> > > *RangeCounting::splitingTable = NULL;
RangeCounting::RangeCounting(const RangeCounting& rhs)://{{{
P(rhs.P), case0_ds(rhs.case0_ds),
case00_pTilde(rhs.case00_pTilde), case00_sublists(rhs.case00_sublists),
case1_sublists(rhs.case1_sublists), case1_pTilde(NULL),
case2_sublists(rhs.case2_sublists), case2_pTilde(NULL),
n(rhs.n), l(rhs.l), pow2A(rhs.pow2A), log2n(rhs.log2n), case00_h(rhs.case00_h){
    if( rhs.case1_pTilde != NULL ){
        case1_pTilde = new RangeCounting( *(rhs.case1_pTilde) );
    }
    if( rhs.case2_pTilde != NULL ){
        case2_pTilde = new RangeCounting( *(rhs.case2_pTilde) );
    }
};//}}}
RangeCounting::RangeCounting(vector<int>& _P)://{{{
n((int)_P.size()), pow2A(1 << A), log2n((int)ceil( log2((double)n) )),
case1_pTilde(NULL), case2_pTilde(NULL){
    if( !(n > 0) ){ return; } // exit when given _P contain no integer

    copy(_P.begin(), _P.end(), back_inserter(P)); // retain given P to refer each value as needed

    // l is # of bits that is needed to represent each integer in given P
    l = 1 + MAX(0, (int)floor( log2( (double)( *max_element(_P.begin(), _P.end()) ) ) ) );
    // int maximum = *max_element(_P.begin(), _P.end());
    // l = 1 + MAX(0, (int)floor( log2( (double)( maximum ) ) ) );
    // printf("n:%d, MAX(P[i])=%d, then l=%d\n", n, maximum, l);
    if( (l <= H) and (n <= pow2A) ){ // case: 0
        // cout << "*** construction: case0 ***" << endl;
        constructInCase0();
    }else if( (l <= H) and (n > pow2A) ){ // case: 00
        // cout << "*** construction: case00 ***" << endl;
        constructInCase00();
    }else if( (H < l) and (l <= L) ){ // case: 1
        // cout << "*** construction: case1 ***" << endl;
        constructInCase1();
    }else if( l > L ){ // case: 2
        // cout << "*** construction: case2 ***" << endl;
        constructInCase2();
    }else{
        fprintf(stderr, "not matched for any case\n");
        exit(1);
    }
};//}}}
RangeCounting::RangeCounting(vector<int>& _P, int _l)://{{{
l(_l), n((int)_P.size()), pow2A(1 << A), log2n((int)ceil( log2((double)n) )),
case1_pTilde(NULL), case2_pTilde(NULL){
    if( !(n > 0) ){ return; } // exit when given _P contain no integer

    copy(_P.begin(), _P.end(), back_inserter(P)); // retain given P to refer each value as needed

    // l is # of bits that is needed to represent each integer in given P
    // in this case, use the bigger one of given _l or required # of bits
    int maximum = *max_element(_P.begin(), _P.end());
    // printf("MAX(%d,%d)\n", l, (1 + (int)floor( log2( (double)( maximum ) ) ) ));
    l = MAX( l, (1 + (int)floor( log2( (double)( maximum ) ) ) ) );
    // printf("n:%d, l:%d\n", n, l);
    if( (H < l) and (l <= L) ){ // case: 1
        // cout << "*** construction: case1 ***" << endl;
        constructInCase1();
    }else if( l > L ){ // case: 2
        // cout << "*** construction: case2 ***" << endl;
        constructInCase2();
    }else{
        fprintf(stderr, "not matched for any case\n");
        exit(1);
    }
};//}}}
void RangeCounting::constructSplitingTable(){//{{{
    if( splitingTable != NULL ){ return; }
    splitingTable = new vector< vector< pair<UINT16, UINT16> > >(
                            (1<<WORD_SIZE), vector< pair<UINT16, UINT16> >(
                                        // to simplify, reserve spaces of (L+1)-elements
                                        (L+1), pair<UINT16, UINT16>(0,0)
                                    )
                        );

    // for(int l = L; l > 0; l--){ // l is # of bits that representing each value
    for(int l = L; l > H; l--){ // l is # of bits that representing each value
        UINT16 head_mask = (1 << (l-1)) - 1; // to delete MSB, the size of this mask is (l-1)
        for(int wi = 0, end_wi = (1 << WORD_SIZE); wi < end_wi; wi++){ // use int to avoid overflow of end_wi
            // cout << "-------------" << wi << "--------------" << endl;
            // showOneWord(l, wi);
            UINT16 P0 = 0, P1 = 0;
            for(int i = 0, g = WORD_SIZE/l, num_p0 = 0, num_p1 = 0; i < g; i++){ // g is # of integers in each word
                UINT16 shifted = wi >> (i*l); // shift to set the current integer at the head
                if( shifted & (1 << (l-1)) ){ // check MSB of the current integer
                    P1 |= ( shifted & head_mask ) << ((num_p1++) * (l-1));
                }else{
                    P0 |= ( shifted & head_mask ) << ((num_p0++) * (l-1));
                }
            }
            // showOneWord(l-1, P0);
            // showOneWord(l-1, P1);
            // cout << "---------------------------" << endl;
            (*splitingTable)[wi][l].first = P0;
            (*splitingTable)[wi][l].second = P1;
        }
    }
}//}}}
int RangeCounting::count16bit(UINT16 v) const{//{{{
    unsigned short count = (v & 0x5555) + ((v >> 1) & 0x5555);
    count = (count & 0x3333) + ((count >> 2) & 0x3333);
    count = (count & 0x0f0f) + ((count >> 4) & 0x0f0f);
    return (count & 0x00ff) + ((count >> 8) & 0x00ff);
}//}}}
void RangeCounting::dividePackedIntegersIntoPow2(int h, PackedIntegers& P, vector<PackedIntegers>& res){//{{{
    if( res.size() < (1 << h) ){
        fprintf(stderr, "short of reserved space. need at least %d elements reserved.", (1 << h));
        exit(1);
    }
    vector<PackedIntegers> p_tree(2 * (1 << h)); // the tree has 2^h leaves, then its size is (2 * 2^h - 1)
    p_tree[1] = P;
    for(int i = 1; i <= h; i++){ // reserve P_0, P_1, ... P_(2^h-1)
        for(int j = 0, pow2i = (1 << i); j < pow2i; j++){
            p_tree[pow2i + j] = PackedIntegers(L-i);
        }
    }
    for(int i = 1, end_i = (1 << h); i < end_i; i++){ // create P_0, P_1, ... P_(2^h-1), recursively
        PackedIntegers pi = p_tree[i]; // divide each word in P_i
        for(int j = 0, end_j = pi.words.size(); j < end_j; j++){
            // mask the current word and count one bit in the masked word
            int num_p1 = count16bit( pi.words[j].bits & pi.counting_mask ),
                num_p0 = pi.words[j].n - num_p1;
            p_tree[(i << 1)].append(num_p0, (*splitingTable)[ pi.words[j].bits ][ pi.l ].first);
            p_tree[(i << 1) + 1].append(num_p1, (*splitingTable)[ pi.words[j].bits ][ pi.l ].second);
        }
    }
    // for(int i = 1; i < p_tree.size(); i++){ // debug
    //     cout << i;
    //     p_tree[i].showWords();
    // }
    for(int i = 0, pow2h = (1 << h); i < pow2h; i++){ // copy to vector to be return
        res[i] = p_tree[pow2h + i];
    }
    // for(int i = 0; i < res.size(); i++){ // debug
    //     printf("res[%d] storing %d-integers\n", i, res[i].n);
    //     res[i].showWords();
    // }
};//}}}
void RangeCounting::constructInCase0(){//{{{
    case0_ds = PackedIntegers(H, P);
    case0_ds.constructDataStructure();
};//}}}
void RangeCounting::constructInCase00(){//{{{
    case00_h = log2n - A;
    case00_sublists = vector< PackedIntegers >(); // create data structures of p_0, p_1, ..., p_(2^h-1)
    int st = 0; // split P into 2^h blocks, each of size 2^A. where h = logn - A
    for(int i=0; st + pow2A < n; st += pow2A){
        vector<int> _P(P.begin() + st, P.begin() + st + pow2A);
        PackedIntegers pi(H, _P);
        pi.constructDataStructure();
        case00_sublists.push_back(pi);
    }
    vector<int> _P(P.begin() + st, P.end()); // for the rest of P
    PackedIntegers pi(H, _P);
    pi.constructDataStructure();
    case00_sublists.push_back(pi);

    vector<int> P_last(0); // for the last block as a dummy
    pi = PackedIntegers(H, P_last);
    pi.constructDataStructure();
    case00_sublists.push_back(pi);
    // printf("create sublists\n");

    // printf("h:%d, H:%d, A:%d, n:%d, ceil(log2n):%d\n", case00_h, H, A, n, log2n);
    case00_pTilde = vector< vector<int> >( // construct P~
                // taking care of the case given a greater x (,y) than n (,P[n-1] respectively),
                // reserve spaces in surplus for P~
                (1 << case00_h) + 1, vector<int>((1 << H) + 1, 0)
            );
    // use a temporary array to calc. all query-answers using DynamicProgramming
    vector< vector<int> > count( (1 << case00_h), vector<int>((1 << H), 0) );
    for(int i = 0, end_i = P.size(); i < end_i; i++){ // count points in each block
        count[ (i / pow2A) ][ P[i] ]++;
    }
    // // debug//{{{
    // for(int i=0, end_i=count.size(); i<end_i; i++){
    //     for(int j=0, end_j=count[0].size(); j<end_j; j++){
    //         cout << count[i][j] << " ";
    //     }
    //     cout << endl;
    // }
    // cout << "------------------------------------" << endl;//}}}
    for(int i=0, end_i=count.size(); i<end_i; i++){ // pile up along with y-axis
        for(int j=1, end_j=count[0].size(); j<end_j; j++){
            count[i][j] += count[i][j - 1];
        }
    }
    // // debug//{{{
    // for(int i=0, end_i=count.size(); i<end_i; i++){
    //     for(int j=0, end_j=count[0].size(); j<end_j; j++){
    //         cout << count[i][j] << " ";
    //     }
    //     cout << endl;
    // }
    // cout << "------------------------------------" << endl;//}}}
    for(int i=1, end_i=case00_pTilde.size(); i<end_i; i++){ // store query-answers in case00_ptilde
        for(int j=1, end_j=case00_pTilde[0].size(); j<end_j; j++){
            case00_pTilde[i][j] = case00_pTilde[i-1][j] + count[i-1][j-1];
        }
    }
    // debug//{{{
    // for(int i=0, end_i=case00_pTilde.size(); i<end_i; i++){
    //     printf("%3d | ", i);
    //     for(int j=0, end_j=case00_pTilde[0].size(); j<end_j; j++){
    //         cout << case00_pTilde[i][j] << " ";
    //     }
    //     cout << endl;
    // }
    // cout << "------------------------------------" << endl;//}}}
    // printf("create P~\n");
};//}}}
void RangeCounting::constructInCase1(){//{{{
    PackedIntegers _P(L, P);
    // divide given P into P_0, P_1, ..., P_(2^H-1)
    constructSplitingTable();
    vector<PackedIntegers> P_i( (1 << H) );
    dividePackedIntegersIntoPow2(H, _P, P_i);
    case1_sublists = vector<RangeCounting>();
    for(int i=0, end_i=P_i.size(); i<end_i; i++){
        vector<int> vec;
        P_i[i].sequence(vec);
        // cout << "vec.size(): " << vec.size() << endl;
        case1_sublists.push_back( RangeCounting(vec) );
        // printf("case1_sublists.push_back( RangeCounting(vec) )\n");
    }

    // create P~
    vector<int> vec(n);
    for(int i=0; i<n; i++){
        vec[i] = (int)( P[i] >> (L - H) );
    }
    case1_pTilde = new RangeCounting(vec);
};//}}}
void RangeCounting::constructInCase2(){//{{{
    vector< vector<int> > P_i( (1 << L) );
    vector<int> P_tilde(n);
    UINT16 mask = (1 << (l - L)) - 1;
    for(int j=0; j<n; j++){
        // divide P[j] by 2^(l-L) so that i is in [0, 2^L) and pi_y is in [0, 2^(l-L))]
        int i = (int)( (UINT16)P[j] ) >> (l - L),
            pi_y = (int)( (UINT16)P[j] ) & mask;
        // printf("x:%d, y:%d, i:%d, pi_y:%d\n", j, P[j], i, pi_y);
        P_i.at(i).push_back(pi_y);
        P_tilde.at(j) = i;
    }

    // create a data structure for each P_i
    case2_sublists = vector<RangeCounting>();
    for(int i=0, end_i=P_i.size(); i<end_i; i++){
        case2_sublists.push_back( RangeCounting( P_i[i], L ) );
    }

    // create P~
    case2_pTilde = new RangeCounting(P_tilde);
    // cout << "P~: ";
    // for(int i=0; i<(*case2_pTilde).P.size(); i++){
    //     cout << (*case2_pTilde).P[i] << " ";
    // }
    // cout << endl;
};//}}}
int RangeCounting::queryCase0(int x){//{{{
    // printf("x:%d, P[%d]:%d\n", x, x, P[x]);
    return case0_ds.rc_query(x);
}//}}}
int RangeCounting::queryCase0(int x, int y){//{{{
    // printf("case0\t x:%d, y:%d\n", x, y);
    return case0_ds.rc_query(x, y);
}//}}}
int RangeCounting::queryCase00(int x){//{{{
    // printf("x:%d, P[%d]:%d\n", x, x, P[x]);
    int i = x / pow2A, j = x % pow2A,
        y = case00_sublists.at(i).get_integer(j);
    return case00_pTilde.at(i).at(y) + case00_sublists.at(i).rc_query(j);
}//}}}
int RangeCounting::queryCase00(int x, int y){//{{{
    // printf("case00\t x:%d, y:%d\n", x, y);
    int i = x / pow2A, j = x % pow2A;
    return case00_pTilde.at(i).at(y) + case00_sublists.at(i).rc_query(j, y);
    // debug
    // int ret1 = case00_pTilde.at(i).at(y);
    // int ret2 = case00_sublists.at(i).rc_query(j, y);
    // return ret1 + ret2;
}//}}}
int RangeCounting::queryCase1(int x){//{{{
    int tilde_y = (*case1_pTilde).P[x],
        pi_x = (*case1_pTilde).query(x, tilde_y+1) - (*case1_pTilde).query(x, tilde_y);
    // printf("queryCase1: x:%d, tilde_y:%d\n", x, tilde_y);
    return (*case1_pTilde).query(x) + case1_sublists[tilde_y].query(pi_x);
}//}}}
int RangeCounting::queryCase1(int x, int y){//{{{
    if( !( y < (1 << L) ) ){ return x; }
    // Note that; supposed H < l <= L in case1. In other words, y is in [0, 2^L)
    int tilde_y = ( y >> (L - H) ), // divide y by 2^(L-H), then tilde_y is in [0, 2^H]
        pi_y = (y & ( (1 << H) - 1 )), // the remainder, that corresponds pi_y, is in [0, 2^(L-H))
        pi_x = (*case1_pTilde).query(x, tilde_y+1) - (*case1_pTilde).query(x, tilde_y);
    // int plus1 = (*case1_pTilde).query(x, tilde_y+1);
    // printf(" case1 plus1: %d\n", plus1);
    // int plus0 = (*case1_pTilde).query(x, tilde_y);
    // printf(" case1 plus0: %d\n", plus0);
    // int pi_x = plus1 - plus0;
    return (*case1_pTilde).query(x, tilde_y) + case1_sublists[tilde_y].query(pi_x, pi_y);
}//}}}
int RangeCounting::queryCase2(int x){//{{{
    // printf("----- step in queryCase2: x:%d\n", x);
    int tilde_y = (*case2_pTilde).P.at(x);

    int pi_x = (*case2_pTilde).query(x, tilde_y+1) - (*case2_pTilde).query(x, tilde_y);
    // debug
    // printf("queryCase2: x:%d, tilde_y:%d\n", x, tilde_y);
    // int plus1 = (*case2_pTilde).query(x, tilde_y+1);
    // printf("case2 plus1: %d\n", plus1);
    // int plus0 = (*case2_pTilde).query(x, tilde_y);
    // printf("case2 plus0: %d\n", plus0);
    // int pi_x = plus1 - plus0;

    return (*case2_pTilde).query(x) + case2_sublists[tilde_y].query(pi_x);
    // debug
    // int ret1 = (*case2_pTilde).query(x);
    // printf("(*case2_pTilde).query(%d):%d\n", x, ret1);
    // int ret2 = case2_sublists[tilde_y].query(pi_x);
    // printf("case2_sublists[%d].query(%d):%d\n", tilde_y, pi_x, ret2);
    // printf("----- step out queryCase2: x:%d\n", x);
    // return ret1+ret2;
}//}}}
int RangeCounting::query(int x){//{{{
    if( n == 0 ){
        // cout << "*** query: case(-1) ***" << endl;
        return 0;
    }else if( (l <= H) and (n <= pow2A) ){ // case: 0
        // cout << "*** query: case0 ***" << endl;
        return queryCase0(x);
    }else if( (l <= H) and (n > pow2A) ){ // case: 00
        // cout << "*** query: case00 ***" << endl;
        return queryCase00(x);
    }else if( (H < l) and (l <= L) ){ // case: 1
        // cout << "*** query: case1 ***" << endl;
        return queryCase1(x);
    }else if( l > L ){ // case: 2
        // cout << "*** query: case2 ***" << endl;
        return queryCase2(x);
    }else{
        fprintf(stderr, "not matched for any case\n");
        exit(1);
    }
    return 0;
};//}}}
int RangeCounting::query(int x, int y){//{{{
    // =========================================================
    // TODO:
    // It is needed that dealing with the cases give (x,y)
    // those values are greater than n-1, P[n-1], respectively
    // =========================================================
    if( n == 0 ){
        // cout << "*** construction: case(-1) ***" << endl;
        return 0;
    }else if( (l <= H) and (n <= pow2A) ){ // case: 0
        // cout << "*** query: case0 ***" << endl;
        return queryCase0(x, y);
    }else if( (l <= H) and (n > pow2A) ){ // case: 00
        // cout << "*** query: case00 ***" << endl;
        return queryCase00(x, y);
    }else if( (H < l) and (l <= L) ){ // case: 1
        // cout << "*** query: case1 ***" << endl;
        return queryCase1(x, y);
    }else{
        fprintf(stderr, "not matched for any case\n");
        exit(1);
    }
    return 0;
};//}}}

/* vim:set foldmethod=marker commentstring=//%s : */
