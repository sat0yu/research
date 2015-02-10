#ifndef OP
#define OP
#include "./order_preserving.h"
#endif

#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<map>
#include<queue>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<ctime>

#define TYPE_TO_BITS(X) (8 * sizeof(X))
#define DIVIDE8(X) ((X) >> 3)
#define MOD8(X) ((X) & 0x07)
#define LSB(X) ((X) & 0x01)
#define MAX(X,Y) ((X > Y) ? X : Y)

using namespace std;

typedef unsigned char block_b; /* DO NOT CHANGE THE SIZE OF B[i] */
typedef unsigned short block_s;
typedef unsigned int block_l;
typedef unsigned char block_p;

typedef unsigned char UCHAR;
typedef unsigned char UINT;
typedef unsigned short USHORT;
typedef unsigned long long ULL;

typedef size_t size_array;
typedef size_t size_bits;

class BitContainer{
private:
    size_array len_B;
    size_bits bits_B;
    vector<block_b> B;
public:
    size_bits n, tail_idx;
    ~BitContainer(){};
    BitContainer();
    BitContainer(int);
    const char access(int);
    block_s set(int, char);
    block_s append(char);
};

BitContainer::BitContainer()://{{{
    n(UB_TEXT_SIZE),
    tail_idx(0),
    bits_B(TYPE_TO_BITS(block_b))
{
    len_B = (size_array)ceil( (int)n / (double)bits_B );
    B.resize( len_B, 0 );
};

BitContainer::BitContainer(int _n):
    n(_n),
    tail_idx(0),
    bits_B(TYPE_TO_BITS(block_b))
{
    len_B = (size_array)ceil( (int)n / (double)bits_B );
    B.resize( len_B, 0 );
};//}}}

const char BitContainer::access(int i){//{{{
    return (B[ DIVIDE8(i) ] & ( 1 << MOD8(i) )) ? '1' : '0';
};//}}}

block_s BitContainer::set(int i, char b){//{{{
    if( i > tail_idx ){ tail_idx = i; }
    if( b - '0' > 0 ){
        return B[ DIVIDE8(i) ] |= 1 << MOD8(i);
    }else{
        return B[ DIVIDE8(i) ] &= ~( 1 << MOD8(i) );
    }
};//}}}

block_s BitContainer::append(char b){//{{{
    if( b - '0' > 0 ){
        B[ DIVIDE8(tail_idx) ] |= 1 << MOD8(tail_idx);
    }else{
        B[ DIVIDE8(tail_idx) ] &= ~( 1 << MOD8(tail_idx) );
    }
    return B[ DIVIDE8(tail_idx++) ];
};//}}}

class BitVector{
private:
    size_array len_S, len_L, len_P, len_B;
    size_bits bits_B, bits_P, bits_S, bits_L, b, s, l;
    vector<block_b> B; /* in which a bit-vector is stored */
    vector<block_s> S; /* S: small blocks */
    vector<block_l> L; /* L: large blocks */
    vector<block_p> P; /* P: a dictionary of # of 1 in each bit-vector */
    void initialize(size_bits);
public:
    size_bits n, _n;
    ~BitVector(){};
    BitVector(const BitVector&);
    BitVector(BitContainer&);
    BitVector(const char*);
    const char access(int);
    int rank1(int);
    int rank0(int);
    int select(int);
    block_s set(int, char);
};

void BitVector::initialize(size_bits N){//{{{
    if( N == 0 ){
        printf("BitVector initialization error\n");
        exit(1);
    }
    _n = N; /* store original bit length for binary-search */
    for(int pow = 1, tmp = _n; /* make 'n' a power of 2 */
            pow <= tmp;
            pow <<= 1)
    { n = ( pow << 1 ); }
    double log2n = log2( (int)n );
    s = (size_bits)MAX(1, ceil(log2n / 2.)); /* s = lg(n)/2 bits in B are covered by S[i] */
    l = (size_bits)MAX(2, (int)(log2n * log2n)); /* l = lg^2(n) bits in B are covered by L[i] */
    while( l % s ){ l++; } /* make l become a multiple of s */
    len_B = (size_array)ceil( (int)n / (double)b ); /* |B[]| */
    len_S = (size_array)ceil( (int)n / (double)s ); /* |S[]| */
    len_L = (size_array)ceil( (int)n / (double)l ); /* |L[]|, for rank1(n), reserve more by 1 */
    len_P = (size_array)(1 << s); /* |P[0, 2^s)| */

    //----- check if each specified size is enough -----
    size_bits required_bits_p, required_bits_s, required_bits_l;
    required_bits_p = (size_bits)ceil( log2((int)s) );
    required_bits_s = (size_bits)ceil( log2((int)l) );
    required_bits_l = (size_bits)ceil( log2((int)n) );
    if( bits_P < required_bits_p ){
        printf("exceed # of bits of P[i], P[i] needs %d bits\n", (int)required_bits_p );
        exit(1);
    }
    if( bits_S < required_bits_s ){
        printf("exceed # of bits of S[i], S[i] needs %d bits\n", (int)required_bits_s );
        exit(1);
    }
    if( bits_L < required_bits_l ){
        printf("exceed # of bits of L[i], L[i] needs %d bits\n", (int)required_bits_l );
        exit(1);
    }

    //----- show status -----
//    printf("n: %d\n", (int)n);
//    printf("# of bits in B, that are covered by L[i]:%d\n", (int)l);
//    printf("# of bits in B, that are covered by S[i]:%d\n", (int)s);
//    printf("len(B): %d, # of bits of B[i]: %d\n", (int)len_B, (int)b);
//    printf("len(S): %d, # of bits of S[i]: %d (%d bits are required)\n", \
//            (int)len_S, (int)bits_S, (int)required_bits_s );
//    printf("len(L): %d, # of bits of L[i]: %d (%d bits are required)\n", \
//            (int)len_L, (int)bits_L, (int)required_bits_l );
//    printf("len(P): %d  # of bits of P[i]: %d (%d bits are required)\n", \
//            (int)len_P, (int)bits_P, (int)required_bits_p);

    //----- reserve mamories for B, L, S, and P -----
    B.resize( len_B, 0 );
    L.resize( len_L, 0 );
    S.resize( len_S, 0 );
    P.resize( len_P, 0 );

    //----- initialize P[0,i) -----
    for(int i=0; i < len_P; ++i){
        P[i] = 0;
        for(int j=i; j > 0; j >>= 1){ /* check bits from right to left */
            if( LSB(j) ){ ++P[i]; }
        }
    }
};//}}}

BitVector::BitVector(const BitVector &bv)://{{{
    bits_B(TYPE_TO_BITS(block_b)),
    b(TYPE_TO_BITS(block_b)),
    bits_P(TYPE_TO_BITS(block_p)),
    bits_S(TYPE_TO_BITS(block_s)),
    bits_L(TYPE_TO_BITS(block_l))
{
    initialize(bv.n);
}

BitVector::BitVector(BitContainer &bc):
    bits_B(TYPE_TO_BITS(block_b)),
    b(TYPE_TO_BITS(block_b)),
    bits_P(TYPE_TO_BITS(block_p)),
    bits_S(TYPE_TO_BITS(block_s)),
    bits_L(TYPE_TO_BITS(block_l))
{
    //----- initialization -----
    initialize( bc.tail_idx );

    //----- create bit vector -----
    L[0] = S[0] = 0;
    for(int i=0, j=0, k=0, tail=bc.tail_idx; i < n; ++i){
        if( !( i % l ) ){
            ++j;
            if( j < len_L ){ L[j] = L[j-1]; };
        }

        if( !( i % s ) ){
            if( !( i % l ) ){ S[k] = 0; } /* if the index of L, j is also incremented */
            ++k;
            if( k < len_S ){ S[k] = S[k-1]; }
        }

        if( i < tail && bc.access(i) - '0' ){
            B[ DIVIDE8(i) ] |= 1 << MOD8(i); /* regist bits from right to left */
            if( j < len_L ){ ++L[j]; }
            if( k < len_S ){ ++S[k]; }
        }
    }
};

BitVector::BitVector(const char* _B):
    bits_B(TYPE_TO_BITS(block_b)),
    b(TYPE_TO_BITS(block_b)),
    bits_P(TYPE_TO_BITS(block_p)),
    bits_S(TYPE_TO_BITS(block_s)),
    bits_L(TYPE_TO_BITS(block_l))
{
    //----- initialization -----
    initialize( (size_bits)strlen(_B) );

    //----- create bit vector -----
    L[0] = S[0] = 0;
    /* i: the index of _B ('i' is NOT the position of a bit),
     * j: the index of L,
     * k: the index of S */
    for(int i=0, j=0, k=0, _strlen=strlen(_B); i < n; ++i){
        if( !( i % l ) ){
            ++j;
            if( j < len_L ){ L[j] = L[j-1]; };
        }

        if( !( i % s ) ){
            if( !( i % l ) ){ S[k] = 0; } /* if the index of L, j is also incremented */
            ++k;
            if( k < len_S ){ S[k] = S[k-1]; }
        }

        if( i < _strlen && _B[i] - '0' ){
            B[ DIVIDE8(i) ] |= 1 << MOD8(i); /* regist bits from right to left */
            if( j < len_L ){ ++L[j]; }
            if( k < len_S ){ ++S[k]; }
        }
    }
};//}}}

const char BitVector::access(int i){//{{{
    return (B[ DIVIDE8(i) ] & ( 1 << MOD8(i) )) ? '1' : '0';
};

int BitVector::rank1(int i){
    // (1) extract bits in B[i], that are covered by S[i]
    ULL bitseq = 0;
    int st = s * (i / s);
    for(int j=0; j<s; j++){
        if( B[ DIVIDE8(j + st) ] & ( 1 << MOD8(j + st) ) ){ bitseq |= (1 << j); }
    }

    // (2) create a mask with the same width of the bits
    //ULL mask = (1 << (i % s)) - 1;
    //ULL mask = 0;
    //for(int j=(i % s)-1; j >= 0; j--){ mask |= (1 << j); }

    return L[ i / l ] + S[ i / s ] + P[ (int)(bitseq & (1 << (i % s)) - 1) ];
};

int BitVector::rank0(int i){
    return i - rank1(i);
}

int BitVector::select(int i){
    int s=0, e=_n, m, r;
    while(s != e){
        m = (s+e)/2;
        r = rank1(m+1);
        if( i < r ){
            e = m;
        }else{
            s = m + 1;
        }
    }

    return s;
}

block_s BitVector::set(int i, char b){
    if( b - '0' > 0 ){
        return B[ DIVIDE8(i) ] |= 1 << MOD8(i);
    }else{
        return B[ DIVIDE8(i) ] &= ~( 1 << MOD8(i) );
    }
};//}}}

class WaveletTree{
private:
    class Node{
    private:
    public:
        int max, th, n;
        BitVector *BV;
        BitContainer BC;
        ~Node();
        Node(int);
        void constructBitVector();
    };

    size_t sigma, n, num_nodes, log2sigma;
    vector<int> dict;
    map<int, int> inverse_dict;
    vector<int> pos_st, pos_en;
    vector<Node> nodes;
    void showTree() const;
    int traverse_on_wavelet(int, char) const;
    bool isLeaf(int) const;
    bool isEmptyNode(int) const;
    int idx2character(int) const;
public:
    ~WaveletTree(){};
    WaveletTree(vector<int>&);
    int access(int) const;
    int rank(int, int) const;
    void rangemink(int, int, int, vector<int>&);
    void rangemink_hash(int, int, int, map<int,int>&);
    void createNatRepKgramVector(vector<int>&, int, natRepKgramVector&);
    int rankLessThan(int, int, int) const;
    int rankLessThan_forany(int, int, int) const;
    void rankLessThanEqual(int, int, int, int*, int*) const;
    int rangefreq(int, int, int, int) const;
    void createCountingCodingKgramVector(vector<int>&, int, countingCodingKgramVector&);
    double createCountingCodingKgramVectorWithSliding(vector<int>&, int, countingCodingKgramVector&);
};

void WaveletTree::Node::constructBitVector(){//{{{
    BV = new BitVector(BC);
};
WaveletTree::Node::~Node(){
    delete BV;
};
WaveletTree::Node::Node(int n):
    max(0),
    th(0),
    BV(NULL),
    BC(n)
{};//}}}

void WaveletTree::showTree() const{//{{{
    cout << endl << "dict[i]: ";
    for(int i=0; i<dict.size(); i++){ cout << dict[i] << ", "; }
    cout << endl;
    for(int d=0; d<log2sigma; d++){
        cout << d << "-generation" << endl;
        for(int j=(1 << d); j<(1 << (d+1)) && !isEmptyNode(j); j++){
            for(int k=d; k<log2sigma; k++){ cout << "\t"; }
            printf("node(%d):", j);
            for(int k=0; k<nodes[j].BV->n; k++){
                cout << nodes[j].BV->access(k);
            }
        }
        cout << endl;
    }
    cout << "illustrated." << endl;
};//}}}

int WaveletTree::traverse_on_wavelet(int n_idx, char b) const{//{{{
    if( b - '0' > 0 ){
        return (n_idx << 1) + 1;
    }else{
        return n_idx << 1;
    }
};//}}}

bool WaveletTree::isLeaf(int n_idx) const{//{{{
    /* internal nodes are stored in first half of nodes[] 
     * while leaves are stored in second half.
     * note: nodes.size() = 2*sigma */
    return (n_idx < sigma) ? false : true;
};//}}}

bool WaveletTree::isEmptyNode(int n_idx) const{//{{{
    return (nodes[n_idx].BC.tail_idx > 0) ? false : true;
};//}}}

int WaveletTree::idx2character(int n_idx) const{//{{{
    /* leaves are stored in latter half of nodes[],
     * and also stored in dict[0:sigma-1] */
    if( n_idx < sigma ){
        fprintf(stderr,
                "nodes[%d] is not a leaf node\n",
                n_idx);
        exit(1);
    }else{
        return dict[n_idx - sigma];
    }
};//}}}

WaveletTree::WaveletTree(vector<int>& _S){//{{{
    n = _S.size();
    if( n > UB_TEXT_SIZE ){
        fprintf(stderr,
                "too large text. (the expected maximal size is %d)\n",
                UB_TEXT_SIZE);
        exit(1);
    }

    //----- sort characters by bucket-sort -----
    vector<bool> bucket(UB_ALPHABET_SIZE, false);
    vector<int>::iterator it_S = _S.begin(), end_it_S = _S.end();
    for(sigma=0; it_S != end_it_S; ++it_S){
        if( *it_S > UB_ALPHABET_SIZE ){
            fprintf(stderr,
                    "too large alphabet size. (the expected maximal size is %d)\n",
                    UB_ALPHABET_SIZE);
            exit(1);
        }
        sigma += bucket[*it_S] ? 0 : 1;
        bucket[*it_S] = true;
    }
    //printf("given characters are sorted.\n");

    //----- create dictionaries that convert character <-> index -----
    for(log2sigma = 1; (1 << log2sigma) < sigma; ++log2sigma){}
    sigma = ( 1 << (log2sigma) ); /* make sigma a power of 2 */
    dict.resize(sigma, 0);
    vector<bool>::iterator it_b = bucket.begin(), end_it_b = bucket.end();
    for(int d_idx=0; it_b != end_it_b; ++it_b){
        if( *it_b ){
            int character = (int)distance(bucket.begin(), it_b);
            inverse_dict[character] = d_idx;
            dict[d_idx] = character;
            d_idx++;
        }
    }
    //printf("dectionaries are created.\n");

    //----- construct wavelet tree -----
    num_nodes = (sigma << 1); /* using 1-origin indices */
    nodes.resize(num_nodes, n);
    for(int i=0, end_i=_S.size(); i<end_i; ++i){ /* regist characters, one by one */
        int path = inverse_dict[ _S[i] ];
        for(int n_idx=1, d=log2sigma-1; d >= 0; d--){
            if( path & ( 1 << d ) ){
                nodes[n_idx].BC.append('1');
                n_idx = (n_idx << 1) + 1;
            }else{
                nodes[n_idx].BC.append('0');
                n_idx = (n_idx << 1);
            }
        }
    }
    //printf("a wavelet tree is constructed.\n");

    //----- convert BitContainer to BitVector -----
    vector<Node>::iterator it_n=nodes.begin(), end_it_n=nodes.end();
    for(; it_n != end_it_n; ++it_n){
        if( (*it_n).BC.tail_idx > 0 ){
            (*it_n).constructBitVector();
        }
    }
    //printf("BitContainers are converted to BitVectors.\n");

    //----- reserve two arrays for tree search -----
    pos_st.resize(num_nodes);
    pos_en.resize(num_nodes);

    //----- show tree -----
    // showTree();
};//}}}

int WaveletTree::access(int i) const{//{{{
    int c=0;
    char b;
    for(
            /* the index of root node is 1 */
            int n_idx=1;
            /* an elem in nodes[], whose index is larger than sigma, is corresspond a leaf*/
            n_idx < sigma;
            n_idx = traverse_on_wavelet(n_idx, b)
        ){
        b = nodes[n_idx].BV->access(i);
        c = (c << 1) | (b - '0');
        if(b - '0' > 0){
            i = nodes[n_idx].BV->rank1(i);
        }else{
            i = nodes[n_idx].BV->rank0(i);
        }
    }
    return dict[c];
};//}}}

int WaveletTree::rank(int c, int i) const{//{{{
    int rank = i, path = inverse_dict.at(c);
    for(int d=log2sigma-1, n_idx=1; d >= 0; d--){
        if( path & ( 1 << d ) ){
            rank = nodes[n_idx].BV->rank1(rank);
            n_idx = (n_idx << 1) + 1;
        }else{
            rank = nodes[n_idx].BV->rank0(rank);
            n_idx = (n_idx << 1);
        }
    }
    return rank;
};//}}}

void WaveletTree::rangemink(int st, int en, int k, vector<int>& res){//{{{
    priority_queue<int, vector<int>, greater<int> > que;
    que.push(1);
    pos_st[1] = st;
    pos_en[1] = en;
    for(int i=0; !que.empty() && i<k;){
        const int n_idx = que.top();
        if( isLeaf(n_idx) ){
            for(int j=0, end_j=(pos_en[n_idx] - pos_st[n_idx]), ch=idx2character(n_idx); j<end_j; ++j){
                res[i++] = ch;
            }
        }else{
            int ost = nodes[n_idx].BV->rank1(pos_st[n_idx]),
                oen = nodes[n_idx].BV->rank1(pos_en[n_idx]),
                zst = pos_st[n_idx] - ost,
                zen = pos_en[n_idx] - oen;
            if( oen - ost ){
                int r_child = (n_idx<<1)+1;
                que.push(r_child);
                pos_st[r_child] = ost;
                pos_en[r_child] = oen;
            }
            if( zen - zst ){
                int l_child = (n_idx<<1);
                que.push(l_child);
                pos_st[l_child] = zst;
                pos_en[l_child] = zen;
            }
        }
        que.pop();
    }
};//}}}

void WaveletTree::rangemink_hash(int st, int en, int k, map<int, int>& hash){//{{{
    priority_queue<int, vector<int>, greater<int> > que;
    que.push(1);
    pos_st[1] = st;
    pos_en[1] = en;
    for(int i=0; !que.empty() && i<k;){
        const int n_idx = que.top();
        if( isLeaf(n_idx) ){
            i = hash[idx2character(n_idx)] = (i + pos_en[n_idx] - pos_st[n_idx]);
        }else{
            int ost = nodes[n_idx].BV->rank1(pos_st[n_idx]),
                oen = nodes[n_idx].BV->rank1(pos_en[n_idx]),
                zst = pos_st[n_idx] - ost,
                zen = pos_en[n_idx] - oen;
            if( oen - ost ){
                int r_child = (n_idx<<1)+1;
                que.push(r_child);
                pos_st[r_child] = ost;
                pos_en[r_child] = oen;
            }
            if( zen - zst ){
                int l_child = (n_idx<<1);
                que.push(l_child);
                pos_st[l_child] = zst;
                pos_en[l_child] = zen;
            }
        }
        que.pop();
    }
};//}}}

void WaveletTree::createNatRepKgramVector(vector<int>& S, int k, natRepKgramVector& res){//{{{
    vector<int> kgram(k);
    map<int, int> hash;
    for(int i=0, end_i=n-k+1; i<end_i; i++){
        hash.clear();
        rangemink_hash(i, i+k, k, hash); /* hashing val -> order */

        for(int j=0; j<k; j++){
            kgram[j] = hash[ S[i+j] ]; /* create an encoded kgram */
        }

        if( res.find(kgram) == res.end() ){ /* regist the kgram */
            res[kgram] = 1;
        }else{
            res[kgram]++;
        }
    }
};//}}}

int WaveletTree::rankLessThan(int c, int st, int en) const{//{{{
    int rank=0, path=inverse_dict.at(c), ost, oen;
    for(int d=log2sigma-1, n_idx=1; st != en && d >= 0; d--){
        ost = nodes[n_idx].BV->rank1(st),
        oen = nodes[n_idx].BV->rank1(en);
        if( path & ( 1 << d ) ){
            rank += en - st - oen + ost;
            st = ost;
            en = oen;
            n_idx <<= 1;
            n_idx++;
        }else{
            st = st - ost;
            en = en - oen;
            n_idx <<= 1;
        }
    }
    return rank;
};//}}}

void WaveletTree::rankLessThanEqual(int c, int st, int en, int* lt, int* eq) const{//{{{
    int rank=0, path=inverse_dict.at(c), ost, oen;
    for(int d=log2sigma-1, n_idx=1; st != en && d >= 0; d--){
        ost = nodes[n_idx].BV->rank1(st),
        oen = nodes[n_idx].BV->rank1(en);
        if( path & ( 1 << d ) ){
            rank += en - st - oen + ost;
            st = ost;
            en = oen;
            n_idx <<= 1;
            n_idx++;
        }else{
            st = st - ost;
            en = en - oen;
            n_idx <<= 1;
        }
    }
    *lt = rank;
    *eq = en - st;
};//}}}

int WaveletTree::rangefreq(int st, int en, int x, int y) const{//{{{
    if( !( (x < y) && (st < en) ) ){ return 0; }
    return rankLessThan(y, st, en) - rankLessThan(x, st, en);
};//}}}

void WaveletTree::createCountingCodingKgramVector(vector<int>& S, int k, countingCodingKgramVector& res){//{{{
    vector<c_code> kgram(k);
    for(int i=0, end_i=n-k+1; i<end_i; i++){
        for(int j=0, lt, eq; j<k; ++j){
            rankLessThanEqual(S[i+j], i, i+j, &lt, &eq);
            kgram[j] = c_code(lt, eq);
        }

        if( res.find(kgram) == res.end() ){ /* regist the kgram */
            res[kgram] = 1;
        }else{
            res[kgram]++;
        }
    }
};//}}}

double WaveletTree::createCountingCodingKgramVectorWithSliding(vector<int>& S, int k, countingCodingKgramVector& res){//{{{
    double insertion_duration = 0.;
    clock_t s_time;

    vector<c_code> kgram(k);
    for( int j=0, lt, eq; j<k; ++j ){ // the first kgram is calcucated naively
        rankLessThanEqual(S[j], 0, j, &lt, &eq);
        kgram[j] = c_code(lt, eq);
    }
    res[kgram] = 1;

    for(int i=1, end_i=n-k+1, lt, eq; i<end_i; i++){
        for(int j=0; j<k-1; ++j){ // utilize the past-head value
            if(S[i-1] < S[i+j]){
                kgram[j+1].first--;
            }else if(S[i-1] == S[i+j]){
                kgram[j+1].second--;
            }
            kgram[j] = kgram[j+1];
        }

        rankLessThanEqual(S[i+k-1], i, i+k-1, &lt, &eq); // the new-tail value is given by WT query
        kgram[k-1] = c_code(lt, eq);

        s_time = clock();
        if( res.find(kgram) == res.end() ){ /* regist the kgram */
            res[kgram] = 1;
        }else{
            res[kgram]++;
        }
        insertion_duration += (clock() - s_time);
    }

    return insertion_duration;
};//}}}

/* vim:set foldmethod=marker commentstring=//%s : */
