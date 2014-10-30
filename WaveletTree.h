#include<iostream>
#include<fstream>
#include<string>
#include<vector>
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

#define UB_ALPHABET_SIZE 100000
#define UB_TEXT_SIZE 1000000

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
    size_bits n;
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
    n = N;
    double log2n = log2( (int)n );
    s = (size_bits)MAX(1, ceil(log2n / 2.)); /* s = lg(n)/2 bits in B are covered by S[i] */
    l = (size_bits)MAX(2, (int)(log2n * log2n)); /* l = lg^2(n) bits in B are covered by L[i] */
    while( l % s ){ l++; } /* make l become a multiple of s */
    len_B = (size_array)ceil( (int)n / (double)b ); /* |B[]| */
    len_S = (size_array)ceil( (int)n / (double)s ); /* |S[]| */
    len_L = (size_array)ceil( (int)n / (double)l ) + 1; /* |L[]|, for rank1(n), reserve more by 1 */
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
    for(int i=0, j=0, k=0; i < n; ++i){
        if( !( i % l ) ){
            ++j;
            if( j < len_L ){ L[j] = L[j-1]; };
        }

        if( !( i % s ) ){
            if( !( i % l ) ){ S[k] = 0; } /* if the index of L, j is also incremented */
            ++k;
            if( k < len_S ){ S[k] = S[k-1]; }
        }

        if( bc.access(i) - '0' ){
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
    for(int i=0, j=0, k=0; i < n; ++i){
        if( !( i % l ) ){
            ++j;
            if( j < len_L ){ L[j] = L[j-1]; };
        }

        if( !( i % s ) ){
            if( !( i % l ) ){ S[k] = 0; } /* if the index of L, j is also incremented */
            ++k;
            if( k < len_S ){ S[k] = S[k-1]; }
        }

        if( _B[i] - '0' ){
            /* |-B[i+1]-|-B[ i ]-|
             * |87654321|87654321|*/
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
    if( i == n ){ return L[ len_L - 1 ]; }

    // (1) extract bits in B[i], that are covered by S[i]
    ULL bitseq = 0;
    int st = s * (i / s), en = st + s;
    for(int j=0; j<s; j++){
        if( B[ DIVIDE8(j + st) ] & ( 1 << MOD8(j + st) ) ){ bitseq |= (1 << j); }
    }

    // (2) create a mask with the same width of the bits
    ULL mask = (1 << (i % s)) - 1;
    //ULL mask = 0;
    //for(int j=(i % s)-1; j >= 0; j--){ mask |= (1 << j); }

    return L[ i / l ] + S[ i / s ] + P[ (int)(bitseq & mask) ];
};

int BitVector::rank0(int i){
    return i - rank1(i);
}

int BitVector::select(int i){
    int s=0, e=n-1, m, r;
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
        Node(int);
        void constructBitVector();
    };

    size_t sigma, n;
    vector<int> dict;
    vector<Node> nodes;
    void showTree();
    pair<int, char> traverse_on_alphabet(int, int);
    int traverse_on_wavelet(int, char);
    bool isLeaf(int);
    bool isEmptyNode(int);
    int idx2character(int);
public:
    ~WaveletTree(){};
    WaveletTree(vector<int>&);
    int access(int);
    void rangemink(int, int, int, vector<int>&);
};

void WaveletTree::Node::constructBitVector(){//{{{
    BV = new BitVector(BC);
};
WaveletTree::Node::Node(int n):
    max(0),
    th(0),
    BV(NULL),
    BC(n)
{};//}}}

void WaveletTree::showTree(){//{{{
    cout << endl << "dict[i]: ";
    for(int i=0; i<dict.size(); i++){ cout << dict[i] << ", "; }
    cout << endl;
    size_bits log2sigma = (size_bits)ceil(log2(sigma));
    for(int d=0; d<log2sigma; d++){
        cout << d << "-generation" << endl;
        for(int i=(1 << d); i<(1 << (d+1)) && !isEmptyNode(i); i++){
            for(int j=d; j<log2sigma; j++){ cout << "\t"; }
            printf("node(%d):", i);
            for(int k=0; k<nodes[i].BV->n; k++){
                cout << nodes[i].BV->access(k);
            }
        }
        cout << endl;
    }
    cout << "illustrated." << endl;
};
pair<int, char> WaveletTree::traverse_on_alphabet(int n_idx, int v){
    if(n_idx >= nodes.size()){
        return pair<int, char>(-1, '0');
    }
    int ch_idx;
    char direc;
    if( v <= nodes[n_idx].th ){
        ch_idx = 2 * n_idx;
        direc = '0';
    }else{
        ch_idx = 2 * n_idx + 1;
        direc = '1';
    }
    return pair<int, char>(ch_idx, direc);
};
int WaveletTree::traverse_on_wavelet(int n_idx, char b){
    if( b - '0' > 0 ){
        return (n_idx << 1) + 1;
    }else{
        return n_idx << 1;
    }
};
bool WaveletTree::isLeaf(int n_idx){
    /* internal nodes are stored in first half of nodes[] 
     * while leaves are stored in second half.
     * note: nodes.size() = 2*sigma */
    return (n_idx < sigma) ? false : true;
};
bool WaveletTree::isEmptyNode(int n_idx){
    return (nodes[n_idx].BC.tail_idx > 0) ? false : true;
};
int WaveletTree::idx2character(int n_idx){
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
    vector<char> bucket(UB_ALPHABET_SIZE, 0);
    vector<int>::iterator it_S = _S.begin(), end_it_S = _S.end();
    for(sigma=0; it_S != end_it_S; ++it_S){
        if( *it_S > UB_ALPHABET_SIZE ){
            fprintf(stderr,
                    "too large alphabet size. (the expected maximal size is %d)\n",
                    UB_ALPHABET_SIZE);
            exit(1);
        }
        sigma += ( bucket[*it_S] > 0 ) ? 0 : 1;
        bucket[*it_S] = 1;
    }
    for(int pow=1,tmp=sigma; pow<tmp; pow<<=1){sigma=(pow<<1);} /* make sigma a power of 2 */
    dict.resize(sigma, 0);
    vector<char>::iterator it_b = bucket.begin(), end_it_b = bucket.end();
    for(int d_idx=0; it_b != end_it_b; ++it_b){
        if( *it_b > 0 ){
            dict[d_idx++] = (int)distance(bucket.begin(), it_b);
        }
    }
    //printf("given characters are sorted.\n");

    //----- construct alphabet tree -----
    nodes.resize(2*sigma, n); /* using 1-origin indices */
    for(int i=0; i<sigma; i++){
        /* later half of nodes[] corresspond to leaves*/
        nodes[i+sigma].max = dict[i];
        nodes[i+sigma].th = dict[i];
    }
    for(int i=sigma-1; i>0; i--){
        /* create non-leaf nodes in a way like merge sort */
        int l_child = nodes[2*i].max, r_child = nodes[2*i+1].max;
        /* if the max of child is 0(empty), then use left child */
        nodes[i].max = (r_child == 0) ? l_child : r_child;
        nodes[i].th = l_child;
    }
    //printf("an alphabet tree is constructed.\n");

    //----- construct wavelet tree -----
    for(int i=0, end_i=_S.size(); i<end_i; ++i){ /* regist characters, one by one */
        for(int n_idx=1;;){ /* the index of root node is 1 */
            pair<int, char> ret = traverse_on_alphabet(n_idx, _S[i]);
            if(ret.first < 0){ break; }
            nodes[n_idx].BC.append(ret.second);
            n_idx = ret.first;
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

    //----- show tree -----
    //showTree();
};//}}}

int WaveletTree::access(int i){//{{{
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

class wt_queue_elem{
public:
    int n_idx, st, en;
    wt_queue_elem(int _n_idx, int _st, int _en):
        n_idx(_n_idx), st(_st), en(_en){}
};
struct rangemaxk_comparison{
public:
    bool operator () (wt_queue_elem* p1, wt_queue_elem* p2){
        return p1->n_idx > p2->n_idx;
    };
};
void WaveletTree::rangemink(int s, int e, int k, vector<int>& res){
    priority_queue<wt_queue_elem*, vector<wt_queue_elem* >, rangemaxk_comparison> que;
    que.push( new wt_queue_elem(1, s, e) );
    for(int i=0; !que.empty() && i<k;){
        const wt_queue_elem* elem = que.top();
        int n_idx   = elem->n_idx,
            st      = elem->st,
            en      = elem->en;
        if( isLeaf(n_idx) ){
            for(int j=0, end_j=en-st, ch=idx2character(n_idx); j<end_j; ++j){
                res[i++] = ch;
            }
        }else{
            int ost = nodes[n_idx].BV->rank1(st),
                oen = nodes[n_idx].BV->rank1(en),
                zst = st - ost,
                zen = en - oen;
            if( oen - ost ){
                que.push( new wt_queue_elem((n_idx<<1)+1, ost, oen) );
            }
            if( zen - zst ){
                que.push( new wt_queue_elem((n_idx<<1), zst, zen) );
            }
        }
        que.pop();
    }
};

/* vim:set foldmethod=marker commentstring=//%s : */
