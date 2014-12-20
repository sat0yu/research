#include<map>
#include<vector>

using namespace std;

#define UB_ALPHABET_SIZE 100000
#define UB_TEXT_SIZE 1000000

typedef map< vector<int>, int> natRepKgramVector;
typedef pair<int, int> c_code;
typedef pair<int, int> sc_code;
typedef map< vector<c_code>, int> countingCodingKgramVector;
typedef map< vector<c_code>, int> suffixCountingCodingKgramVector;

void naive_natRepKgramVector(vector<int>&, int, natRepKgramVector&);
void naive_countingCodingKgramVector(vector<int>&, int, countingCodingKgramVector&);
void naive_countingCodingKgramVectorWithSliding(vector<int>&, int, countingCodingKgramVector&);
void kgramVector_SuffixCountingCoding(vector<int>&, int, suffixCountingCodingKgramVector&);
void kgramVector_SuffixCountingCodingAndWindowSliding(vector<int>&, int, suffixCountingCodingKgramVector&);

void naive_natRepKgramVector(vector<int>& S, int k, natRepKgramVector& res){//{{{
    vector<int> substring(k), kgram(k);
    map<int, int> hash;
    for(int i=0, end_i=S.size()-k+1; i<end_i; i++){
        for(int j=0; j<k; j++){ /* slice substring */
            substring[j] = S[i+j];
        }

        stable_sort(substring.begin(), substring.end()); /* merge sort */

        hash.clear();
        for(int j=0; j<k; j++){
            hash[ substring[j] ] = j+1; /* hashing val -> order */
        }

        for(int j=0; j<k; j++){
            kgram[j] = hash[ S[i+j] ]; /* create an encoded kgram */
        }

        if( res.find(kgram) == res.end() ){ /* regist the kgram */
            res[kgram] = 1;
        }else{
            res[kgram]++;
        }
    }
}//}}}

void naive_countingCodingKgramVector(vector<int>& S, int k, countingCodingKgramVector& res){//{{{
    vector<c_code> kgram(k);
    for(int i=0, end_i=S.size()-k+1; i<end_i; i++){
        for(int j=0; j<k; j++){ /* for each a substring of length k */
            int lt=0, eq=0;
            for(int l=j-1; l>=0; l--){ /* suffix counting coding */
                if(S[i+l] == S[i+j]){
                    eq++;
                }else if(S[i+l] < S[i+j]){
                    lt++;
                }
            }
            kgram[j] = c_code(lt, eq);
        }

        if( res.find(kgram) == res.end() ){ /* regist the kgram */
            res[kgram] = 1;
        }else{
            res[kgram]++;
        }
    }
}//}}}

void naive_countingCodingKgramVectorWithSliding(vector<int>& S, int k, countingCodingKgramVector& res){//{{{
    vector<c_code> kgram(k);
    for( int j=0; j<k; ++j ){ // the first kgram is calcucated naively
        int lt=0, eq=0;
        for(int l=j-1; l>=0; l--){ /* counting coding */
            if(S[l] < S[j]){
                lt++;
            }else if(S[l] == S[j]){
                eq++;
            }
        }
        kgram[j] = c_code(lt, eq);
    }
    res[kgram] = 1;

    for(int i=1, end_i=S.size()-k+1; i<end_i; i++){
        for(int j=0; j<k-1; j++){
            if(S[i-1] < S[i+j]){ // utilize the past-head value
                kgram[j+1].first--;
            }else if(S[i-1] == S[i+j]){
                kgram[j+1].second--;
            }
            kgram[j] = kgram[j+1];
        }

        int lt=0, eq=0; // the new-tail value is calcucated naively
        for(int l=k-2; l>=0; l--){ /* counting coding */
            if(S[i+l] < S[i+k-1]){
                lt++;
            }else if(S[i+l] == S[i+k-1]){
                eq++;
            }
        }
        kgram[k-1] = c_code(lt, eq);

        if( res.find(kgram) == res.end() ){ /* regist the kgram */
            res[kgram] = 1;
        }else{
            res[kgram]++;
        }
    }
}//}}}

void kgramVector_SuffixCountingCoding(vector<int>& S, int k, suffixCountingCodingKgramVector& res){//{{{
    vector<sc_code> kgram(k);
    for(int i=0, end_i=S.size()-k+1; i<end_i; i++){
        for(int j=0; j<k; j++){ /* for each a substring of length k */
            int lt=0, eq=0;
            for(int l=j+1; l<k; l++){ /* suffix counting coding */
                if(S[i+j] > S[i+l]){
                    lt++;
                }else if(S[i+j] == S[i+l]){
                    eq++;
                }
            }
            kgram[j] = sc_code(lt, eq);
        }

        if( res.find(kgram) == res.end() ){ /* regist the kgram */
            res[kgram] = 1;
        }else{
            res[kgram]++;
        }
    }
}//}}}

void kgramVector_SuffixCountingCodingAndWindowSliding(vector<int>& S, int k, suffixCountingCodingKgramVector& res){//{{{
    vector<sc_code> kgram(k);
    for( int j=0; j<k; ++j ){ // the first kgram code(S[0:k-1]) is calcucated naively
        int lt=0, eq=0;
        for(int l=j+1; l<k; l++){ /* suffix counting coding */
            if(S[j] > S[l]){
                lt++;
            }else if(S[j] == S[l]){
                eq++;
            }
        }
        kgram[j] = sc_code(lt, eq);
    }
    kgram[k-1] = sc_code(0, 0); // the tail of a kgram using SuffixCountingCoding is always (0,0)
    res[kgram] = 1;

    for(int i=1, end_i=S.size()-k+1; i<end_i; i++){ // calc. code(S[i:i+k-1]) for each i
        for(int j=0, tail_idx=i+k-1; j<k-1; j++){
            kgram[j] = kgram[j+1];
            if(S[i+j] > S[tail_idx]){ // utilize the tail value
                kgram[j].first++;
            }else if(S[i+j] == S[tail_idx]){
                kgram[j].second++;
            }
        }
        kgram[k-1] = sc_code(0, 0); // the tail of a kgram using SuffixCountingCoding is always (0,0)

        if( res.find(kgram) == res.end() ){ /* regist the kgram */
            res[kgram] = 1;
        }else{
            res[kgram]++;
        }
    }
}//}}}

/* vim:set foldmethod=marker commentstring=//%s : */
