#ifndef OP
#define OP
#include "./order_preserving.h"
#endif

#include<iostream>
#include<string>
#include<vector>
#include<map>
#include<algorithm>
#include<stdlib.h>
#include<string.h>
#include<ctime>

#include "./WaveletTree_double.h"

using namespace std;

void test_for_kgram_kernel(int, int);
int naive_natRepKgram(vector<double>&, vector<double>&, int);

int main(){
    srand(time(0));
    int textsize = 10000, max_k=10;

    for(int N_i=1; N_i<textsize; N_i<<=1){
        int n = N_i + (rand() % N_i); /* create random sequence */
        int m = N_i + (rand() % N_i);
        vector<double> S(n), T(m);

        for(int i=0; i<n; i++){ S[i] = (float)rand() / (float)(RAND_MAX); }
        for(int i=0; i<m; i++){ T[i] = (float)rand() / (float)(RAND_MAX); }

        for(int k=1; k<max_k; k++){
            int naive = naive_natRepKgram(S, T, k);
            int wt = op_kgram(S, T, k);
            printf("S=%d,T=%d,k=%d:\tnaive=%d, wt=%d\n", n, m, k, naive, wt);
            if(naive != wt){
                exit(1);
            }
        }
    }
}

int naive_natRepKgram(vector<double>& S, vector<double>& T, int k){
    vector<double> substring(k);
    vector<int> kgram(k);
    map<double, int> hash;
    natRepKgramVector vec;
    int ret=0;
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

        if( vec.find(kgram) == vec.end() ){ /* regist the kgram */
            vec[kgram] = 1;
        }else{
            vec[kgram]++;
        }
    }

    for(int i=0, end_i=T.size()-k+1; i<end_i; i++){
        for(int j=0; j<k; j++){ /* slice substring */
            substring[j] = T[i+j];
        }

        stable_sort(substring.begin(), substring.end()); /* merge sort */

        hash.clear();
        for(int j=0; j<k; j++){
            hash[ substring[j] ] = j+1; /* hashing val -> order */
        }

        for(int j=0; j<k; j++){
            kgram[j] = hash[ T[i+j] ]; /* create an encoded kgram */
        }

        if( vec.find(kgram) != vec.end() ){ /* regist the kgram */
            ret += vec[kgram];
        }
    }

    return ret;
}

/* vim:set foldmethod=marker commentstring=//%s : */
