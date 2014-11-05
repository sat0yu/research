#ifndef OP
#define OP
#include "./order_preserving.h"
#endif

#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<map>
#include<algorithm>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<ctime>

#include "./WaveletTree_double.h"

using namespace std;

int test_for_bitvector(int);
int test_for_wavelettree(int, int);
void test_for_kgram_kernel(int, int);
int naive_natRepKgram(vector<double>&, vector<double>&, int);
int wavelettree_natRepKgram(vector<double>&, vector<double>&, int);
void naive_natRepKgramVector(vector<int>&, int, natRepKgramVector&);

int main(){
    srand(time(0));
    int bv_textsize = 10000;
    test_for_bitvector(bv_textsize);
    int wt_textsize = 10000, k=10;
    test_for_wavelettree(wt_textsize, k);
    test_for_kgram_kernel(wt_textsize, 10);
}

int test_for_bitvector(int N){//{{{
    for(int N_i=1; N_i<N; N_i<<=1){
        int n = N_i + (rand() % N_i);

        printf("\na test in the condition n=%d starts;\n", n);

        char *B = (char*)malloc((n+1)*sizeof(char));
        for(int i=0; i<n; i++){
            B[i] = (rand() % 2) ? '1' : '0';
        }
        B[n] = '\0';

        clock_t s_time, e_time;
        double duration;
        bool result;

        s_time = clock();
        // <construst an instance>
        BitVector bv = BitVector(B);
        // </construst an instance>
        e_time = clock();
        duration = (double)(e_time - s_time) / (double)CLOCKS_PER_SEC;
        printf("construction: OK\t %f [s]\n", duration);

        s_time = clock();
        // <a test for access>
        result = true;
        for(int i=0; i<strlen(B); ++i){
            if( bv.access(i) != B[i] ){
                printf("bv[%d]=%c, B[%d]=%c\n", i, bv.access(i), i, B[i]);
                result = false;
            }
        }
        if(!result){ exit(1); }
        // </a test for access>
        e_time = clock();
        duration = (double)(e_time - s_time) / (double)CLOCKS_PER_SEC;
        printf("Acsess(i) test: OK\t %f [s]\n", duration);

        s_time = clock();
        // <a test for rank>
        result = true;
        int naive=0;
        for(int i=0; i<=strlen(B); ++i){
            if( i > 0 ){ naive += ( B[i-1] - '0' ) ? 1 : 0; }
            int r = bv.rank1(i);
            if( r != naive ){
                printf("rank1(B,%d)=%d (counted naively:rank1(B,%d)=%d)\n", i, r, i, naive);
                result = false;
            }
        }
        if(!result){ exit(1); }
        // <a test for rank>
        e_time = clock();
        duration = (double)(e_time - s_time) / (double)CLOCKS_PER_SEC;
        printf("Rank(i) test: OK\t %f [s]\n", duration);

        s_time = clock();
        // <a test for select>
        result = true;
        for(int i=0, k=0, s=0; i<strlen(B); ++i){
            if(B[i] - '0'){
                if( i != (s = bv.select(k)) ){
                    printf("select(B,%d)=%d (counted naively:select(B,%d)=%d)\n", k, s, k, i);
                    result = false;
                }
                k++;
            }
        }
        if(!result){ printf("B:%s\n", B); exit(1); }
        if(!result){ exit(1); }
        // </a test for select>
        e_time = clock();
        duration = (double)(e_time - s_time) / (double)CLOCKS_PER_SEC;
        printf("Select(i) test: OK\t %f [s]\n", duration);

        s_time = clock();
        // <a test for set>
        result = true;
        for(int i=0; i<strlen(B); ++i){
            int ir = rand();
            char b = ((ir*ir) & 1) + '0';
            bv.set(i, b);
            if( bv.access(i) != b ){
                printf("bv[%d]=%c, r=%c\n", i, bv.access(i), b);
                result = false;
            };
        }
        if(!result){ exit(1); }
        // </a test for set>
        e_time = clock();
        duration = (double)(e_time - s_time) / (double)CLOCKS_PER_SEC;
        printf("Set(i,b) test: OK\t %f [s]\n", duration);

        free(B);
    }

    return 0;
};//}}}

int test_for_wavelettree(int length, int k){//{{{
    for(int _j=1; _j<length; _j<<=1){
        int j = _j + (rand() % _j);

        vector<double> S(j);
        for(int l=0; l<j; l++){ S[l] = (float)rand() / (float)(RAND_MAX); }
        printf("---------------------------------------------------");
        printf("\na test in the condition |T|=%d, starts;\n", j);

        clock_t s_time, e_time;
        double duration;
        bool result;

        s_time = clock();
        // <construst an instance>
        WaveletTree wt = WaveletTree(S);
        // </construst an instance>
        e_time = clock();
        duration = (double)(e_time - s_time) / (double)CLOCKS_PER_SEC;
        printf("construction: OK \t %f [s]\n", duration);

        s_time = clock();
        // <a test for access>
        result = true;
        for(int i=0, end_i=S.size(); i<end_i; ++i){
            if( S[i] != wt.access(i) ){
                printf("S[%d]=%lf, WT.access[%d]=%lf\n", i, S[i], i, wt.access(i));
                result = false;
            }
        }
        if(!result){ exit(1); }
        // </a test for access>
        e_time = clock();
        duration = (double)(e_time - s_time) / (double)CLOCKS_PER_SEC;
        printf("Acsess(i) test: OK\t %.10lf [s]\n", duration / S.size());
    }

    return 0;
};//}}}

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

void test_for_kgram_kernel(int length, int max_k){//{{{
    for(int N_i=1; N_i<length; N_i<<=1){
        int n = N_i + (rand() % N_i); /* create random sequence */
        int m = N_i + (rand() % N_i);
        vector<double> S(n), T(m);

        for(int i=0; i<n; i++){ S[i] = (float)rand() / (float)(RAND_MAX); }
        for(int i=0; i<m; i++){ T[i] = (float)rand() / (float)(RAND_MAX); }

        for(int k=1; k<max_k; k++){
            int naive = naive_natRepKgram(S, T, k);
            int wt = wavelettree_natRepKgram(S, T, k);
            printf("S=%d,T=%d,k=%d:\tnaive=%d, wt=%d\n", n, m, k, naive, wt);
            if(naive != wt){
                exit(1);
            }
        }
    }
};//}}}

int naive_natRepKgram(vector<double>& S, vector<double>& T, int k){//{{{
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
}//}}}

int wavelettree_natRepKgram(vector<double>& S, vector<double>& T, int k){//{{{
    int ret=0, n = S.size(), m=T.size();

    vector<double> X; /* create a new vector concatenating S and T*/
    copy(S.begin(), S.end(), back_inserter(X));
    X.insert(X.end(), T.begin(), T.end());
    WaveletTree wt(X);
    //printf("S:%d, T:%d, X:%d, k:%d\n", n, m, X.size(), k);

    natRepKgramVector vec; /* use 'map' to represent feature vector */
    wt.createNatRepKgramVector(0, n, k, vec); /* create the feature vector of S */

    map<double, int> hash;
    vector<int> kgram(k);
    for(int i=n, end_i=n+m-k+1; i<end_i; i++){
        hash.clear();
        wt.rangemink_hash(i, i+k, k, hash);
        //printf("i:%d, i+k:%d, hash.size():%d\n", i, i+k, hash.size());
        for(int j=0; j<k; j++){
            kgram[j] = hash[ wt.access(i+j) ]; /* create an encoded kgram */
        }

        if( vec.find(kgram) != vec.end() ){ /* if the kgram exists, then add that count */
            ret += vec[kgram];
        }
    }

    return ret;
}//}}}
/* vim:set foldmethod=marker commentstring=//%s : */
