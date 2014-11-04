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

#include "./WaveletTree.h"

using namespace std;

int test_for_bitvector(int);
int test_for_wavelettree(int, int, int);
int comparison_kgram_vector_construct(int, int, int);
void test_for_kgram_kernel(int, int, int);
int naive_natRepKgram(vector<double>&, vector<double>&, int);
void naive_natRepKgramVector(vector<int>&, int, natRepKgramVector&);

int main(){
    srand(time(0));
    int bv_textsize = 10000;
    test_for_bitvector(bv_textsize);
    int wt_textsize = 1000, wt_alphabetsize = 1000, k=10;
    test_for_wavelettree(wt_textsize, wt_alphabetsize, k);
    comparison_kgram_vector_construct(wt_textsize, wt_alphabetsize, k);
    test_for_kgram_kernel(wt_textsize, wt_alphabetsize, 100);
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

int test_for_wavelettree(int length, int range, int k){//{{{
    for(int _i=2; _i<range; _i<<=1){
        int i = _i + (rand() % _i);
        for(int _j=1; _j<length; _j<<=1){
            int j = _j + (rand() % _j);

            vector<int> S(j);
            for(int l=0; l<j; l++){ S[l] = rand() % i; }
            printf("\na test in the condition |T|=%d, |Σ|=%d starts;\n", j, i);

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
                    printf("S[%d]=%d, WT.access[%d]=%d\n", i, S[i], i, wt.access(i));
                    result = false;
                }
            }
            if(!result){ exit(1); }
            // </a test for access>
            e_time = clock();
            duration = (double)(e_time - s_time) / (double)CLOCKS_PER_SEC;
            printf("Acsess(i) test: OK\t %.10lf [s]\n", duration / S.size());

            // <a test for rangemaxk>
            result = true;
            duration = 0.;
            for(int s=0, end_s=S.size()-k; s<=end_s; ++s){
                vector<int> naive(k);
                for(int x=s, end_x=s+k, y=0; x<end_x;){ naive[y++] = S[x++]; }
                sort(naive.begin(), naive.end());

                vector<int> code(k);
                s_time = clock();
                wt.rangemink(s, s+k, k, code);
                duration += (clock() - s_time);

                if( naive != code ){
                    printf("naive coding(%d, %d, %d): ", s, s+k, k);
                    for(int j=0, j_end=naive.size(); j<j_end; j++){
                        cout << naive[j] << " ";
                    }
                    cout << endl;
                    printf("WT.rangemaxk(%d, %d, %d): ", s, s+k, k);
                    for(int j=0, j_end=code.size(); j<j_end; j++){
                        cout << code[j] << " ";
                    }
                    cout << endl;
                    result = false;
                }
            }
            if(!result){ exit(1); }
            // </a test for rangemaxk>
            duration = (double)(duration) / (double)CLOCKS_PER_SEC;
            printf("Rangemaxk(s,e,%d) test: OK\t %.10lf [s]\n", k, duration);
        }
    }

    return 0;
};//}}}

int comparison_kgram_vector_construct(int length, int range, int k){//{{{
    for(int _i=2; _i<range; _i<<=1){
        int i = _i + (rand() % _i);
        for(int _j=1; _j<length; _j<<=1){
            int j = _j + (rand() % _j);

            vector<int> S(j);
            for(int l=0; l<j; l++){ S[l] = rand() % i; }

            clock_t s_time, e_time, naive_time;
            double duration, naive_duration;
            WaveletTree wt(S);

            // <a test for kgram using natural rep.>
            duration = naive_duration = 0.;
            natRepKgramVector naive_nat_vec, wt_nat_vec;

            naive_time = clock();
            naive_natRepKgramVector(S, k, naive_nat_vec);
            naive_duration += (clock() - s_time);

            s_time = clock();
            wt.createNatRepKgramVector(k, wt_nat_vec);
            duration += (clock() - s_time);

            if( naive_nat_vec != wt_nat_vec ){
                printf("error: something worse happen.\n");
                exit(1);
            }
            // </a test for kgram using natural rep.>
            duration = (double)(duration) / (double)CLOCKS_PER_SEC;
            naive_duration = (double)(naive_duration) / (double)CLOCKS_PER_SEC;
            printf("natRepKgramVector(S,%d,vec) test: OK\t WT:%.10lf [s], naive:%.10lf\n", k, duration, naive_duration);
        }
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

void test_for_kgram_kernel(int length, int alphabetsize, int max_k){
    for(int N_i=1; N_i<length; N_i<<=1){
        int n = N_i + (rand() % N_i); /* create random sequence */
        int m = N_i + (rand() % N_i);
        vector<double> S(n), T(m);

        for(int i=0; i<n; i++){ S[i] = (float)rand() / (float)(RAND_MAX); }
//        for(int i=0; i<n; i++){ cout << S[i] << " "; }
//        cout << endl;
        for(int i=0; i<m; i++){ T[i] = (float)rand() / (float)(RAND_MAX); }
//        for(int i=0; i<m; i++){ cout << T[i] << " "; }
//        cout << endl;

        for(int k=1; k<max_k; k++){
            cout << k << ":" << naive_natRepKgram(S, T, k) << endl;
        }
    }
};

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
/* vim:set foldmethod=marker commentstring=//%s : */
