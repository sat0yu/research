#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<ctime>

#include "./WaveletTree.h"

using namespace std;

int test_for_bitvector(int);
int test_for_wavelettree(int, int);

int main(){
    srand(time(0));
    int bv_textsize = 100000;
    test_for_bitvector(bv_textsize);
    int wt_textsize = 10000, wt_alphabetsize = 100000;
    test_for_wavelettree(wt_textsize, wt_alphabetsize);
}

int test_for_bitvector(int N){//{{{
    for(int N_i=1; N_i<N; N_i<<=1){
        int n = N_i + (rand() % N_i);

        printf("\na test in the condition n=%d starts;\n", n);

        char *B = (char*)malloc((n+1)*sizeof(char));
        for(int i=0; i<n; i++){
            B[i] = (rand() % 2) ? '1' : '0';
        }
        B[N_i] = '\0';

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

int test_for_wavelettree(int length, int range){//{{{
    for(int _i=2; _i<range; _i<<=1){
        int i = _i + (rand() % _i);
        for(int _j=1; _j<length; _j<<=1){
            int j = _j + (rand() % _j);

            vector<int> S(j);
            for(int k=0; k<j; k++){ S[k] = rand() % i; }
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

            s_time = clock();
            // <a test for rangemaxk>
            result = true;
            for(int k=1, end_k=S.size(); k<=end_k; k++){
                for(int s=0, end_s=S.size()-k; s<=end_s; ++s){
                    vector<int> naive(k);
                    /* sort a substring of length k using bucket-sort */
                    vector<char> bucket(UB_ALPHABET_SIZE, 0);
                    for(int j=0; j<k; j++){
                        if( j > UB_ALPHABET_SIZE ){
                            fprintf(stderr,
                                    "too large alphabet size. (the expected maximal size is %d)\n",
                                    UB_ALPHABET_SIZE);
                            exit(1);
                        }
                        bucket[ S[s+j] ] = 1;
                    }
                    vector<char>::iterator it = bucket.begin(), end_it = bucket.end();
                    for(int d_idx=0; it != end_it; ++it){
                        if( *it > 0 ){
                            naive[d_idx++] = (int)distance(bucket.begin(), it);
                        }
                    }

                    /* compare each encoded substrings */
                    vector<int> code = wt.rangemaxk(s,s+k,k);
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
            }
            if(!result){ exit(1); }
            // </a test for rangemaxk>
            e_time = clock();
            duration = (double)(e_time - s_time) / (double)CLOCKS_PER_SEC;
            printf("Rangemaxk(s,e,k) test: OK\t %.10lf [s]\n", duration);
        }
    }

    return 0;
};//}}}

/* vim:set foldmethod=marker commentstring=//%s : */
