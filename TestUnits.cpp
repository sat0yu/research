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
int test_for_wavelettree(int, int);
void naive_natRepKgramVector(vector<int>&, int, natRepKgramVector&);

int main(){
    srand(0);
    int wt_textsize = 10000, loop=100;
    clock_t s_time;
    test_for_wavelettree(wt_textsize, loop);
    double duration = (double)(clock() - s_time) / (double)CLOCKS_PER_SEC;
    printf("T=10000, loop=10; %.10lf\n", duration);
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

int test_for_wavelettree(int length, int loop){//{{{
    int array[] = {2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000};
    vector<int> sigma_list(array, array+(sizeof(array)/sizeof(int)));
    vector<int>::iterator it = sigma_list.begin(), end_it = sigma_list.end();
    for(; it!=end_it; it++){
        int sigma=*it;

        vector<int> S(length);
        for(int k=0; k<length; k++){ S[k] = rand() % sigma; }
        printf("\na test in the condition |T|=%d, |Σ|=%d starts;\n", length, sigma);

        clock_t s_time, e_time;
        double duration;
        bool result;

        for(int l=0; l<loop; l++){
            s_time = clock();
            // <construst an instance>
            WaveletTree wt(S);
            // </construst an instance>
            duration += (clock() - s_time);
        }
        duration = (double)(duration) / (double)CLOCKS_PER_SEC;
        printf("|T|=%d, |Σ|=%d; construction: OK \t %.10lf [s]\n", length, sigma, duration / loop);

        WaveletTree wt(S);
        for(int k=1; k<length; k<<=1){
            clock_t naive_time;
            double naive_duration = 0.;
            result = true;
            duration = 0.;
            for(int l=0; l<loop; l++){
                natRepKgramVector vec_naive, vec_wt;

                naive_time = clock();
                naive_natRepKgramVector(S, k, vec_naive);
                naive_duration += (clock() - s_time);

                s_time = clock();
                wt.createNatRepKgramVector(k, vec_wt);
                duration += (clock() - s_time);

                if( vec_naive != vec_wt ){
                    printf("error: something worse happen.\n");
                    result = false;
                }
                if(!result){ exit(1); }
                // </a test for kgram>
            }
            duration = (double)(duration) / (double)CLOCKS_PER_SEC;
            naive_duration = (double)(naive_duration) / (double)CLOCKS_PER_SEC;
            printf("|T|=%d, |Σ|=%d, k=%d; \t WT:%.10lf [s], naive:%.10lf [s]\n", length, sigma, k, duration/loop, naive_duration/loop);
            cout << flush;
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

/* vim:set foldmethod=marker commentstring=//%s : */
