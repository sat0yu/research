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
#include "./RangeCounting.h"

using namespace std;


int test_for_bitvector(int);
int test_for_wavelettree(int, int, int);
int test_for_rangecounting(int);
int comparison_kgram_vector_construct(int, int, int);
int test_for_SuffixCountingCoding(int);

int main(){
    srand(time(0));

    int bv_textsize = 10000;
    test_for_bitvector(bv_textsize);
    int wt_textsize = 1000, wt_alphabetsize = 1000, k=10;
    test_for_wavelettree(wt_textsize, wt_alphabetsize, k);
    int rc_textsize = 1000;
    test_for_rangecounting(rc_textsize);
    int scc_textsize = 2000;
    test_for_SuffixCountingCoding(scc_textsize);

    int textsize = 1000, range = 1000, max_k = 100;
    comparison_kgram_vector_construct(textsize, range, max_k);

    printf("\nAll Tests are passed. Acceptance\n");
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

int test_for_rangecounting(int length){//{{{
    for(int _i=2; _i<length; _i<<=1){
        int i = _i + (rand() % _i);

        vector<int> S(i);
        for(int j=0; j<i; j++){ S[j] = rand() % i; }
        printf("\na test in the condition |T|=%d starts;\n", i);

        clock_t s_time;
        double construct_duration, query_duration;
        bool result;

        s_time = clock();
        // <construst an instance>
        RangeCounting rc(S);
        // </construst an instance>
        construct_duration = clock() - s_time;
        construct_duration = (double)(construct_duration) / (double)CLOCKS_PER_SEC;
        printf("construction: OK \t %f [s]\n", construct_duration);

        // <a test for rangecounting>
        result = true;
        query_duration = 0.;
        for(int j=0; j<i; ++j){
            int naive=0;
            for(int k=0; k<j; k++){
                if( S[k] < S[j] ){ naive++; }
            }

            s_time = clock();
            int rc_query = rc.query(j);
            query_duration += (clock() - s_time);
            if( naive != rc_query ){
                printf("RangeCount query(%d): naive: %d, RangeCounting:%d\n", j, naive, rc_query);
                result = false;
            }
        }
        if(!result){ exit(1); }
        // </a test for rangecounting>
        query_duration = (double)(query_duration) / (double)CLOCKS_PER_SEC;
        printf("RC.query test: OK\t %.10lf [s]\n", query_duration / i);
    }
    return 0;
};//}}}

int test_for_SuffixCountingCoding(int length){//{{{
    for(int i=1; i<length; i<<=1){
        for(int k=1; k<i; k<<=1){
            vector<int> S(i);
            for(int j=0; j<i; j++){ S[j] = rand(); }
            printf("|T|=%d, k=%d\n", i, k);

            clock_t s_time;
            double iterate_duration = 0.,
                   slide_duration = 0.;

            suffixCountingCodingKgramVector iterate_scc_vec, slide_scc_vec;

            s_time = clock();
            kgramVector_SuffixCountingCoding(S, k, iterate_scc_vec);
            iterate_duration += (clock() - s_time);

            s_time = clock();
            kgramVector_SuffixCountingCodingAndWindowSliding(S, k, slide_scc_vec);
            slide_duration += (clock() - s_time);

            if( iterate_scc_vec != slide_scc_vec ){//{{{
                printf("error: something worse happen.\n");

                suffixCountingCodingKgramVector::iterator
                    n_it=iterate_scc_vec.begin(),
                    end_n_id=iterate_scc_vec.end();
                cout << "iteration coding" << endl;
                for(; n_it != end_n_id; n_it++){
                    vector<c_code> kgram = n_it->first;
                    vector<c_code>::iterator it=kgram.begin(), end_it=kgram.end();
                    for(; it!=end_it; it++){
                        printf("(%d, %d) ", it->first, it->second);
                    }
                    cout << endl;
                }

                suffixCountingCodingKgramVector::iterator
                    w_it=slide_scc_vec.begin(),
                    end_w_id=slide_scc_vec.end();
                cout << "sliding coding" << endl;
                for(; w_it != end_w_id; w_it++){
                    vector<c_code> kgram = w_it->first;
                    vector<c_code>::iterator it=kgram.begin(), end_it=kgram.end();
                    for(; it!=end_it; it++){
                        printf("(%d, %d) ", it->first, it->second);
                    }
                    cout << endl;
                }
                exit(1);
            }//}}}

            // </a test for kgram using suffix counting rep.>
            iterate_duration /= (double)CLOCKS_PER_SEC;
            slide_duration /= (double)CLOCKS_PER_SEC;
            printf("suffixCountingCodingKgramVector(S,%d,vec) test: OK\n "
                    "iterate:%.10lf, slide:%.10lf\n", k, iterate_duration, slide_duration);
            cout << flush;
        }
    }
    return 0;
};//}}}

int comparison_kgram_vector_construct(int length, int range, int k){//{{{
    for(int _i=2; _i<range; _i<<=1){
        int i = _i + (rand() % _i);
        for(int _j=1; _j<length; _j<<=1){
            int j = _j + (rand() % _j);
            if( k > j ){ continue; }

            vector<int> S(j);
            for(int l=0; l<j; l++){ S[l] = rand() % i; }
            printf("|T|=%d, sigma=%d, k=%d\n", j, i, k);

            clock_t s_time;
            double wt_construct_duration = 0.,
                   wt_sliding_duration = 0.,
                   naive_nat_duration = 0.,
                   naive_iterate_duration = 0.,
                   naive_sliding_duration = 0.,
                   naive_sliding_wo_oracle_duration = 0.,
                   scc_iterate_duration = 0.,
                   scc_sliding_duration = 0.,
                   rc_construct_duration = 0.,
                   rc_sliding_duration = 0.;


            s_time = clock();
            // <construst an instance>
            WaveletTree wt(S);
            // </construst an instance>
            wt_construct_duration += (clock() - s_time);

            s_time = clock();
            // <construst an instance>
            RangeCounting rc(S);
            // </construst an instance>
            rc_construct_duration += (clock() - s_time);

            natRepKgramVector naive_nat_vec;
            countingCodingKgramVector naive_iterate_cc_vec,
                                      naive_slide_cc_vec,
                                      naive_slide_wo_oracle_cc_vec,
                                      wt_slide_cc_vec,
                                      rc_slide_cc_vec;
            suffixCountingCodingKgramVector scc_iterate_vec,
                                            scc_slide_vec;

            s_time = clock();
            naive_natRepKgramVector(S, k, naive_nat_vec);
            naive_nat_duration += (clock() - s_time);

            s_time = clock();
            naive_countingCodingKgramVector(S, k, naive_iterate_cc_vec);
            naive_iterate_duration += (clock() - s_time);

            s_time = clock();
            naive_countingCodingKgramVectorWithSliding(S, k, naive_slide_cc_vec);
            naive_sliding_duration += (clock() - s_time);

            s_time = clock();
            kgramVector_CountingCodingAndWindowSlidingWithoutCharacterOracle(S, k, naive_slide_wo_oracle_cc_vec);
            naive_sliding_wo_oracle_duration += (clock() - s_time);

            s_time = clock();
            kgramVector_SuffixCountingCoding(S, k, scc_iterate_vec);
            scc_iterate_duration += (clock() - s_time);

            s_time = clock();
            kgramVector_SuffixCountingCodingAndWindowSliding(S, k, scc_slide_vec);
            scc_sliding_duration += (clock() - s_time);

            s_time = clock();
            wt.createCountingCodingKgramVectorWithSliding(S, k, wt_slide_cc_vec);
            wt_sliding_duration += (clock() - s_time);

            s_time = clock();
            rc.createCountingCodingKgramVectorWithSliding(S, k, rc_slide_cc_vec);
            rc_sliding_duration += (clock() - s_time);

            if( naive_iterate_cc_vec != naive_slide_cc_vec ){//{{{
                printf("error: something worse happen.\n");

                countingCodingKgramVector::iterator n_it=naive_iterate_cc_vec.begin(),
                                                    end_n_id=naive_iterate_cc_vec.end();
                cout << "naive coding" << endl;
                for(; n_it != end_n_id; n_it++){
                    vector<c_code> kgram = n_it->first;
                    vector<c_code>::iterator it=kgram.begin(), end_it=kgram.end();
                    for(; it!=end_it; it++){
                        printf("(%d, %d) ", it->first, it->second);
                    }
                    cout << endl;
                }

                countingCodingKgramVector::iterator w_it=naive_slide_cc_vec.begin(), end_w_id=naive_slide_cc_vec.end();
                cout << "naive sliding coding" << endl;
                for(; w_it != end_w_id; w_it++){
                    vector<c_code> kgram = w_it->first;
                    vector<c_code>::iterator it=kgram.begin(), end_it=kgram.end();
                    for(; it!=end_it; it++){
                        printf("(%d, %d) ", it->first, it->second);
                    }
                    cout << endl;
                }
                exit(1);
            }//}}}

            if( naive_iterate_cc_vec != naive_slide_wo_oracle_cc_vec ){//{{{
                printf("error: something worse happen.\n");

                countingCodingKgramVector::iterator n_it=naive_iterate_cc_vec.begin(),
                                                    end_n_id=naive_iterate_cc_vec.end();
                cout << "naive coding" << endl;
                for(; n_it != end_n_id; n_it++){
                    vector<c_code> kgram = n_it->first;
                    vector<c_code>::iterator it=kgram.begin(), end_it=kgram.end();
                    for(; it!=end_it; it++){
                        printf("(%d, %d) ", it->first, it->second);
                    }
                    cout << endl;
                }

                countingCodingKgramVector::iterator w_it=naive_slide_wo_oracle_cc_vec.begin(), \
                                                 end_w_id=naive_slide_wo_oracle_cc_vec.end();
                cout << "naive sliding coding without oracle" << endl;
                for(; w_it != end_w_id; w_it++){
                    vector<c_code> kgram = w_it->first;
                    vector<c_code>::iterator it=kgram.begin(), end_it=kgram.end();
                    for(; it!=end_it; it++){
                        printf("(%d, %d) ", it->first, it->second);
                    }
                    cout << endl;
                }
                exit(1);
            }//}}}

            if( naive_iterate_cc_vec != wt_slide_cc_vec ){//{{{
                printf("error: something worse happen.\n");

                countingCodingKgramVector::iterator n_it=naive_iterate_cc_vec.begin(),
                                                    end_n_id=naive_iterate_cc_vec.end();
                cout << "naive coding" << endl;
                for(; n_it != end_n_id; n_it++){
                    vector<c_code> kgram = n_it->first;
                    vector<c_code>::iterator it=kgram.begin(), end_it=kgram.end();
                    for(; it!=end_it; it++){
                        printf("(%d, %d) ", it->first, it->second);
                    }
                    cout << endl;
                }

                countingCodingKgramVector::iterator w_it=wt_slide_cc_vec.begin(), end_w_id=wt_slide_cc_vec.end();
                cout << "wavelet sliding coding" << endl;
                for(; w_it != end_w_id; w_it++){
                    vector<c_code> kgram = w_it->first;
                    vector<c_code>::iterator it=kgram.begin(), end_it=kgram.end();
                    for(; it!=end_it; it++){
                        printf("(%d, %d) ", it->first, it->second);
                    }
                    cout << endl;
                }
                exit(1);
            }
//}}}

            if( naive_iterate_cc_vec != rc_slide_cc_vec ){//{{{
                printf("error: something worse happen.\n");

                countingCodingKgramVector::iterator n_it=naive_iterate_cc_vec.begin(),
                                                    end_n_id=naive_iterate_cc_vec.end();
                cout << "naive coding" << endl;
                for(; n_it != end_n_id; n_it++){
                    vector<c_code> kgram = n_it->first;
                    vector<c_code>::iterator it=kgram.begin(), end_it=kgram.end();
                    for(; it!=end_it; it++){
                        printf("(%d, %d) ", it->first, it->second);
                    }
                    cout << endl;
                }

                countingCodingKgramVector::iterator w_it=rc_slide_cc_vec.begin(), end_w_id=rc_slide_cc_vec.end();
                cout << "counting sliding coding" << endl;
                for(; w_it != end_w_id; w_it++){
                    vector<c_code> kgram = w_it->first;
                    vector<c_code>::iterator it=kgram.begin(), end_it=kgram.end();
                    for(; it!=end_it; it++){
                        printf("(%d, %d) ", it->first, it->second);
                    }
                    cout << endl;
                }
                exit(1);
            }
//}}}

            if( scc_iterate_vec != scc_slide_vec ){//{{{
                printf("error: something worse happen.\n");

                suffixCountingCodingKgramVector::iterator
                    n_it=scc_iterate_vec.begin(),
                    end_n_id=scc_iterate_vec.end();
                cout << "iteration coding" << endl;
                for(; n_it != end_n_id; n_it++){
                    vector<c_code> kgram = n_it->first;
                    vector<c_code>::iterator it=kgram.begin(), end_it=kgram.end();
                    for(; it!=end_it; it++){
                        printf("(%d, %d) ", it->first, it->second);
                    }
                    cout << endl;
                }

                suffixCountingCodingKgramVector::iterator
                    w_it=scc_slide_vec.begin(),
                    end_w_id=scc_slide_vec.end();
                cout << "sliding coding" << endl;
                for(; w_it != end_w_id; w_it++){
                    vector<c_code> kgram = w_it->first;
                    vector<c_code>::iterator it=kgram.begin(), end_it=kgram.end();
                    for(; it!=end_it; it++){
                        printf("(%d, %d) ", it->first, it->second);
                    }
                    cout << endl;
                }
                exit(1);
            }//}}}

            // </a test for kgram using counting rep.>
            wt_construct_duration /= (double)CLOCKS_PER_SEC;
            rc_construct_duration /= (double)CLOCKS_PER_SEC;
            wt_sliding_duration /= (double)CLOCKS_PER_SEC;
            naive_nat_duration /= (double)CLOCKS_PER_SEC;
            naive_iterate_duration /= (double)CLOCKS_PER_SEC;
            naive_sliding_duration /= (double)CLOCKS_PER_SEC;
            naive_sliding_wo_oracle_duration /= (double)CLOCKS_PER_SEC;
            scc_iterate_duration /= (double)CLOCKS_PER_SEC;
            scc_sliding_duration /= (double)CLOCKS_PER_SEC;
            rc_sliding_duration /= (double)CLOCKS_PER_SEC;
            printf("countingCodingKgramVector(S,%d,vec) test: OK\n "
                    "WT(construction):%.10lf, WT(slide):%.10lf, "
                    "naive(nat):%.10lf, naive(iterate):%.10lf, naive(slide):%.10lf, naive(slide w/ oracle):%.10lf, "
                    "scc(iterate):%.10lf, scc(slide):%.10lf, "
                    "RC(construction):%.10lf, RC(slide):%.10lf \n",
                    k, wt_construct_duration, wt_sliding_duration,
                    naive_nat_duration, naive_iterate_duration, naive_sliding_duration, naive_sliding_wo_oracle_duration,
                    scc_iterate_duration, scc_sliding_duration,
                    rc_construct_duration, rc_sliding_duration);
            cout << flush;
        }
    }

    return 0;
};//}}}

/* vim:set foldmethod=marker commentstring=//%s : */
