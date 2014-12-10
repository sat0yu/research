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

void naive_natRepKgramVector(vector<int>&, int, natRepKgramVector&);
void naive_rangeCountingKgramVector(vector<int>&, int, rangeCountingKgramVector&);
void naive_rangeCountingKgramVectorWithSliding(vector<int>&, int, rangeCountingKgramVector&);
int comparison_kgram_vector_construct(const int*, int, const int*, int, const int*, int, int);

int main(){
    srand(0);
    const int length_list[] = {100000};
    const int sigma_list[] = {\
              100, 200, 300, 400, 500, 600, 700, 800, 900,\
        1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,\
        2000,2100,2200,2300,2400,2500,2600,2700,2800,2900,\
        3000,3100,3200,3300,3400,3500,3600,3700,3800,3900,\
        4000,4100,4200,4300,4400,4500,4600,4700,4800,4900,\
        5000,5100,5200,5300,5400,5500,5600,5700,5800,5900,\
        6000,6100,6200,6300,6400,6500,6600,6700,6800,6900,\
        7000,7100,7200,7300,7400,7500,7600,7700,7800,7900,\
        8000,8100,8200,8300,8400,8500,8600,8700,8800,8900,\
        9000,9100,9200,9300,9400,9500,9600,9700,9800,9900,\
        10000};
    const int k_list[] = {10,100,1000};
    // const int length_list[] = {10,100,1000,10000};
    // const int sigma_list[] = {10,100,1000,10000};
    // const int k_list[] = {10,100,1000,10000};
    int length_list_size = sizeof(length_list) / sizeof(length_list[0]),
        sigma_list_size = sizeof(sigma_list) / sizeof(sigma_list[0]),
        k_list_size = sizeof(k_list) / sizeof(k_list[0]);
    comparison_kgram_vector_construct(
            length_list, length_list_size,
            sigma_list, sigma_list_size,
            k_list, k_list_size, 100);
}

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

void naive_rangeCountingKgramVector(vector<int>& S, int k, rangeCountingKgramVector& res){//{{{
    vector<rc_code> kgram(k);
    for(int i=0, end_i=S.size()-k+1; i<end_i; i++){
        for(int j=0; j<k; j++){ /* for each a substring of length k */
            int lt=0, eq=0;
            for(int l=j-1; l>=0; l--){ /* range counting */
                if(S[i+l] == S[i+j]){
                    eq++;
                }else if(S[i+l] < S[i+j]){
                    lt++;
                }
            }
            kgram[j] = rc_code(lt, eq);
        }

        if( res.find(kgram) == res.end() ){ /* regist the kgram */
            res[kgram] = 1;
        }else{
            res[kgram]++;
        }
    }
}//}}}

void naive_rangeCountingKgramVectorWithSliding(vector<int>& S, int k, rangeCountingKgramVector& res){//{{{
    vector<rc_code> kgram(k);
    for( int j=0; j<k; ++j ){ // the first kgram is calcucated naively
        int lt=0, eq=0;
        for(int l=j-1; l>=0; l--){ /* range counting */
            if(S[l] < S[j]){
                lt++;
            }else if(S[l] == S[j]){
                eq++;
            }
        }
        kgram[j] = rc_code(lt, eq);
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
        for(int l=k-2; l>=0; l--){ /* range counting */
            if(S[i+l] < S[i+k-1]){
                lt++;
            }else if(S[i+l] == S[i+k-1]){
                eq++;
            }
        }
        kgram[k-1] = rc_code(lt, eq);

        if( res.find(kgram) == res.end() ){ /* regist the kgram */
            res[kgram] = 1;
        }else{
            res[kgram]++;
        }
    }
}//}}}

int comparison_kgram_vector_construct(//{{{
        const int* length_list, int length_list_size,
        const int* sigma_list, int sigma_list_size,
        const int* k_list, int k_list_size, int loop){
    for(int _i=0; _i<length_list_size; _i++){
        for(int _j=0; _j<sigma_list_size; _j++){
            for(int _k=0; _k<k_list_size; _k++){
                // i:text size, j:sigma size, k:k-parameter
                int i = length_list[_i],
                    j = sigma_list[_j],
                    k = k_list[_k];
                if( (k > i) or (j > i) ){ continue; }

                clock_t s_time;
                double wt_construct_duration = 0.,
                       wt_nat_duration = 0.,
                       wt_iterate_duration = 0.,
                       wt_sliding_duration = 0.,
                       naive_nat_duration = 0.,
                       naive_iterate_duration = 0.,
                       naive_sliding_duration = 0.,
                       rc_construct_duration = 0.,
                       rc_iterate_duration = 0.,
                       rc_sliding_duration = 0.;


                for(int l=0; l<loop; l++){
                    vector<int> S(i);
                    for(int m=0; m<i; m++){ S[m] = rand() % j; }

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

                    // s_time = clock();
                    // naive_natRepKgramVector(S, k, naive_nat_vec);
                    // naive_nat_duration += (clock() - s_time);

                    // s_time = clock();
                    // naive_rangeCountingKgramVector(S, k, naive_iterate_rc_vec);
                    // naive_iterate_duration += (clock() - s_time);

                    {
                        s_time = clock();
                        rangeCountingKgramVector naive_slide_rc_vec;
                        naive_rangeCountingKgramVectorWithSliding(S, k, naive_slide_rc_vec);
                        naive_sliding_duration += (clock() - s_time);
                    }

                    // s_time = clock();
                    // wt.createNatRepKgramVector(S, k, wt_nat_vec);
                    // wt_nat_duration += (clock() - s_time);

                    // s_time = clock();
                    // wt.createRangeCountingKgramVector(S, k, wt_iterate_rc_vec);
                    // wt_iterate_duration += (clock() - s_time);

                    {
                        s_time = clock();
                        rangeCountingKgramVector wt_slide_rc_vec;
                        wt.createRangeCountingKgramVectorWithSliding(S, k, wt_slide_rc_vec);
                        wt_sliding_duration += (clock() - s_time);
                    }

                    // s_time = clock();
                    // rc.createRangeCountingKgramVector(S, k, rc_iterate_rc_vec);
                    // rc_iterate_duration += (clock() - s_time);

                    {
                        s_time = clock();
                        rangeCountingKgramVector rc_slide_rc_vec;
                        rc.createRangeCountingKgramVectorWithSliding(S, k, rc_slide_rc_vec);
                        rc_sliding_duration += (clock() - s_time);
                    }
                }
                wt_construct_duration /= (double)CLOCKS_PER_SEC * loop;
                rc_construct_duration /= (double)CLOCKS_PER_SEC;
                wt_nat_duration /= (double)CLOCKS_PER_SEC * loop;
                wt_iterate_duration /= (double)CLOCKS_PER_SEC * loop;
                wt_sliding_duration /= (double)CLOCKS_PER_SEC * loop;
                naive_nat_duration /= (double)CLOCKS_PER_SEC * loop;
                naive_iterate_duration /= (double)CLOCKS_PER_SEC * loop;
                naive_sliding_duration /= (double)CLOCKS_PER_SEC * loop;
                rc_iterate_duration /= (double)CLOCKS_PER_SEC * loop;
                rc_sliding_duration /= (double)CLOCKS_PER_SEC * loop;
                // printf("(|T|=%d, sigma=%d, k=%d) WT(construction):%.10lf, RC(construction):%.10lf, "
                //         "WT(nat):%.10lf, WT(iterate):%.10lf, WT(slide):%.10lf, "
                //         "naive(nat):%.10lf, naive(iterate):%.10lf, naive(slide):%.10lf, "
                //         "RC(iterate):%.10lf, RC(slide):%.10lf \n",
                //         i, j, k, wt_construct_duration, rc_construct_duration,
                //         wt_nat_duration, wt_iterate_duration, wt_sliding_duration,
                //         naive_nat_duration, naive_iterate_duration, naive_sliding_duration,
                //         rc_iterate_duration, rc_sliding_duration);
                printf("%d,%d,%d,%.10lf,%.10lf,%.10lf,%.10lf,%.10lf,%.10lf,%.10lf,%.10lf,%.10lf,%.10lf \n",
                        i, j, k, wt_construct_duration, rc_construct_duration,
                        wt_nat_duration, wt_iterate_duration, wt_sliding_duration,
                        naive_nat_duration, naive_iterate_duration, naive_sliding_duration,
                        rc_iterate_duration, rc_sliding_duration);
                cout << flush;
            }
        }
    }

    return 0;
};//}}}

/* vim:set foldmethod=marker commentstring=//%s : */
