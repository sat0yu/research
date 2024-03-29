#ifndef OP
#define OP
#include "./order_preserving.h"
#endif

#include<fstream>
#include<string>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<ctime>

#include "./WaveletTree.h"
#include "./RangeCounting.h"

using namespace std;

int comparison_kgram_vector_construct(const int*, int, const int*, int, const int*, int, int, double);

int main(){
    srand(0);
    const int length_list[] = {1000,10000,100000};
    const int sigma_list[] = {100,1000,10000};
    const int k_list[] = {\
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
    // const int length_list[] = {10,100,1000,10000};
    // const int sigma_list[] = {10,100,1000,10000};
    // const int k_list[] = {10,100,1000,10000};
    int length_list_size = sizeof(length_list) / sizeof(length_list[0]),
        sigma_list_size = sizeof(sigma_list) / sizeof(sigma_list[0]),
        k_list_size = sizeof(k_list) / sizeof(k_list[0]);

    clock_t s_time = clock();
    comparison_kgram_vector_construct(
            length_list, length_list_size,
            sigma_list, sigma_list_size,
            k_list, k_list_size, 1000, 10.0);
    double all_duration = (clock() - s_time) / (double)CLOCKS_PER_SEC;
    printf("the time to complete the all processing: %.10lf\n", all_duration);
}

int comparison_kgram_vector_construct(//{{{
        const int* length_list, int length_list_size,
        const int* sigma_list, int sigma_list_size,
        const int* k_list, int k_list_size, int loop, double threshold){

    for(int _i=0; _i<length_list_size; _i++){
        for(int _j=0; _j<sigma_list_size; _j++){
            for(int _k=0; _k<k_list_size; _k++){
                // i:text size, j:sigma size, k:k-parameter
                int i = length_list[_i],
                    j = sigma_list[_j],
                    k = k_list[_k];
                if( (k > i) or (j > i) ){ continue; }

                bool nat_iterate_flg = false;
                clock_t s_time;
                double wt_construct_duration = 0.,
                       rc_construct_duration = 0.;
                double wt_sliding_duration = 0.,
                       nat_iterate_duration = 0.,
                       nat_sliding_duration = 0.,
                       naive_rc_sliding_duration = 0.,
                       naive_sliding_wo_oracle_duration = 0.,
                       rc_sliding_duration = 0.,
                       scc_sliding_duration = 0.,
                       sab_sliding_duration = 0.,
                       temp = 0.;

                // check the processing time takes more than "threshold" sec, only once
                if( !nat_iterate_flg ){
                    vector<int> S(i);
                    for(int m=0; m<i; m++){ S[m] = rand() % j; }
                    natRepKgramVector nat_iterate_vec;
                    s_time = clock();
                    naive_natRepKgramVector(S, k, nat_iterate_vec);
                    temp = (clock() - s_time) / (double)CLOCKS_PER_SEC;
                    // printf("naive_nat_flg: %lf\n",temp);
                    if( temp < threshold ){ nat_iterate_flg = true; }
                }

                for(int l=0; l<loop; l++){
                    vector<int> S(i);
                    for(int m=0; m<i; m++){ S[m] = rand() % j; }

                    s_time = clock();
                    WaveletTree wt(S);
                    wt_construct_duration += (clock() - s_time);

                    s_time = clock();
                    RangeCounting rc(S);
                    rc_construct_duration += (clock() - s_time);

                    if( nat_iterate_flg ){
                        natRepKgramVector nat_iterate_vec;
                        s_time = clock();
                        naive_natRepKgramVector(S, k, nat_iterate_vec);
                        nat_iterate_duration += (clock() - s_time);
                    }else{
                        nat_iterate_duration = -(double)CLOCKS_PER_SEC * loop;
                    }

                    {
                        natRepKgramVector nat_slide_vec;
                        s_time = clock();
                        kgramVector_NaturalRepresentationAndWindowSliding(S, k, nat_slide_vec);
                        nat_sliding_duration += (clock() - s_time);
                    }


                    {
                        countingCodingKgramVector naive_slide_rc_vec;
                        s_time = clock();
                        naive_countingCodingKgramVectorWithSliding(S, k, naive_slide_rc_vec);
                        naive_rc_sliding_duration += (clock() - s_time);
                    }

                    {
                        countingCodingKgramVector wt_slide_rc_vec;
                        s_time = clock();
                        wt.createCountingCodingKgramVectorWithSliding(S, k, wt_slide_rc_vec);
                        wt_sliding_duration += (clock() - s_time);
                    }

                    {
                        countingCodingKgramVector rc_slide_rc_vec;
                        s_time = clock();
                        rc.createCountingCodingKgramVectorWithSliding(S, k, rc_slide_rc_vec);
                        rc_sliding_duration += (clock() - s_time);
                    }

                    {
                        countingCodingKgramVector naive_slide_wo_oracle_cc_vec;
                        s_time = clock();
                        kgramVector_CountingCodingAndWindowSlidingWithoutCharacterOracle(S, k, naive_slide_wo_oracle_cc_vec);
                        naive_sliding_wo_oracle_duration += (clock() - s_time);
                    }

                    {
                        s_time = clock();
                        suffixCountingCodingKgramVector scc_slide_vec;
                        kgramVector_SuffixCountingCodingAndWindowSliding(S, k, scc_slide_vec);
                        scc_sliding_duration += (clock() - s_time);
                    }

                    {
                        suffixAplhaBetaCodingKgramVector sab_slide_vec;
                        s_time = clock();
                        kgramVector_SuffixAlphaBetaCodingWithWindowSliding(S, k, sab_slide_vec);
                        sab_sliding_duration += (clock() - s_time);
                    }

                }
                nat_iterate_duration /= (double)CLOCKS_PER_SEC * loop;
                nat_sliding_duration /= (double)CLOCKS_PER_SEC * loop;
                wt_construct_duration /= (double)CLOCKS_PER_SEC * loop;
                wt_sliding_duration /= (double)CLOCKS_PER_SEC * loop;
                rc_construct_duration /= (double)CLOCKS_PER_SEC * loop;
                rc_sliding_duration /= (double)CLOCKS_PER_SEC * loop;
                naive_rc_sliding_duration /= (double)CLOCKS_PER_SEC * loop;
                naive_sliding_wo_oracle_duration /= (double)CLOCKS_PER_SEC * loop;
                scc_sliding_duration /= (double)CLOCKS_PER_SEC * loop;
                sab_sliding_duration /= (double)CLOCKS_PER_SEC * loop;
                printf("%d,%d,%d,%.10lf,%.10lf,%.10lf,%.10lf,%.10lf,%.10lf,%.10lf,%.10lf,%.10lf,%.10lf,%.10lf,%.10lf \n",
                        i, j, k,
                        wt_construct_duration,
                        wt_sliding_duration,
                        wt_sliding_duration + wt_construct_duration,
                        rc_construct_duration,
                        rc_sliding_duration,
                        rc_sliding_duration + rc_construct_duration,
                        nat_iterate_duration,
                        nat_sliding_duration,
                        naive_rc_sliding_duration,
                        naive_sliding_wo_oracle_duration,
                        scc_sliding_duration,
                        sab_sliding_duration
                    );
                cout << flush;
            }
        }
    }

    return 0;
};//}}}

/* vim:set foldmethod=marker commentstring=//%s : */
