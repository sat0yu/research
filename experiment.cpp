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

int comparison_kgram_vector_construct(const int*, int, const int*, int, const int*, int, int);

int main(){
    srand(0);
    const int length_list[] = {100000};
    const int sigma_list[] = {100,1000,10000};
    const int k_list[] = {\
        2,3,4,5,6,7,8,9,\
        10,20,30,40,50,60,70,80,90,\
        100,200,300,400,500,600,700,800,900,\
        1000};
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

int comparison_kgram_vector_construct(//{{{
        const int* length_list, int length_list_size,
        const int* sigma_list, int sigma_list_size,
        const int* k_list, int k_list_size, int loop){

    double threshold = 10.0;
    bool naive_nat_flg = true,
         naive_iterate_flg = true,
         rc_iterate_flg = true;
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
                       rc_sliding_duration = 0.,
                       temp = 0.;

                // check the processing time takes more than "threshold" sec, only once
                if( !naive_nat_flg ){
                    vector<int> S(i);
                    for(int m=0; m<i; m++){ S[m] = rand() % j; }
                    natRepKgramVector naive_nat_vec;
                    s_time = clock();
                    naive_natRepKgramVector(S, k, naive_nat_vec);
                    temp = (clock() - s_time) / (double)CLOCKS_PER_SEC;
                    // printf("naive_nat_flg: %lf\n",temp);
                    if( temp < threshold ){ naive_nat_flg = true; }
                }
                // check the processing time takes more than "threshold" sec, only once
                if( !naive_iterate_flg ){
                    vector<int> S(i);
                    for(int m=0; m<i; m++){ S[m] = rand() % j; }
                    countingCodingKgramVector naive_iterate_rc_vec;
                    s_time = clock();
                    naive_countingCodingKgramVector(S, k, naive_iterate_rc_vec);
                    temp = (clock() - s_time) / (double)CLOCKS_PER_SEC;
                    // printf("naive_iterate_flg: %lf\n",temp);
                    if( temp < threshold ){ naive_iterate_flg = true; }
                }
                // check the processing time takes more than "threshold" sec, only once
                if( !rc_iterate_flg ){
                    vector<int> S(i);
                    for(int m=0; m<i; m++){ S[m] = rand() % j; }
                    RangeCounting rc(S);
                    countingCodingKgramVector rc_iterate_rc_vec;
                    s_time = clock();
                    rc.createCountingCodingKgramVector(S, k, rc_iterate_rc_vec);
                    temp = (clock() - s_time) / (double)CLOCKS_PER_SEC;
                    // printf("rc_iterate_flg: %lf\n",temp);
                    if( temp < threshold ){ rc_iterate_flg = true; }
                }

                for(int l=0; l<loop; l++){
                    vector<int> S(i);
                    for(int m=0; m<i; m++){ S[m] = rand() % j; }

                    // s_time = clock();
                    // <construst an instance>
                    // WaveletTree wt(S);
                    // </construst an instance>
                    // wt_construct_duration += (clock() - s_time);

                    s_time = clock();
                    // <construst an instance>
                    RangeCounting rc(S);
                    // </construst an instance>
                    rc_construct_duration += (clock() - s_time);

                    // natRepKgramVector naive_nat_vec, wt_nat_vec;
                    // countingCodingKgramVector naive_iterate_rc_vec,
                    //                          naive_slide_rc_vec,
                    //                          wt_iterate_rc_vec,
                    //                          wt_slide_rc_vec,
                    //                          rc_iterate_rc_vec,
                    //                          rc_slide_rc_vec;

                    if( naive_nat_flg ){
                        s_time = clock();
                        natRepKgramVector naive_nat_vec;
                        naive_natRepKgramVector(S, k, naive_nat_vec);
                        naive_nat_duration += (clock() - s_time);
                    }else{
                        naive_nat_duration = -(double)CLOCKS_PER_SEC * loop;
                    }

                    if( naive_iterate_flg ){
                        s_time = clock();
                        countingCodingKgramVector naive_iterate_rc_vec;
                        naive_countingCodingKgramVector(S, k, naive_iterate_rc_vec);
                        naive_iterate_duration += (clock() - s_time);
                    }else{
                        naive_iterate_duration = -(double)CLOCKS_PER_SEC * loop;
                    }

                    {
                        s_time = clock();
                        countingCodingKgramVector naive_slide_rc_vec;
                        naive_countingCodingKgramVectorWithSliding(S, k, naive_slide_rc_vec);
                        naive_sliding_duration += (clock() - s_time);
                    }

                    // s_time = clock();
                    // wt.createNatRepKgramVector(S, k, wt_nat_vec);
                    // wt_nat_duration += (clock() - s_time);
                    //
                    // s_time = clock();
                    // wt.createCountingCodingKgramVector(S, k, wt_iterate_rc_vec);
                    // wt_iterate_duration += (clock() - s_time);
                    //
                    // s_time = clock();
                    // wt.createCountingCodingKgramVectorWithSliding(S, k, wt_slide_rc_vec);
                    // wt_sliding_duration += (clock() - s_time);

                    if( rc_iterate_flg ){
                        s_time = clock();
                        countingCodingKgramVector rc_iterate_rc_vec;
                        rc.createCountingCodingKgramVector(S, k, rc_iterate_rc_vec);
                        rc_iterate_duration += (clock() - s_time);
                    }else{
                        rc_iterate_duration = -(double)CLOCKS_PER_SEC * loop;
                    }


                    // s_time = clock();
                    // rc.createCountingCodingKgramVectorWithSliding(S, k, rc_slide_rc_vec);
                    // rc_sliding_duration += (clock() - s_time);
                }
                wt_construct_duration /= (double)CLOCKS_PER_SEC * loop;
                rc_construct_duration /= (double)CLOCKS_PER_SEC * loop;
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

                if( naive_nat_duration > threshold ){ naive_nat_flg = false; }
                if( naive_iterate_duration > threshold ){ naive_iterate_flg = false; }
                if( rc_iterate_duration > threshold ){ rc_iterate_flg = false; }
            }
        }
    }

    return 0;
};//}}}

/* vim:set foldmethod=marker commentstring=//%s : */
