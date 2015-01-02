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

int comparison_character_oracle(const int*, int, const int*, int, int);

int main(){
    srand(0);
    const int length_list[] = {
        10,20,50,
        100,200,500,
        1000,2000,5000,
        10000,20000,50000,
        100000,200000,500000,
        1000000,2000000,5000000,
        10000000 };
    const int sigma_list[] = {10,100,1000,10000};
    int length_list_size = sizeof(length_list) / sizeof(length_list[0]),
        sigma_list_size = sizeof(sigma_list) / sizeof(sigma_list[0]);

    clock_t s_time = clock();
    comparison_character_oracle(
            length_list, length_list_size,
            sigma_list, sigma_list_size, 1000);
    double all_duration = (clock() - s_time) / (double)CLOCKS_PER_SEC;
    printf("the time to complete the all processing: %.10lf\n", all_duration);
}

int comparison_character_oracle(//{{{
        const int* length_list, int length_list_size,
        const int* sigma_list, int sigma_list_size, int loop){

    for(int _i=0; _i<length_list_size; _i++){
        for(int _j=0; _j<sigma_list_size; _j++){
            // i:text size, j:sigma size
            int i = length_list[_i],
                j = sigma_list[_j];

            clock_t s_time;
            double naive_duration = 0.,
                   wt_construct_duration = 0.,
                   rc_construct_duration = 0.,
                   wt_oracle_duration = 0.,
                   rc_oracle_duration = 0.;

            for(int l=0; l<loop; l++){
                vector<int> S(i);
                for(int m=0; m<i; m++){ S[m] = rand() % j; }

                s_time = clock();
                WaveletTree wt(S);
                wt_construct_duration += (clock() - s_time);

                s_time = clock();
                RangeCounting rc(S);
                rc_construct_duration += (clock() - s_time);

                const int st = 0;
                for(int en=1; en<i; en++){
                    int naive_lt = 0, naive_eq = 0,
                        wt_lt = 0, wt_eq = 0,
                        rc_lt = 0, rc_eq = 0;

                    // naive
                    s_time = clock();
                    for(int l=st; l<en; l++){
                        if(S[l] == S[en]){
                            naive_eq++;
                        }else if(S[l] < S[en]){
                            naive_lt++;
                        }
                    }
                    naive_duration += (clock() - s_time);

                    // wavelet
                    s_time = clock();
                    wt.rankLessThanEqual(S[en], st, en, &wt_lt, &wt_eq);
                    wt_oracle_duration += (clock() - s_time);

                    // Chan's data-structure
                    s_time = clock();
                    int ij_Sij = rc.query(en),
                        i_Sij  = rc.query(st, S[en]);
                    rc_lt = ij_Sij - i_Sij;
                    rc_eq = rc.query(en, S[en]+1) - rc.query(st, S[en]+1) - ij_Sij + i_Sij;
                    rc_oracle_duration += (clock() - s_time);

                    // if( (naive_lt != wt_lt) or (naive_eq != wt_eq) ){
                    //     printf("st:%d.en%d; wt fault: naive->(%d,%d), wt->(%d,%d)\n", st, en, naive_lt,naive_eq,wt_lt,wt_eq);
                    //     exit(1);
                    // }
                    // if( (naive_lt != rc_lt) or (naive_eq != rc_eq) ){
                    //     printf("st:%d.en%d; rc fault: naive->(%d,%d), rc->(%d,%d)\n", st, en, naive_lt,naive_eq,rc_lt,rc_eq);
                    //     exit(1);
                    // }

                }
            }
            naive_duration /= (double)CLOCKS_PER_SEC * loop;
            wt_construct_duration /= (double)CLOCKS_PER_SEC * loop;
            wt_oracle_duration /= (double)CLOCKS_PER_SEC * loop;
            rc_construct_duration /= (double)CLOCKS_PER_SEC * loop;
            rc_oracle_duration /= (double)CLOCKS_PER_SEC * loop;
            printf("%d,%d,%.10lf,%.10lf,%.10lf,%.10lf,%.10lf,%.10lf,%.10lf,%.10lf,%.10lf,%.10lf,%.10lf,%.10lf \n",
                    i, j,
                    naive_duration,
                    naive_duration / i,
                    wt_construct_duration,
                    wt_oracle_duration,
                    wt_oracle_duration / i,
                    wt_construct_duration + wt_oracle_duration,
                    (wt_construct_duration / i) + wt_oracle_duration,
                    rc_construct_duration,
                    rc_oracle_duration,
                    rc_oracle_duration / i,
                    rc_construct_duration + rc_oracle_duration,
                    (rc_construct_duration / i) + rc_oracle_duration
                );
            cout << flush;
        }
    }

    return 0;
};//}}}

/* vim:set foldmethod=marker commentstring=//%s : */
