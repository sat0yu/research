#ifndef OP
#define OP
#include "./order_preserving.h"
#endif

#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<math.h>
#include<ctime>

using namespace std;

int count_inversions(vector< pair<int, bool> >& S){
    return 0;
};

class RangeCounting{
private:
    static int const w = 32; // w: 32 bits
    static int const L = 6; // L: ceil( sqrt(w) )
    static int const A = 5; // A: (1 + epsilon) * ceil( log(w) )
    static int const H = 5; // H: epsilon * ceil( log(w) )
    int n, l;
public:
    ~RangeCounting(){};
    RangeCounting(vector<int>&);
    int query(int, int);
};

RangeCounting::RangeCounting(vector<int>& S): n( (int)S.size() ){
};

int RangeCounting::query(int x, int y){
    if( l == 0 ){ // initialization
        l = (int)ceil(log(n));
        vector<int> bucket(UB_ALPHABET_SIZE);
    }

    if( (l <= H) and (n <= (1 << A)) ){ // case: 0'
    }else if( (l <= H) and (n > (1 << A)) ){ // case: 0''
    }else if( (H < l) and (l <= L) ){ // case: 1
    }else if( l > L ){ // case: 2
    }else{
        fprintf(stderr, "not matched for any case\n");
        exit(1);
    }

    return 0;
}
/* vim:set foldmethod=marker commentstring=//%s : */
