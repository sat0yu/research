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

class RangeCounting{
private:
public:
    ~RangeCounting();
    RangeCounting(vector<int>&);
    int query(int, int) const;
};

RangeCounting::RangeCounting(vector<int>& S){
};

int RangeCounting::query(int x, int y) const{
    return 0;
}
/* vim:set foldmethod=marker commentstring=//%s : */
