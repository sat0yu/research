#include<map>
#include<vector>

using namespace std;

#define UB_ALPHABET_SIZE 100000
#define UB_TEXT_SIZE 1000000

typedef map< vector<int>, int> natRepKgramVector;
typedef pair<int, int> rc_code;
typedef map< vector<rc_code>, int> rangeCountingKgramVector;
