#ifndef  __HUNGARIAN4_H_
#define  __HUNGARIAN4_H_


#include <vector>
#include <algorithm>
using namespace std;
typedef vector<double> VD;
typedef vector<VD> VVD;
typedef vector<int> VI;

double MinCostMatching(const VVD& cost, VI& Lmate, VI& Rmate);




#endif // ! __HUNGARIAN4_H_
