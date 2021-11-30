#ifndef __MATCHINFO_H_
#define __MATCHINFO_H_
#include <unordered_map>
using namespace std;
typedef struct Matchinfo {
	int dim;
	double** cost;
	double* u;
	double* v;
	int* colsol;
	int* rowsol;
	int* index;
	//unordered_map<int, int> index;
} Mi;

#endif // !__MATCHINFO_H_
