#ifndef __MWM_H_
#define __MWM_H_
#include <vector>

using namespace std;

class MWM
{
private:
	int MAXN;
	int dim;
	double* u, * v;
	int* rowsol, * colsol;
	double** costmatrix;
public:
	MWM(int maxn);
	double getmincost(vector<vector<double>>& G);
	double getmaxweight(vector<vector<double>>& G);
	double getmaxweight_by_hungarian(vector<vector<double>>& G);
	double getmaxweight_by_hungarian2(vector<vector<double>>& G);
	double getmaxweight_by_hungarian3(vector<vector<double>>& G);
};

#endif // !__MWM_H_