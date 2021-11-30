#ifndef _DYMATCH_H_
#define _DYMATCH_H_

#include <vector>
#include "matchinfo.h"
#include "lap.h"
using namespace std;

//mfΪ��һ�α����ƥ����Ϣ��mc���������˴μ���Ľ��common����һ����ȵĹ�ͬ�ڵ��ţ�diffΪ��ͬ�Ľڵ���


double dymatch_col(int f_dim,int dim, int idx, double** cost, int* rowsol, int* colsol, double* u, double* v); 
double hungarian(int dim, double** cost, int* rowsol, int* colsol, double* u, double* v);
bool dfs(int i, int dim,double** cost, int* rowsol, int* colsol, double* u, double* v);
bool bfs(int r, int dim, double** cost, int* rowsol, int* colsol, double* u, double* v);

double hungarian2(double** cost, int dimx);
bool dfs2(int r, double** cost);


#endif // !_DYMATCH_H_