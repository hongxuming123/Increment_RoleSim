#ifndef _DYMATCH_H_
#define _DYMATCH_H_

#include <vector>
#include "matchinfo.h"
#include "lap.h"
using namespace std;

//mf为上一次保存的匹配信息，mc用来保留此次计算的结果common与上一次相比的共同节点编号，diff为不同的节点编号


double dymatch_col(int f_dim,int dim, int idx, double** cost, int* rowsol, int* colsol, double* u, double* v); 
double hungarian(int dim, double** cost, int* rowsol, int* colsol, double* u, double* v);
bool dfs(int i, int dim,double** cost, int* rowsol, int* colsol, double* u, double* v);
bool bfs(int r, int dim, double** cost, int* rowsol, int* colsol, double* u, double* v);

double hungarian2(double** cost, int dimx);
bool dfs2(int r, double** cost);


#endif // !_DYMATCH_H_