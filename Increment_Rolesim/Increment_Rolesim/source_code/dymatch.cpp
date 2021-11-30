#include "dymatch.h"
#include <queue>
#include <set>
#include <iostream>
#include <ctime>
using namespace std;
static const int MAXN = 15000;
static bool visrow[MAXN];
static bool viscol[MAXN];
static double slack[MAXN];
static double oo = 0x3f3f3f3f;
static double eps = 1e-6;

static double rowdul[MAXN], coldul[MAXN];
static int match[MAXN];
int dim;
static int q[MAXN];
static int fau[MAXN];
static int fav[MAXN];
int head, tail;

double dymatch_col(int f_dim,int c_dim, int idx, double** cost,int* rowsol,int* colsol, double* u, double* v){
	for (int i = idx; i < f_dim; ++i) {
		double minval=oo;
		for (int j = 0; j < f_dim; ++j) {
			minval=min(minval,cost[j][i]-u[j]);
		}
		v[i]=minval;
	}
	set<int> matched;
	for (int i = 0; i < idx; ++i) {
		matched.insert(colsol[i]);
	}
	for (int i = 0; i < f_dim; ++i) {
		if (matched.find(i) != matched.end()) continue;
		for (int j = 0; j < f_dim; ++j)slack[j] = oo;
		while(1){
			for(int j=0;j< f_dim;++j) visrow[j]=0;
			for(int j=0;j< f_dim;++j) viscol[j]=0;
			if (dfs(i, f_dim, cost, rowsol, colsol, u, v)) break;
			double theta=oo;
			for(int j=0;j< f_dim;++j){
				if(!viscol[j]){
					theta=min(theta,slack[j]);
				}
			}
			for (int j = 0; j < f_dim; ++j) {
				visrow[j] ? u[j] += theta : u[j] -= theta; //Adjust the dual variable
				viscol[j] ? v[j] -= theta : v[j] += theta;
			}
		}
	}
	for (int k = f_dim; k < c_dim; ++k) {
		double minval = oo;
		for (int i = 0; i < k; ++i) {
			minval = min(minval, cost[i][k] - u[i]);
		}
		v[k] = min(minval, cost[k][k]);
		minval = oo;
		for (int i = 0; i <= k; ++i) {
			minval = min(minval, cost[k][i] - v[i]);
		}
		u[k] = minval;
		for (int i = 0; i <= c_dim; ++i) slack[i] = oo;
		while (1) {
			for (int j = 0; j <= c_dim; ++j) visrow[j] = 0;
			for (int j = 0; j <= c_dim; ++j) viscol[j] = 0;

			if (dfs(k, k+1, cost, rowsol, colsol, u, v)) break;
			double theta = oo;
			for (int i = 0; i <= k; ++i) {
				if (!viscol[i]) {
					theta = min(theta, slack[i]);
				}
			}
			for (int j = 0; j <= k; ++j) {
				visrow[j] ? u[j] += theta : u[j] -= theta; //Adjust the dual variable
				viscol[j] ? v[j] -= theta : v[j] += theta;
			}
		}
	}
	double costval = 0;
	for (int i = 0; i < c_dim; ++i) {
		//if (colsol[i] != -1) {
			//costval += cost[colsol[i]][i];
		//}
		costval += (u[i] + v[i]);
	}
	return costval;
}

double hungarian(int dimc, double** cost, int* rowsol, int* colsol, double* u, double* v)
{
	for (int i = 0; i < dimc; ++i) {
		//rowdul[i] = 0;
		u[i] = 0;
	}
	for (int i = 0; i < dimc; ++i) {
		double minval = oo;
		for (int j = 0; j < dimc; ++j) {
			minval = min(minval, cost[j][i]);
		}
		v[i] = minval;
	}
	//memset(match, -1, sizeof match);
	for (int i = 0; i < dimc; ++i) colsol[i] = -1;
	for (int r = 0; r < dimc; ++r) {
		//if (rowsol[r] == -1) continue;
		//memset(slack, oo, sizeof slack);
		for (int i = 0; i < dimc; ++i) slack[i] = oo;
		while (1) {
			//memset(visrow, 0, sizeof visrow);
			//memset(viscol, 0, sizeof viscol);
			for (int i = 0; i < dimc; ++i) visrow[i] = 0;
			for (int i = 0; i < dimc; ++i) viscol[i] = 0;
			if (dfs(r, dimc, cost, rowsol, colsol, u, v)) break;
			double theta = oo;
			for (int i = 0; i < dimc; ++i) {
				if (!viscol[i]) {
					theta = min(theta, slack[i]);
				}
			}
			for (int i = 0; i < dimc; ++i) {
				visrow[i] ? u[i] += theta : u[i] -= theta;
				viscol[i] ? v[i] -= theta : v[i] += theta;
				/*
				if (visrow[i]) {
					rowdul[i] += theta;
				}
				else {
					rowdul[i] -= theta;
				}
				if (viscol[i]) {
					coldul[i] -= theta;
				}
				else {
					coldul[i] += theta;
				}*/
			}
		}
	}
	double costval = 0;

	for (int i = 0; i < dimc; ++i) {
		/*if (match[i] != -1) {
			costval += cost[match[i]][i];
		}*/
		costval += (u[i] + v[i]);
	}
	/*for (int i = 0; i < dim; ++i) {
		mc.match[i] = match[i];
		mc.u[i] = rowdul[i];
		mc.v[i] = coldul[i];
	}*/
	return costval;
}

bool dfs(int r, int dimx, double** cost, int* rowsol, int* colsol, double* u, double* v) {
	visrow[r] = true;
	for (int i = 0; i < dimx; ++i) {
		if (viscol[i]) continue;
		double d = 0.5 * fabs(cost[r][i] - u[r] - v[i]);
		if (d <= eps) {
			viscol[i] = true;
			if (colsol[i] == -1 || dfs(colsol[i], dimx, cost, rowsol, colsol, u, v)) {
				colsol[i] = r, rowsol[r] = i;
				return true;
			}
		}
		else {
			slack[i] = min(slack[i], d);
		}
	}
	return false;
}

double hungarian2(double** cost, int dimx)
{
	dim = dimx;
	for (int i = 0; i < dim; ++i) {
		rowdul[i] = 0;
	}
	for (int i = 0; i < dim; ++i) {
		double minval = oo;
		for (int j = 0; j < dim; ++j) {
			minval = min(minval, cost[j][i]);
		}
		coldul[i] = minval;
	}
	for (int i = 0; i < dimx; ++i) match[i] = -1;
	for (int i = 0; i < dim; ++i) {
		for (int i = 0; i < dim; ++i) slack[i] = oo;
		//memset(slack, oo, sizeof slack);
		while (1) {
			for (int i = 0; i < dimx; ++i) visrow[i] = 0;
			for (int i = 0; i < dimx; ++i) viscol[i] = 0;
			if (dfs2(i, cost)) break;
			double theta = oo;
			for (int i = 0; i < dim; ++i) {
				if (!viscol[i]) {
					theta = min(theta, slack[i]);
				}
			}
			for (int i = 0; i < dim; ++i) {
				visrow[i] ? rowdul[i] += theta : rowdul[i] -= theta;
				viscol[i] ? coldul[i] -= theta : coldul[i] += theta;
			}
		}
	}
	double costval = 0;
	for (int i = 0; i < dim; ++i) {
			costval += (rowdul[i] + coldul[i]);
	}
	return costval;
}

bool dfs2(int r, double** cost) {
	visrow[r] = true;
	for (int i = 0; i < dim; ++i) {
		if (viscol[i]) continue;
		double d = 0.5 * fabs(cost[r][i] - rowdul[r] - coldul[i]);
		if (d <= eps) {
			viscol[i] = true;
			if (match[i] == -1 || dfs2(match[i], cost)) {
				match[i] = r;
				return true;
			}
		}
		else {
			slack[i] = min(slack[i], d);
		}
	}
	return false;
}

bool bfs(int r, int dim, double** cost, int* rowsol, int* colsol, double* u, double* v) {
	head = tail = 0;
	q[head] = r;
	int x;
	double d;
	for (int i = 0; i < dim; ++i) fau[i] = fav[i] = -1;
	while (head <= tail) {
		x = q[head++];
		visrow[x] = true;
		for (int i = 0; i < dim; ++i) {
			if (viscol[i]) continue;
			d = 0.5 * fabs(cost[x][i] - u[x] - v[i]);
			if (d < eps) {
				viscol[i] = true;
				fau[i] = x;
				if (colsol[i] == -1) {//iÎ´Æ¥Åä
					int u = i;
					int fv;
					while (u != -1) {
						fv = fau[u];
						colsol[u] = fv;
						u = fav[fv];
					}
					return true;
				}
				else {
					q[++tail] = colsol[i];
					fav[colsol[i]] = i;
				}
			}
			else {
				slack[i] = min(slack[i], d);
			}
		}
	}
	return false;
}