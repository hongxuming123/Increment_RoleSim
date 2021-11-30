#include "MWM.h"
#include "lap.h"
#include "dymatch.h"
#include "hungarian1.h"
#include "hungarian3.h"

MWM::MWM(int maxn) {
	this->dim = 0;
	MAXN = maxn;
	costmatrix = new double* [MAXN];
	for (int i = 0; i < MAXN; ++i) {
		costmatrix[i] = new double[MAXN];
	}
	u = new double[MAXN];
	v = new double[MAXN];
	colsol = new int[MAXN];
	rowsol = new int[MAXN];
}


double MWM::getmincost(vector<vector<double>>& G) {
	int dim = G.size();
	for (int i = 0; i < dim; ++i) {
		for (int j = 0; j < dim; ++j) {
			costmatrix[i][j] = G[i][j];
		}
	}
	double cost;
	cost = lap(dim, costmatrix, rowsol, colsol, u, v);
	return cost;
}

double MWM::getmaxweight(vector<vector<double>>& G) {
	int dim = G.size();
	for (int i = 0; i < dim; ++i) {
		for (int j = 0; j < dim; ++j) {
			costmatrix[i][j] = 1 - G[i][j];
		}
	}
	double  weight;
	weight = dim - lap(dim, costmatrix, rowsol, colsol, u, v);
	return weight;
}

double MWM::getmaxweight_by_hungarian(vector<vector<double>>& G) {
	int dim = G.size();
	for (int i = 0; i < dim; ++i) {
		for (int j = 0; j < dim; ++j) {
			costmatrix[i][j] = 1 - G[i][j];
		}
	}
	double  weight;
	weight = dim - hungarian2(costmatrix,dim);
	return weight;
}

double MWM::getmaxweight_by_hungarian2(vector<vector<double>>& G)
{
	int dim = G.size();
	vector<int> rm, cm;
	for (int i = 0; i < dim; ++i) {
		for (int j = 0; j < dim; ++j) {
			G[i][j] = 1 - G[i][j];
		}
	}
	return dim - MinCostMatching(G, rm, cm);
}

double MWM::getmaxweight_by_hungarian3(vector<vector<double>>& G)
{
	KM km(G);
	km.compute();
	return km.maxWeight();
}
