#include "rolesim.h"
#include "lap.h"
#include "MWM.h"
#include <algorithm>
using namespace std;


vector<vector<double>> OriginalRoleSim(vector<vector<int>>& G, double beta, int k)
{//G为有向图 G[u][v]==1 表示 v是u的入邻居
	int n = G.size();
	vector<vector<double>> H1(n, vector<double>(n, 1.0));
	vector<vector<double>> H2(n, vector<double>(n, 1.0));
	MWM wei(n);
	/*
	vector<vector<int>> I(n, vector<int>());
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (i == j) continue;
			if (G[i][j] == 1) {
				I[i].push_back(j);
			}
		}
	}*/
	for (int i = 0; i < k; ++i) {
		if (i % 2 == 0) {
			IterateRoleSim(H1, H2, G, wei, beta);
		}
		else {
			IterateRoleSim(H2, H1, G, wei, beta);
		}
	}
	if (k % 2 == 1) {
		return H2;
	}
	return H1;
}

void IterateRoleSim(vector<vector<double>>& Ha, vector<vector<double>>& Hb, vector<vector<int>>& nei, MWM& wei, double beta)
{//Ha --> Hb
	int n = Ha.size();
	for (int u = 0; u < n; ++u) {
		for (int v = 0; v < n; ++v) {
			//if (u == v) {
				//Hb[u][v] = 1;
				//continue;
			//}
			int Nu = nei[u].size(), Nv = nei[v].size();
			int maxn = max(Nu, Nv);
			if (maxn == 0) {
				Hb[u][v] = beta;
				continue;
			}
			vector<vector<double>> weight(maxn, vector<double>(maxn, 0));
			for (int i = 0; i < nei[u].size(); ++i) {
				for (int j = 0; j < nei[v].size(); ++j) {
					weight[i][j] = Ha[nei[u][i]][nei[v][j]];
				}
			}
			double maxweight = wei.getmaxweight_by_hungarian3(weight);
			//double maxweight = wei.getmaxweight(weight);
			Hb[u][v] = (1 - beta) * maxweight / maxn + beta;
		}
	}
}