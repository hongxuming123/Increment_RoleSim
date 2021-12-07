#include <iostream>
#include "dymatch.h"
#include "lap.h"
#include <vector>
#include "increament_rolesim.h"
#include <stdio.h>
#include <ctime>
#include <fstream>
#include "rolesim.h"
#include <algorithm>
#include "getsharedpath.h"
using namespace std;

vector<vector<int>> Load_Graph(int& n, int& m) {
	int nodes, edges;
	ofstream ofs;
	FILE* fn1;
	fopen_s(&fn1, "D:\\datasets\\CA-GrQc.txt", "r");
	//fopen_s(&fn1, "D:\\datasets\\soc-sign-bitcoinalpha.txt", "r");
	//fopen_s(&fn1, "D:\\datasets\\soc-sign-bitcoinotc.txt", "r");
	//fopen_s(&fn1, "D:\\datasets\\facebook_combined.txt", "r");
	//fopen_s(&fn1, "D:\\datasets\\CA-HepTh.txt", "r");

	fscanf_s(fn1, "%d\t%d\n", &nodes, &edges);

	vector<vector<int>> G(nodes, vector<int>());
	int i, j;
	while (fscanf_s(fn1, "%d\t%d\n", &i, &j) != EOF) {
		G[j].push_back(i);
	}
	cout << "# of nodes:  " << nodes << endl;
	cout << "# of edges:  " << edges << endl;
	cout << "degree:      " << edges / nodes << endl;
	n = nodes;
	m = edges;
	fclose(fn1);
	return G;
}

int main()
{
	srand(time(nullptr));
	int n, m;
	vector<vector<int>> g = Load_Graph(n, m);

	vector<vector<double>> S1, S2;
	clock_t start, end;
	
	count_pruning_edges(g);

	cout << "Incremental time = ";
	start = clock();
	S1 = Incremental_Rolesim(g, 0.8, 5);  //(1 - beta) * weight / max(Nu,Nv) + beta;
	end = clock();
	cout << (end - start) / 1000.0 << endl << endl;

	cout << "original time = ";	
	start = clock();
	S2 = OriginalRoleSim(g, 0.8, 5);
	end = clock();
	cout << (end - start) / 1000.0 << endl << endl;
	
	double err = 0.0;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			err += fabs(S1[i][j] - S2[i][j]);
	cout << "accuracy:  " << err  << endl;

	return 0;
}
