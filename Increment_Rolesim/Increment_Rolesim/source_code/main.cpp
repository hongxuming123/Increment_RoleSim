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
using namespace std;

double getrand0_1() {
	double k = rand() % 433512;
	return k / (double)433512;
}

vector < vector<double>> getG(int n) {
	vector<vector<double>> G(n, vector<double>(n, 0));
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			G[i][j] = getrand0_1();
		}
	}
	return G;
}

void printview(vector<vector<double>>& Rsim) {
	cout.precision(2);
	int n = Rsim.size();
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			//cout << Rsim[i][j] << " ";
			printf("%4.2f ", Rsim[i][j]);
		}
		cout << endl << endl;
	}
	cout << endl;
}

vector<vector<int>> Get_Random_Graph(int n) {
	int tot = 0;
	vector<vector<int>> G(n, vector<int>(n, 0));
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			if (i == j) continue;
			//G[i][j] = rand() % 2;
			//if (G[i][j])++tot;
			int maxval = rand() % (10*n);
			if (maxval < 13)
				G[i][j] = 1;
			else
				G[i][j] = 0;
			if (G[i][j]) ++tot;
		}
	}
	cout << "edges = " << tot << endl;
	cout << "nodes = " << n << endl;
	cout << "degree = " << tot / n << endl;
	return G;
}

vector<vector<int>> Load_Graph(int& n, int&m ) {
	int maxval = 0;

	int nodes, edges;
	ofstream ofs;
	FILE* fn1;
	//fopen_s(&fn1,"D:\\datasets\\CA-GrQc.txt", "r");
	//fopen_s(&fn1, "D:\\datasets\\soc-sign-bitcoinalpha.txt", "r");
	fopen_s(&fn1, "D:\\datasets\\soc-sign-bitcoinotc.txt", "r");
	//fopen_s(&fn1, "D:\\datasets\\facebook_combined.txt", "r");


	fscanf_s(fn1, "%d\t%d\n", &nodes, &edges);

	vector<vector<int>> G(nodes, vector<int>());
	int i, j;


	while (fscanf_s(fn1, "%d\t%d\n", &i, &j) != EOF) {
		G[j].push_back(i);
		//maxval = max(maxval, max(i, j));
	}
	cout << "# of nodes:  " << nodes << endl;
	cout << "# of edges:  " << edges << endl;
	cout << "degree:      " << edges/nodes << endl;
	cout << "maxval = " << maxval << endl;

	n = nodes;
	m = edges;

	fclose(fn1);
	return G;
}

int main()
{
	srand(time(nullptr));
	int n;
	int m;
	vector<vector<int>> gg = {
		//a,b,c,d,e,f,g
		{},//a
		{0,2,7},//b
		{},//c
		{1,2,7},//d
		{},//e
		{0,1,2,3,4},//f
		{0,1,2,3,4,8},//g
		{6},
		{1,2,3,4,5}
	};
	vector<vector<int>> g = Load_Graph(n, m);
	vector<vector<double>> S1, S2;
	clock_t start, end;
	cout << endl << endl;

	
	cout << "Incremental time = \n";
	start = clock();
	S1 = Incremental_Rolesim(g, 0.4, 5);  //(1 - beta) * weight / max(Nu,Nv) + beta;
	end = clock();
	cout << (end - start) / 1000.0 << endl;
	
	
	/*cout << "original time = ";	
	start = clock();
	S2 = OriginalRoleSim(g, 0.4, 5);
	end = clock();
	cout << (end - start) / 1000.0 << endl;*/
	

	
	/*double err = 0.0;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			err += fabs(S1[i][j] - S2[i][j]);
	
	cout << "accuracy:  " << err << endl;*/
	

	return 0;
}