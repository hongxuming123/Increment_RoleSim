#pragma once
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

class KM {

public:
	//	KM(float* data, int m, int n);
	KM(float* data, int m, int n) {
		init(data, m, n);
	}
	KM(vector<vector<double>>& costmatrix) {
		init(costmatrix);
	}

	~KM();

	int N;
	int front;
	int back;
	int* matchX;
	int* matchY;
	float* weights;

	void init(float* data, int m, int n);
	void init(vector<vector<double>>& costmatrix);
	void del();
	void compute();
	float maxWeight() {
		return max_w;
	}
	vector<int> getMatch(bool front2back = true);

private:

	float max_w;
	float* flagX;
	float* flagY;
	char* usedX;
	char* usedY;

	void constructMatrix(float* data, int m, int n);
	void constructMatrix(vector<vector<double>>& costmatrix);
	bool dfs(int v);
};