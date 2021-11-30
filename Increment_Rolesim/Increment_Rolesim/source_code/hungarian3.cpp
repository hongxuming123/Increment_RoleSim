#include "hungarian3.h"

#include <sstream>
using namespace std;

//float* weights;
//float* flagX;
//float* flagY;
//char* usedX;
//char* usedY;
//int* matchX;
//int* matchY;

void KM::init(float* data, int m, int n) {
	N = max(m, n);
	front = m;
	back = n;
	weights = new float[N * N];
	flagX = new float[N];
	flagY = new float[N];
	usedX = new char[N];
	usedY = new char[N];
	matchX = new int[N];
	matchY = new int[N];
	max_w = 0.0;

	constructMatrix(data, m, n);
}

void KM::init(vector<vector<double>>& costmatrix) {
	N = costmatrix.size();
	front = N;
	back = N;
	weights = new float[N * N];
	flagX = new float[N];
	flagY = new float[N];
	usedX = new char[N];
	usedY = new char[N];
	matchX = new int[N];
	matchY = new int[N];
	max_w = 0.0;

	constructMatrix(costmatrix);
}


void KM::del() {
	delete[] weights;
	delete[] flagX;
	delete[] flagY;
	delete[] usedX;
	delete[] usedY;
	delete[] matchX;
	delete[] matchY;
}

KM::~KM() {
	// The destruction will not be called in tp_free(Python)
	del();
	//cout << "Destruct KM" << endl;
}

void KM::constructMatrix(float* data, int m, int n) {
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			weights[i * N + j] = data[i * n + j];
		}
		for (int j = n; j < N; ++j) {
			weights[i * N + j] = 0.0;
		}
	}
	for (int i = m; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			weights[i * N + j] = 0.0;
		}
	}

	for (int i = 0; i < N; ++i) {
		flagX[i] = -1e7;
		flagY[i] = 0;
		for (int j = 0; j < N; ++j) {
			flagX[i] = max(flagX[i], weights[i * N + j]);
		}
		usedX[i] = 0;
		usedY[i] = 0;
		matchX[i] = -1;
		matchY[i] = -1;
	}
}

void KM::constructMatrix(vector<vector<double>>& costmatrix) {
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			weights[i * N + j] = (float)costmatrix[i][j];
		}
	}
	for (int i = 0; i < N; ++i) {
		flagX[i] = -1e7;
		flagY[i] = 0;
		for (int j = 0; j < N; ++j) {
			flagX[i] = max(flagX[i], weights[i * N + j]);
		}
		usedX[i] = 0;
		usedY[i] = 0;
		matchX[i] = -1;
		matchY[i] = -1;
	}
}



bool KM::dfs(int v) {
	usedX[v] = 1;
	float* data = weights + v * N;
	for (int i = 0; i < N; ++i) {
		float k = flagX[v] + flagY[i];
		bool eq = (k - data[i] < 1e-6 && data[i] - k < 1e-6);
		if (!usedY[i] && eq) {
			usedY[i] = 1;
			if (matchY[i] == -1 || dfs(matchY[i])) {
				matchY[i] = v;
				matchX[v] = i;
				return true;
			}
		}
	}
	return false;
}

void KM::compute() {
	for (int i = 0; i < N; ++i) {
		while (true) {
			for (int j = 0; j < N; ++j) {
				usedX[j] = 0;
				usedY[j] = 0;
			}
			if (dfs(i))
				break;

			float d = 1e7;
			for (int j = 0; j < N; ++j) {
				if (!usedX[j])
					continue;
				for (int k = 0; k < N; ++k) {
					if (!usedY[k]) {
						d = min(d, flagX[j] + flagY[k] - weights[j * N + k]);
					}
				}
			}
			if (d == 0) {
				max_w = -1;
				cout << "No max weight match" << endl;
				return;
			}
			for (int j = 0; j < N; ++j) {
				if (usedX[j])
					flagX[j] -= d;
				if (usedY[j])
					flagY[j] += d;
			}
		}
	}
	for (int i = 0; i < N; ++i) {
		int j = matchX[i];
		if (j >= 0) {
			max_w += weights[i * N + j];
		}
	}
}

vector<int> KM::getMatch(bool front2back) {
	vector<int> matches;
	if (front2back) {
		for (int i = 0; i < front; ++i)
			matches.push_back(matchX[i]);
	}
	else {
		for (int i = 0; i < back; ++i)
			matches.push_back(matchY[i]);
	}
	return matches;
}

ostream& operator<<(ostream& os, KM& km) {
	for (int i = 0; i < km.N; ++i) {
		for (int j = 0; j < km.N; ++j) {
			os << km.weights[i * km.N + j];
			if (j < km.N - 1)
				os << ',';
			else
				os << endl;
		}
	}
	return os;
}