#ifndef __ROLESIM1_H_
#define __ROLESIM1_H_

#include <vector>
#include "MWM.h"
using namespace std;

vector<vector<double>> OriginalRoleSim(vector<vector<int>>& G, double beta, int k);
void IterateRoleSim(vector<vector<double>>& Ha, vector<vector<double>>& Hb, vector<vector<int>>& nei, MWM& wei, double beta);


#endif // !__ROLESIM_H_
