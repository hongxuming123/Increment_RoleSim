#ifndef __ROLESIM_H_
#define __ROLESIM_H_

#include "matchinfo.h"
#include <vector>
#include "pathtree.h"
#include <unordered_map>
using namespace std;

vector<vector<double>> Incremental_Rolesim(vector<vector<int>>& G, double beta, int k);

void Iterative(pathnode* root, vector<vector<int>>& nei, vector<vector<double>>& Ha, vector<vector<double>>& Hb, double beta);
void dfs(int u, pathnode* root, vector<vector<int>>& nei, unordered_map<int, Mi*>& forlapjv, vector<vector<double>>& Ha, vector<vector<double>>& Hb, double beta);
void freeupspace(pathnode* root);

void dfs2(pathnode* root);

int create(Mi* mfptr, Mi* mcptr, pathnode* curent_node, int uid, vector<vector<int>>& nei, vector<vector<double>>& Ha);
int create_v2(Mi* mfptr, Mi* mcptr, pathnode* curent_node, int uid, vector<vector<int>>& nei, vector<vector<double>>& Ha);
int create_first_level(Mi* mcptr,pathnode* curent_node,int uid,vector<vector<int>>& nei, vector<vector<double>>& Ha);
void dfs_init_v_common_index(pathnode* root);


#endif // !__ROLESIM_H_
