#ifndef __GETSHAREDPATH_H_
#define __GETSHAREDPATH_H_

#include <vector>
#include "approximate_steiner_tree.h"
using namespace std;


vector<vector<int>> generategraph(vector<vector<int>>& sets, vector<vector<int>>& aftersets);


vector<vector<int>> mapreduce(vector<vector<int>>& sets, vector<vector<int>>& aftersets);
int distince(vector<int>& set1, vector<int>& set2);
vector<int> distince2(vector<int>& set1, vector<int>& set2);
void purning(vector<vector<int>>& graph, vector<vector<int>>& W, int n);

//pathnode* getsharepath(vector<vector<int>>& G, vector<vector<int>>& nei);
pathnode* getsharepath(vector<vector<int>>& G);



#endif // __GETSHAREDPATH_H_
