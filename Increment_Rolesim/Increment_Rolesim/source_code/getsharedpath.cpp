#include "getsharedpath.h"
#include <string>
#include <ctime>
#include <iostream>
using namespace std;
#define oo 10000
pathnode* getsharepath(vector<vector<int>>& G) {
    int n = G.size();
    vector<vector<int>> sets;   //sets表示各个节点的邻居节点集合以及空集
    sets.push_back({});
    for (int i = 0; i < n; ++i) {
        sets.push_back(G[i]);
    }

    vector<vector<int>> aftermapsets;
    
    clock_t start, end;
    
    vector<vector<int>> W = generategraph(sets, aftermapsets);
    purning(aftermapsets, W, n + 1);

    vector<int> ids;
    for (int i = 0; i <= n; ++i)
        ids.push_back(i);
    start = clock();
    pathnode* root = approximate_steiner_tree(W, ids, aftermapsets);
    end = clock();
    cout << "get steiner tree: ";
    cout << (end - start) / 1000.0 << endl;

    return root;
}

vector<vector<int>> generategraph(vector<vector<int>>& sets,  vector<vector<int>>& aftersets)
{
    //int n = sets.size();
    clock_t start, end;
    start = clock();
    vector<vector<int>> newsets = mapreduce(sets, aftersets);
    end = clock();
    cout << "get intersection: ";
    cout << (end - start) / 1000.0 << endl;

    int n = newsets.size();
    vector<vector<int>> G(n, vector<int>(n, 0));

    start = clock();
    for (int i = 0; i < n; ++i) {
        G[i][i] = oo;
        for (int j = i + 1; j < n; ++j) {
            vector<int> d = distince2(newsets[i], newsets[j]);
            G[i][j] = d[0], G[j][i] = d[1];
            //G[i][j] = G[j][i] = distince(newsets[i], newsets[j]);
        }
    }

    end = clock();
    cout << "get intersection w: ";
    cout << (end - start) / 1000.0 << endl;

    return G;
}

vector<vector<int>> mapreduce(vector<vector<int>>& sets, vector<vector<int>>& aftersets) //可能会导致重复出现初始集合的元素
{
    vector<vector<int>> newsets;
    unordered_map<int, string> mapped;
    for (int i = 0; i < sets.size(); ++i) {
        string tag = std::to_string(i);
        tag += '-';
        for (int j = 0; j < sets[i].size(); ++j) {
            mapped[sets[i][j]] += tag;
        }
    }
    unordered_map<string, vector<int>> remapped; //reverse
    unordered_map<int, string>::iterator it;
    for (it = mapped.begin(); it != mapped.end(); ++it) {
        if (it->second.size() != 0) {
            remapped[it->second].push_back(it->first);
        }
    }
    unordered_map<string, vector<int>>::iterator it2;
    for (it2 = remapped.begin(); it2 != remapped.end(); ++it2) {
        string s = it2->first;
        int cnt = 0;
        for (int i = 0; i < s.size(); ++i) {
            if (s[i] == '-') {
                ++cnt;
            }
        }
        if (cnt > 1 && cnt < 2)
            newsets.push_back(it2->second);
    }
    /*
    int lens = sets.size(), lenns = newsets.size();
    int k = (lenns + lens) / lens;
    int cnt = 0, cnt2 = 0;
    vector<vector<int>> resset;
    for (int i = 0; i < (lenns + lens); ++i) {
        if (i % k == 0 && cnt < lens) {
            resset.push_back(sets[cnt]);
            ++cnt;
        }
        else {
            resset.push_back(newsets[cnt2]);
            ++cnt2;
        }
    }*/
    vector<vector<int>> resset = sets;
    for (int i = 0; i < newsets.size(); ++i) {
        resset.push_back(newsets[i]);
    }
    aftersets = resset;
    return resset;
}

int distince(vector<int>& set1, vector<int>& set2)
{
    int n = set1.size();
    int m = set2.size();

    unordered_map<int, int> nums;
    for (int i = 0; i < n; ++i) {
        ++nums[set1[i]];
    }

    for (int i = 0; i < m; ++i) {
        ++nums[set2[i]];
    }

    int sum = 0;
    unordered_map<int, int>::iterator it;
    for (it = nums.begin(); it != nums.end(); ++it) {
        if (it->second == 2) {
            ++sum;
        }
    }
    return (n - sum) + (m - sum);
}

vector<int> distince2(vector<int>& set1, vector<int>& set2)
{
    int n = set1.size();
    int m = set2.size();
    unordered_map<int, int> nums;
    for (int i = 0; i < n; ++i) {
        ++nums[set1[i]];
    }
    for (int i = 0; i < m; ++i) {
        ++nums[set2[i]];
    }
    int sum = 0;
    unordered_map<int, int>::iterator it;
    for (it = nums.begin(); it != nums.end(); ++it) {
        if (it->second == 2) {
            ++sum;
        }
    }
    vector<int> ans(2, 0);
    if (sum < n && n <= min(2 * sum, m)) {
        ans[0] = n + m - 2 * sum;
    }
    else if (sum == n) {
        ans[0] = m - n;
    }
    else {
        ans[0] = oo;
    }

    if (sum < m && m <= min(2 * sum, n)) {
        ans[1] = n + m - 2 * sum;
    }
    else if (sum == m) {
        ans[1] = n - m;
    }
    else {
        ans[1] = oo;
    }
    return ans;
}

void purning(vector<vector<int>>& graph, vector<vector<int>>& W, int L)
{
    int n = graph.size();
    //Purning Rule 2
    //int R1 = 0, R2 = 0, R3 = 0;
    for (int j = 0; j < L; ++j) {
        for (int i = 1; i < n; ++i) {
            if (W[i][j] == oo)continue;
            if (W[i][j] <= W[0][j]) {
                W[0][j] = oo;
                //++R2;
                break;
            }
        }
    }

    //Purning Rule 1
    for (int i = 1; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (W[i][j] == oo) continue;
            if (W[i][j] >= graph[j].size() && graph[i].size() > 0) {
                W[i][j] = oo;
                //++R1;
            }
        }
    }

    //Purning Rule 3
    vector<int> minval(n, oo);
    for (int j = 0; j < n; ++j) {
        for (int i = L; i < n; ++i) {
            if (W[i][j] != oo && W[i][j] < minval[j]) {
                minval[j] = W[i][j];
            }
        }
    }
    for (int j = 0; j < n; ++j) {
        for (int i = L; i < n; ++i) {
            if (W[i][j] == oo) continue;
            if (minval[j] < W[i][j]) W[i][j] = oo;
        }
    }
    /*cout << "R1 = " << R1 << endl;
    cout << "R2 = " << R2 << endl;
    cout << "R3 = " << R3 << endl;*/
}

