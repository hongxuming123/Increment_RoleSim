#ifndef __PATHTREE_H_
#define __PATHTREE_H_

#include <vector>
#include <set>
using namespace std;
typedef struct pathnode {
	vector<pathnode*> nexts;
	int val; //ָʾ��˹̹�����еĽڵ��ţ���ͼ�еı��ƫ��Ϊ1 ������˹̹���� 1-->4 == 5 ��ʾͼ��0-->3��ת�ƴ�����5
	pathnode* fa;

	set<int> commom;
	set<int> diff;

	int* nodes_seq;
	int* index; //common elem index,get fahter's common elem index in cost matrix
	int common_size;
	int diff_size;

	pathnode(int v, pathnode* f) :val(v), fa(f) {}
}pathnode;

#endif // !__PATHTREE_H_
