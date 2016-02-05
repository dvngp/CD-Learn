#ifndef BE_H_
#define BE_H_

#include <vector>
#include "GM.h"
struct BE
{
	Double1 pe;
	vector<vector<Function*> > bucketfuncs;
	vector<vector<vector<Function*> > > buckets;
	
	vector<int> order;
	double currPRVal;
	double oldPRVal;
	double diffPR;
	double addPR;
	int ibound;
	//self-joins
	vector<vector<vector<bool> > > dirtyBitArrays;
	BE(){

	}
	void init(vector<Variable*>& variables, vector<Function*>& functions, vector<int>& order);
	void initandstore(vector<Variable*>& variables, vector<Function*>& functions, vector<int>& order);
	//bool BESample(vector<Variable*> &variables,vector<Function*> &functions);
	void updateBuckets(int flippedpred,bool isSelfJoin = false);
	void revertBuckets(bool change);
	//void setBuckets(vector<vector<int> > changedentries,vector<vector<double> > vals,double PR);
	void propagateall();
	vector<vector<vector<int> > > hashtable;
	vector<vector<vector<int> > > changedindexes;
	vector<vector<vector<Variable*> > > margVarList;

};
#endif
