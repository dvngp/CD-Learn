#include "BE.h"

void createFunction(Function& function, Function& out)
{
	out.variables().clear();
	out.table().clear();
	if(function.variables().empty() && function.table().empty())
	{
		out.table()=vector<Double>(1);
		out.table()[0]=Double(1.0);
		return;
	}
	//Remove all variables that have been assignend
	for(int i=0;i<function.variables().size();i++)
	{
		if(function.variables()[i]->value()==INVALID_VALUE)
			out.variables().push_back(function.variables()[i]);
	}
	if(out.variables().empty())
	{
		out.table()=vector<Double> (1);
		out.table()[0]=function.table()[Variable::getAddress(function.variables())];
		return;
	}
	if(out.variables().size()==function.variables().size())
	{
		out.table()=function.table();
	}
	else
	{
		int domain_size=Variable::getDomainSize(out.variables());
		out.table()=vector<Double>(domain_size);
		int num_variables=out.variables().size();
		/***
		int g[num_variables];
		int a[num_variables];
		int f[num_variables+1];
		int o[num_variables];
		int m[num_variables];
		*/
		vector<int> g(num_variables);
		vector<int> a(num_variables);
		vector<int> f(num_variables+1);
		vector<int> o(num_variables);
		vector<int> m(num_variables);

		int multiplier=1;
		for(int i=0;i<num_variables;i++)
		{
			a[i]=0;
			g[i]=multiplier;
			f[i]=i;
			o[i]=1;
			m[i]=out.variables()[i]->domain_size();
			multiplier*=m[i];
		}
		f[num_variables]=num_variables;
		//**int h[num_variables];
		vector<int> h(num_variables);
		int func_address=0;
		multiplier=1;
		int k=0;
		for(int i=0;i<function.variables().size();i++)
		{
			if(function.variables()[i]->value()!=INVALID_VALUE)
				func_address+=(multiplier*function.variables()[i]->value());
			else
			{
				assert(k<num_variables);
				h[k++]=multiplier;
			}
			multiplier*=function.variables()[i]->domain_size();
		}
		int address=0;

		while(1)
		{
			out.table()[address]=function.table()[func_address];
			int j=f[0];
			f[0]=0;

			if(j==num_variables) break;
			int old_aj=a[j];
			a[j]=a[j]+o[j];
			if(a[j]==0 || a[j]==(m[j]-1))
			{
				o[j]=-o[j];
				f[j]=f[j+1];
				f[j+1]=j+1;
			}
			address-=(g[j]*old_aj);
			address+=(g[j]*a[j]);
			func_address-=(h[j]*old_aj);
			func_address+=(h[j]*a[j]);
		}
	}
}

Double1 get_norm_const(Function& function)
{
	Double1 norm_const=0.0;
	for(int i=0;i<function.table().size();i++){
		norm_const+=function.table()[i];
	}
	/*for(int i=0;i<function.table().size();i++){
			function.table()[i]/=norm_const;
		}
		*/
	return norm_const;
}

void BE::init(std::vector<Variable*> &variables, std::vector<Function*> &functions, std::vector<int> &order)
{
	//pe=Double(1.0);
	//vector<vector<Function*> > buckets (order.size());
	vector<vector<Function*> > tmp_buckets(order.size());
	this->order = order;
	vector<int> var_in_pos(order.size());
	for(int i=0;i<var_in_pos.size();i++) {
		var_in_pos[order[i]]=i;
	}

	// First put the functions in the proper buckets
	for(int i=0;i<functions.size();i++)
	{
		int pos=order.size();
		//LogFunction* function=new LogFunction(*functions[i]);
		if(functions[i]->variables().empty())
		{
			pe*=functions[i]->table()[0];

			continue;
		}
		//Function* function=new Function();
		//createFunction(*functions[i],*function);
		Function* function = functions[i];
		for(int j=0;j<functions[i]->variables().size();j++)
		{
			if(var_in_pos[functions[i]->variables()[j]->id()] < pos)
				pos=var_in_pos[functions[i]->variables()[j]->id()];
		}
		assert(pos!=(int)order.size());
		tmp_buckets[pos].push_back(function);
	}
	int prval;
	for(int i=0;i<tmp_buckets.size();i++)
	{
		if(tmp_buckets[i].empty())
			continue;

		vector<Variable*> bucket_variables;
		for(int j=0;j<tmp_buckets[i].size();j++)
		{
			do_set_union(bucket_variables,tmp_buckets[i][j]->variables(),bucket_variables,less_than_comparator_variable);
		}
		//cout<<bucket_variables.size()<<" "<<flush;
		//cout<<"buck-vars.size ="<<bucket_variables.size()<<endl;
		//cerr<<"Processing bucket "<<i<<" out of "<<buckets.size()<<endl;


		// Compute variables required for marginalization
		//cerr<<bucket_variables.size()<<endl;
		vector<Variable*> bucket_variable;
		bucket_variable.push_back(variables[order[i]]);
		vector<Variable*> marg_variables;
		do_set_difference(bucket_variables,bucket_variable,marg_variables,less_than_comparator_variable);
		
		Function* function= new Function();
		
		Function::multiplyAndMarginalize(marg_variables,tmp_buckets[i],*function,false);
		
		if(i==tmp_buckets.size()-1) {
			currPRVal = get_norm_const(*function);
		}
		if(function->variables().empty())
		{
			assert((int)function->table().size()==1);
			//pe*=function->table()[0];
			//currPRVal = function->table()[0];
			//delete(function);
			continue;
		}
		
		//prval = get_norm_const(*function);
		//pe*=prval;
		//Put the function in the appropriate bucket
		int pos=order.size();
		//function->print(cout);
		//assert(!function->log_table.empty());

		for(int j=0;j<function->variables().size();j++)
		{
			if(var_in_pos[function->variables()[j]->id()] < pos)
				pos=var_in_pos[function->variables()[j]->id()];
		}
		assert(pos!=(int)order.size());
		assert(pos > i);
		tmp_buckets[pos].push_back(function);
		
	}
	

	/*
	for(int i=0;i<buckets.size();i++){
		for(int j=0;j<buckets[i].size();j++){
			if (buckets[i][j]!=NULL){
				delete(buckets[i][j]);
			}
		}
	}
	buckets.clear();
	*/
	//BESample(variables,functions,order,buckets);
}

void BE::initandstore(std::vector<Variable*> &variables, std::vector<Function*> &functions, std::vector<int> &order)
{
	//pe=Double(1.0);
	//vector<vector<Function*> > buckets (order.size());
	buckets.resize(order.size());
	hashtable.resize(order.size());
	dirtyBitArrays.resize(order.size());
	bucketfuncs.resize(order.size());
	margVarList.resize(buckets.size());
	this->order = order;
	changedindexes.resize(order.size());
	vector<int> var_in_pos(order.size());
	for(int i=0;i<var_in_pos.size();i++) {
		var_in_pos[order[i]]=i;
	}

	// First put the functions in the proper buckets
	for(int i=0;i<functions.size();i++)
	{
		int pos=order.size();
		//LogFunction* function=new LogFunction(*functions[i]);
		if(functions[i]->variables().empty())
		{
			pe*=functions[i]->table()[0];

			continue;
		}
		//Function* function=new Function();
		//createFunction(*functions[i],*function);
		Function* function = functions[i];
		for(int j=0;j<functions[i]->variables().size();j++)
		{
			if(var_in_pos[functions[i]->variables()[j]->id()] < pos)
				pos=var_in_pos[functions[i]->variables()[j]->id()];
		}
		assert(pos!=(int)order.size());
		//check if exceeds ibound
		if(buckets[pos].size()==0) {
			buckets[pos].resize(1);
		}
		vector<Variable*> bucket_variables;
		for(int j=0;j<buckets[pos][buckets[pos].size()-1].size();j++)
		{
			do_set_union(bucket_variables,buckets[pos][buckets[pos].size()-1][j]->variables(),bucket_variables,less_than_comparator_variable);
		}
		if(bucket_variables.size() > ibound) {
			vector<Function*> tmpfuncs;
			tmpfuncs.push_back(function);
			buckets[pos].push_back(tmpfuncs);
		} else {
			buckets[pos][buckets[pos].size()-1].push_back(function);
		}
	}

	int prval;
	currPRVal = 1;
	for(int i=0;i<buckets.size();i++)
	{
		if(buckets[i].empty())
			continue;
		for(int j=0;j<buckets[i].size();j++)
		{
			vector<Variable*> bucket_variables;
			for(int k=0;k<buckets[i][j].size();k++) {
				do_set_union(bucket_variables,buckets[i][j][k]->variables(),bucket_variables,less_than_comparator_variable);
			}
			vector<Variable*> bucket_variable;
			bucket_variable.push_back(variables[order[i]]);
			vector<Variable*> marg_variables;
			do_set_difference(bucket_variables,bucket_variable,marg_variables,less_than_comparator_variable);
			Function* function= new Function();
			bucketfuncs[i].push_back(function);
			vector<int> htable;
			hashtable[i].push_back(htable);
			//Function::multiplyAndMarginalize(marg_variables,buckets[i],*function,false);
			Function::multiplyAndMarginalizeNew(marg_variables,buckets[i][j],*function,hashtable[i][hashtable[i].size()-1]);
			if(i==buckets.size()-1) {
				currPRVal *= get_norm_const(*function);
			}
			//cout<<bucket_variables.size()<<" "<<flush;
			//cout<<"buck-vars.size ="<<bucket_variables.size()<<endl;
			//cerr<<"Processing bucket "<<i<<" out of "<<buckets.size()<<endl;
			if(function->variables().empty())
			{
				assert((int)function->table().size()==1);
				//pe*=function->table()[0];
				//currPRVal = function->table()[0];
				//delete(function);
				continue;
			}
		
			//prval = get_norm_const(*function);
			//pe*=prval;
			//Put the function in the appropriate bucket
			int pos=order.size();
			//function->print(cout);
			//assert(!function->log_table.empty());

			for(int k=0;k<function->variables().size();k++)
			{
				if(var_in_pos[function->variables()[k]->id()] < pos)
					pos=var_in_pos[function->variables()[k]->id()];
			}
			assert(pos!=(int)order.size());
			assert(pos > i);
			int indextoadd=-1;
			vector<Variable*> new_vars;
			for(int k=0;k<buckets[pos].size();k++)
			{
				vector<Variable*> bucket_variables;
				for(int m=0;m<buckets[pos][k].size();m++) {
					do_set_union(bucket_variables,buckets[pos][k][m]->variables(),bucket_variables,less_than_comparator_variable);
				}
				do_set_union(bucket_variables,function->variables(),bucket_variables,less_than_comparator_variable);
				if(bucket_variables.size() <= ibound) {
					indextoadd = k;
					new_vars = bucket_variables;
					break;
				}
			}
			if(indextoadd==-1) {
				vector<Function*> tmp_f;
				tmp_f.push_back(function);
				buckets[pos].push_back(tmp_f);
			} else {
				buckets[pos][indextoadd].push_back(function);
			}
		}
		
	}

	for(int i=0;i<buckets.size();i++)
	{
		changedindexes[i].resize(buckets[i].size());
		vector<Variable*> bucket_variable;
		bucket_variable.push_back(variables[order[i]]);
		for(int j=0;j<buckets[i].size();j++) {
			vector<Variable*> bucket_variables;
			for(int k=0;k<buckets[i][j].size();k++) {
				do_set_union(bucket_variables,buckets[i][j][k]->variables(),bucket_variables,less_than_comparator_variable);
			}
			vector<Variable*> marg_variables;
			do_set_difference(bucket_variables,bucket_variable,marg_variables,less_than_comparator_variable);
			margVarList[i].push_back(marg_variables);
		}
	}

	/*
	for(int i=0;i<buckets.size();i++){
		for(int j=0;j<buckets[i].size();j++){
			if (buckets[i][j]!=NULL){
				delete(buckets[i][j]);
			}
		}
	}
	buckets.clear();
	*/
	//BESample(variables,functions,order,buckets);
}

void BE::updateBuckets(int flippedpred, bool isSelfJoin) {
	oldPRVal = currPRVal;
	//changedentries.resize(buckets.size());
	//newvals.resize(buckets.size());
	bool startpropagation = false;
	for(int i=0;i<buckets.size();i++)
	{
		if(buckets[i].empty())
			continue;
		if(!startpropagation) {
			for(int j=0;j<buckets[i].size();j++) {
				for(int k=0;k<buckets[i][j].size();k++) {
					if(buckets[i][j][k]->predid==flippedpred) {
						startpropagation = true;
						break;
					}
				}
				if(startpropagation)
					break;
			}
		}
		if(startpropagation) {
			if(isSelfJoin) {
				if(dirtyBitArrays[i].size()==0) {
					dirtyBitArrays[i].resize(hashtable[i].size());
					for(int j=0;j<hashtable[i].size();j++) {
						dirtyBitArrays[i][j] = vector<bool>(hashtable[i][j].size());
					}
				}
				for(int j=0;j<buckets[i].size();j++) {
					Function::projectMultiplyAndMarginalizeSJ(margVarList[i][j],buckets[i][j],*bucketfuncs[i][j],
						hashtable[i][j],changedindexes[i][j],dirtyBitArrays[i][j]);
				}
			} else {
				for(int j=0;j<buckets[i].size();j++) {
					Function::projectMultiplyAndMarginalize(margVarList[i][j],buckets[i][j],*bucketfuncs[i][j],
					hashtable[i][j],changedindexes[i][j]);
				}
			}
				//changedentries[i].push_back(changedindexes[i][j]);
				//newvals[i].push_back(bucketfuncs[i]->table()[changedindexes[i][j]]);
			if(i==buckets.size()-1) {
				diffPR = 1;
				addPR = 1;
				for(int j=0;j<changedindexes[i].size();j++) {
					double tmp_sub=0;
					double tmp_add=0;
					for(int k=0;k<changedindexes[i][j].size();k++) {
						tmp_sub += hashtable[i][j][changedindexes[i][j][k]];
						tmp_add += bucketfuncs[i][j]->table()[changedindexes[i][j][k]];
					}
					diffPR *= tmp_sub;
					addPR *= tmp_add;
				}
				currPRVal = currPRVal - diffPR + addPR;
			}
		}
	}
}

/*
//update without ever needing revert back
void BE::updateBuckets(int flippedpred,bool isSelfJoin) {
	oldPRVal = currPRVal;
	bool startpropagation=false;
	for(int i=0;i<buckets.size();i++)
	{
		if(buckets[i].empty())
			continue;
		if(!startpropagation) {
			for(int j=0;j<buckets[i].size();j++) {
				if(buckets[i][j]->predid==flippedpred) {
					startpropagation = true;
					break;
				}
			}
		}
		if(startpropagation) {
			if(isSelfJoin) {
				if(dirtyBitArrays[i].size()==0) {
					dirtyBitArrays[i] = vector<bool>(hashtable[i].size());
				}
				Function::projectMultiplyAndMarginalizeSJ(margVarList[i],buckets[i],*bucketfuncs[i],
					hashtable[i],changedindexes[i],dirtyBitArrays[i]);
			} else {
				Function::projectMultiplyAndMarginalize(margVarList[i],buckets[i],*bucketfuncs[i],hashtable[i],changedindexes[i]);
			}
			diffPR = 0;
			addPR = 0;
			if(i==buckets.size()-1) {
				for(int j=0;j<changedindexes[i].size();j++) {
					diffPR += hashtable[i][changedindexes[i][j]];
					addPR += bucketfuncs[i]->table()[changedindexes[i][j]];
				}
			}
			if(i==buckets.size()-1) {
				//currPRVal = get_norm_const(*bucketfuncs[i]);
				currPRVal = currPRVal - diffPR + addPR;
			}
		}
	}
	for(int i=0;i<buckets.size();i++) {
		for(int j=0;j<changedindexes[i].size();j++) {
			hashtable[i][changedindexes[i][j]] = -1;
		}
		changedindexes[i].clear();
	}
}


void BE::setBuckets(vector<vector<int> > changedentries,vector<vector<double> > vals,double PR) {
	currPRVal = PR;
	for(int i=0;i<buckets.size();i++)
	{
		if(buckets[i].empty())
			continue;
		for(int j=0;j<changedentries[i].size();j++) {
			bucketfuncs[i]->table()[changedentries[i][j]] = vals[i][j];
		}
	}
	double d1 = get_norm_const(*bucketfuncs[bucketfuncs.size()-1]);
}
*/
void BE::revertBuckets(bool change) {
	if(change)
		currPRVal = oldPRVal;
	//exchange hashtable values
	for(int i=0;i<buckets.size();i++) {
		for(int j=0;j<changedindexes[i].size();j++) {
			for(int k=0;k<changedindexes[i][j].size();k++) {
				if(change)
					bucketfuncs[i][j]->table()[changedindexes[i][j][k]] = hashtable[i][j][changedindexes[i][j][k]];
				hashtable[i][j][changedindexes[i][j][k]] = -1;
			}
			changedindexes[i][j].clear();
		}
	}
}


void BE::propagateall() {
	for (int i = 0; i<buckets.size(); i++)
	{
		if (buckets[i].empty())
			continue;
		for (int j = 0; j<buckets[i].size(); j++) {
			Function::multiplyAndMarginalizeStored(margVarList[i][j], buckets[i][j], *bucketfuncs[i][j]);
		}
		if (i == buckets.size() - 1) {
			currPRVal = get_norm_const(*bucketfuncs[i][bucketfuncs[i].size()-1]);
		}
	}
}