#include "Function.h"
//#include "myRandom.h"

#include <cassert>
void Function::reduceDomains() {
	int new_num_values = Variable::getDomainSize(variables_);
	vector<Double> new_table(new_num_values);
	for (int i = 0; i < new_num_values; i++) {
		Variable::setAddress(variables_, i);
		int address = 0;
		int multiplier = 1;
		for (int j = 0; j < variables_.size(); j++) {
			address += (multiplier
					* variables_[j]->mapping[variables_[j]->addr_value()]);
			multiplier *= (int) variables_[j]->old_domain.size();
		}
		new_table[i] = table_[address];
	}
	table_ = new_table;
}
void Function::removeEvidence() {

	vector<Variable*> other_variables;
	for (int i = 0; i < variables_.size(); i++)
		if (variables_[i]->value() == INVALID_VALUE)
			other_variables.push_back(variables_[i]);

	int other_num_values = Variable::getDomainSize(other_variables);
	vector<Double> new_table(other_num_values);
	for (int j = 0; j < other_num_values; j++) {
		Variable::setAddress(other_variables, j);
		int entry = Variable::getAddress(variables_);
		new_table[j] = table()[entry];
	}
	variables_ = other_variables;
	table_ = new_table;
}
void Function::normalize() {
	Double norm_const = 0.0;
	for (int i = 0; i < table_.size(); i++) {
		norm_const += table_[i];
	}
	for (int i = 0; i < table_.size(); i++) {
		table_[i] /= norm_const;
	}
}


void Function::multiplyAndMarginalize(vector<Variable*>& marg_variables_, vector<Function*>& functions, Function& new_func,bool to_normalize)
{
	if (functions.empty()) return;
	vector<Variable*> variables;
	vector<Variable*> marg_variables;
	for(int i=0;i<functions.size();i++)
	{
		do_set_union(variables,functions[i]->variables(),variables,less_than_comparator_variable);
	}
	do_set_intersection(variables,marg_variables_,marg_variables,less_than_comparator_variable);
	//if (variables.empty()) return;
	new_func.variables()=marg_variables;
	new_func.table()=vector<Double>(Variable::getDomainSize(marg_variables),0);
	//cout<<"num-marg-vars = "<<marg_variables.size()<<" ";
	//if(variables.empty()){
		//new_func.table()=vector<Double>(1);
		//return;
	//}

	int num_variables=variables.size();
	int num_functions=functions.size();
	int num_marg_variables=marg_variables.size();

	// Compute gray index for all variables and functions
	vector<vector<pair<int,int> > > gray_index(num_variables);
	//**int old_temp_value[num_variables];
	vector<int> old_temp_value(num_variables);
	for(int i=0;i<num_variables;i++)
	{
		old_temp_value[i]=variables[i]->temp_value;
		variables[i]->temp_value=i;
	}
	for(int i=0;i<num_functions;i++)
	{
		int multiplier=1;
		for(int j=0;j<functions[i]->variables().size();j++) {
			gray_index[functions[i]->variables()[j]->temp_value].push_back(pair<int,int>(i,multiplier));
			multiplier*=functions[i]->variables()[j]->domain_size();
		}
	}

	//Initialize the data structure for gray code
	//***
	//int a[num_variables];
	//int f[num_variables+1];
	//int o[num_variables+1];
	//int m[num_variables];
	
	vector<int> a(num_variables);
	vector<int> f(num_variables+1);
	vector<int> o(num_variables+1);
	vector<int> m(num_variables);

	for(int i=0;i<num_variables;i++)
	{
		a[i]=0;
		f[i]=i;
		o[i]=1;
		m[i]=variables[i]->domain_size();
	}
	f[num_variables]=num_variables;
	//**int func_address[num_functions];
	vector<int> func_address(num_functions);
	int marg_address=0;
	Double mult(1.0);
	for(int i=0;i<num_functions;i++)
	{
		if(!functions[i]->table().empty()) {
			if(functions[i]->table()[0] < 0)
				continue;
			mult*=functions[i]->table()[0];
		}
		func_address[i]=0;
	}
	int domain_size=Variable::getDomainSize(variables);
	//**int gray_marg_index[num_variables];
	vector<int> gray_marg_index(num_variables);
	for(int i=0;i<num_variables;i++)
		gray_marg_index[i]=0;
	int multiplier=1;
	for(int i=0;i<new_func.variables().size();i++)
	{
		gray_marg_index[new_func.variables()[i]->temp_value]=multiplier;
		multiplier*=new_func.variables()[i]->domain_size();
	}
	for(int i=0;i<domain_size;i++)
	{
		new_func.table()[marg_address]+=mult;
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
		//cout<<"j="<<j<<",oldaj="<<old_aj<<endl;
		//cout<<"A ="<<a[0]<<","<<a[1]<<","<<a[2]<<endl;
		//cout<<"O ="<<o[0]<<","<<o[1]<<","<<o[2]<<endl;
		//cout<<"F ="<<f[0]<<","<<f[1]<<","<<f[2]<<endl;
		//cout<<j<<endl;
		//cout<<marg_address<<endl;
		for(int k=0;k<gray_index[j].size();k++)
		{
			int index=gray_index[j][k].first;
			int addr_multiplier=gray_index[j][k].second;
			func_address[index]-=addr_multiplier*old_aj;
			func_address[index]+=addr_multiplier*a[j];
		}
		mult=Double(1.0);
		for(int k=0;k<num_functions;k++)
		{
			if(func_address[k] > functions[k]->table().size()-1) {
				continue;
			}
			if(functions[k]->table()[func_address[k]] < 0) 
				continue;
			mult*=functions[k]->table()[func_address[k]];
		}
		//End Hack
		if(gray_marg_index[j]>0)
		{
			marg_address-=gray_marg_index[j]*old_aj;
			marg_address+=gray_marg_index[j]*a[j];
		}
		//cout<<func_address[0]<<":"<<func_address[1]<<":"<<marg_address<<endl;
		//cout<<marg_address<<":"<<gray_marg_index[j]*(a[j]-old_aj)<<endl;
	}
	if (to_normalize){
		new_func.normalize();
	}
}

void Function::multiplyAndMarginalizeStored(vector<Variable*>& marg_variables_, vector<Function*>& functions, Function& new_func)
{
	if (functions.empty()) return;
	vector<Variable*> variables;
	vector<Variable*> marg_variables;
	for (int i = 0; i<functions.size(); i++)
	{
		do_set_union(variables, functions[i]->variables(), variables, less_than_comparator_variable);
	}
	do_set_intersection(variables, marg_variables_, marg_variables, less_than_comparator_variable);
	//if (variables.empty()) return;
	//new_func.variables() = marg_variables;
	//new_func.table() = vector<Double>(Variable::getDomainSize(marg_variables), 0);
	//cout<<"num-marg-vars = "<<marg_variables.size()<<" ";
	//if(variables.empty()){
	//new_func.table()=vector<Double>(1);
	//return;
	//}
	for (int i = 0; i < new_func.table().size(); i++) {
		new_func.table()[i] = 0;
	}
	int num_variables = variables.size();
	int num_functions = functions.size();
	int num_marg_variables = marg_variables.size();

	// Compute gray index for all variables and functions
	vector<vector<pair<int, int> > > gray_index(num_variables);
	//**int old_temp_value[num_variables];
	vector<int> old_temp_value(num_variables);
	for (int i = 0; i<num_variables; i++)
	{
		old_temp_value[i] = variables[i]->temp_value;
		variables[i]->temp_value = i;
	}
	for (int i = 0; i<num_functions; i++)
	{
		int multiplier = 1;
		for (int j = 0; j<functions[i]->variables().size(); j++) {
			gray_index[functions[i]->variables()[j]->temp_value].push_back(pair<int, int>(i, multiplier));
			multiplier *= functions[i]->variables()[j]->domain_size();
		}
	}

	//Initialize the data structure for gray code
	//***
	//int a[num_variables];
	//int f[num_variables+1];
	//int o[num_variables+1];
	//int m[num_variables];

	vector<int> a(num_variables);
	vector<int> f(num_variables + 1);
	vector<int> o(num_variables + 1);
	vector<int> m(num_variables);

	for (int i = 0; i<num_variables; i++)
	{
		a[i] = 0;
		f[i] = i;
		o[i] = 1;
		m[i] = variables[i]->domain_size();
	}
	f[num_variables] = num_variables;
	//**int func_address[num_functions];
	vector<int> func_address(num_functions);
	int marg_address = 0;
	Double mult(1.0);
	for (int i = 0; i<num_functions; i++)
	{
		if (!functions[i]->table().empty()) {
			if (functions[i]->table()[0] < 0)
				continue;
			mult *= functions[i]->table()[0];
		}
		func_address[i] = 0;
	}
	int domain_size = Variable::getDomainSize(variables);
	//**int gray_marg_index[num_variables];
	vector<int> gray_marg_index(num_variables);
	for (int i = 0; i<num_variables; i++)
		gray_marg_index[i] = 0;
	int multiplier = 1;
	for (int i = 0; i<new_func.variables().size(); i++)
	{
		gray_marg_index[new_func.variables()[i]->temp_value] = multiplier;
		multiplier *= new_func.variables()[i]->domain_size();
	}
	for (int i = 0; i<domain_size; i++)
	{
		new_func.table()[marg_address] += mult;
		int j = f[0];
		f[0] = 0;
		if (j == num_variables) break;
		int old_aj = a[j];
		a[j] = a[j] + o[j];
		if (a[j] == 0 || a[j] == (m[j] - 1))
		{
			o[j] = -o[j];
			f[j] = f[j + 1];
			f[j + 1] = j + 1;
		}
		//cout<<"j="<<j<<",oldaj="<<old_aj<<endl;
		//cout<<"A ="<<a[0]<<","<<a[1]<<","<<a[2]<<endl;
		//cout<<"O ="<<o[0]<<","<<o[1]<<","<<o[2]<<endl;
		//cout<<"F ="<<f[0]<<","<<f[1]<<","<<f[2]<<endl;
		//cout<<j<<endl;
		//cout<<marg_address<<endl;
		for (int k = 0; k<gray_index[j].size(); k++)
		{
			int index = gray_index[j][k].first;
			int addr_multiplier = gray_index[j][k].second;
			func_address[index] -= addr_multiplier*old_aj;
			func_address[index] += addr_multiplier*a[j];
		}
		mult = Double(1.0);
		for (int k = 0; k<num_functions; k++)
		{
			if (func_address[k] > functions[k]->table().size() - 1) {
				continue;
			}
			if (functions[k]->table()[func_address[k]] < 0)
				continue;
			mult *= functions[k]->table()[func_address[k]];
		}
		//End Hack
		if (gray_marg_index[j]>0)
		{
			marg_address -= gray_marg_index[j] * old_aj;
			marg_address += gray_marg_index[j] * a[j];
		}
		//cout<<func_address[0]<<":"<<func_address[1]<<":"<<marg_address<<endl;
		//cout<<marg_address<<":"<<gray_marg_index[j]*(a[j]-old_aj)<<endl;
	}
}

/*
void Function::multiplyAndMarginalize(vector<Variable*>& marg_variables_, vector<Function*>& functions, Function& new_func,bool to_normalize)
{
	if (functions.empty()) return;
	vector<Variable*> variables;
	vector<Variable*> marg_variables;
	for(int i=0;i<functions.size();i++)
	{
		do_set_union(variables,functions[i]->variables(),variables,less_than_comparator_variable);
	}
	do_set_intersection(variables,marg_variables_,marg_variables,less_than_comparator_variable);
	//if (variables.empty()) return;
	new_func.variables()=marg_variables;
	new_func.table()=vector<Double>(Variable::getDomainSize(marg_variables),0);
	//cout<<"num-marg-vars = "<<marg_variables.size()<<" ";
	//if(variables.empty()){
		//new_func.table()=vector<Double>(1);
		//return;
	//}

	int num_variables=variables.size();
	int num_functions=functions.size();
	int num_marg_variables=marg_variables.size();

	// Compute gray index for all variables and functions
	vector<vector<pair<int,int> > > gray_index(num_variables);
	//**int old_temp_value[num_variables];
	vector<int> old_temp_value(num_variables);
	for(int i=0;i<num_variables;i++)
	{
		old_temp_value[i]=variables[i]->temp_value;
		variables[i]->temp_value=i;
	}
	for(int i=0;i<num_functions;i++)
	{
		int multiplier=1;
		for(int j=0;j<functions[i]->variables().size();j++) {
			gray_index[functions[i]->variables()[j]->temp_value].push_back(pair<int,int>(i,multiplier));
			if(functions[i]->variables()[j]->id()==0)
				multiplier*=functions[i]->variables()[j]->domain_size();
			else
				multiplier*=functions[i]->variables()[j]->domain_size();
		}
	}

	//Initialize the data structure for gray code
	//***
	//int a[num_variables];
	//int f[num_variables+1];
	//int o[num_variables+1];
	//int m[num_variables];
	
	vector<int> a(num_variables);
	vector<int> f(num_variables+1);
	vector<int> o(num_variables+1);
	vector<int> m(num_variables);

	for(int i=0;i<num_variables;i++)
	{
		a[i]=0;
		f[i]=i;
		o[i]=1;
		if(variables[i]->id()==0)
			m[i]=variables[i]->domain_size();
		else
			m[i]=variables[i]->domain_size();
		
	}
	f[num_variables]=num_variables;
	//**int func_address[num_functions];
	vector<int> func_address(num_functions);
	int marg_address=0;
	Double mult(1.0);
	for(int i=0;i<num_functions;i++)
	{
		if(!functions[i]->table().empty())
			mult*=functions[i]->table()[0];
		func_address[i]=0;
	}
	int domain_size=Variable::getDomainSize(variables);
	//**int gray_marg_index[num_variables];
	vector<int> gray_marg_index(num_variables);
	for(int i=0;i<num_variables;i++)
		gray_marg_index[i]=0;
	int multiplier=1;
	for(int i=0;i<new_func.variables().size();i++)
	{
		gray_marg_index[new_func.variables()[i]->temp_value]=multiplier;
		if(new_func.variables()[i]->id()==0)
			multiplier*=new_func.variables()[i]->domain_size();
		else
			multiplier*=new_func.variables()[i]->domain_size();
	}
	domain_size = 9;
	//for(int i=0;i<domain_size;i++)
	int i = 1;
	while(true)
	{
		new_func.table()[marg_address]+=mult;
		int j=f[0];
		f[0]=0;
		if(j==num_variables) break;
		int old_aj=a[j];
		if(m[j]!=1){
		a[j]=a[j]+o[j];
		if(a[j]==0 || a[j]==(m[j]-1))
		{
			o[j]=-o[j];
			f[j]=f[j+1];
			f[j+1]=j+1;
		}
		}
		else{
		//if(a[j]==0)
		{
			f[0] = j+1;
			continue;	
		}
		}

		//cout<<"j="<<j<<",oldaj="<<old_aj<<endl;
		cout<<"j="<<j<<" A ="<<a[0]<<","<<a[1]<<","<<a[2]<<endl;
		//cout<<"O ="<<o[0]<<","<<o[1]<<","<<o[2]<<endl;
		cout<<"F ="<<f[0]<<","<<f[1]<<","<<f[2]<<endl;
		//cout<<j<<endl;
		for(int k=0;k<gray_index[j].size();k++)
		{
			int index=gray_index[j][k].first;
			int addr_multiplier=gray_index[j][k].second;
			func_address[index]-=addr_multiplier*old_aj;
			func_address[index]+=addr_multiplier*a[j];
		}
		mult=Double(1.0);
		for(int k=0;k<num_functions;k++)
		{
			if(func_address[k] > functions[k]->table().size()-1) {
				continue;
			}
			mult*=functions[k]->table()[func_address[k]];
		}
		//End Hack
		if(gray_marg_index[j]>0)
		{
			marg_address-=gray_marg_index[j]*old_aj;
			marg_address+=gray_marg_index[j]*a[j];
		}
		//cout<<func_address[0]<<" "<<func_address[1]<<" "<<marg_address<<endl;
		i++;
		if(i>=domain_size)
			break;
		//cout<<marg_address<<":"<<gray_marg_index[j]*(a[j]-old_aj)<<endl;
	}
	if (to_normalize){
		new_func.normalize();
	}

}
*/

void Function::multiplyAndMarginalizeNew(vector<Variable*>& marg_variables_, vector<Function*>& functions, 
	Function& new_func,vector<int>& hashtable)
{
	if (functions.empty()) return;
	vector<Variable*> variables;
	vector<Variable*> marg_variables;
	for(int i=0;i<functions.size();i++)
	{
		do_set_union(variables,functions[i]->variables(),variables,less_than_comparator_variable);
	}
	do_set_intersection(variables,marg_variables_,marg_variables,less_than_comparator_variable);
	//if (variables.empty()) return;
	new_func.variables()=marg_variables;
	int dsize = Variable::getDomainSize(marg_variables);
	new_func.table()=vector<Double>(dsize,0);
	hashtable = vector<int>(dsize,-1);
	//cout<<"num-marg-vars = "<<marg_variables.size()<<" ";
	//if(variables.empty()){
		//new_func.table()=vector<Double>(1);
		//return;
	//}

	int num_variables=variables.size();
	int num_functions=functions.size();
	int num_marg_variables=marg_variables.size();

	// Compute gray index for all variables and functions
	vector<vector<pair<int,int> > > gray_index(num_variables);
	//**int old_temp_value[num_variables];
	vector<int> old_temp_value(num_variables);
	for(int i=0;i<num_variables;i++)
	{
		old_temp_value[i]=variables[i]->temp_value;
		variables[i]->temp_value=i;
	}
	for(int i=0;i<num_functions;i++)
	{
		int multiplier=1;
		for(int j=0;j<functions[i]->variables().size();j++) {
			gray_index[functions[i]->variables()[j]->temp_value].push_back(pair<int,int>(i,multiplier));
			multiplier*=functions[i]->variables()[j]->domain_size();
		}
	}

	//Initialize the data structure for gray code
	//***
	//int a[num_variables];
	//int f[num_variables+1];
	//int o[num_variables+1];
	//int m[num_variables];
	
	vector<int> a(num_variables);
	vector<int> f(num_variables+1);
	vector<int> o(num_variables+1);
	vector<int> m(num_variables);

	for(int i=0;i<num_variables;i++)
	{
		a[i]=0;
		f[i]=i;
		o[i]=1;
		m[i]=variables[i]->domain_size();
	}
	f[num_variables]=num_variables;
	//**int func_address[num_functions];
	vector<int> func_address(num_functions);
	int marg_address=0;
	Double mult(1.0);
	for(int i=0;i<num_functions;i++)
	{
		if(!functions[i]->table().empty())
			mult*=functions[i]->table()[0];
		func_address[i]=0;
	}
	int domain_size=Variable::getDomainSize(variables);
	//**int gray_marg_index[num_variables];
	vector<int> gray_marg_index(num_variables);
	for(int i=0;i<num_variables;i++)
		gray_marg_index[i]=0;
	int multiplier=1;
	for(int i=0;i<new_func.variables().size();i++)
	{
		gray_marg_index[new_func.variables()[i]->temp_value]=multiplier;
		multiplier*=new_func.variables()[i]->domain_size();
	}
	for(int i=0;i<domain_size;i++)
	{
		new_func.table()[marg_address]+=mult;
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
		for(int k=0;k<gray_index[j].size();k++)
		{
			int index=gray_index[j][k].first;
			int addr_multiplier=gray_index[j][k].second;
			func_address[index]-=addr_multiplier*old_aj;
			func_address[index]+=addr_multiplier*a[j];
		}
		mult=Double(1.0);
		for(int k=0;k<num_functions;k++)
		{
			if(func_address[k] > functions[k]->table().size()-1) {
				continue;
			}
			mult*=functions[k]->table()[func_address[k]];
		}
		//End Hack
		if(gray_marg_index[j]>0)
		{
			marg_address-=gray_marg_index[j]*old_aj;
			marg_address+=gray_marg_index[j]*a[j];
		}
	}
}

void Function::projectMultiplyAndMarginalize(vector<Variable*>& marg_variables_, vector<Function*>& functions, 
Function& new_func,vector<int>& hashtable,vector<int>& margindexes)
{
	if (functions.empty()) return;
	vector<Variable*> variables;
	vector<Variable*> marg_variables;
	for(int i=0;i<functions.size();i++)
	{
		do_set_union(variables,functions[i]->variables(),variables,less_than_comparator_variable);
	}
	do_set_intersection(variables,marg_variables_,marg_variables,less_than_comparator_variable);
	//if (variables.empty()) return;
	//if(new_func.variables().size()==0) {
	//	new_func.variables()=marg_variables;
	//	new_func.table()=vector<Double>(Variable::getDomainSize(marg_variables),0);
	//}
	//cout<<"num-marg-vars = "<<marg_variables.size()<<" ";
	//if(variables.empty()){
		//new_func.table()=vector<Double>(1);
		//return;
	//}
	//for(int i=0;i<marg_variables.size();i++) {
	//	vector<int>::iterator it = find(vars.begin(),vars.end(),marg_variables[i]->id());
	//}

	int num_variables=variables.size();
	int num_functions=functions.size();
	int num_marg_variables=marg_variables.size();

	// Compute gray index for all variables and functions
	vector<vector<pair<int,int> > > gray_index(num_variables);
	//**int old_temp_value[num_variables];
	vector<int> old_temp_value(num_variables);
	for(int i=0;i<num_variables;i++)
	{
		old_temp_value[i]=variables[i]->temp_value;
		variables[i]->temp_value=i;
	}
	vector<int> offsets(num_functions);
	for(int i=0;i<num_functions;i++)
	{
		int multiplier=1;
		int offset = 0;
		for(int j=0;j<functions[i]->variables().size();j++) {
			gray_index[functions[i]->variables()[j]->temp_value].push_back(pair<int,int>(i,multiplier));

			if(find(marg_variables.begin(),marg_variables.end(),functions[i]->variables()[j])!=marg_variables.end()) {
				int val = functions[i]->variables()[j]->value();
				if(val!=INVALID_VALUE) {
					offset += val*multiplier;
				}
			}
			multiplier*=functions[i]->variables()[j]->domain_size();
		}
		offsets[i] = offset;
	}

	//Initialize the data structure for gray code
	//***
	//int a[num_variables];
	//int f[num_variables+1];
	//int o[num_variables+1];
	//int m[num_variables];
	
	vector<int> a(num_variables);
	vector<int> f(num_variables+1);
	vector<int> o(num_variables+1);
	vector<int> m(num_variables);
	vector<int>  mapvars(num_variables);
	int effnumvars = 0;
	int domain_size = 1;
	for(int i=0;i<num_variables;i++)
	{
		//constant check
		if(find(marg_variables.begin(),marg_variables.end(),variables[i])!=marg_variables.end()) {
			if(variables[i]->value()!=INVALID_VALUE) {
			continue;
			}
		}
		a[effnumvars]=0;
		f[effnumvars]=effnumvars;
		o[effnumvars]=1;
		m[effnumvars]=variables[i]->domain_size();
		domain_size *= variables[i]->domain_size();
		mapvars[effnumvars] = i;
		effnumvars++;
	}
	f[effnumvars]=effnumvars+1;
	//**int func_address[num_functions];
	vector<int> func_address(num_functions);
	Double mult(1.0);
	for(int i=0;i<num_functions;i++)
	{
		//if(!functions[i]->table().empty())
			//mult*=functions[i]->table()[0];
		mult*=functions[i]->table()[offsets[i]];
		//func_address[i]=0;
		func_address[i]=offsets[i];
	}
	//int domain_size=Variable::getDomainSize(variables);
	//**int gray_marg_index[num_variables];
	vector<int> gray_marg_index(num_variables);
	for(int i=0;i<num_variables;i++)
		gray_marg_index[i]=0;
	int multiplier=1;
	int margoffset = 0;
	for(int i=0;i<new_func.variables().size();i++)
	{
		gray_marg_index[new_func.variables()[i]->temp_value]=multiplier;
		int val = new_func.variables()[i]->value();
		if(val!=INVALID_VALUE) {
			margoffset += val*multiplier;
		}
		multiplier*=new_func.variables()[i]->domain_size();
	}
	int marg_address=margoffset;
	//for(int i=0;i<domain_size;i++)
	int icnt = 1;
	//f.resize(3);
	//f[0] = 0;
	//f[1] = 1;
	//f[2] = 2;
	//m.resize(2);
	//m[0] = 3;
	//m[1] = 3;
	//vector<int>  mapvars(2);
	//mapvars[0]  = 0;
	//mapvars[1] = 1;
	//cout<<offsets[0]<<" "<<offsets[1]<<" "<<margoffset<<endl;
	while(true)
	{
		if(hashtable[marg_address]==-1) {
			hashtable[marg_address] = 1;
			hashtable[marg_address] = new_func.table()[marg_address];
			new_func.table()[marg_address] = 0;
			margindexes.push_back(marg_address);
		}
		new_func.table()[marg_address]+=mult;
		int j=f[0];
		f[0]=0;
		if(j>=effnumvars) break;
		int old_aj=a[j];
		a[j]=a[j]+o[j];
		if(a[j]==0 || a[j]==(m[j]-1))
		{
			o[j]=-o[j];
			f[j]=f[j+1];
			f[j+1]=j+1;
		}

		//cout<<"j="<<j<<",oldaj="<<old_aj<<endl;
		//cout<<"j="<<j<<" A = "<<a[0]<<","<<a[1]<<","<<a[2]<<endl;
		//cout<<" A = 0 "<<a[0]<<" "<<a[1]<<endl;
		//cout<<"O ="<<o[0]<<","<<o[1]<<","<<o[2]<<endl;
		//cout<<"F ="<<f[0]<<","<<f[1]<<","<<f[2]<<endl;
		//cout<<j<<endl;
		for(int k=0;k<gray_index[mapvars[j]].size();k++)
		{
			int index=gray_index[mapvars[j]][k].first;
			int addr_multiplier=gray_index[mapvars[j]][k].second;
			func_address[index]-=addr_multiplier*old_aj;
			func_address[index]+=addr_multiplier*a[j];
		}
		mult=Double(1.0);
		for(int k=0;k<num_functions;k++)
		{
			if(func_address[k] > functions[k]->table().size()-1) {
				continue;
			}
			mult*=functions[k]->table()[func_address[k]];
		}
		//End Hack
		if(gray_marg_index[mapvars[j]]>0)
		{
			marg_address-=gray_marg_index[mapvars[j]]*old_aj;
			marg_address+=gray_marg_index[mapvars[j]]*a[j];
		}
		//cout<<func_address[0]+offsets[0]<<" "<<func_address[1]+offsets[1]<<" "<<marg_address+margoffset<<endl;
		icnt++;
		if(icnt>domain_size)
			break;
		//cout<<marg_address<<":"<<gray_marg_index[j]*(a[j]-old_aj)<<endl;
	}
}


void Function::projectMultiplyAndMarginalizeSJ(vector<Variable*>& marg_variables_, vector<Function*>& functions, 
Function& new_func,vector<int>& hashtable,vector<int>& margindexes,vector<bool>& dirtyBitArray)
{
	if (functions.empty()) return;

	vector<Variable*> variables;
	vector<Variable*> marg_variables;
	for(int i=0;i<functions.size();i++)
	{
		do_set_union(variables,functions[i]->variables(),variables,less_than_comparator_variable);
	}
	do_set_intersection(variables,marg_variables_,marg_variables,less_than_comparator_variable);
	int num_variables=variables.size();
	int num_functions=functions.size();
	int num_marg_variables=marg_variables.size();

	//**int old_temp_value[num_variables];
	vector<int> old_temp_value(num_variables);
	for(int i=0;i<num_variables;i++)
	{
		old_temp_value[i]=variables[i]->temp_value;
		variables[i]->temp_value=i;
	}

	//CHECK for self-joined vars
	vector<int> sj_index;
	for(int i=0;i<marg_variables.size();i++) {
		if(marg_variables[i]->sj_values.size() > 0) {
			marg_variables[i]->value() = INVALID_VALUE;
			sj_index.push_back(i);
		}
	}
	int lastidx = 0;
	for(int sj=0;sj<sj_index.size();sj++) {
		for(int s = 0;s<marg_variables[sj]->sj_values.size();s++) {
			marg_variables[sj_index[sj]]->value() = marg_variables[sj_index[sj]]->sj_values[s];
		// Compute gray index for all variables and functions
		vector<vector<pair<int,int> > > gray_index(num_variables);
		vector<int> offsets(num_functions);
		for(int i=0;i<num_functions;i++)
		{
			int multiplier=1;
			int offset = 0;
			for(int j=0;j<functions[i]->variables().size();j++) {
				gray_index[functions[i]->variables()[j]->temp_value].push_back(pair<int,int>(i,multiplier));

				if(find(marg_variables.begin(),marg_variables.end(),functions[i]->variables()[j])!=marg_variables.end()) {
					int val = functions[i]->variables()[j]->value();
					if(val!=INVALID_VALUE) {
						offset += val*multiplier;
					}
				}
				multiplier*=functions[i]->variables()[j]->domain_size();
			}
			offsets[i] = offset;
		}

		vector<int> a(num_variables);
		vector<int> f(num_variables+1);
		vector<int> o(num_variables+1);
		vector<int> m(num_variables);
		vector<int>  mapvars(num_variables);
		int effnumvars = 0;
		int domain_size = 1;
		for(int i=0;i<num_variables;i++)
		{
			//constant check
			if(find(marg_variables.begin(),marg_variables.end(),variables[i])!=marg_variables.end()) {
				if(variables[i]->value()!=INVALID_VALUE) {
				continue;
				}
			}
			a[effnumvars]=0;
			f[effnumvars]=effnumvars;
			o[effnumvars]=1;
			m[effnumvars]=variables[i]->domain_size();
			domain_size *= variables[i]->domain_size();
			mapvars[effnumvars] = i;
			effnumvars++;
		}
		f[effnumvars]=effnumvars+1;
		//**int func_address[num_functions];
		vector<int> func_address(num_functions);
		Double mult(1.0);
		for(int i=0;i<num_functions;i++)
		{
			mult*=functions[i]->table()[offsets[i]];
			func_address[i]=offsets[i];
		}
		//int domain_size=Variable::getDomainSize(variables);
		//**int gray_marg_index[num_variables];
		vector<int> gray_marg_index(num_variables);
		for(int i=0;i<num_variables;i++)
			gray_marg_index[i]=0;
		int multiplier=1;
		int margoffset = 0;
		for(int i=0;i<new_func.variables().size();i++)
		{
			gray_marg_index[new_func.variables()[i]->temp_value]=multiplier;
			int val = new_func.variables()[i]->value();
			if(val!=INVALID_VALUE) {
				margoffset += val*multiplier;
			}
			multiplier*=new_func.variables()[i]->domain_size();
		}
		int marg_address=margoffset;
		//for(int i=0;i<domain_size;i++)
		int icnt = 1;
		while(true)
		{
			if(hashtable[marg_address]==-1) {
				hashtable[marg_address] = new_func.table()[marg_address];
				new_func.table()[marg_address] = 0;
				margindexes.push_back(marg_address);
			}
			if(!dirtyBitArray[marg_address]) {
				new_func.table()[marg_address]+=mult;
			}
			int j=f[0];
			f[0]=0;
			if(j>=effnumvars) break;
			int old_aj=a[j];
			a[j]=a[j]+o[j];
			if(a[j]==0 || a[j]==(m[j]-1))
			{
				o[j]=-o[j];
				f[j]=f[j+1];
				f[j+1]=j+1;
			}
			for(int k=0;k<gray_index[mapvars[j]].size();k++)
			{
				int index=gray_index[mapvars[j]][k].first;
				int addr_multiplier=gray_index[mapvars[j]][k].second;
				func_address[index]-=addr_multiplier*old_aj;
				func_address[index]+=addr_multiplier*a[j];
			}
			mult=Double(1.0);
			for(int k=0;k<num_functions;k++)
			{
				if(func_address[k] > functions[k]->table().size()-1) {
					continue;
				}
				mult*=functions[k]->table()[func_address[k]];
			}
			//End Hack
			if(gray_marg_index[mapvars[j]]>0)
			{
				marg_address-=gray_marg_index[mapvars[j]]*old_aj;
				marg_address+=gray_marg_index[mapvars[j]]*a[j];
			}
			icnt++;
			if(icnt>domain_size)
				break;
			//cout<<marg_address<<":"<<gray_marg_index[j]*(a[j]-old_aj)<<endl;
		}
		for(int j=lastidx;j<margindexes.size();j++) {
			dirtyBitArray[margindexes[j]] = true;
		}
		lastidx = margindexes.size();
		marg_variables[sj_index[sj]]->value() = INVALID_VALUE;
	}
	}

	for(int i=0;i<margindexes.size();i++) {
		dirtyBitArray[margindexes[i]] = false;
	}
}

/*
//for every grounding of Self join marg --Wrong
void Function::projectMultiplyAndMarginalizeSJ(vector<Variable*>& marg_variables_, vector<Function*>& functions, 
Function& new_func,vector<int>& hashtable,vector<int>& margindexes)
{
	if (functions.empty()) return;
	vector<Variable*> variables;
	vector<Variable*> marg_variables;
	for(int i=0;i<functions.size();i++)
	{
		do_set_union(variables,functions[i]->variables(),variables,less_than_comparator_variable);
	}
	do_set_intersection(variables,marg_variables_,marg_variables,less_than_comparator_variable);
	int num_variables=variables.size();
	int num_functions=functions.size();
	int num_marg_variables=marg_variables.size();

	//**int old_temp_value[num_variables];
	vector<int> old_temp_value(num_variables);
	for(int i=0;i<num_variables;i++)
	{
		old_temp_value[i]=variables[i]->temp_value;
		variables[i]->temp_value=i;
	}

	//count for every combination of the self joined marg variable values
	int nxtchange = marg_variables.size()-1;
	int totalsize = 1;
	for(int i=0;i<marg_variables.size();i++) {
		totalsize *= marg_variables[i]->tmp_domain.size();
	}
	vector<int> inds(marg_variables.size());
	for(int t=0;t<totalsize;t++) {
		for(int i=0;i<marg_variables.size();i++) {
			marg_variables[i]->value() = marg_variables[i]->tmp_domain[inds[i]];
			cout<<marg_variables[i]->tmp_domain[inds[i]]<<" ";
		}
		cout<<endl;
		// Compute gray index for all variables and functions
		vector<vector<pair<int,int> > > gray_index(num_variables);
		vector<int> offsets(num_functions);
		for(int i=0;i<num_functions;i++)
		{
			int multiplier=1;
			int offset = 0;
			for(int j=0;j<functions[i]->variables().size();j++) {
				gray_index[functions[i]->variables()[j]->temp_value].push_back(pair<int,int>(i,multiplier));

				if(find(marg_variables.begin(),marg_variables.end(),functions[i]->variables()[j])!=marg_variables.end()) {
					int val = functions[i]->variables()[j]->value();
					if(val!=INVALID_VALUE) {
						offset += val*multiplier;
					}
				}
				multiplier*=functions[i]->variables()[j]->domain_size();
			}
			offsets[i] = offset;
		}

		vector<int> a(num_variables);
		vector<int> f(num_variables+1);
		vector<int> o(num_variables+1);
		vector<int> m(num_variables);
		vector<int>  mapvars(num_variables);
		int effnumvars = 0;
		int domain_size = 1;
		for(int i=0;i<num_variables;i++)
		{
			//constant check
			if(find(marg_variables.begin(),marg_variables.end(),variables[i])!=marg_variables.end()) {
				if(variables[i]->value()!=INVALID_VALUE) {
				continue;
				}
			}
			a[effnumvars]=0;
			f[effnumvars]=effnumvars;
			o[effnumvars]=1;
			m[effnumvars]=variables[i]->domain_size();
			domain_size *= variables[i]->domain_size();
			mapvars[effnumvars] = i;
			effnumvars++;
		}
		f[effnumvars]=effnumvars+1;
		//**int func_address[num_functions];
		vector<int> func_address(num_functions);
		Double mult(1.0);
		for(int i=0;i<num_functions;i++)
		{
			mult*=functions[i]->table()[offsets[i]];
			func_address[i]=offsets[i];
		}
		//int domain_size=Variable::getDomainSize(variables);
		//**int gray_marg_index[num_variables];
		vector<int> gray_marg_index(num_variables);
		for(int i=0;i<num_variables;i++)
			gray_marg_index[i]=0;
		int multiplier=1;
		int margoffset = 0;
		for(int i=0;i<new_func.variables().size();i++)
		{
			gray_marg_index[new_func.variables()[i]->temp_value]=multiplier;
			int val = new_func.variables()[i]->value();
			if(val!=INVALID_VALUE) {
				margoffset += val*multiplier;
			}
			multiplier*=new_func.variables()[i]->domain_size();
		}
		int marg_address=margoffset;
		//for(int i=0;i<domain_size;i++)
		int icnt = 1;
		while(true)
		{
			if(hashtable[marg_address]==-1) {
				hashtable[marg_address] = new_func.table()[marg_address];
				new_func.table()[marg_address] = 0;
				margindexes.push_back(marg_address);
			}
			new_func.table()[marg_address]+=mult;
			int j=f[0];
			f[0]=0;
			if(j>=effnumvars) break;
			int old_aj=a[j];
			a[j]=a[j]+o[j];
			if(a[j]==0 || a[j]==(m[j]-1))
			{
				o[j]=-o[j];
				f[j]=f[j+1];
				f[j+1]=j+1;
			}
			for(int k=0;k<gray_index[mapvars[j]].size();k++)
			{
				int index=gray_index[mapvars[j]][k].first;
				int addr_multiplier=gray_index[mapvars[j]][k].second;
				func_address[index]-=addr_multiplier*old_aj;
				func_address[index]+=addr_multiplier*a[j];
			}
			mult=Double(1.0);
			for(int k=0;k<num_functions;k++)
			{
				if(func_address[k] > functions[k]->table().size()-1) {
					continue;
				}
				mult*=functions[k]->table()[func_address[k]];
			}
			//End Hack
			if(gray_marg_index[mapvars[j]]>0)
			{
				marg_address-=gray_marg_index[mapvars[j]]*old_aj;
				marg_address+=gray_marg_index[mapvars[j]]*a[j];
			}
			icnt++;
			if(icnt>domain_size)
				break;
			//cout<<marg_address<<":"<<gray_marg_index[j]*(a[j]-old_aj)<<endl;
		}
		inds[nxtchange]++;
		if(inds[nxtchange]==marg_variables[nxtchange]->tmp_domain.size()) {
			inds[nxtchange] = 0;
			for(int j=nxtchange-1;j>=0;j--) {
				inds[j]++;
				if(inds[j] > marg_variables[j]->tmp_domain.size()-1) {
					inds[j] = 0;
				} else {
					break;
				}
			}
		}
	}
}
*/
void Function::product(Function& function)
{
	if(function.table().empty() || function.variables().empty())
		return;
	vector<Variable*> new_variables;
	do_set_union(variables(),function.variables(),new_variables,less_than_comparator_variable);
	int num_values=Variable::getDomainSize(new_variables);
	if(new_variables.size()==variables().size())
	{
		for(int i=0;i<num_values;i++)
		{
			Variable::setAddress(variables(),i);
			int func_entry=Variable::getAddress(function.variables());
			table()[i]*=function.table()[func_entry];
		}
	}
	else
	{
		vector<Double> old_table;
		old_table=table_;
		table_=vector<Double> (num_values);
		for(int i=0;i<num_values;i++)
		{
			Variable::setAddress(new_variables,i);
			int entry=Variable::getAddress(variables());
			int func_entry=Variable::getAddress(function.variables());
			table()[i]=Double(function.table()[func_entry]*old_table[entry]);
		}
		variables_=new_variables;
	}
}

/*
void Function::multiplyAndMarginalize(vector<Variable*>& marg_variables_,
		vector<Function*>& functions, Function& out_function,
		bool to_normalize) {

	vector<Variable*> variables;
	vector<Variable*> marg_variables;
	int num_functions = functions.size();
	if (num_functions == 0)
		return;
	for (int i = 0; i < num_functions; i++) {
		do_set_union(variables, functions[i]->variables(), variables,
				less_than_comparator_variable);
	}
	int num_variables = variables.size();
	if (num_variables == 0)
		return;
	do_set_intersection(variables, marg_variables_, marg_variables,
			less_than_comparator_variable);
	int num_marg_variables = marg_variables.size();
	//if(num_marg_variables==0)
	//	return;

	// 6 arrays for graycoding using Knuth's algorithm
	vector<vector<pair<int, int> > > g(num_variables);
	int c[num_functions];
	int a[num_variables];
	int f[num_variables + 1];
	int o[num_variables];
	int m[num_variables];
	int t[num_variables];
	Double mult(1.0);
	int address = 0;

	// Init variables for graycoding
	for (int i = 0; i < num_variables; i++) {
		a[i] = 0;
		f[i] = i;
		o[i] = 1;
		t[i] = 0;
		m[i] = variables[i]->domain_size();
		variables[i]->temp_value = i;
	}
	for (int i = 0; i < num_functions; i++)
		c[i] = 0;
	f[num_variables] = num_variables;
	for (int i = 0; i < num_functions; i++) {
		int multiplier = 1;
		assert(functions[i]!=NULL);
		for (int j = 0; j < functions[i]->variables().size(); j++) {
			g[functions[i]->variables()[j]->temp_value].push_back(
					pair<int, int>(i, multiplier));
			multiplier *= functions[i]->variables()[j]->domain_size();
		}
		if (!functions[i]->table().empty()) {
			mult *= functions[i]->table()[0];
		}
	}
	int multiplier = 1;
	//cout<<"mult here\n";
	for (int i = 0; i < num_marg_variables; i++) {
		t[marg_variables[i]->temp_value] = multiplier;
		multiplier *= marg_variables[i]->domain_size();
	}
	//cout<<"mult initing log function\n";
	//Gray  code algorithm
	//Initialize LogFunction
	out_function.variables() = marg_variables;
	out_function.table() = vector<Double>(
			Variable::getDomainSize(marg_variables),0.0);

	//cout<<"Log function inited\n";
	while (1) {
		//cout<<address<<endl;
		// Step 1: Visit
		out_function.table()[address] += mult;
		// Step 2: Choose j
		int j = f[0];
		f[0] = 0;

		if (j == num_variables)
			break;
		int old_aj = a[j];
		a[j] = a[j] + o[j];
		if (a[j] == 0 || a[j] == (m[j] - 1)) {
			o[j] = -o[j];
			f[j] = f[j + 1];
			f[j + 1] = j + 1;
		}
		if (mult == 0.0){
			for (int i = 0; i < g[j].size(); i++) {
				int index = g[j][i].first;
				int multiplier = g[j][i].second;
				mult /= functions[index]->table()[c[index]];
				c[index] -= multiplier * old_aj;
				c[index] += multiplier * a[j];
				mult *= functions[index]->table()[c[index]];
			}
		}
		else {
			for (int i = 0; i < g[j].size(); i++) {
				int index = g[j][i].first;
				int multiplier = g[j][i].second;
				c[index] -= multiplier * old_aj;
				c[index] += multiplier * a[j];
			}
			mult = Double(1.0);
			for (int i = 0; i < num_functions; i++)
				mult *= functions[i]->table()[c[i]];
		}
		if (t[j] > 0) {
			address -= t[j] * old_aj;
			address += t[j] * a[j];
		}
	}
	if (to_normalize)
		out_function.normalize();
}
*/


void Function::initApproxInfer(vector<Variable*>& marg_variables_, vector<Function*>& functions, 
	Function& new_func,vector<int>& hashtable,vector<double>& alt_message)
{
	if (functions.empty()) return;
	vector<Variable*> variables;
	vector<Variable*> marg_variables;
	for(int i=0;i<functions.size();i++)
	{
		do_set_union(variables,functions[i]->variables(),variables,less_than_comparator_variable);
	}
	do_set_intersection(variables,marg_variables_,marg_variables,less_than_comparator_variable);

	int num_variables=variables.size();
	int num_functions=functions.size();
	int num_marg_variables=marg_variables.size();

	new_func.variables()=marg_variables;
	int dsize = Variable::getDomainSize(marg_variables);
	new_func.table()=vector<Double>(dsize,0);
	alt_message=vector<Double>(dsize,0);
	hashtable = vector<int>(dsize,-1);
}

//void Function::approxMultiplyAndMarginalize(vector<Variable*>& marg_variables_, vector<Function*>& functions, 
//	Function& new_func,vector<int>& hashtable,vector<int>& margindexes,vector<double>& alt_table,vector<int> ftable)
void Function::approxMultiplyAndMarginalize(vector<Variable*>& marg_variables_, vector<Function*>& functions,
	Function& new_func, vector<int>& hashtable, vector<int>& margindexes)
{
	if (functions.empty()) return;
	vector<Variable*> variables;
	vector<Variable*> marg_variables;
	for(int i=0;i<functions.size();i++)
	{
		do_set_union(variables,functions[i]->variables(),variables,less_than_comparator_variable);
	}
	do_set_intersection(variables,marg_variables_,marg_variables,less_than_comparator_variable);

	int num_variables=variables.size();
	int num_functions=functions.size();
	int num_marg_variables=marg_variables.size();

	// Compute gray index for all variables and functions
	vector<vector<pair<int,int> > > gray_index(num_variables);
	//**int old_temp_value[num_variables];
	vector<int> old_temp_value(num_variables);
	for(int i=0;i<num_variables;i++)
	{
		old_temp_value[i]=variables[i]->temp_value;
		variables[i]->temp_value=i;
	}
	vector<int> offsets(num_functions);

	for(int i=0;i<num_functions;i++)
	{
		int multiplier=1;
		int offset = 0;
		for(int j=0;j<functions[i]->variables().size();j++) {
			gray_index[functions[i]->variables()[j]->temp_value].push_back(pair<int,int>(i,multiplier));
			int val = functions[i]->variables()[j]->value();
			if(val!=INVALID_VALUE) {
				offset += val*multiplier;
			}
			multiplier*=functions[i]->variables()[j]->domain_size();
		}
		offsets[i] = offset;
	}

	vector<int> a(num_variables);
	vector<int> f(num_variables+1);
	vector<int> o(num_variables+1);
	vector<int> m(num_variables);
	vector<int>  mapvars(num_variables);
	int effnumvars = 0;
	int domain_size = 1;
	for(int i=0;i<num_variables;i++)
	{
		//constant check
		/*if(find(marg_variables.begin(),marg_variables.end(),variables[i])!=marg_variables.end()) {
			if(variables[i]->value()!=INVALID_VALUE) {
				continue;
			}
		}*/
		if(variables[i]->value()!=INVALID_VALUE) {
			continue;
		}
		a[effnumvars]=0;
		f[effnumvars]=effnumvars;
		o[effnumvars]=1;
		m[effnumvars]=variables[i]->domain_size();
		domain_size *= variables[i]->domain_size();
		mapvars[effnumvars] = i;
		effnumvars++;
	}
	f[effnumvars]=effnumvars+1;
	//**int func_address[num_functions];
	vector<int> func_address(num_functions);
	Double mult(1.0);
	Double mult1(1.0);
	for(int i=0;i<num_functions;i++)
	{
		//if(!functions[i]->table().empty())
			//mult*=functions[i]->table()[0];
		mult*=functions[i]->table()[offsets[i]];
		double negml = functions[i]->table()[offsets[i]];
		//if(functions[i]->id()!=-1) {
		//	if(ftable[functions[i]->id()]==offsets[i]) {
		//		negml = 1 - functions[i]->table()[offsets[i]];
		//	}
		//}
		mult1 *= negml;
		//func_address[i]=0;
		func_address[i]=offsets[i];
	}
	//int domain_size=Variable::getDomainSize(variables);
	//**int gray_marg_index[num_variables];
	vector<int> gray_marg_index(num_variables);
	for(int i=0;i<num_variables;i++)
		gray_marg_index[i]=0;
	int multiplier=1;
	int margoffset = 0;
	for(int i=0;i<new_func.variables().size();i++)
	{
		gray_marg_index[new_func.variables()[i]->temp_value]=multiplier;
		int val = new_func.variables()[i]->value();
		if(val!=INVALID_VALUE) {
			margoffset += val*multiplier;
		}
		multiplier*=new_func.variables()[i]->domain_size();
	}
	int marg_address=margoffset;
	//for(int i=0;i<domain_size;i++)
	int icnt = 1;
	while(true)
	{
		if(hashtable[marg_address]==-1) {
			hashtable[marg_address] = 1;
			hashtable[marg_address] = new_func.table()[marg_address];
			new_func.table()[marg_address] = 0;
			//alt_table[marg_address] = 0;
			margindexes.push_back(marg_address);
		}
		new_func.table()[marg_address]+=mult;
		//alt_table[marg_address]+=mult1;
		int j=f[0];
		f[0]=0;
		if(j>=effnumvars) break;
		int old_aj=a[j];
		a[j]=a[j]+o[j];
		if(a[j]==0 || a[j]==(m[j]-1))
		{
			o[j]=-o[j];
			f[j]=f[j+1];
			f[j+1]=j+1;
		}

		//cout<<"j="<<j<<",oldaj="<<old_aj<<endl;
		//cout<<"j="<<j<<" A = "<<a[0]<<","<<a[1]<<","<<a[2]<<endl;
		//cout<<" A = 0 "<<a[0]<<" "<<a[1]<<endl;
		//cout<<"O ="<<o[0]<<","<<o[1]<<","<<o[2]<<endl;
		//cout<<"F ="<<f[0]<<","<<f[1]<<","<<f[2]<<endl;
		//cout<<j<<endl;
		for(int k=0;k<gray_index[mapvars[j]].size();k++)
		{
			int index=gray_index[mapvars[j]][k].first;
			int addr_multiplier=gray_index[mapvars[j]][k].second;
			func_address[index]-=addr_multiplier*old_aj;
			func_address[index]+=addr_multiplier*a[j];
		}
		mult=Double(1.0);
		mult1=Double(1.0);
		for(int k=0;k<num_functions;k++)
		{
			if(func_address[k] > functions[k]->table().size()-1) {
				continue;
			}
			mult*=functions[k]->table()[func_address[k]];
			double negml = functions[k]->table()[func_address[k]];
			//if(functions[k]->id()!=-1) {
			//	if(ftable[functions[k]->id()]==func_address[k]) {
			//		negml = 1 - functions[k]->table()[func_address[k]];
			//	}
			//}
			mult1 *= negml;
		}
		//End Hack
		if(gray_marg_index[mapvars[j]]>0)
		{
			marg_address-=gray_marg_index[mapvars[j]]*old_aj;
			marg_address+=gray_marg_index[mapvars[j]]*a[j];
		}
		//cout<<func_address[0]+offsets[0]<<" "<<func_address[1]+offsets[1]<<" "<<marg_address+margoffset<<endl;
		icnt++;
		if(icnt>domain_size)
			break;
		//cout<<marg_address<<":"<<gray_marg_index[j]*(a[j]-old_aj)<<endl;
	}
}
