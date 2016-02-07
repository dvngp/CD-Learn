#ifndef FUNCTION_H_
#define FUNCTION_H_

#include "Globals.h"
#include "Variable.h"



using namespace std;
struct Function
{
protected:
	int id_;
	vector<Variable*> variables_;
	vector<Double> table_;
	
public:
	Function():id_(INVALID_VALUE){}
	Function(const int& id): id_(id),variables_(vector<Variable*>()),table_(vector<Double>())
	{
	}
	Function(const int& id, const vector<Variable*>& variables): id_(id),variables_(variables), table_(vector<Double>(Variable::getDomainSize(variables)))
	{
	}
	virtual ~Function(){ }
	int& id() { return id_;}
	int id() const {return id_;}
	vector<Variable*>& variables()  { return variables_;}
	vector<Double>& table() { return table_;}
	void normalize();
	void print(ostream& out=cout){ 
		for (int i = 0; i < table_.size(); i++) {
			out << table_[i] << " ";
		}
		cout << endl;
	}
	virtual void removeEvidence();
	virtual void reduceDomains();
	static void multiplyAndMarginalize( vector<Variable*>& marg_variables, vector<Function*>& functions, Function& function, bool to_normalize=true);
	void product(Function& function);
	static void multiplyAndMarginalizeNew(vector<Variable*>& marg_variables_, vector<Function*>& functions, 
	Function& new_func,vector<int>& hashtable);
	static void projectMultiplyAndMarginalize(vector<Variable*>& marg_variables_, vector<Function*>& functions, 
	Function& new_func,vector<int>& hashtable,vector<int>& margindexes);
	static void projectMultiplyAndMarginalizeSJ(vector<Variable*>& marg_variables_, vector<Function*>& functions, 
	Function& new_func,vector<int>& hashtable,vector<int>& margindexes,vector<bool>& dirtyBitArray);
	static void multiplyAndMarginalizeStored(vector<Variable*>& marg_variables, vector<Function*>& functions, Function& function);
	
	int predid;
	vector<int> relationalOrder;
	//MAP inference
	vector<bool> isEvid;
	//approx
	static void initApproxInfer(vector<Variable*>& marg_variables_, vector<Function*>& functions, 
	Function& new_func,vector<int>& hashtable,vector<double>& alt_message);
	static void approxMultiplyAndMarginalize(vector<Variable*>& marg_variables_, vector<Function*>& functions, 
	Function& new_func,vector<int>& hashtable,vector<int>& margindexes);
	vector<int> addressTable;
};


#endif
