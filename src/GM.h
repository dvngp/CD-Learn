#ifndef GM_H_
#define GM_H_


#include "Globals.h"

#include "Variable.h"
#include "Function.h"
#include "CPT.h"
using namespace std;



// A graphical model has Random Variables and Functions
struct GM {
	NETWORK_TYPE type;
	NETWORK_MODE mode;

	// Variable, functions and copy of variables
	vector<Variable*> variables;
	vector<Function*> functions;
	vector<Variable*> copy_of_variables;
	vector<vector<int> > csp_to_sat_variables;
	// multiplicative constant, useful for computing partition functions
	Double mult_factor;
	//
	GM():
		mult_factor(1.0), type(BAYES), mode(POSITIVE){
	}

	void readUAI08(const char* infile);
	void setEvidenceBeliefsUAI08(vector<int>& evidence);
	void reduceDomains();

	// Orderings information
	void getMinFillOrdering(vector<int>& order,
			vector<set<int> >& clusters, double& estimate, int& max_cluster_size);
	void getLexOrdering(vector<int>& order, vector<set<int> >& clusters,
			double& estimate);
	// Cutset algorithms
	void rearrangeOrdering_randomized(std::vector<int> &order, std::vector<set<
				int> > &clusters, std::vector<int> &new_order, double& limit);
	// printing marginals
	void printMarginalsUAI10(vector<vector<Double1> >& marginals, ostream& out);
	void printMarginalsUAI10(vector<Function >& marginals, ostream& out);
	void funcReduce();
};

#endif
