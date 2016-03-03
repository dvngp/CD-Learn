#ifndef MLN_H_
#define MLN_H_

#include <vector>
#include <iostream>
#include <set>
#include <map>
#include <hash_map>
#include <algorithm>
#include <fstream>
#include <string>
#include<iomanip>
using namespace std;
#include "stringconversionutils.h"
#include "GM.h"
#include "BE.h"
//#include "JT.h"
//#include "VEC.h"
struct Term
{
	int domsize;
	int groundval;    
};

struct PredicateSymbol
{
	int id;
	string symbol;
	vector<int> domsizes;
	vector<int> cumulativesz;
	void setcsz() {
		cumulativesz.resize(domsizes.size());
		for(int i=0;i<domsizes.size();i++) {
			int mlt = 1;
			for(int j=i+1;j<domsizes.size();j++) {
				mlt *= domsizes[j];
			}
			cumulativesz[i] = mlt;
		}
	}
};

struct Atom
{
	int predid;
	//terms may be shared across atoms
	vector<Term*> terms;
};

struct WClause
{
	vector<Atom*> atoms;
	vector<bool> sign;
	double weight;
	double numgroundings;
	bool isHard;
};

enum InferenceType{
	GIBBS = 0,
	WS,
	AGIBBS
};
struct MLN
{
	vector<PredicateSymbol*> predicates;
	vector<WClause*> clauses;
	vector<vector<int> > assignments;
	vector<vector<int> > numsamples;
	vector<vector<int> > mapassignment;
	vector<vector<int> > nevids;
	vector<vector<double> > probabilities;
	vector<vector<int> > assignmentsevid;
	vector<double> avgcounts;
	int itersperupdate;
	vector<GM> models;
	vector<vector<double> > softevidences;
	vector<vector<bool> > isEvidence;
	vector<bool> isQuery;
	//vector<JT*> jtrees;
	vector<BE*> btrees;
	vector<double> currentcounts;
	vector<double> datacounts;
	vector<double> totalgroundings;
	vector<vector<double> > currentestimates;
	//vector<int> qrypreds;
	vector<vector<int> > constraints;
	int burnin;
	double bestmapcost;
	double gprob;
	InferenceType iType;
	double alpha;
	double lambda;
	int samplecnt;

	vector<vector<int> > allsamplecounts;
	vector<double> weights;
	vector<double> averageweights;
	vector<double> oldweights;
	vector<double> oldgradient;
	double cglambda;
	double momentum;
	double dotprod(vector<double> a,vector<double> b) {
		double dp = 0;
		for(int i=0;i<a.size();i++)
			dp += a[i]*b[i];
		return dp;
	}

	void getpredterms(string s,string& symbol,vector<string>& terms) {
		if(s[0]=='!') {
			s = s.substr(1);
		}
		int ind = s.find("(");
		string s1 = s.substr(0,ind);
		string s2 = s.substr(ind+1);
		ind = s2.find(")");
		string s3 = s2.substr(0,ind);
		LStringConversionUtils::tokenize(s3,terms,",");
		symbol = s1;
	}

	void getsoftevid(string s,string& symbol,vector<string>& terms,double& wt) {
		int ind = s.find(" ");
		string ws = s.substr(0,ind);
		wt = LStringConversionUtils::toDouble(ws);
		s = s.substr(ind+1);
		ind = s.find("(");
		string s1 = s.substr(0,ind);
		string s2 = s.substr(ind+1);
		ind = s2.find(")");
		string s3 = s2.substr(0,ind);
		LStringConversionUtils::tokenize(s3,terms,",");
		symbol = s1;
	}

	int getaddr(int id,vector<int> terms) {
		int addr = 0;
		for(int i=0;i<terms.size();i++) {
			addr += terms[i]*predicates[id]->cumulativesz[i];
		}
		return addr;
	}

	int getaddr(string symbol,vector<int> terms) {
		int id = 0;
		for(int i=0;i<predicates.size();i++) {
			if(predicates[i]->symbol==symbol) {
				id = predicates[i]->id;
				break;
			}
		}
		return getaddr(id,terms);
	}

	void getgrounding(int predid, int addr, vector<int>& grd) {
		grd.resize(predicates[predid]->domsizes.size());
		for(int idx=grd.size()-1;idx>=0;idx--) {
			grd[idx] = addr%predicates[predid]->domsizes[idx];
			addr = addr/predicates[predid]->domsizes[idx];
		}
	}
	int getFunctionIndex(int predid, int addr) {
		vector<int> grd(predicates[predid]->domsizes.size());
		int idx = grd.size()-1;
		for(int i=0;i<predicates[predid]->domsizes.size();i++) {
			grd[idx--] = addr%predicates[predid]->domsizes[i];
			addr = addr/predicates[predid]->domsizes[i];
		}
		int outaddr = 0;
		int didx = predicates[predid]->domsizes.size()-1;
		int mult = 1;
		for(int i=0;i<grd.size();i++) {
			outaddr += grd[i]*mult;
			mult *= predicates[predid]->domsizes[didx--];
		}
		return outaddr;
	}
	void init(string mlnfile,string evidfile,string qryfile,double alpha_=0.0001) {
		burnin = 100;
		gprob = 0.25;
		//itersperupdate = 100;
		alpha = alpha_;
		lambda = 0.001;
		samplecnt = 0;
		ifstream infile(mlnfile.c_str());
		string predline;
		getline(infile,predline);
		vector<string> toks;
		LStringConversionUtils::tokenize(predline,toks," ");
		predicates.resize(toks.size());
		for(int i=0;i<toks.size();i++) {
			vector<string> toks1;
			LStringConversionUtils::tokenize(toks[i],toks1,":");
			PredicateSymbol* ps = new PredicateSymbol();
			ps->id = i;
			ps->symbol = toks1[0];
			for(int j=1;j<toks1.size();j++) {
				int dsize = LStringConversionUtils::toInt(toks1[j]);
				ps->domsizes.push_back(dsize);
			}
			predicates[i] = ps;
		}
		if(iType==GIBBS){
			isQuery.resize(predicates.size());
			ifstream qfile(qryfile.c_str());
			while(qfile.good()) {
				string ln;
				getline(qfile,ln);
                if(ln.size()<1)
					continue;
				for(int jj=0;jj<predicates.size();jj++) {
					if(predicates[jj]->symbol.compare(ln)==0) {
						//qrypreds.push_back(jj);
						isQuery[jj] = true;
                        break;
					}
				}
			}
            qfile.close();
		}
		for(int i=0;i<predicates.size();i++) {
			predicates[i]->setcsz();
		}
		while(infile.good()) {
			string ln;
			getline(infile,ln);
            //cout<<ln<<endl;
			if(ln.size()<2) {
				continue;
			}
			WClause* wc = new WClause();
			double ngroundings = 1;
			vector<string> toks1;
			LStringConversionUtils::tokenize(ln,toks1,":");
			wc->weight = LStringConversionUtils::toDouble(toks1[0]);
			if (wc->weight >= 10) {
				wc->isHard = true;
			}
			else {
				wc->isHard = false;
			}
			vector<string> toks2;
			LStringConversionUtils::tokenize(toks1[1],toks2," v ");
			vector<Term*> vartoterm;
			vector<string> keys;
			for(int j=0;j<toks2.size();j++) {
                //cout<<toks2[j]<<endl;
				if(toks2[j][0]=='!') {
					wc->sign.push_back(false);
				}
				else {
					wc->sign.push_back(true);
				}
				string symbol;
				vector<string> terms;
				getpredterms(toks2[j],symbol,terms);
				int atid = 0;
				vector<int> atdomsizes;
				for(int k=0;k<predicates.size();k++) {
					if(predicates[k]->symbol==symbol) {
						atid = predicates[k]->id;
						atdomsizes = predicates[k]->domsizes;
						break;
					}
				}
				Atom* atom = new Atom();
				for(int k=0;k<terms.size();k++) {
					vector<string>::iterator tp = find(keys.begin(),keys.end(),terms[k]);
					if(tp == keys.end()) {
						Term* term = new Term();
						term->groundval = -1;
						term->domsize = atdomsizes[k];
						ngroundings *= term->domsize; 
						atom->terms.push_back(term);
						//vartoterm.insert(pair<string,Term*>(terms[k],term));
						keys.push_back(terms[k]);
						vartoterm.push_back(term);
					} else {
						Term* term = vartoterm[tp-keys.begin()];
						atom->terms.push_back(term);
					}
				}
				atom->predid = atid;
				wc->atoms.push_back(atom);
				wc->numgroundings = ngroundings;
			}
			clauses.push_back(wc);
		}
		infile.close();
		assignments.resize(predicates.size());
		numsamples.resize(predicates.size());
		assignmentsevid.resize(predicates.size());
		isEvidence.resize(predicates.size());
#ifdef _SOFTEVID_
		softevidences.resize(predicates.size());
#endif
		for(int i=0;i<predicates.size();i++) {
			int dsize = 1;
			for(int j=0;j<predicates[i]->domsizes.size();j++) {
				dsize *= predicates[i]->domsizes[j];
			}
			if(iType==GIBBS) {
				assignments[i].resize(dsize);
				numsamples[i].resize(dsize);
				assignmentsevid[i].resize(dsize);
				isEvidence[i].resize(dsize);
#ifdef _SOFTEVID_
				softevidences[i].resize(dsize);
#endif
			}
			for(int j=0;j<assignments[i].size();j++) {
				double r = (double)rand()/RAND_MAX;
				if(r<0.5) {
					assignments[i][j] = 0;
				} else {
					assignments[i][j] = 1;
				}
			}
		}
#ifdef _SOFTEVID_
		readSoftEvids(evidfile);
#else
		readHardEvids(evidfile);
#endif
	}

	void readSoftEvids(string evidfile) {
		ifstream infile1(evidfile.c_str());
		while (infile1.good()) {
			string ln;
			getline(infile1, ln);
			if (ln.size()<2) {
				continue;
			}
			string symbol;
			vector<string> terms;
			double wt;
			getsoftevid(ln, symbol, terms, wt);
			vector<int> termsint;
			LStringConversionUtils::toIntArr(terms, termsint);
			int atid = -1;
			for (int k = 0; k<predicates.size(); k++) {
				if (predicates[k]->symbol.compare(symbol) == 0) {
					atid = predicates[k]->id;
					break;
				}
			}
			if (atid == -1)
				continue;
			int idx = getaddr(atid, termsint);
			softevidences[atid][idx] = wt;
			if (wt == 1) {
				isEvidence[atid][idx] = true;
			}
		}
		infile1.close();
	}

	void readHardEvids(string evidfile) {
		ifstream infile1(evidfile.c_str());
		while (infile1.good()) {
			string ln;
			getline(infile1, ln);
			if (ln.size()<2) {
				continue;
			}
			string symbol;
			vector<string> terms;
			getpredterms(ln, symbol, terms);
			vector<int> termsint;
			LStringConversionUtils::toIntArr(terms, termsint);
			int atid = -1;
			for (int k = 0; k<predicates.size(); k++) {
				if (predicates[k]->symbol.compare(symbol) == 0) {
					atid = predicates[k]->id;
					break;
				}
			}
			if (atid == -1)
				continue;
			int idx = getaddr(atid, termsint);
			if (ln[0] == '!') {
				assignments[atid][idx] = 0;
				assignmentsevid[atid][idx] = -1;
			}
			else {
				assignments[atid][idx] = 1;
				assignmentsevid[atid][idx] = -2;
			}
			isEvidence[atid][idx] = true;
		}
		infile1.close();
		nevids.resize(predicates.size());
		for (int i = 0; i<assignmentsevid.size(); i++) {
			for (int j = 0; j<assignmentsevid[i].size(); j++) {
				if (assignmentsevid[i][j] == 0) {
					nevids[i].push_back(j);
				}
			}
		}
	}


	void readConstraints(string fname) {
		ifstream infile(fname.c_str());
		while (infile.good()) {
			string ln;
			getline(infile, ln);
			if (ln.size()<2)
				continue;
			vector<string> tokens;
			LStringConversionUtils::tokenize(ln, tokens, ",");
			vector<int> tmp;
			for (int i = 0; i < tokens.size(); i++) {
				for (int j = 0; j < predicates.size(); j++) {
					if (predicates[j]->symbol.compare(tokens[i]) == 0) {
						tmp.push_back(j);
						break;
					}
				}
			}
			constraints.push_back(tmp);
		}
	}

		void countInData(int bound) {
			for(int i=0;i<models.size();i++) {
				double estimate;
				vector<int> bw_order;
				vector<set<int> > bw_clusters; 
				int curr_max_cluster_size = models[i].variables.size();
				models[i].getMinFillOrdering(bw_order, bw_clusters,estimate, curr_max_cluster_size);		
				BE* be = new BE();
				be->ibound = bound;
				be->initandstore(models[i].variables, models[i].functions, bw_order);
				datacounts.push_back(totalgroundings[i] - be->currPRVal);
				delete be;
			}
		}


	void initinfer() {
		probabilities.resize(predicates.size());
		currentestimates.resize(predicates.size());
		for(int i=0;i<predicates.size();i++) {
			//bool isquery = false;
			//if(find(qrypreds.begin(),qrypreds.end(),i)!=qrypreds.end()) {
			if(isQuery[i]) {
				//isquery=true;
				//nevids[i].clear();
				//nevids[i].resize(assignments[i].size());
				probabilities[i].resize(assignments[i].size());
				//currentestimates[i]=vector<double>(assignments[i].size(),0.5);

//#ifdef _SOFTEVID_
//				
//				for(int j=0;j<probabilities[i].size();j++) {
//					double p = softevidences[i][j];
//					probabilities[i][j] = p;
					//currentestimates[i].push_back(p);
//				}
//#else
//				for (int j = 0; j<probabilities[i].size(); j++) {
//					probabilities[i][j] = 0;
					//learning only
					/*if(isquery) {
					nevids[i][j] = j;
					}*/

					//currentestimates[i].push_back(p);
//				}
//#endif
			}
		}

	}

	void compileLearning(int bound) {
		for (int i = 0; i<assignmentsevid.size(); i++) {
			for (int j = 0; j<assignmentsevid[i].size(); j++) {
				if (assignmentsevid[i][j] == 0) {
					//closed world counting
					assignments[i][j] = 0;
				}
			}
		}
		for (int i = 0; i<clauses.size(); i++) {
			createModels(i);
			cout << "created model " << i << endl;
		}
		//cout<<"Exact counting..."<<endl;
		countInData(bound);
		//cout<<"Init Learning..."<<endl;
		initinfer();
		for (int i = 0; i<clauses.size(); i++) {
			initfunctions(i);
			cout << "Init model " << i << endl;
		}
		for (int i = 0; i < nevids.size();i++) {
			//ignore evidence on all qry atoms
			//All qry atoms are NE
			//ALL non-qry atoms are closed world
			nevids[i].clear();
			if (isQuery[i]) {
				for (int j = 0; j < assignments[i].size(); j++) {
					nevids[i].push_back(j);
					double r = rand() / (double)RAND_MAX;
					if (r < 0.5) {
						assignments[i][j] = 0;
					} 
					else {
						assignments[i][j] = 1;
					}
				}
			}
		}
		for (int i = 0; i < clauses.size(); i++) {
			initfunctions(i);
			cout << "Set model " << i << endl;
		}
	}

	void compileall(int bound) {
		for(int i=0;i<clauses.size();i++) {
			createModels(i);
			cout<<"created model "<<i<<endl;
		}
		//cout<<"Exact counting..."<<endl;
		countInData(bound);
		//cout<<"Init Learning..."<<endl;
		initinfer();
		for(int i=0;i<clauses.size();i++) {
			initfunctions(i);
			cout<<"Init model "<<i<<endl;
		}
	}

	void createModels(int clausenum) {
		double totgrnds = 1;
		vector<Term*> uniqueterms;
		vector<int> varidx;
		int numfuncs = 0;
		vector<int> constantatoms;
		for(int i=0;i<clauses[clausenum]->atoms.size();i++) {
			numfuncs++;
			bool isconstant = true;
			for(int k=0;k<clauses[clausenum]->atoms[i]->terms.size();k++) {
				vector<Term*>::iterator it = find(uniqueterms.begin(),uniqueterms.end(),clauses[clausenum]->atoms[i]->terms[k]);
				if(it == uniqueterms.end()) {
					uniqueterms.push_back(clauses[clausenum]->atoms[i]->terms[k]);
					varidx.push_back(uniqueterms.size()-1);
					totgrnds *= clauses[clausenum]->atoms[i]->terms[k]->domsize;
				} else {
					varidx.push_back(it-uniqueterms.begin());
				}
			}
			if(isconstant) {
				constantatoms.push_back(i);
			}
		}

		GM gm;
		gm.mult_factor = 1.0;
		gm.type = MARKOV;
		gm.variables = vector<Variable*>(uniqueterms.size());
		for(int i=0;i<uniqueterms.size();i++) {
			int domsize = 1;
			if(uniqueterms[i]->groundval==-1) {
				domsize = uniqueterms[i]->domsize;
			}
			vector<int> domain(domsize);
			for (int j = 0; j < domsize; j++) {
				domain[j] = j;
			}
			gm.variables[i] = new Variable(i, domain);
		}
		gm.functions = vector<Function*>(numfuncs);
		int fid = 0;
		int vidx = 0;
		for(int i=0;i<clauses[clausenum]->atoms.size();i++) {
			vector<Variable*> scope;
			int num_probabilities = 1;
			for(int k=0;k<clauses[clausenum]->atoms[i]->terms.size();k++) {
				scope.push_back(gm.variables[varidx[vidx]]);
				num_probabilities *= gm.variables[varidx[vidx]]->domain_size();
				vidx++;
			}
			gm.functions[fid] = new Function(fid, scope);
			gm.functions[fid]->predid = clauses[clausenum]->atoms[i]->predid;
			vector<Variable*> origorder;
			for(int j=0;j<gm.functions[fid]->variables().size();j++) {
				origorder.push_back(gm.functions[fid]->variables()[j]);
			}
			Function* cpt = gm.functions[fid];
			sort(cpt->variables().begin(), cpt->variables().end(),
					less_than_comparator_variable);
			for(int j=0;j<cpt->variables().size();j++) {
				for(int k=0;k<origorder.size();k++) {
					if(cpt->variables()[j]==origorder[k]){
						cpt->relationalOrder.push_back(k);
						break;
					}
				}
			}
			cpt->table() = vector<Double>(num_probabilities);
			if(iType==WS) {
				cpt->isEvid = vector<bool>(num_probabilities);
				cpt->addressTable = vector<int>(num_probabilities);
			}
			for(int j=0;j<num_probabilities;j++) {
				Variable::setAddressVIB(scope, j);
				int entry = Variable::getAddress(cpt->variables());
				vector<int> tgrd(clauses[clausenum]->atoms[i]->terms.size());
				for(int k=0;k<clauses[clausenum]->atoms[i]->terms.size();k++) {
					if(clauses[clausenum]->atoms[i]->terms[k]->groundval!=-1) {
						tgrd[k] = clauses[clausenum]->atoms[i]->terms[k]->groundval;
					} else {
						tgrd[k] = scope[k]->addr_value();
					}
				}
				int addr = getaddr(clauses[clausenum]->atoms[i]->predid,tgrd);
				int val = assignments[clauses[clausenum]->atoms[i]->predid][addr];
				if(val < 0) {
					cpt->table()[entry] = val;
				} else {
					if(clauses[clausenum]->sign[i]) {
						cpt->table()[entry] = 1 - val;
					} else {
						cpt->table()[entry] = val;
					}
				}
			}
			fid++;
		}
		models.push_back(gm);
		totalgroundings.push_back(clauses[clausenum]->numgroundings);
	}

	

	void initfunctions(int clausenum) {
		double totgrnds = 1;
		vector<Term*> uniqueterms;
		vector<int> varidx;
		int numfuncs = 0;
		vector<int> constantatoms;
		for(int i=0;i<clauses[clausenum]->atoms.size();i++) {
			numfuncs++;
			bool isconstant = true;
			for(int k=0;k<clauses[clausenum]->atoms[i]->terms.size();k++) {
				vector<Term*>::iterator it = find(uniqueterms.begin(),uniqueterms.end(),clauses[clausenum]->atoms[i]->terms[k]);
				if(it == uniqueterms.end()) {
					uniqueterms.push_back(clauses[clausenum]->atoms[i]->terms[k]);
					varidx.push_back(uniqueterms.size()-1);
					totgrnds *= clauses[clausenum]->atoms[i]->terms[k]->domsize;
				} else {
					varidx.push_back(it-uniqueterms.begin());
				}
			}
			if(isconstant) {
				constantatoms.push_back(i);
			}
		}

		int fid = 0;
		int vidx = 0;
		for(int i=0;i<clauses[clausenum]->atoms.size();i++) {
			vector<Variable*> scope;
			int num_probabilities = 1;
			for(int k=0;k<clauses[clausenum]->atoms[i]->terms.size();k++) {
				scope.push_back(models[clausenum].variables[varidx[vidx]]);
				num_probabilities *= models[clausenum].variables[varidx[vidx]]->domain_size();
				vidx++;
			}
			Function* cpt = models[clausenum].functions[fid];
			sort(cpt->variables().begin(), cpt->variables().end(),
					less_than_comparator_variable);
			for(int j=0;j<num_probabilities;j++) {
				Variable::setAddressVIB(scope, j);
				int entry = Variable::getAddress(cpt->variables());
				vector<int> tgrd(clauses[clausenum]->atoms[i]->terms.size());
				for(int k=0;k<clauses[clausenum]->atoms[i]->terms.size();k++) {
					if(clauses[clausenum]->atoms[i]->terms[k]->groundval!=-1) {
						tgrd[k] = clauses[clausenum]->atoms[i]->terms[k]->groundval;
					} else {
						tgrd[k] = scope[k]->addr_value();
					}
				}
				int addr = getaddr(clauses[clausenum]->atoms[i]->predid,tgrd);
				int val = assignments[clauses[clausenum]->atoms[i]->predid][addr];
				if(clauses[clausenum]->sign[i]) {
					cpt->table()[entry] = 1 - val;
				} else {
					cpt->table()[entry] = val;
				}
			}
			fid++;
		}
	}
	/*
	void countInData() {
		for(int i=0;i<models.size();i++) {
			double estimate;
			vector<int> bw_order;
			vector<set<int> > bw_clusters; 
			int curr_max_cluster_size = models[i].variables.size();
			models[i].getMinFillOrdering(bw_order, bw_clusters,estimate, curr_max_cluster_size);		
			BE* be = new BE();
			be->init(models[i].variables, models[i].functions, bw_order);
			datacounts.push_back(totalgroundings[i] - be->currPRVal);
			delete be;
		}
	}
	*/
	void initGibbs(int bound) {
		for(int i=0;i<models.size();i++) {
			double estimate;
			vector<int> bw_order;
			vector<set<int> > bw_clusters; 
			int curr_max_cluster_size = models[i].variables.size();
			models[i].getMinFillOrdering(bw_order, bw_clusters,estimate, curr_max_cluster_size);		
			BE* be = new BE();
			be->ibound = bound;
			be->initandstore(models[i].variables, models[i].functions, bw_order);
			btrees.push_back(be);
			currentcounts.push_back(totalgroundings[i] - be->currPRVal);
		}
	}


	//Gibbs
	void gibbsflip(int flippedPred,int nidx,int iter,bool doinfer=false) {
		//int ftableidx = getFunctionIndex(pidx,nidx);
		double negwt = 0;
		double poswt = 0;
		vector<int> grounding;
		int flippedassg = 1 - assignments[flippedPred][nidx];
		getgrounding(flippedPred,nidx,grounding);
		vector<vector<int> > flippedPList(models.size());
		vector<vector<int> > flippedEList(models.size());
		for(int btreeidx = 0;btreeidx<models.size();btreeidx++) {
			vector<int> fliptables;
			vector<int> flipentries;
			vector<Variable*> nf;
			for(int j=0;j<models[btreeidx].functions.size();j++) {
				if(flippedPred == models[btreeidx].functions[j]->predid) {
					fliptables.push_back(j);
					//compute the entry to flip
					for(int k=0;k<models[btreeidx].functions[j]->relationalOrder.size();k++){
						int v_id = models[btreeidx].functions[j]->relationalOrder[k];
						models[btreeidx].functions[j]->variables()[v_id]->addr_value() = grounding[k];
						models[btreeidx].functions[j]->variables()[v_id]->value() = grounding[k];
						if(find(models[btreeidx].functions[j]->variables()[v_id]->sj_values.begin(),
							models[btreeidx].functions[j]->variables()[v_id]->sj_values.end(),grounding[k])==
							models[btreeidx].functions[j]->variables()[v_id]->sj_values.end()) {
								models[btreeidx].functions[j]->variables()[v_id]->sj_values.push_back(grounding[k]);
						}
					}
					int v_entry = Variable::getAddress(models[btreeidx].functions[j]->variables());
					models[btreeidx].functions[j]->table()[v_entry] = 1 - models[btreeidx].functions[j]->table()[v_entry];
					flipentries.push_back(v_entry);
					if(nf.size()==0) {
						do_set_union(models[btreeidx].functions[j]->variables(),nf,nf,less_than_comparator_variable);
					} else {
						do_set_intersection(models[btreeidx].functions[j]->variables(),nf,nf,less_than_comparator_variable);
					}
				}
			}
			if(fliptables.size()==0) {
				continue;
			}
			//current state's weight
			if(assignments[flippedPred][nidx]==0) {
				negwt += clauses[btreeidx]->weight*(totalgroundings[btreeidx]-btrees[btreeidx]->currPRVal);
			} else {
				poswt += clauses[btreeidx]->weight*(totalgroundings[btreeidx]-btrees[btreeidx]->currPRVal);
			}
			//count for the new state after flip
			//non-selfjoins
			if(fliptables.size() == 1) {
				btrees[btreeidx]->updateBuckets(flippedPred);
			} else {
				if(nf.size()==0) {
					for(int j=0;j<models[btreeidx].variables.size();j++) {
						models[btreeidx].variables[j]->value() = INVALID_VALUE;
					}
					btrees[btreeidx]->updateBuckets(flippedPred);
				} else {
					btrees[btreeidx]->updateBuckets(flippedPred,true);
				}
			}

			double cnt = btrees[btreeidx]->currPRVal;
			//cout<<"wt="<<cnt<<endl;
			if(flippedassg==1) {
				poswt += clauses[btreeidx]->weight*(totalgroundings[btreeidx] - cnt);
			} else {
				negwt += clauses[btreeidx]->weight*(totalgroundings[btreeidx] - cnt);
			}
			//reset variables
			/*for(int j=0;j<fliptables.size();j++){
				int fidx = fliptables[j];
				for(int k=0;k<models[btreeidx].functions[fidx]->variables().size();k++) {
					models[btreeidx].functions[fidx]->variables()[k]->value() = INVALID_VALUE;
					models[btreeidx].functions[fidx]->variables()[k]->sj_values.clear();
					models[btreeidx].functions[fidx]->variables()[k]->tmp_domain.clear();
				}
			}*/
			for(int j=0;j<models[btreeidx].variables.size();j++) {
				models[btreeidx].variables[j]->value() = INVALID_VALUE;
				models[btreeidx].variables[j]->sj_values.clear();
			}
			//store state
			flippedPList[btreeidx] = fliptables;
			flippedEList[btreeidx] = flipentries;
		}
		//sample
		int new_value;
		int old_assg = assignments[flippedPred][nidx];
		double z;
		//poswt = poswt + softevidences[flippedPred][nidx];
		if (poswt > negwt) {
			z = log(1 + exp(negwt - poswt)) + poswt;
		} else {
			z = negwt + log(1 + exp(poswt - negwt));
		}
		//double prob = negwt - z;
		//cout<<exp(prob);
		double prob = exp(poswt - z);
		double r = rand()/(double)RAND_MAX;
		//if(log(r) < prob) {
		if (r < prob) {
			new_value = 1;
			assignments[flippedPred][nidx] = 1;
		} else {
			new_value = 0;
			assignments[flippedPred][nidx] = 0;
		}
		if(doinfer && iter > burnin &&  isQuery[flippedPred]) {
			probabilities[flippedPred][nidx] += prob;
			numsamples[flippedPred][nidx]++;
			//for(int q=0;q<qrypreds.size();q++) {
			//	int p = qrypreds[q];
				//if(flippedPred==p) {
					//currentestimates[p][nidx] = exp(poswt-z);
					//probabilities[flippedPred][nidx] += prob;
					//numsamples[flippedPred][nidx]++;
				//}
				//for(int n=0;n<probabilities[p].size();n++) {
				//	probabilities[p][n] += currentestimates[p][n];
				//}			
			//}
		}		
		bool revert=false;
		if(assignments[flippedPred][nidx]==old_assg) {
			revert = true;
		}
		//Revert the junction tree state and all the flipped variables
		for(int i=0;i<flippedPList.size();i++) {
			if(flippedPList[i].size()==0) {
				continue;
			}
			btrees[i]->revertBuckets(revert);
			if(revert) {
				for(int j=0;j<flippedPList[i].size();j++) {
					int idx = flippedPList[i][j];
					int e_idx = flippedEList[i][j];
					models[i].functions[idx]->table()[e_idx] = 1 - models[i].functions[idx]->table()[e_idx];
				}
			}
		}
	}


	void clamp(int flippedPred, int nidx) {
		//int ftableidx = getFunctionIndex(pidx,nidx);
		vector<int> grounding;
		int flippedassg = 1 - assignments[flippedPred][nidx];
		getgrounding(flippedPred, nidx, grounding);
		for (int btreeidx = 0; btreeidx<models.size(); btreeidx++) {
			vector<int> fliptables;
			vector<int> flipentries;
			vector<Variable*> nf;
			for (int j = 0; j<models[btreeidx].functions.size(); j++) {
				if (flippedPred == models[btreeidx].functions[j]->predid) {
					fliptables.push_back(j);
					//compute the entry to flip
					for (int k = 0; k<models[btreeidx].functions[j]->relationalOrder.size(); k++) {
						int v_id = models[btreeidx].functions[j]->relationalOrder[k];
						models[btreeidx].functions[j]->variables()[v_id]->addr_value() = grounding[k];
						models[btreeidx].functions[j]->variables()[v_id]->value() = grounding[k];
						if (find(models[btreeidx].functions[j]->variables()[v_id]->sj_values.begin(),
							models[btreeidx].functions[j]->variables()[v_id]->sj_values.end(), grounding[k]) ==
							models[btreeidx].functions[j]->variables()[v_id]->sj_values.end()) {
							models[btreeidx].functions[j]->variables()[v_id]->sj_values.push_back(grounding[k]);
						}
					}
					int v_entry = Variable::getAddress(models[btreeidx].functions[j]->variables());
					models[btreeidx].functions[j]->table()[v_entry] = 1 - models[btreeidx].functions[j]->table()[v_entry];
					flipentries.push_back(v_entry);
					if (nf.size() == 0) {
						do_set_union(models[btreeidx].functions[j]->variables(), nf, nf, less_than_comparator_variable);
					}
					else {
						do_set_intersection(models[btreeidx].functions[j]->variables(), nf, nf, less_than_comparator_variable);
					}
				}
			}
			if (fliptables.size() == 0) {
				continue;
			}
			//count for the new state after flip
			//non-selfjoins
			if (fliptables.size() == 1) {
				btrees[btreeidx]->updateBuckets(flippedPred);
			}
			else {
				if (nf.size() == 0) {
					for (int j = 0; j<models[btreeidx].variables.size(); j++) {
						models[btreeidx].variables[j]->value() = INVALID_VALUE;
					}
					btrees[btreeidx]->updateBuckets(flippedPred);
				}
				else {
					btrees[btreeidx]->updateBuckets(flippedPred, true);
				}
			}

			for (int j = 0; j<models[btreeidx].variables.size(); j++) {
				models[btreeidx].variables[j]->value() = INVALID_VALUE;
				models[btreeidx].variables[j]->sj_values.clear();
			}
		}
		assignments[flippedPred][nidx] = 1;
		if (isQuery[flippedPred]) {
			probabilities[flippedPred][nidx] += 1;
			numsamples[flippedPred][nidx]++;
		}
	}


	void block(vector<int> flippedPreds, int nidx, int iter, bool doinfer=false) {
		//int ftableidx = getFunctionIndex(pidx,nidx);
		vector<int> grounding;
		//int flippedassg = 1 - assignments[flippedPred][nidx];
		vector<double> probs;
		for (int t = 0; t < flippedPreds.size(); t++) {
			double prob = 0;
			//bool first = true;
			//if (t > 0)
				//first = false;
			vector<int> topropagate;
			for (int t1 = 0; t1 < flippedPreds.size(); t1++) {
				int flippedassg = 1;
				if (t == t1) {
					flippedassg = 0;
				}
				int flippedPred = flippedPreds[t1];

				getgrounding(flippedPred, nidx, grounding);
				for (int btreeidx = 0; btreeidx < models.size(); btreeidx++) {
					bool found = false;
					for (int j = 0; j < models[btreeidx].functions.size(); j++) {
						if (flippedPred == models[btreeidx].functions[j]->predid) {
							//compute the entry to flip
							for (int k = 0; k < models[btreeidx].functions[j]->relationalOrder.size(); k++) {
								int v_id = models[btreeidx].functions[j]->relationalOrder[k];
								models[btreeidx].functions[j]->variables()[v_id]->addr_value() = grounding[k];
								models[btreeidx].functions[j]->variables()[v_id]->value() = grounding[k];
							}
							int v_entry = Variable::getAddress(models[btreeidx].functions[j]->variables());
							if (!isEvidence[flippedPred][nidx]) {
								models[btreeidx].functions[j]->table()[v_entry] = flippedassg;
								found = true;
							}
						}
					}
					if (found) {
						if (find(topropagate.begin(), topropagate.end(), btreeidx) == topropagate.end()) {
							topropagate.push_back(btreeidx);
						}
					}
				}
			}
			for (int j = 0; j < topropagate.size(); j++) {
				int btreeidx = topropagate[j];
				int d1 = btrees[btreeidx]->currPRVal;
				btrees[btreeidx]->propagateall();
				if (btrees[btreeidx]->currPRVal > 1) {
					int dum = 1;
				}
				prob = prob + (totalgroundings[btreeidx] - btrees[btreeidx]->currPRVal)*clauses[btreeidx]->weight;
			}
			probs.push_back(prob);
		}
		//cout << "********" << endl;
		//cout << assignments[0][0] << " " << assignments[1][0] << " " << assignments[2][0] << " " << assignments[3][0] << endl;
		//cout << flippedPreds[0] << " " << flippedPreds[1] << endl;
		//cout << probs[0] << " " << probs[1] << endl;
		//cout << "********" << endl;
		double z = probs[0];
		for (int i = 1; i < probs.size(); i++) {
			if (z > probs[i]) {
				z = log(1 + exp(probs[i] - z)) + z;
			}
			else {
				z = probs[i] + log(1 + exp(z - probs[i]));
			}
		}
		int sampled = -1;
		double r = rand() / (double)RAND_MAX;
		double cm = probs[0] - z;
		if (r < exp(cm)) {
			sampled = 0;
		}
		else {
			for (int i = 1; i < probs.size(); i++) {
				double x = probs[i] - z;

				if (cm > x) {
					cm = log(1 + exp(x - cm)) + cm;
				}
				else {
					cm = x + log(1 + exp(cm - x));
				}
				if (r < exp(cm)) {
					sampled = i;
					break;
				}
			}
		}
		if (sampled == -1) {
			sampled = probs.size() - 1;
		}
		//cout << "sampled="<< sampled << endl;
		//cout << flippedPreds[0]<<" "<< flippedPreds[1]<<" "<<exp(probs[0] - z) << " " << exp(probs[1] - z) << endl;

		for (int t1 = 0; t1 < flippedPreds.size(); t1++) {
			int flippedassg = 1;
			if (t1 == sampled) {
				flippedassg = 0;
			}
			int flippedPred = flippedPreds[t1];
			if (isEvidence[flippedPred][nidx])
				continue;
			if (iter > burnin && doinfer && isQuery[flippedPred]) {
				numsamples[flippedPred][nidx]++;
				probabilities[flippedPred][nidx] += exp(probs[t1] - z);
			}
			getgrounding(flippedPred, nidx, grounding);
			for (int btreeidx = 0; btreeidx < models.size(); btreeidx++) {
				bool found = false;
				for (int j = 0; j < models[btreeidx].functions.size(); j++) {
					if (flippedPred == models[btreeidx].functions[j]->predid) {
						found = true;
						//compute the entry to flip
						for (int k = 0; k < models[btreeidx].functions[j]->relationalOrder.size(); k++) {
							int v_id = models[btreeidx].functions[j]->relationalOrder[k];
							models[btreeidx].functions[j]->variables()[v_id]->addr_value() = grounding[k];
							models[btreeidx].functions[j]->variables()[v_id]->value() = grounding[k];
						}
						int v_entry = Variable::getAddress(models[btreeidx].functions[j]->variables());
						models[btreeidx].functions[j]->table()[v_entry] = flippedassg;
					}
				}
				if (found) {
					btrees[btreeidx]->propagateall();
				}
				for (int j = 0; j<models[btreeidx].variables.size(); j++) {
					models[btreeidx].variables[j]->value() = INVALID_VALUE;
					models[btreeidx].variables[j]->sj_values.clear();
				}
				assignments[flippedPred][nidx] = 1-flippedassg;
			}
		}
	}

	//LEARNING USING CD
	void sample(string mlnfile,string evidfile,string qryfile,int niters,int bound,
		string outmln,string csfile) {
		/*
			assignments[0][0] = 1;assignments[0][1] = 0;assignments[0][2] = 0;
			assignments[1][0] = 1;assignments[1][1] = 0;assignments[1][2] = 1;
			assignments[1][3] = 1;assignments[1][4] = 0;assignments[1][5] = 0;
			assignments[1][6] = 1;assignments[1][7] = 0;assignments[1][8] = 1;
			*/
		string outconv = "convergence.dat";
		iType = GIBBS;
		itersperupdate = 10;
		init(mlnfile,evidfile,qryfile);
		alpha = 0.001;
		readConstraints(csfile);
		cout<<"done init..."<<endl;
		compileLearning(bound);
		ofstream wfile("wlog");

		cout<<"Init Gibbs..."<<endl;
		initGibbs(bound);
		cout<<"Starting Learning Iters..."<<endl;
		//cout<<"Done compiling.."<<endl;
		int iters = 0;
		avgcounts.resize(models.size());
		averageweights.resize(models.size());
		oldweights.resize(models.size());
		weights.resize(models.size());
		oldgradient.resize(models.size());
		vector<double> delta_pred(models.size());
		vector<double> sqavgcounts(models.size());
		double lltreshold = 0.0001;
		momentum = 0;
		cglambda = 100;
		//init weights
                //cout<<"Initial Weights"<<endl;
		for(int j=0;j<avgcounts.size();j++) {
			/*clauses[j]->weight = alpha*(datacounts[j]/totalgroundings[j]);
			double r = (double)rand()/RAND_MAX;
			if(datacounts[j]/totalgroundings[j] < 0.1)
				clauses[j]->weight = -1*alpha*r;
			else
				clauses[j]->weight = alpha*r;
                        */
            //cout<<clauses[j]->weight<<" ";
			weights[j] = clauses[j]->weight;
			oldweights[j] = clauses[j]->weight;
		}
                cout<<endl;
                cout<<"Learning rate="<<alpha<<endl;
		cout<<"Samples per update="<<itersperupdate<<endl;
		time_t start;
		time(&start);
		int t=1;
		bool done = false;
		ofstream fs1(outconv.c_str());
                int effiters = 0;
		while(!done) {
		//for(int i=1;i<niters;i++) {
                        //cout<<predicates.size()<<endl;
			int pidx = rand()%predicates.size();
			if(nevids[pidx].size()==0)
				continue;
			int nidx = rand()%nevids[pidx].size();
			int blocked = -1;
			for (int i = 0; i < constraints.size(); i++) {
				for (int j = 0; j < constraints[i].size(); j++) {
					if (constraints[i][j] == pidx) {
						blocked = i;
						break;
					}
				}
				if (blocked != -1)
					break;
			}
			if (blocked == -1)
				gibbsflip(pidx, nevids[pidx][nidx], t);
			else
				block(constraints[blocked], nidx, t);
			//cout<<"here2"<<endl;
			//int thin = 1;
			//if(i%thin==0){
			samplecnt++;
			vector<int> tmpcounts(avgcounts.size());
			for(int j=0;j<avgcounts.size();j++) {
				avgcounts[j] += (totalgroundings[j]-btrees[j]->currPRVal);
				tmpcounts[j] = (totalgroundings[j]-btrees[j]->currPRVal);
				sqavgcounts[j] += (totalgroundings[j]-btrees[j]->currPRVal)*(totalgroundings[j]-btrees[j]->currPRVal);
				//clauses[j]->weight += (datacounts[j]-totalgroundings[j]+btrees[j]->currPRVal)*alpha - lambda*clauses[j]->weight;
				//avgcounts[j] = 0;
				if (j >4000) {
					cout << datacounts[j] << "," << totalgroundings[j] - btrees[j]->currPRVal << endl;
				}
				//cout<<clauses[j]->weight<<",";
			}
			allsamplecounts.push_back(tmpcounts);
			//}
			if(samplecnt==itersperupdate) {
           		effiters++;
				//cout<<"Running DN"<<endl;
				int numWeights=avgcounts.size();
				vector<double> dw(numWeights);
				vector<double> gradient(numWeights);
				double realdist = 1.0;
      			double preddist = 1.0;
				bool backtrack = false;
      			if (effiters > 1)
				{
					vector<double> dist(numWeights, 0.0);
					for (int ii = 0; ii < numWeights; ii++)
					  dist[ii] = weights[ii] - oldweights[ii];
					// Predicted change is quadratic approximation
					vector<double> avgPred(numWeights);
					for (int ii = 0; ii < numWeights; ii++)
					  avgPred[ii] = oldgradient[ii] + delta_pred[ii]/2.0;
					preddist = dotprod(avgPred, dist);
					// Real change is lower bound on actual change
					realdist = dotprod(gradient, dist);
					// Adjust lambda using technique of (Fletcher, 1987)
					double delta = realdist/preddist;

					if (preddist == 0)
					  cglambda /= 4;

					if (preddist != 0 && delta > 0.75)
					  cglambda /= 2;
					// if (delta < 0.25)   // Update schedule from (Fletcher, 1987)
					if (delta < 0.0)       // Gentler update schedule, to handle noise
					{
					  if (cglambda * 4 > cglambda)
						cglambda = cglambda;
					  else
						cglambda *= 4;
					  //cout<<"Backtracking!!"<<endl;
					}
				}
			    cout<<"Iteration "<<effiters<<endl;
				//ofstream tmpfs("currentweights.txt");
				//DN
				bool allzeros = true;
				for(int v1=0;v1<numWeights;v1++) {
					double varn = 0;
					//for(int v2=0;v2<allsamplecounts.size();v2++) {
					 //cout<<"SQ="<<sqavgcounts[v1]<<" "<<avgcounts[v1]<<" "<<samplecnt<<endl;
					  varn = sqavgcounts[v1]/(double)samplecnt - (avgcounts[v1]/(double)samplecnt)*(avgcounts[v1]/(double)samplecnt);
					//}
			  		double grad = ((avgcounts[v1]/itersperupdate)-datacounts[v1])/(double)datacounts[v1];
                    if(abs(grad)<0.0000001)
                    	grad = 0;
                    //cout<<"varn="<<varn<<endl;
					if(varn==0) {
						//done = true;
						dw[v1]=0;
						//cout<<sqavgcounts[v1]<<" "<<avgcounts[v1]<<endl;
					} else {
					allzeros = false;
                    dw[v1]=-1*grad/varn;
					}
			 	 	gradient[v1]=grad;
					//cout<<dw[v1]<<" "<<gradient[v1]<<endl;
				}
				if(allzeros) {
					cout << "allzeros" << endl;
					done=true;
					break;
				}
				//Hessian vector product
				double sumVN = 0;
                vector<double> sumN(avgcounts.size());
                vector<double> sumNiVN(avgcounts.size());
		    	for (int s = 0; s < allsamplecounts.size(); s++) {
		      		// Compute v * n
		      		double vn = 0;
		      		for (int ii = 0; ii < allsamplecounts[s].size(); ii++)
						vn += dw[ii] * allsamplecounts[s][ii];
		      		// Tally v*n, n_i, and n_i v*n
		      		sumVN += vn;
		      		for (int ii = 0; ii < allsamplecounts[s].size(); ii++)
		      		{
						sumN[ii]    += allsamplecounts[s][ii];
						sumNiVN[ii] += allsamplecounts[s][ii] * vn;
		      		}
		    	}
				//cout<<"here4"<<endl;
				vector<double> Hd(avgcounts.size());
		    	for (int clauseno = 0; clauseno < avgcounts.size(); clauseno++) {
		      		double E_vn = sumVN/itersperupdate;
		      		double E_ni = sumN[clauseno]/itersperupdate;
		      		double E_nivn = sumNiVN[clauseno]/itersperupdate;
		      		Hd[clauseno] = E_nivn - E_ni * E_vn;
					//cout<<"E="<<sumVN<<" "<<sumN[clauseno]<<" "<<sumNiVN[clauseno]<<endl;
					//cout<<"hd="<<Hd[clauseno]<<" ";
		    	}
				//cout<<endl;
		    //cout<<"here5"<<endl;
		    // Compute step length using trust region approach
		    double dHd = dotprod(dw, Hd);
		    double dd = dotprod(dw, dw);
		    double dg = dotprod(gradient, dw);
		    double alpha = -dg/(dHd + cglambda * dd);
			vector<double> wchange(avgcounts.size());
//			cout<<"alpha="<<alpha<<endl;
			//cout<<dHd<<" "<<dd<<endl;
			for (int w = 0; w < numWeights; w++) {
				if (clauses[w]->isHard)
					continue;
				wchange[w] = dw[w] * alpha + (weights[w] - oldweights[w]) * momentum;
				//cout<<"wchange="<<wchange[w]<<endl;
			}
			//cout<<"here6"<<endl;
		    // Convergence criteria for 2nd order methods:
		    // Stop when the maximum predicted improvement in log likelihood
		    // is very small.
        	double maxchange = -1*dotprod(gradient, wchange);
			if(maxchange<lltreshold){
				done=true;
				break;
			}
       		for (int w = 0; w < numWeights; w++) {
				if (clauses[w]->isHard) {
					wfile << "10 ";
					continue;
				}
				oldweights[w] = weights[w];
				oldgradient[w] = gradient[w];
				weights[w] += wchange[w];
				clauses[w]->weight = weights[w];
				averageweights[w] = ((effiters - 1) * averageweights[w] + weights[w]) / (effiters);
				avgcounts[w] = 0;
				sqavgcounts[w] = 0;
				wfile << averageweights[w] << " ";
			}
			allsamplecounts.clear();
			wfile<<endl;
			wfile.flush();
			{
				double sumVN = 0;
                vector<double> sumN(avgcounts.size());
                vector<double> sumNiVN(avgcounts.size());
		    	for (int s = 0; s < allsamplecounts.size(); s++) {
		      		// Compute v * n
		      		double vn = 0;
		      		for (int ii = 0; ii < allsamplecounts[s].size(); ii++)
						vn += wchange[ii] * allsamplecounts[s][ii];
		      		// Tally v*n, n_i, and n_i v*n
		      		sumVN += vn;
		      		for (int ii = 0; ii < allsamplecounts[s].size(); ii++)
		      		{
						sumN[ii]    += allsamplecounts[s][ii];
						sumNiVN[ii] += allsamplecounts[s][ii] * vn;
		      		}
		    	}
		    	for (int clauseno = 0; clauseno < avgcounts.size(); clauseno++) {
		      		double E_vn = sumVN/itersperupdate;
		      		double E_ni = sumN[clauseno]/itersperupdate;
		      		double E_nivn = sumNiVN[clauseno]/itersperupdate;
		      		delta_pred[clauseno] = E_nivn - E_ni * E_vn;
		    	}
				//cout<<"hereend"<<endl;
			}
			/*	
			for(int j=0;j<avgcounts.size();j++) {
				double orig = clauses[j]->weight;
				clauses[j]->weight -= (((avgcounts[j]/itersperupdate)-datacounts[j])/datacounts[j])*alpha;
				if(clauses[j]->weight > 50) {
					clauses[j]->weight = 50;
				} else if(clauses[j]->weight < -50) {
					clauses[j]->weight = -50;
				}
				double chg = clauses[j]->weight;
				double perchg = (chg-orig)/orig;
				if(abs(perchg)>maxchange) {
					maxchange = abs(perchg);
				}
				cout<<clauses[j]->weight<<" ";
			
				avgcounts[j] = 0;
			}
			*/
			//tmpfs<<endl;
            //cout<<endl;
				samplecnt = 0;
			//tmpfs.close();
			//if(maxchange<0.001) {
			//	done=true;
			//}
			//if(effiters==1000)
				//exit(0);
			}
			if(done) {
				time_t tm;
				time(&tm);
				cout<<"Convergence time="<<difftime(tm,start)<<endl;
				fs1<<difftime(tm,start)<<endl;
				break;

			}
			//cout<<endl;
			time_t currtime;
			time(&currtime);
			int cp = difftime(currtime,start);
			if(cp > niters) {
				cout<<"Convergence time=-1"<<endl;
				fs1<<"-1"<<endl;
				break;
			}
			t++;
		}
		ifstream infile(mlnfile.c_str());
		ofstream ofile(outmln.c_str());
		string ln1;
		getline(infile,ln1);
		ofile << ln1<<endl;
		int cid = 0;
		while(infile.good()) {
			string ln;
			getline(infile,ln);
			if(ln.size()<2) {
				continue;
			}
			vector<string> toks1;
			LStringConversionUtils::tokenize(ln,toks1,":");
			ofile<<averageweights[cid]<<":"<<toks1[1]<<endl;
			cid++;
		}
		infile.close();
		ofile.close();
		fs1.close();
		wfile.close();
	}

	void testing(string mlnfile,string evidfile,string qryfile,int niters,int bound,string outfile) {
		/*
			assignments[0][0] = 1;assignments[0][1] = 0;assignments[0][2] = 0;
			assignments[1][0] = 1;assignments[1][1] = 0;assignments[1][2] = 1;
			assignments[1][3] = 1;assignments[1][4] = 0;assignments[1][5] = 0;
			assignments[1][6] = 1;assignments[1][7] = 0;assignments[1][8] = 1;
			*/
		iType = GIBBS;
		init(mlnfile,evidfile,qryfile);
                vector<vector<int> > counts;
                ifstream qfile(qryfile.c_str());
                while(!qfile.eof()){
                string s;
                qfile>>s;
                cout<<s<<endl;
                if(s.size()<2)
                    continue;
                string s1 = s+".map";
                ifstream ofs(s1.c_str());
                vector<int> tmpcounts;
                while(!ofs.eof()){
                  int i;
                  ofs >> i;
                  cout<<i<<" ";
                  tmpcounts.push_back(i);
                }
                //cout<<endl;
                ofs.close();
                counts.push_back(tmpcounts);
                }
                qfile.close();
		cout<<"done init..."<<endl;
		for(int i=0;i<clauses.size();i++) {
			createModels(i);
		}
		cout<<"Created Models..."<<endl;
		initinfer();
		for(int i=0;i<clauses.size();i++) {
			initfunctions(i);
		}
		cout<<"Init Gibbs..."<<endl;
		initGibbs(bound);

		//cout<<"Done compiling.."<<endl;
		//int i = 0;
		int iters = 0;
		avgcounts.resize(models.size());
		time_t start;
		time(&start);
		//while(true) {
		for(int i=0;i<niters;i++) {
			int pidx = rand()%predicates.size();
			if(nevids[pidx].size()==0)
				continue;
			int nidx = rand()%nevids[pidx].size();
			cout<<"Iter="<<i<<endl;
			gibbsflip(pidx,nevids[pidx][nidx],i,true);
			//time_t currtime;
			//time(&currtime);
			//int cp = difftime(currtime,start);
			//if(cp > niters) {
				//break;
			//}
			//i++;
		}
		double cll=0;
		double cnt=0;
		/*for(int q=0;q<qrypreds.size();q++) {
			int i = qrypreds[q];
			for(int j1=0;j1<nevids[i].size();j1++) {
				int j = nevids[i][j1];
				probabilities[i][j] = (probabilities[i][j])/(double)(niters-burnin-1);
				if(probabilities[i][j]==1)
					cll += log(0.99)*counts[i][j];
					//cll += log(0.99);
				else if(probabilities[i][j]==0)
					cll += log(0.0001)*counts[i][j];
					//cll += log(0.0001);
				else
					cll += log(probabilities[i][j])*counts[i][j];
					//cll += log(probabilities[i][j]);
				cnt += counts[i][j];
				//cnt++;
			}
		}
		ofstream tmpf(outfile.c_str());
		tmpf<<cll/cnt<<endl;
		tmpf.close();
		*/
	}


	void doSampling(string mlnfile,string evidfile,string qryfile,int niters, int bound, string outfile, string csfile, string logfilename) {
		iType = GIBBS;
		init(mlnfile,evidfile,qryfile);
		readConstraints(csfile);
		cout<<"done init..."<<endl;
		for(int i=0;i<clauses.size();i++) {
			createModels(i);
		}
		cout<<"Created Models..."<<endl;
		initinfer();
		cout<<"Inference Init..."<<endl;
		for(int i=0;i<clauses.size();i++) {
			initfunctions(i);
		}
		cout<<"Init Gibbs..."<<endl;
		//set bound 10
		initGibbs(bound);
		int fileid = 0;
		cout<<"starting Gibbs..."<<endl;
		//cout<<"Done compiling.."<<endl;
		//int i = 0;
		int iters = 0;
		avgcounts.resize(models.size());
		time_t start;
		time(&start);
		time_t lastlog;
		time(&lastlog);
		//cout << softevidences[0][0] << endl;
		//while(true) {
#ifdef _LOGGING_
		ofstream logfile(logfilename.c_str());
		vector<vector<int> > sampleddata;
		for (int i = 0; i < predicates.size(); i++) {
			vector<int> tmp;
			sampleddata.push_back(tmp);
		}
		int numpoints = 100;
		srand(192230943);
		for (int i = 0; i < predicates.size(); i++) {
			if (!isQuery[i])
				continue;
#ifdef _SOFTEVID_
			for (int j = 0; j < numpoints; j++) {
				int n1 = rand() % softevidences[i].size();
				sampleddata[i].push_back(n1);
			}
#else
			if (nevids[i].size() == 0)
				continue;
			for (int j = 0; j < numpoints; j++) {
				int n1 = rand() % nevids[i].size();
				sampleddata[i].push_back(nevids[i][n1]);
			}
#endif
		}
#endif
		int r = rand();
		srand(r + time(NULL));
		int totnevids = 0;
#ifdef _SOFTEVID_
		for (int i = 0; i < softevidences.size(); i++) {
			totnevids += softevidences[i].size();
		}
#else
		for (int i = 0; i < predicates.size(); i++) {
			totnevids = totnevids + nevids[i].size();
		}
#endif
		int maxiters = totnevids * 10000;
		int t = 0;
		while(1) {
		//for(int t=0;t<maxiters;t++) {
			int pidx=0;
			int nidx = 0;
			pidx = rand()%predicates.size();
#ifdef _SOFTEVID_
			nidx = rand()%assignments[pidx].size();
			double r1 = rand() / (double)(RAND_MAX);
			if(r1 < softevidences[pidx][nidx]) {
				if (assignments[pidx][nidx] != 1)
					clamp(pidx, nidx);
				continue;
			}
#else
			if (nevids[pidx].size() == 0)
				continue;
			int n1 = rand() % nevids[pidx].size();
			nidx = nevids[pidx][n1];
#endif
			int blocked = -1;
			for (int i = 0; i < constraints.size(); i++) {
				for (int j = 0; j < constraints[i].size(); j++) {
					if (constraints[i][j] == pidx) {
						blocked = i;
						break;
					}
				}
				if (blocked != -1)
					break;
			}
			if (blocked == -1)
				gibbsflip(pidx, nidx, t, true);
			else
				block(constraints[blocked], nidx,t, true);
#ifdef _LOGGING_
			time_t ntime;
			time(&ntime);
			int ndiff = difftime(ntime, lastlog);
			if (t >= burnin && ndiff > 100) {
				lastlog = ntime;
				/*logfile << t - burnin << " ";
				for (int t1 = 0; t1 < sampleddata.size(); t1++) {
					for (int t2 = 0; t2 < sampleddata[t1].size(); t2++) {
						int sx = sampleddata[t1][t2];
						double p = 0;
						if (numsamples[t1][sx] > 0)
							p = probabilities[t1][sx] / (numsamples[t1][sx]);
						logfile << p << " ";
					}
				}
				logfile << endl;
				logfile.flush();
				*/
				stringstream sh;
				sh << fileid;
				string sname = sh.str() + "-" + outfile;
				writeOutfile(sname);
				fileid++;
			}
#endif
			//if (assignments[0][0] == assignments[2][0] || assignments[1][0] == assignments[3][0]) {
			//	int dum = 1;
			//}
			//cout<<t<<" ";
			//update transitions
			/*if (t > burnin && t%aval == 0) {
				for (int i = 0; i < probabilities.size(); i++) {
					for (int j = 0; j < probabilities[i].size(); j++) {
						vector<int> grd;
						getgrounding(i, j, grd);
						softevidences[i][j] = probabilities[i][j] / (numsamples[i][j]+1);
					}
				}
			}*/
			time_t currtime;
			time(&currtime);
			int cp = difftime(currtime, start);
			if (cp > niters) {
				cout << "Timeout" << endl;
				break;
			}
			t++;
		}
#ifdef _LOGGING_
		logfile.close();
#endif
		//writeOutfile(outfile);
		}

	void writeOutfile(string outfile)
	{
		ofstream out(outfile.c_str());
#ifdef _SOFTEVID_
		for (int i = 0; i < probabilities.size(); i++) {
			if (!isQuery[i])
				continue;
			for (int j = 0; j < probabilities[i].size(); j++) {
				vector<int> grd;
				getgrounding(i, j, grd);
				if (numsamples[i][j] == 0) {
					out << softevidences[i][j] << " " << predicates[i]->symbol << "(";
				}
				else {
					double p = probabilities[i][j] / (numsamples[i][j]);
					out << p << " " << predicates[i]->symbol << "(";
				}
				for (int k = 0; k < grd.size(); k++) {
					out << grd[k];
					if (k != grd.size() - 1)
						out << ",";
				}
				out << ")\n";
			}
#else

		for (int i = 0; i < probabilities.size(); i++) {
			if (!isQuery[i])
				continue;
			for (int j = 0; j < probabilities[i].size(); j++) {
				vector<int> grd;
				getgrounding(i, j, grd);
				if (numsamples[i][j] == 0) {
					if (assignmentsevid[i][j] == -1) {
						out << "0 " << predicates[i]->symbol << "(";
					}
					else if (assignmentsevid[i][j] == -2) {
						out << "1 " << predicates[i]->symbol << "(";
					}
					else {
						out << "0.5 " << predicates[i]->symbol << "(";
					}
				}
				else {
					double p = probabilities[i][j] / (numsamples[i][j]);
					out << p << " " << predicates[i]->symbol << "(";
				}
				for (int k = 0; k < grd.size(); k++) {
					out << grd[k];
					if (k != grd.size() - 1)
						out << ",";
				}
				out << ")\n";
			}
#endif
		}
		out.close();
	}

};
#endif /* MLN_H_ */
