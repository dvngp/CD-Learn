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
	vector<vector<int> > mapassignment;
	vector<vector<int> > nevids;
	vector<vector<double> > probabilities;
	vector<vector<int> > assignmentsevid;
	vector<double> avgcounts;
	int itersperupdate;
	vector<GM> models;
	//vector<JT*> jtrees;
	vector<BE*> btrees;
	vector<double> currentcounts;
	vector<double> datacounts;
	vector<double> totalgroundings;
	vector<vector<double> > currentestimates;
	vector<int> qrypreds;
	int burnin;
	double bestmapcost;
	double gprob;
	InferenceType iType;
	double alpha;
	double lambda;
	int samplecnt;
	
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
			ifstream qfile(qryfile.c_str());
			while(qfile.good()) {
				string ln;
				getline(qfile,ln);
                                if(ln.size()<2)
                                   continue;
				for(int jj=0;jj<predicates.size();jj++) {
					if(predicates[jj]->symbol.compare(ln)==0) {
						qrypreds.push_back(jj);
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
                        cout<<ln<<endl;
			if(ln.size()<2) {
				continue;
			}
			WClause* wc = new WClause();
			double ngroundings = 1;
			vector<string> toks1;
			LStringConversionUtils::tokenize(ln,toks1,":");
			wc->weight = LStringConversionUtils::toDouble(toks1[0]);
			vector<string> toks2;
			LStringConversionUtils::tokenize(toks1[1],toks2," v ");
			vector<Term*> vartoterm;
			vector<string> keys;
			for(int j=0;j<toks2.size();j++) {
                                cout<<toks2[j]<<endl;
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
		assignmentsevid.resize(predicates.size());
		for(int i=0;i<predicates.size();i++) {
			int dsize = 1;
			for(int j=0;j<predicates[i]->domsizes.size();j++) {
				dsize *= predicates[i]->domsizes[j];
			}
			assignments[i].resize(dsize);
			if(iType==GIBBS) {
				assignmentsevid[i].resize(dsize);
			}
			for(int j=0;j<assignments[i].size();j++) {
				/*double r = (double)rand()/RAND_MAX;
				if(r<0.5) {
					assignments[i][j] = 0;
				} else {
					assignments[i][j] = 1;
				}*/
				assignments[i][j] = -1;
			}
		}

		ifstream infile1(evidfile.c_str());
		while(infile1.good()) {
			string ln;
			getline(infile1,ln);
			if(ln.size()<2) {
				continue;
			}
			string symbol;
			vector<string> terms;
			getpredterms(ln,symbol,terms);
			vector<int> termsint;
			LStringConversionUtils::toIntArr(terms,termsint);
			int atid=-1;
			for(int k=0;k<predicates.size();k++) {
				if(predicates[k]->symbol.compare(symbol)==0) {
					atid = predicates[k]->id;
					break;
				}
			}
			if(atid==-1)
				continue;
			int idx = getaddr(atid,termsint);
			if(ln[0]=='!') {
				assignments[atid][idx] = 0;
				assignmentsevid[atid][idx] = -1;
			} else {
				assignments[atid][idx] = 1;
				assignmentsevid[atid][idx] = -2;
			}
		}
		infile1.close();
		nevids.resize(predicates.size());
		for(int i=0;i<assignmentsevid.size();i++) {
			for(int j=0;j<assignmentsevid[i].size();j++) {
				if(assignmentsevid[i][j]==0) {
					nevids[i].push_back(j);
				}
			}
		}
	}

		void countInData(int bound) {
			for(int i=0;i<assignmentsevid.size();i++) {
				for(int j=0;j<assignmentsevid[i].size();j++) {
					if(assignmentsevid[i][j]==0) {
						//closed world counting
						assignments[i][j] = 0;
					}
				}
			}
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
			bool isquery = false;
			if(find(qrypreds.begin(),qrypreds.end(),i)!=qrypreds.end()) {
				isquery=true;
				//nevids[i].clear();
				//nevids[i].resize(assignments[i].size());
				probabilities[i].resize(assignments[i].size());
				currentestimates[i]=vector<double>(assignments[i].size(),0.5);
			}
			for(int j=0;j<assignments[i].size();j++) {
				if(assignmentsevid[i][j]==0) {
					double r = (double)rand()/RAND_MAX;
					if(r<0.5) {
						assignments[i][j] = 0;
					} else {
						assignments[i][j] = 1;
					}
				} else if(isquery){
					probabilities[i][j] = assignmentsevid[i][j];
				}
				
				//learning only
				/*if(isquery) {
					nevids[i][j] = j;
				}*/
				
			}
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
		if (poswt > negwt) {
			z = log(1 + exp(negwt - poswt)) + poswt;
		} else {
			z = negwt + log(1 + exp(poswt - negwt));
		}
		double prob = negwt - z;
		double r = rand()/(double)RAND_MAX;
		if(log(r) < prob) {
			new_value = 0;
			assignments[flippedPred][nidx] = 0;
		} else {
			new_value = 1;
			assignments[flippedPred][nidx] = 1;
		}
		if(doinfer && iter > burnin) {
			for(int q=0;q<qrypreds.size();q++) {
				int p = qrypreds[q];
				if(flippedPred==p)
					currentestimates[p][nidx] = exp(poswt-z);
				for(int n=0;n<nevids[p].size();n++) {
					probabilities[p][nevids[p][n]] += currentestimates[p][nevids[p][n]];
					/*if(assignments[p][nevids[p][n]]==1) {
						probabilities[p][nevids[p][n]] = probabilities[p][nevids[p][n]] + 1;
					}*/
				}
			}
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





	void sample(string mlnfile,string evidfile,string qryfile,int niters,int bound,int updatefreq,string outmln,string outconv,double lrate) {
		/*
			assignments[0][0] = 1;assignments[0][1] = 0;assignments[0][2] = 0;
			assignments[1][0] = 1;assignments[1][1] = 0;assignments[1][2] = 1;
			assignments[1][3] = 1;assignments[1][4] = 0;assignments[1][5] = 0;
			assignments[1][6] = 1;assignments[1][7] = 0;assignments[1][8] = 1;
			*/
		iType = GIBBS;
		itersperupdate = updatefreq;
		init(mlnfile,evidfile,qryfile,lrate);
		cout<<"done init..."<<endl;
		compileall(bound);
		cout<<"Init Gibbs..."<<endl;
		initGibbs(bound);
		cout<<"Starting Learning Iters..."<<endl;
		//cout<<"Done compiling.."<<endl;
		int iters = 0;
		avgcounts.resize(models.size());
		//init weights
                cout<<"Initial Weights"<<endl;
		for(int j=0;j<avgcounts.size();j++) {
			/*clauses[j]->weight = alpha*(datacounts[j]/totalgroundings[j]);
			double r = (double)rand()/RAND_MAX;
			if(datacounts[j]/totalgroundings[j] < 0.1)
				clauses[j]->weight = -1*alpha*r;
			else
				clauses[j]->weight = alpha*r;
                        */
                        cout<<clauses[j]->weight<<" ";
		}
                cout<<endl;
                cout<<"Learning rate="<<alpha<<endl;
		cout<<"Samples per update="<<itersperupdate<<endl;
		time_t start;
		time(&start);
		int i=1;
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
			gibbsflip(pidx,nevids[pidx][nidx],i);
			int thin = 1;
			if(i%thin==0){
			samplecnt++;
			for(int j=0;j<avgcounts.size();j++) {
				avgcounts[j] += (totalgroundings[j]-btrees[j]->currPRVal);
				//clauses[j]->weight += (datacounts[j]-totalgroundings[j]+btrees[j]->currPRVal)*alpha - lambda*clauses[j]->weight;
				//avgcounts[j] = 0;
				//cout<<datacounts[j]<<","<<totalgroundings[j]-btrees[j]->currPRVal<<endl;
				//cout<<clauses[j]->weight<<",";
			}
			}
			if(samplecnt==itersperupdate) {
			        cout<<"Iteration "<<effiters<<endl;
				double maxchange=0;
			//ofstream tmpfs("currentweights.txt");
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
			//tmpfs<<endl;
                        cout<<endl;
			samplecnt = 0;
			//tmpfs.close();
			if(maxchange<0.001) {
				done=true;
			}
                        effiters++;
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
			i++;
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
			ofile<<clauses[cid]->weight<<":"<<toks1[1]<<endl;
			cid++;
		}
		infile.close();
		ofile.close();
		fs1.close();
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
		for(int q=0;q<qrypreds.size();q++) {
			int i = qrypreds[q];
			for(int j1=0;j1<nevids[i].size();j1++) {
				int j = nevids[i][j1];
				probabilities[i][j] = (probabilities[i][j])/(double)(niters-burnin-1);
				if(probabilities[i][j]==1)
					cll += log(0.99)*counts[i][j];
				else if(probabilities[i][j]==0)
					cll += log(0.0001)*counts[i][j];
				else
					cll += log(probabilities[i][j])*counts[i][j];
				cnt += counts[i][j];
			}
		}
		ofstream tmpf(outfile.c_str());
		tmpf<<cll/cnt<<endl;
		tmpf.close();
	}


	void doSampling(string mlnfile,string evidfile,string qryfile,int niters, string outfile) {
		iType = GIBBS;
		init(mlnfile,evidfile,qryfile);
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
		initGibbs(10);

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
		}
			ofstream out(outfile.c_str());
			for(int q=0;q<qrypreds.size();q++) {
				int i = qrypreds[q];
				for(int j=0;j<probabilities[i].size();j++) {
					vector<int> grd;
					getgrounding(i,j,grd);
					if(probabilities[i][j]==-1) {
						probabilities[i][j] = 0;
					} else if(probabilities[i][j]==-2)
					{
						probabilities[i][j] = 1;

					} else {
						probabilities[i][j] = probabilities[i][j]/(double)(niters-burnin-1);
					}
					out<<predicates[i]->symbol<<":";
					for(int k=0;k<grd.size();k++) {
					    out<<grd[k]<<":";
					}
					out<<probabilities[i][j]<<endl;
				}
				
			}
			out.close();
		}
};



#endif /* MLN_H_ */
