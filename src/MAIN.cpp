#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;
#include <stdio.h>
#include <time.h>
#include "mln.h"

void rec(vector<vector<int> > tmpx) {
	vector<int> seq(tmpx.size());
	int nxtchange = tmpx.size()-1;
	int totalsize = 81;
	vector<int> inds(tmpx.size());
	for(int i=0;i<totalsize;i++) {
		for(int j=0;j<tmpx.size();j++) {
			cout<<tmpx[j][inds[j]]<<" ";
		}
		cout<<endl;
		inds[nxtchange]++;
		if(inds[nxtchange]==tmpx[nxtchange].size()) {
			inds[nxtchange] = 0;
			for(int j=nxtchange-1;j>=0;j--) {
				inds[j]++;
				if(inds[j] > tmpx[j].size()-1) {
					inds[j] = 0;
				} else {
					break;
				}
			}
		}

	}
}

int main(int argc, char* argv[])
{
	//CD learning
	if(argc < 9) {
		cout<<"Usage::./cdlearn mlnfile evidfile queryfile totaliters ibound LEARN outputmlnfile constraintsfile"<<endl;
		//cout<<"Usage::./cdlearn mlnfile evidfile queryfile totaliters ibound TEST resultfilename"<<endl;
		cout<<"Usage::./cdlearn mlnfile evidfile queryfile totaliters ibound MAR resultfilename constraints"<<endl;
		//cout<<"queryfile:Required for sampling, for MAP -1 otherwise"<<endl;
		//cout<<"gprob:random step probability for MAP, -1 otherwise"<<endl;
		//cout<<"TYPE:GIBBS/WS"<<endl;
		return 0;
	}
	//vector<vector<int> > tmpx(4);
	//tmpx[0].push_back(0);tmpx[0].push_back(1);tmpx[0].push_back(2);
	//tmpx[1].push_back(0);tmpx[1].push_back(1);tmpx[1].push_back(2);
	//tmpx[2].push_back(0);tmpx[2].push_back(1);tmpx[2].push_back(2);
	//tmpx[3].push_back(0);tmpx[3].push_back(1);tmpx[3].push_back(2);
	//rec(tmpx);
	srand(time(NULL));
	string mlnfile(argv[1]);
	string evidfile(argv[2]);
	string qryfile(argv[3]);
	stringstream st(argv[4]);
	int itrs;
	st >> itrs;
	stringstream st1(argv[5]);
	int bound;
	st1 >> bound;
	string option(argv[6]);

	MLN mln;
	if(option.compare("LEARN")==0) {
		//int updatefreq = 1;
		string outmln(argv[7]);
		string csfile(argv[8]);
		//stringstream st2(argv[9]);
		//int updatefreq;
		//st2 >> updatefreq;
		//stringstream st3(argv[10]);
		//double lrate;
		//st3 >> lrate;
		mln.sample(mlnfile,evidfile,qryfile,itrs,bound,outmln,csfile);
	}
	else if(option.compare("TEST")==0){
		string outresfile(argv[7]);
		mln.testing(mlnfile,evidfile,qryfile,itrs,bound,outresfile);
	}
	else if(option.compare("MAR")==0){
		//stringstream st1(argv[7]);
		//int aval;
		//st1 >> aval;
		string outresfile(argv[7]);
		string csfile(argv[8]);
		string logfilename(argv[9]);
		mln.doSampling(mlnfile,evidfile,qryfile,itrs,bound,outresfile,csfile,logfilename);
	}

	//string outresfile(argv[7]);
	//mln.doSampling(mlnfile,evidfile,qryfile,itrs,outresfile);	
}
