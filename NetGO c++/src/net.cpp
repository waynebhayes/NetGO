//============================================================================
// Name        : net.cpp
// Author      : me
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <unordered_set>
#include <assert.h>
#include "net.hpp"

NetGO::NetGO(){
	this->DRACONIAN = true;
}

void NetGO::newCluster(std::string alignFile){
	this->pC.clear();
	this->CA.clear();
	std::ifstream alignFileOpen;
	alignFileOpen.open(alignFile);
	this->get_pC_CA(alignFileOpen);
	alignFileOpen.close();
}

void NetGO::newg2g(std::string g2gFile){
	this->GOp.clear();
	this->pGO.clear();
	this->GOfreq.clear();
	std::ifstream g2gFileOpen;
	g2gFileOpen.open(g2gFile);
	this->get_pGO_GOp(g2gFileOpen);
	g2gFileOpen.close();
}

std::map<std::string, int> NetGO::SetIntersect(std::map<std::string, int> T1, std::map<std::string, int> T2){
	std::map<std::string, int> out;
	if(T1.size()>T2.size()){
		for(auto gpair : T2){
			std::string g = gpair.first;
			if(T1.find(g)!=T1.end()){
				out[g]=1;
			}
		}
	}
	else{
		for(auto gpair : T1){
			std::string g = gpair.first;
			if(T2.find(g)!=T2.end()){
				out[g]=1;
			}
		}
	}
	return out;
}

float NetGO::K_g(std::string g){
	if(this->GOp.find(g)!=this->GOp.end()){
		float out = (float)1/(float)(this->GOp[g].size());
		return out;
	}else{
		return 0;
	}
}

float NetGO::K_gset(std::map<std::string, int> T){
	float out = 0;
	for(auto gpair : T){
		std::string g = gpair.first;
		out = out + this->K_g(g);
	}
	return out;
}

float NetGO::K_p(std::string p){
	if(this->pGO.find(p)!=this->pGO.end()){
		return this->K_gset(this->pGO[p]);
	}else{
		return 0;
	}
}

//K_A2 was not implemented
float NetGO::K_AC(std::map<int, std::vector<std::string>> C){
	float out = 0;
	for(auto clpair : C){
		int cl = clpair.first;
		float K_C = 0.0;
		if(this->VERBOSE){
			std::cout<<"cluster "<<cl<<" num of proteins "<<C[cl].size()<<"\n";
		}
		//std::cout<<"cl: "<<cl<<"\n";
		std::map<std::string, int> M;
		std::map<std::string, int> T;
		float numClusterFields = C[cl].size();
		if(numClusterFields>1){
			std::string u;
			/*
			if(u!="_" && u!="NA"){
				M[u]++;
				for(auto gpair : this->pGO[u]){
					std::string g = gpair.first;
					T[g]++;
				}
			}
			*/
			for(int i = 0; i<numClusterFields; i++){
				u=C[cl][i];
				assert(u!="_"&&u!="NA"&&u!="-");
				//std::cout << "\t int in cl: " << u << "\n";
                if(i==0){
                    u=C[cl][0];
                    M[u]+=1;
                    for(auto gpair : this->pGO[u]){
                    	std::string g = gpair.first;
                        T[g]+=1;
                    }
                }
                else{
					M[u]++;
					if(this->DRACONIAN){
						std::unordered_set<std::string> Tremove;
						for(auto gpair : T){
							std::string g = gpair.first;
							if(this->pGO[u].find(g)==this->pGO[u].end()){
								Tremove.insert(g);
							}
						}
						for(std::string g : Tremove){
							T.erase(g);
						}
						//std::cout << "\t\t\t" << T.size() << "\n";
					}else{
						for(auto gpair : T){
							std::string g = gpair.first;
							T[g]++;
						}
					}
				}
        		if(this->VERBOSE){
        			std::cout<<"\t"<<u<<" ("<<T.size()<<" GOs { ";
        			for(auto gpair: T){
        				std::string g = gpair.first;
        				std::cout<<g.substr(1,g.length())<<"("<<this->GOfreq[g]<<") ";
        			}
        			std::cout << "}, K("<<u<<")="<<this->K_gset(T) <<")\n";
        		}
			}
			//std::cout<<"\tchecks:"<<T.size()<<" "<<M.size()<<" "<<M[u]<<"\n";
			if(T.size()>0 and (M.size()>1 or (M.size()==1 and M[u]>1))){
				if(this->DRACONIAN){
					K_C=K_C+this->K_gset(T);
				}else{
					for(auto gpair : T){
						std::string g = gpair.first;
						if(T[g]>1){
							K_C=K_C+T[g]*this->K_g(g)/numClusterFields;
						}
					}
				}
			}
		}
		if(this->VERBOSE){
			std::cout<<"\tClusterCommonGOs {";
			for(auto gpair: T){
				std::string g = gpair.first;
				std::cout<<g.substr(1,g.length());
			}
			std::cout << "}, K_C=" << K_C << "\n";
		}
		out = out + K_C;
	}
	return out;
}

void NetGO::get_pC_CA(std::ifstream& alignFile){
	int lineCount = 0;
	std::string cluster;
	std::string protein;
	while(getline(alignFile, cluster)){
	    size_t delimPos1 = 0;
	    size_t delimPos2 = cluster.find("\t");
	    std::vector<std::string> newvector;
	    this->CA[lineCount]=newvector;
	    while(delimPos2!=std::string::npos){
	        protein = cluster.substr(delimPos1, delimPos2-delimPos1);
	        if(protein!="_"&&protein!="NA"&&protein!="-"){
				this->encountered_proteins.insert(protein);
				delimPos1 = delimPos2+1;
				delimPos2 = cluster.find("\t", delimPos1);
				//set pC
				if(this->pC.find(protein)==pC.end()){
					std::map<int,int> newDict;
					pC[protein]=newDict;
				}
				if(this->pC[protein].find(lineCount)==pC[protein].end()){
					pC[protein][lineCount]=1;
				}else{
					pC[protein][lineCount]++;
				}
				//set CA
				CA[lineCount].push_back(protein);
	        }
	    }
	    //get last elemt for CA
	    size_t end = cluster.length()-1;
	    protein = cluster.substr(delimPos1, end);
	    this->encountered_proteins.insert(protein);
	    CA[lineCount].push_back(protein);
	    //get last element for pC
        if(this->pC.find(protein)==pC.end()){
        	std::map<int,int> newDict;
        	pC[protein]=newDict;
        }
		if(this->pC[protein].find(lineCount)==pC[protein].end()){
			pC[protein][lineCount]=1;
		}else{
			pC[protein][lineCount]++;
		}
		//next
		lineCount++;
	}
}

void NetGO::get_pGO_GOp(std::ifstream& g2gFile){
	int lineCount = 0;
	std::string cluster;
	std::string protein;
	while(getline(g2gFile, cluster)){
	    size_t delimPos1 = 0;
	    size_t delimPos2 = cluster.find("\t");
	    size_t delimPos3 = cluster.find("\t", delimPos2+1);
	    size_t delimPos4 = cluster.find("\t", delimPos3+1);
	    std::string protein;
	    std::string GOterm;
	    protein = cluster.substr(delimPos2+1, delimPos3-delimPos2-1);
	    GOterm = cluster.substr(delimPos3, delimPos4-delimPos3);
	    //fill GOfreq
	    if(this->GOfreq.find(GOterm)!=this->GOfreq.end()){
	    	this->GOfreq[GOterm]++;
	    }else{
	    	this->GOfreq[GOterm]=1;
	    }
	    if(this->encountered_proteins.find(protein)!= encountered_proteins.end()){
	    	//fill pGO
	    	if(this->pGO.find(protein)==this->pGO.end()){
	    		std::map<std::string, int> newmap;
	    		this->pGO[protein]=newmap;
	    	}
	    	this->pGO[protein][GOterm]=1;
	    	//fill GOp
	    	if(this->GOp.find(GOterm)==this->GOp.end()){
	    		std::map<std::string, int> newmap;
	    		this->GOp[GOterm] = newmap;
	    	}
	    	if(this->GOp[GOterm].find(protein)==this->GOp[GOterm].end()){
	    		this->GOp[GOterm][protein]=1;
	    	}else{
	    		this->GOp[GOterm][protein]++;
	    	}
	    }
	}
}

std::vector<std::string> CA;
std::map<std::string, std::map<int, int>> pC;
std::map<std::string, std::map<std::string, int>> pGO;
std::map<std::string, std::map<std::string, int>> GOp;

int main(int argc, char** argv) {
//int main() {
	NetGO s;
	std::string nextword;
	std::string USAGE="USAGE: $0 [-verbose] [-L] gene2goFile alignFile[s]\n-L: 'Lenient'. The default behavior is what we call 'Dracanion' in the paper, which\ninsists that a GO term must annotate every protein in a cluster for it to count.\nThe Lenient option gives a GO term a weight per-cluster that is scaled by the\nnumber of proteins it annotates (so long as it's more than 1).\nalignFile: each line consists of a cluster of any number of proteins; proteins can\nappear in more than one cluster. Note it must use the same protein naming convention\nas your gene2go file.\ngene2goFile: a standard-format gene2go file downloaded from the GO consortium's\nwebsite. For now we only use columns 2 (protein name) and 3 (GO term). We ignore\nthe species, which is a simplification since about 20% of proteins appear in two\nspecies, but only 2% appear in more than 2.\n";
	std::string clusterFile;
	std::string g2gFile;
	std::vector<std::string> clusterfilenames;
	int g2gindex = 1;
	int clusterindex = 2;
	if((std::string) argv[1]=="-L"){
		s.DRACONIAN=false;
		g2gindex=2;
		clusterindex=3;
		if((std::string) argv[2]=="-verbose"){
			s.VERBOSE=true;
			g2gindex=3;
			clusterindex=4;
		}
	}
	if((std::string) argv[1]=="-verbose"){
			s.VERBOSE=true;
			g2gindex=2;
			clusterindex=3;
	}
	//add support for rregex name
	std::string g2gfilename = argv[g2gindex];
	for(int i = clusterindex; i<argc; i++){
		clusterfilenames.push_back(argv[i]);
	}
	std::cout << g2gfilename << "\n";
	std::cout << clusterfilenames[0] << "\n";
	for(std::string clusterfilename : clusterfilenames){
		//s.newCluster("/home/me/myWS/alex/netGO/alignFiles/1-to-1-homologs-RNorvegicus-AThaliana.tsv");
		//s.newg2g("/home/me/myWS/alex/netGO/gene2go");
		s.newCluster(clusterfilename);
		s.newg2g(g2gfilename);
		float sumGOp;
		int count;
		for(auto ppair : s.pC){
			std::string p = ppair.first;
			if(s.pGO.find(p)==s.pGO.end()){
				std::map<std::string, int> newvector;
				s.pGO[p] = newvector;
				count++;
			}
			for(auto gpair : s.pGO[p]){
				std::string g = gpair.first;
				sumGOp=sumGOp+float(s.pC[p].size())*float(s.K_g(g));
			}
		}
		float know = s.K_AC(s.CA);
		std::cout << "numClus " << s.CA.size();
		std::cout << " numP " << s.pGO.size();
		std::cout << " numGO " << sumGOp;
		std::cout << " GOcorpus " << s.GOfreq.size();
		std::cout << " K(A) " << know;
		std::cout << " score " << float(know)/float(sumGOp) << "\n";
		//std::cout<<s.CA[3][0]<<"\n";
		//std::cout<<s.pC["361107"][3];
	}
}
