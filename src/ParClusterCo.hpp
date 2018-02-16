
#ifndef CLUSTERCO_H_
#define CLUSTERCO_H_

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <vector>
#include <string>
#include <map>
#include <time.h>
#include "utility.hpp"
#include "Set.hpp"

#define LINESIZE 5000
#define ERROR -1

using namespace std;

typedef  unsigned long int  ULI;
typedef  unsigned long int  Vertex;
typedef  unsigned long int  Vsize;


class ClusterCo
{
private:
				
	ULI     N;				// number of nodes
	double  avgcc;          // Average CC
	double  *ccArray;
	Set     *nlist;
	ULI     *name;
	map<ULI, Vertex> vindex;
	int me, numProc;
	
public:
	
	ClusterCo();
	~ClusterCo();
	void initPar(int npro, int self);		
	double NodeCC(Vertex v);
	void AllNodeCC();	
	void WriteCC();
	double AvgCC(); 
	void ReadGraph(const char *fname);
	void catOut(const char *oFileName);
	ULI getN() {return N;}
	double getSumCC() {return avgcc;}
};	

//Constructor
	
ClusterCo::ClusterCo()	
{ 	
	N=0;					
	avgcc = 0;
	nlist=NULL;
	name=NULL;
	ccArray=NULL;
}


// Destructor

ClusterCo::~ClusterCo()
{ 
	if (ccArray) delete []ccArray;
	if(nlist){
		for (ULI i=0; i<vindex.size();i++)
			nlist[i].destroy();
		delete []nlist;
	}
	if (name) delete []name;
	vindex.clear();
	
}


//Initialize parallel computing parameters

void ClusterCo::initPar(int npro, int self)
{
	numProc=npro;
	me=self;

}


double ClusterCo::AvgCC() {
	 
	return avgcc;
	
}


// compute cluster coefficient of node v in graph g

double ClusterCo::NodeCC(Vertex v)
{
	Vertex w = vindex[v];
 	Vsize  degree = nlist[w].Size();
	if (degree <= 1) 
		return 0;
		
	Vertex u;
	Set A(degree);
	int ccount = 0;

	
	
	for (Vsize i=0; i<degree; i++) {
		u = nlist[w][i];
		A.intersect(nlist[w], nlist[vindex[u]]);
		ccount += A.Size();
	}
	
	A.destroy();
	double cc = (double) ccount / (degree * (degree-1));
	return cc;	
}

// AllNodeCC: calls NodeCC to compute CC of all the nodes

void ClusterCo::AllNodeCC()
{
	Vertex v;
	
	int quota, start;
	quota= (int) ceil((double) N/numProc);
	start= me*quota;
	
	if(start>=N)
		return;
		
	ULI ncount=vindex.size();
	for (Vsize i=0; i<ncount; i++)
		nlist[i].sort();
	
	
	ccArray= new double[quota];
	
	for (v=0; v<quota && (start+v)<N; v++) {
		ccArray[v] = NodeCC(name[v]);
		
		avgcc += ccArray[v];
	}
	//avgcc /= N;
	
}


// write data to a file 
void ClusterCo::WriteCC()
{
	int quota, start;
	quota= (int) ceil((double) N/numProc);
	start= me*quota;
	
	if(start>=N)
		return;
	
	char filename[512];
	int NumChar=sprintf(filename,"./tempdir9xcc/ccOut%d.txt",me);
	
	ofstream outfile(filename, ios::out);
	
	for (ULI v=0; v<quota && (v+start)<N; v++) {
		outfile << name[v] << "\t\t"<< ccArray[v]<< endl; 	
	}	
	
	outfile.close();
	
}


// ReadGraph(): Read the input graph and store necessary data into memory

void ClusterCo::ReadGraph(const char *fname) 
{
	ULI     u, v, temp,wt;					// vertices
	Vsize   i, k, tmpdeg, dummy;
	Vertex  vidx;					// vertex index
	
	ifstream ifp(fname);
	
	if(!ifp.is_open()) {
		cout <<"Proc# "<< me<< " Cannot open input graph file: " << fname << endl;
		exit(1);	
	}

	//cout <<"Proc# "<< me<< " Reading graph file: " << fname << " ... "<<endl; 
	//cout.flush();
	
	
	ifp >> N;			// read first line -- the number of vertex
	
	int quota, start;
	quota= (int) ceil((double) N/numProc);
	
	start= me*quota;
	if(start>=N)
		return;
		
	// Skip 'start' number of entries
	
	for (i=0;i<start;i++){
		ifp >> u >> tmpdeg;
		for(k=0; k<tmpdeg; k++)
      		ifp >> v >> wt >> temp;
	
	}
		
    //nlist = new Set[N];
	name = new ULI[quota];
	
	//Fill the map with vertices which need to further explore for inserting neighbors
	vidx = 0;
	
	for (i=0; i<quota && (i+start)<N; i++) {
	
    	ifp >> u >> tmpdeg;						// read nodes and its degree
		name[i]=u;
		if (!vindex.count(u)) {					// if new vertex, map it 
			
			vindex[u] = vidx;
			vidx++;
		}
		
		
		for(k=0; k<tmpdeg; k++) {
      		ifp >> v >> wt >> temp;					// read the adjacent nodes

			if (!vindex.count(v)) {
				
				vindex[v] = vidx;
				vidx++;
			}
						
		}	
		
	}
	
	// Now fill the neighbor list
	ULI ncount=vindex.size();
	nlist = new Set[ncount];
	ifp.seekg(0, ios::beg);
	ULI iter=0;
	ifp >> N;
	
	while(iter<ncount){
		
		ifp >> u >> tmpdeg;	
		if (vindex.count(u)) {		
			iter++;
			nlist[vindex[u]].init(4);
			for(k=0; k<tmpdeg; k++) {
				ifp >> v >> wt >> temp;
				//if(vindex.count(v))
				nlist[vindex[u]].Dinsert(v);
			}	
		}
		
		else{
			for(k=0; k<tmpdeg; k++) {
				ifp >> v >> wt >> temp;
				
			}
				
		}
				
	}
	
	//vindex.clear();	
	ifp.close();
}


//Construct the output file

void ClusterCo::catOut(const char *oFileName)
{
	
	//cout<<"Concatenating intermediate files..."<<endl;
	
	ofstream oFile;
	oFile.open(oFileName);
	if(!oFile.is_open())
	{
		cout<<"Output file creation failed!"<<endl;
		exit(1);
	}
	
	
	char CatCommand[1024];
	strcpy(CatCommand, "cat ./tempdir9xcc/*>>");

	int k, l;
	
	l = strlen(CatCommand);
	
	for (k=0; k<strlen(oFileName); k++)
		CatCommand[l+k] = oFileName[k];
	CatCommand[l+k] = 0;
	
	system(CatCommand);          

	oFile.close();
	cout<<"\nFinal output file created."<<endl;
	

}

 

#endif /* #ifndef CLUSTERCO_H_ */


