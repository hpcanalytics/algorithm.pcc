 
#include "ParClusterCo.hpp"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <mpi.h>


using namespace std;

/****ComputeCC: invoke routine for computing and printing CC***********/

int ComputeCC(const char *gphfile_name, const char *ofile_name, int numProc, int me)
{
	
	ClusterCo cc;
	cc.initPar(numProc, me);
	
	if (me==0)
		system("mkdir tempdir9xcc");
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(me==0) cout<<"Reading graph file..."<<endl;
	cc.ReadGraph(gphfile_name);
			
	if (me==0) cout<< "Computing exact ClusterCo ... "<<endl;
	cc.AllNodeCC();
	cc.WriteCC();
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	double ReducedSum=0, sum=cc.getSumCC();
	MPI_Reduce(&sum, &ReducedSum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		
	if (me==0){
		
		cc.catOut(ofile_name);
		
		cout<< "Avg. CC (global CC): "<<  ReducedSum/ cc.getN()<<endl;
		
		system("rm -r tempdir9xcc");
		
	}
	
	return 0;
}


int main(int argc, char **argv)
{
	
	int numProc, me;
	unsigned int numBlocks=0;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numProc);
	MPI_Comm_rank(MPI_COMM_WORLD, &me);

	
	if (argc!=3 ) {
		cout << "Usage: "<<argv[0]<<" <graph_file_name> <output-file name>" << endl;
		exit (1);
	}
	
	ComputeCC(argv[1], argv[2], numProc, me);
	
	
	MPI_Finalize();
	
	return 0;
}
