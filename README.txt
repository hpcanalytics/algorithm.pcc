
Compile:
make parcc

Run:
Please see the qsub file in this directory

mpiexec -f $PBS_NODEFILE parcc <input graph> <output file name>

input graph should be in galib graph file format.
