#include "mpi.h"
//#include "uberlist.h"

#define MAXPROCESSES 8
#define NXPROB      20                 /* x dimension of problem grid */
#define NYPROB      20                 /* y dimension of problem grid */
#define STEPS       100                  /* number of time steps */
#define UTAG		1				   /* Tag for message from upper neighbour */
#define LTAG		2				   /* Tag for message from left neighbour */
#define RTAG		3				   /* Tag for message from right neighbour */
#define DTAG		4				   /* Tag for message from down neighbour */

int argumentCheck(int,int,int,int,int);
int determineRowCols(int,int*,int*,int*);
int initialiseArray(float*** my_array,int ndims,int coords2D[2],int number_of_processes,int dim_size,int*,int*);
int allocateArray(float****,int,int);
void printInTurns(float**,int,int,int,MPI_Comm,char*,int coords2D[2]);
float **alloc2d(int, int);
int updateWithNewSides(float**,float**,int);
int updateOldArray(float**,float**,int);
int sendMyUpSide(float** my_array,int dim_size,int neighbourU,MPI_Comm comm2D,MPI_Request* request);
int sendMyDownSide(float** my_array,int dim_size,int neighbourD,MPI_Comm comm2D,MPI_Request* request);
int sendMyRightSide(float** my_array,int dim_size,int neighbourR,MPI_Comm comm2D,MPI_Request* request);
int sendMyLeftSide(float** my_array,int dim_size,int neighbourL,MPI_Comm comm2D,MPI_Request* request);
int receiveUpSide(float** my_array,int dim_size,int neighbourU,MPI_Comm comm2D,MPI_Request* request);
int receiveDownSide(float** my_array,int dim_size,int neighbourD,MPI_Comm comm2D,MPI_Request* request);
int receiveRightSide(float** my_array,int dim_size,int neighbourR,MPI_Comm comm2D,MPI_Request* request);
int receiveLeftSide(float** my_array,int dim_size,int neighbourL,MPI_Comm comm2D,MPI_Request* request);
int performCalculations(float**** my_array,int dim_size,int neighbourU,int neighbourD,int neighbourE,int neighbourW,MPI_Comm comm2D,int number_of_processes);
int destroyArray(float**** my_array,int dim_size,int Iam);
int locateNeighbours(MPI_Comm,int*,int*,int*,int*);
