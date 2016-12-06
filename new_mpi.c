#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "utilities.h"

void main(int argc, char *argv[])
{

	int nrow, mcol, i, j, lastrow, number_of_processes, MASTER=0, neighbourU,neighbourD,neighbourE,neighbourW,offsetx,offsety;
	int Iam, id2D, colID, ndims,dim_size=0;
	int coords1D[1], coords2D[2], dims[2], aij[1], alocal[3];
	int belongs[2], periods[2], reorder;

	float ***my_array; /* The array partition the process will use */
	
	clock_t tic,toc;

	MPI_Comm comm2D, commcol, commrow;

/*=====================================MAIN CODE=========================================================================*/
	
	MPI_Init(&argc, &argv); /* Starts MPI processes ... */
	MPI_Comm_rank(MPI_COMM_WORLD, &Iam); /* MPI process id */
	MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);   /* get number of processes */
	
	argumentCheck(number_of_processes,NXPROB,NYPROB,Iam,MAXPROCESSES); /* Check if array and number of processes are compatible */

	if(Iam == 0) /*Start clock!*/
		tic = clock();

	determineRowCols(number_of_processes,&nrow,&mcol,&ndims); /* Determine how to distribute the processes */

	/* create cartesian topology for processes */
	if(ndims == 2){

		if(Iam == 0){
			printf ("Starting mpi_heat2D with %d processes.\n", number_of_processes);

			/* Initialize grid */
			printf("Grid size: X= %d  Y= %d  Time steps= %d\n",NXPROB,NYPROB,STEPS);
		}

		periods[0] = 0; periods[1] = 0; reorder = 1;
		dims[0] = nrow; /* number of rows */
		dims[1] = mcol; /* number of columns */

		//if(Iam == 0){
			//printf("nrow: %d - ncol: %d\n",nrow,mcol);
		//}

		MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &comm2D); /* Create a 2D topology */
																				/* Both rows and columns are not periodic */
																				/* MPI is allowed to reorder the process ranks for efficiency */
																				/* To create the topology, the old Communicator: MPI_COMM_WORLD was used */
	
		MPI_Barrier(MPI_COMM_WORLD);

		MPI_Comm_rank(comm2D, &id2D); // Get ID of this process in the new communicator
		MPI_Cart_coords(comm2D, id2D, ndims, coords2D); /* In the new communicator, use my ID and number of dimensions and retrieve my coordinates */

		//printf("OriginalID: %d - New ID: %d - 2Dcoord: (%d,%d)\n",Iam, id2D, coords2D[0], coords2D[1]);

		if(Iam == 0){
			printf("Created Topology...\n\n");
		}

		dim_size = (int) sqrt((double) ((NXPROB * NYPROB) / number_of_processes)); /* How many squares per dimension belong to me */
		allocateArray(&my_array,dim_size,Iam); /* Allocate a square big enough for my partition and 1 neighbour side on each side of my square */
		initialiseArray(my_array,ndims,coords2D,number_of_processes,dim_size,&offsetx,&offsety); /* Initialiase my area */

		//printInTurns(*my_array,NXPROB,NYPROB,dim_size,comm2D,"initial_mpi.dat",coords2D); /* Print my area to output file */


		/* Catch north,south,west,east neighbours ID if they exist */

		locateNeighbours(comm2D,&neighbourU,&neighbourD,&neighbourW,&neighbourE);

		MPI_Barrier(MPI_COMM_WORLD);

		if(Iam == 0){
			printf("Everyone located their neighbours!\n");
		}

/*========================================END OF INITIALISATIONS=======================================*/

		if(Iam == 0){
			printf("Initialised array! Let's get to work!\n\n");
		}

		performCalculations(&my_array,dim_size,neighbourU,neighbourD,neighbourE,neighbourW,comm2D,number_of_processes); /*The important part */

		//printInTurns

		MPI_Barrier(MPI_COMM_WORLD);
	}

	destroyArray(&my_array,dim_size,Iam); /*Free resources */

	if(Iam == 0){
      	toc = clock(); 
	  	printf("\nElapsed: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
		printf("\nFinished everything! Bye-bye!\n");
	}

	MPI_Finalize();

}











