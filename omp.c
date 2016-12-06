#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifdef _OPENMP
	#include <omp.h>
#else
	int omp_get_thread_num(void) { return 0; }
	int omp_get_num_threads(void) { return 1; }
#endif

#define NXPROB      20                 /* x dimension of problem grid */
#define NYPROB      20                 /* y dimension of problem grid */
#define STEPS       100                /* number of time steps */
#define MAXWORKER   8                  /* maximum number of worker tasks */

#define MINWORKER   3                  /* minimum number of worker tasks */
#define BEGIN       1                  /* message tag */
#define LTAG        2                  /* message tag */
#define RTAG        3                  /* message tag */
#define NONE        0                  /* indicates no neighbor */
#define DONE        4                  /* message tag */
#define MASTER      0                  /* taskid of first process */


void pack_info(int*, int*, int*, int*, int*, MPI_Datatype*);
void pack_results(int*, int*, float*, MPI_Datatype*);

struct Parms { 
  float cx;
  float cy;
} parms = {0.1, 0.1};

int main (int argc, char *argv[]) {
	void inidat(), prtdat(), update();

	float  u[2][NXPROB][NYPROB], f, temp[NXPROB][NYPROB]; 	/* array for grid */
	int	taskid,                     	/* this task's unique id */
		numworkers,                		/* number of worker processes */
		numtasks,                  		/* number of tasks */
		averow,rows,offset,extra,   	/* for sending rows of data */
		dest, source,              		/* to - from for message send-receive */
		left,right, ble=0,     		    /* neighbor tasks */
		msgtype,                   		/* for message types */
		rc,start,end,              		/* misc */
		i,ix,iy,iz,it,j,k,calculated_flag=0;          /* loop variables */
	MPI_Status status;

	typedef struct {
		int offset;
		int rows;
		int left;
		int right;
		int ble;
	} Box;
	Box a_box, c_box;					/* a box and client_box */
	MPI_Datatype package;


/* First, find out my taskid and how many tasks are running */
   MPI_Init(&argc,&argv);  
   MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
   MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
   numworkers = numtasks-1;
   MPI_Request r[numtasks];

   pack_info(&offset, &rows, &left, &right, &ble, &package);   /*create new types */

   if (taskid == MASTER) { 
	   clock_t tic = clock();
      /************************* master code *******************************/
      /* Check if numworkers is within range - quit if not */

      if ((numworkers > MAXWORKER) || (numworkers < MINWORKER)) {
         printf("ERROR: the number of tasks must be between %d and %d.\n",
                 MINWORKER+1,MAXWORKER+1);
         printf("Quitting...\n");
         MPI_Abort(MPI_COMM_WORLD, rc);
         exit(1);
      }
      printf ("Starting mpi_heat2D with %d worker tasks.\n", numworkers);

      /* Initialize grid */
      printf("Grid size: X= %d  Y= %d  Time steps= %d\n",NXPROB,NYPROB,STEPS);
      printf("Initializing grid and writing initial.dat file...\n");
      inidat(NXPROB, NYPROB, u);
      prtdat(NXPROB, NYPROB, u, "initial.dat");

      /* Distribute work to workers.  Must first figure out how many rows to */
      /* send and what to do with extra rows.  */
      averow = NXPROB/numworkers;
      extra = NXPROB%numworkers;
      offset = 0;
      for (i=1; i<=numworkers; i++) {
         rows = (i <= extra) ? averow+1 : averow; 
         /* Tell each worker who its neighbors are, since they must exchange */
         /* data with each other. */  
         if (i == 1) 
            left = NONE;
         else
            left = i - 1;
         if (i == numworkers)
            right = NONE;
         else
            right = i + 1;
         /*  Now send startup information to each worker  */
         dest = i; 
		 
/*	$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */

		 a_box.offset = offset;
	     a_box.rows = rows;
		 a_box.left = left; 
	     a_box.right = right; 
		MPI_Isend(&a_box, 1, package, dest, BEGIN, MPI_COMM_WORLD, &r[dest]);
		MPI_Wait(&r[dest], &status);	
		MPI_Isend(&u[0][offset][0], NYPROB*rows, MPI_FLOAT, dest, BEGIN, MPI_COMM_WORLD, &r[dest]);
		MPI_Wait(&r[dest], &status);
/*  _____________________________________________________________________________________________________ */

         offset = offset + rows; 
      }
      /* Now wait for results from all worker tasks */
      for (i=1; i<=numworkers; i++) { 
         source = i;
         msgtype = DONE;

/*	$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */

		  MPI_Irecv(&c_box, 1, package, source, msgtype, MPI_COMM_WORLD, &r[source]);
		  MPI_Wait(&r[source], &status);
		  offset = c_box.offset;
		  rows = c_box.rows;
		  MPI_Irecv(&u[0][offset][0], rows*NYPROB, MPI_FLOAT, source, msgtype, MPI_COMM_WORLD, &r[source]);
		  MPI_Wait(&r[source], &status);
/*  _____________________________________________________________________________________________________ */

      }
      /* Write final output, call X graph and finalize MPI */
      printf("Writing final.dat file...\n");
      prtdat(NXPROB, NYPROB, &u[0][0][0], "final_omp.dat");
	  	clock_t toc = clock();
		printf("Elapsed: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
   }   /* End of master code */

   /************************* workers code **********************************/
   if (taskid != MASTER) { 
      /* Initialize everything - including the borders - to zero */
      for (iz=0; iz<2; iz++)
        for (ix=0; ix<NXPROB; ix++) 
            for (iy=0; iy<NYPROB; iy++) 
               u[iz][ix][iy] = 0.0;

      /* Receive my offset, rows, neighbors and grid partition from master */
      source = MASTER;
      msgtype = BEGIN;

/*	$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */

	  MPI_Irecv(&c_box, 1, package, source, msgtype, MPI_COMM_WORLD, &r[source]);
	  MPI_Wait(&r[source], &status);
	  offset = c_box.offset;
	  rows = c_box.rows;
	  left = c_box.left;
	  right = c_box.right;
	  MPI_Irecv(&u[0][offset][0], NYPROB*rows, MPI_FLOAT, source, msgtype, MPI_COMM_WORLD, &r[source]);
	  MPI_Wait(&r[source], &status);
/*  _____________________________________________________________________________________________________ */

      /* Determine border elements.  Need to consider first and last columns. */
      /* Obviously, row 0 can't exchange with row 0-1.  Likewise, the last */
      /* row can't exchange with last+1.  */
      if (offset==0) 
         start=1;
      else 
         start=offset;
      if ((offset+rows)==NXPROB) 
         end=start+rows-2;
      else 
         end = start+rows-1;
      /* Begin doing STEPS iterations.  Must communicate border rows with */
      /* neighbors.  If I have the first or last grid row, then I only need */
      /*  to  communicate with one neighbor  */
      printf("Task %d received work. Beginning time steps...\n",taskid);
      iz = 0;
	  
      for (it = 1; it <= 1; it++) {
		 calculated_flag = 0;
         if (left != NONE) {		
            MPI_Isend(&u[iz][offset][0], NYPROB, MPI_FLOAT, left,
                     RTAG, MPI_COMM_WORLD, &r[left]);
            source = left;
            msgtype = LTAG;
            MPI_Irecv(&u[iz][offset-1][0], NYPROB, MPI_FLOAT, source,
                      msgtype, MPI_COMM_WORLD, &r[source]);
			if(calculated_flag == 0){
				update(start,end,NYPROB,&u[iz][0][0],&temp[0][0],NULL);					
				calculated_flag = 1;
			}
			MPI_Wait(&r[source], &status);
         }
         if (right != NONE) {
            MPI_Isend(&u[iz][offset+rows-1][0], NYPROB, MPI_FLOAT, right,
                      LTAG, MPI_COMM_WORLD, &r[right]);
            source = right;
            msgtype = RTAG;
            MPI_Irecv(&u[iz][offset+rows][0], NYPROB, MPI_FLOAT, source, msgtype,
                      MPI_COMM_WORLD, &r[source]);
			if(calculated_flag == 0){
			  update(start,end,NYPROB,&u[iz][0][0],&temp[0][0],NULL);	
			  calculated_flag = 1;
			}
			MPI_Wait(&r[source], &status);
         }
         /* Now call update to update the value of grid points */
         update(start,end,NYPROB,&u[iz][0][0],&u[1-iz][0][0],&temp[0][0]);
         iz = 1 - iz;
      } 

      /* Finally, send my portion of final results back to master */
/*	$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ */

		 MPI_Isend(&c_box, 1, package, MASTER, DONE, MPI_COMM_WORLD, &r[MASTER]);
		 MPI_Wait(&r[MASTER], &status);
		 MPI_Isend(&u[iz][offset][0], rows*NYPROB, MPI_FLOAT, MASTER, DONE, MPI_COMM_WORLD, &r[MASTER]);
		 MPI_Wait(&r[MASTER], &status);
/*  _____________________________________________________________________________________________________ */
   }
   MPI_Type_free(&package);
   MPI_Finalize();
}


/**************************************************************************
 *  subroutine update
 ****************************************************************************/

void update(int start, int end, int ny, float *u1, float *u2, float *temp) {
    int ix, iy;
	int zone_start = start, zone_end = end, zone_lines;
	
	if (zone_start == 1) zone_start = 0;
	if (zone_end == NXPROB-2) zone_end = NXPROB - 1;
	zone_lines = zone_end - zone_start;
	if(zone_start == 0) zone_start = 1;
	if(zone_end == NXPROB-1) zone_end = NXPROB - 2;

	int zone_columns = NYPROB;
   #pragma omp parallel 
   {
		int tid, total, thread_lines, thread_start, thread_end;

		tid = omp_get_thread_num();				
		total = omp_get_num_threads();

	if(temp == NULL){
		thread_lines = zone_lines/total;					/* dexetai beltiwsh */		
		thread_start = zone_start + tid*thread_lines;
		thread_end = thread_start + thread_lines - 1;			/* final line that this thread changes */
		if (thread_start == 0) thread_start = 1;
		if (thread_end == NXPROB-1) thread_end = NXPROB - 2;

		if(zone_start == 0 && zone_end == 7)
	    for (ix = thread_start; ix <= thread_end ; ix++) 
			for (iy = 1; iy <= NYPROB-2; iy++) {
				 *(u2+ix*ny+iy) = *(u1+ix*ny+iy) +
		          parms.cy * (*(u1+ix*ny+iy+1) +
		         *(u1+ix*ny+iy-1) - 
		          2.0 * *(u1+ix*ny+iy)); }
		#pragma omp barrier
	}
	else{
		int thread_columns;

		thread_columns = zone_columns/total;
		thread_start = tid*thread_columns;
		thread_end = thread_start + thread_columns - 1;
		if (tid == 0) thread_start = 1;
		if (tid == total-1) thread_end = NYPROB-2;
		for (iy=thread_start ; iy<=thread_end ; iy++)
			for (ix=1 ; ix<=NXPROB-2 ; ix++) {
				*(u2+ix*ny+iy) = *(temp+ix*ny+iy)  + 
		          parms.cx * (*(u1+(ix+1)*ny+iy) +
		          *(u1+(ix-1)*ny+iy) - 
		          2.0 * *(u1+ix*ny+iy)); }
	}	
   }
}

/*****************************************************************************
 *  subroutine inidat
 *****************************************************************************/
void inidat(int nx, int ny, float *u) {
int ix, iy;

for (ix = 0; ix <= nx-1; ix++) 
  for (iy = 0; iy <= ny-1; iy++)
     *(u+ix*ny+iy) = (float)(ix * (nx - ix - 1) * iy * (ny - iy - 1));
}

/**************************************************************************
 * subroutine prtdat
 **************************************************************************/
void prtdat(int nx, int ny, float *u1, char *fnam) {
int ix, iy;
FILE *fp;

fp = fopen(fnam, "w");
for (iy = ny-1; iy >= 0; iy--) {
  for (ix = 0; ix <= nx-1; ix++) {
    fprintf(fp, "%6.1f", *(u1+ix*ny+iy));
    if (ix != nx-1) 
      fprintf(fp, " ");
    else
      fprintf(fp, "\n");
    }
  }
fclose(fp);
}

/*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  $ 							OUR CODE
  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/

void pack_info(int *offset, int *rows, int *left, int *right, int *ble, MPI_Datatype *package) {
	//int cells = (*rows)*NYPROB;
	int block_lengths[5] = {1, 1, 1, 1, 1};
	MPI_Datatype types[5] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
	MPI_Aint a_addr, b_addr, c_addr, d_addr, e_addr;
	MPI_Aint displacements[5] = {0};
	
	MPI_Get_address(offset, &a_addr);
	MPI_Get_address(rows, &b_addr);
	MPI_Get_address(left, &c_addr);
	MPI_Get_address(right, &d_addr);
	MPI_Get_address(ble, &e_addr);
	displacements[1] = b_addr - a_addr;
	displacements[2] = c_addr - a_addr;
	displacements[3] = d_addr - a_addr;
	displacements[4] = e_addr - a_addr;
	MPI_Type_create_struct(5, block_lengths, displacements, types, package);
	MPI_Type_commit(package); 
} 
