#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define NXPROB      1000                 /* x dimension of problem grid */
#define NYPROB      1000                /* y dimension of problem grid */
#define STEPS       100                /* number of time steps */

void prtdat(int nx, int ny, float *u1, char *fnam);
void update(int start, int end, int ny, float **u1, float **u2);
void inidat(int nx, int ny, float*** u);

struct Parms { 
  float cx;
  float cy;
} parms = {0.1, 0.1};


int main (int argc, char *argv[]) {

	float*** u;        /* array for grid */
	int	taskid,                     /* this task's unique id */
	numworkers,                 /* number of worker processes */
	numtasks,                   /* number of tasks */
	averow,rows,offset,extra,   /* for sending rows of data */
	dest, source,               /* to - from for message send-receive */
	left,right,       			/* neighbor tasks */
	msgtype,                    /* for message types */
	rc,start,end,               /* misc */
	i,j,ix,iy,iz,it;              /* loop variables */


      printf ("Starting mpi_heat2D (serial).\n");

	  clock_t tic = clock();
      /* Initialize grid */

	  u = malloc(sizeof(float**) * 2);
	  if(u == NULL) printf("Malloc failed!\n");
	  for(i=0; i<2; i++){
		u[i] = malloc(sizeof(float*) * NXPROB);
		if(u[i] == NULL) printf("Malloc failed2!\n");
	  }

	  for(i=0; i<2; i++)
		for(j=0; j<NXPROB; j++){
			u[i][j] = malloc(sizeof(float) *NYPROB);
			if(u[i][j] == NULL) printf("Malloc failed3!\n");
		}

      printf("Grid size: X= %d  Y= %d  Time steps= %d\n",NXPROB,NYPROB,STEPS);
      printf("Initializing grid and writing initial.dat file...\n");
      inidat(NXPROB, NYPROB, u);
 		prtdat(NXPROB, NYPROB, &u[0][0][0], "initial_serial.dat");

      /* Begin doing STEPS iterations.  Must communicate border rows with */
      /* neighbors.  If I have the first or last grid row, then I only need */
      /*  to  communicate with one neighbor  */

      iz = 0;
      for (it = 1; it <= STEPS; it++) {
         /* Now call update to update the value of grid points */
        update(start,end,NYPROB,u[iz],u[1-iz]);
         iz = 1 - iz;
      }

	if(u != NULL){
		for(i=0; i<2; i++){
			if(*((u)+i) != NULL){
				for(j=0; j<NXPROB; j++){
					if( *(*((u)+i) + j) != NULL)
						free(*(*((u)+i) + j));
				}
				free(*((u)+i));
			}
		}
		free(u);
	}

	   /* Write final output, call X graph and finalize MPI */
      printf("Writing final.dat file and generating graph...\n");
      prtdat(NXPROB, NYPROB, &u[iz][0][0], "final_serial.dat");
      printf("Click on MORE button to view initial/final states.\n");
      printf("Click on EXIT button to quit program.\n");
      clock_t toc = clock();
	  printf("Elapsed: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
}


/**************************************************************************
 *  subroutine update
 ****************************************************************************/
void update(int start, int end, int ny, float **u1, float **u2)
{
   int ix, iy;
   for (ix = start; ix <= end; ix++) 
      for (iy = 1; iy <= ny-2; iy++) 
         u2[ix][iy] = u1[ix][iy]  + 
                          parms.cx * (u1[ix+1][iy] +
                          u1[ix-1][iy] - 
                          2.0 * u1[ix][iy]) +
                          parms.cy * (u1[ix][iy+1] +
                         u1[ix][iy-1] - 
                          2.0 * u1[ix][iy]);
}

/*****************************************************************************
 *  subroutine inidat
 *****************************************************************************/
void inidat(int nx, int ny, float*** u) {
int ix, iy;

for (ix = 0; ix <= nx-1; ix++) 
  for (iy = 0; iy <= ny-1; iy++)
     *(*((*(u))+ix) + iy) = (float)(ix * (nx - ix - 1) * iy * (ny - iy - 1));
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
