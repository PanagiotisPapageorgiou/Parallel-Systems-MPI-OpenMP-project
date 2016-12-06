#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utilities.h"

int argumentCheck(int number_of_processes,int nx,int ny,int Iam,int max_processes){ /*Checks validity of arguments */

	if((number_of_processes <= 2) || (number_of_processes > max_processes) || (((nx * ny) % number_of_processes) != 0) || (number_of_processes % 2 != 0)){

		if(Iam == 0){
			printf("Error: Number of Processes can only be 4,6 or 8!\n");
			printf("NOTE: Check also that you Grid Dimensions are appropriate!\n");
		}
		MPI_Finalize();
		exit(1);

	}

}

int determineRowCols(int number_of_processes,int* nrow,int* mcol,int* ndims){ /* Determine Rows and Cols of array partition this process will use*/

	if(number_of_processes == 2){
		*mcol = 0;
		*nrow = 2;
		*ndims = 1;
	}
	else if(number_of_processes == 4){
		*mcol = 2;
		*nrow = 2;
		*ndims = 2;
	}
	else if(number_of_processes == 6){
		*mcol = 2;
		*nrow = 3;	
		*ndims = 2;
	}
	else if(number_of_processes == 8){
		*mcol = 2;
		*nrow = 4;
		*ndims = 2;
	}

}

int allocateArray(float**** my_array,int dim_size,int Iam){ /* allocate memory for array partition - allocating in Heap in order to allow for large arrays*/

	int i,j;


	*my_array = malloc(sizeof(float**) * 2);
	
	if(*my_array == NULL){
		printf("ID: %d - Error: allocating memory for array partition!\n",Iam);	
		MPI_Finalize();
		exit(1); 
	}

	for(i=0; i<2; i++){
		*((*my_array)+i) = malloc(sizeof(float*)*(dim_size+2));
		if(*((*my_array)+i) == NULL){
			printf("ID: %d - Error: allocating memory for array partition!\n",Iam);	
			MPI_Finalize();
			exit(1); 
		}
	}
	for(i=0; i<2; i++){
		for(j=0; j<dim_size+2; j++){
			*(*((*my_array)+i)+j) = malloc(sizeof(float) * (dim_size+2));
			if(*(*((*my_array)+i)+j) == NULL){
				printf("ID: %d - Error: allocating memory for array partition!\n",Iam);	
				MPI_Finalize();
				exit(1); 					
			}
		}
	}

}

int initialiseArray(float*** my_array,int ndims,int coords2D[2],int num_of_processes,int dim_size,int* xptr,int* yptr){ 

/* Initialise partition (based on initialisation of original MPI_2DHeat.c */

	int i,j,k,offsetx,offsety,ix,iy;
	
	if((num_of_processes == 4) || (num_of_processes == 6) || (num_of_processes == 8)){

		if(coords2D[0] == 0){
			offsetx=0;
			if(coords2D[1] == 0){
				offsety=0;
			}
			else if(coords2D[1] == 1){
				offsety=dim_size;
			}
			else{
				offsety=2*dim_size;
			}
		}
		else if(coords2D[0] == 1){
			offsetx=dim_size;
			if(coords2D[1] == 0){
				offsety=0;
			}
			else if(coords2D[1] == 1){
				offsety=dim_size;
			}
			else{
				offsety = 2*dim_size;
			}
		}
		else if(coords2D[0] == 2){
			offsetx=2*dim_size;
			if(coords2D[1] == 0){
				offsety=0;
			}
			else if(coords2D[1] == 1){
				offsety=dim_size;
			}
			else{
				offsety = 2*dim_size;
			}
		}

		for(k=0; k<2; k++){
			for(i=0; i<dim_size+2; i++)
				for(j=0; j<dim_size+2; j++){
					if(((i >= 1) && (i <= dim_size)) && ((j >= 1) && (j <= dim_size)) && (k == 0)){
						ix = offsetx + i - 1;
						iy = offsety + j - 1;
						my_array[k][i][j] = (float) (ix * (dim_size - ix - 1) * iy * (dim_size - iy - 1));
					}
					else{
						my_array[k][i][j] = 0.0;
					}
				}
		}

		*xptr = offsetx;
		*yptr = offsety;
	}
	else{
		printf("initialiseArray: Not ready for this number of processes!\n");

		for(i=0; i<2; i++){
			for(j=0; j<dim_size+2; j++){
				if(my_array[i][j] != NULL) free(my_array[i][j]);
			}
			if(my_array[i] != NULL) free(my_array[i]);
		}
		if(my_array[i] != NULL) free(my_array);
		
		MPI_Finalize();
	
		exit(1);
	}

}

int locateNeighbours(MPI_Comm comm2D,int* neighbourU,int* neighbourD,int* neighbourW,int* neighbourE){

		int index = 0;
		int displ = 1;

		MPI_Cart_shift(comm2D, index, displ, neighbourU, neighbourD); /* Allows me to catch immediate neighbours up and below me */

		index = 1;
		MPI_Cart_shift(comm2D, index, displ, neighbourW, neighbourE); /* Allows me to catch immediate neighbours left and right of me */

}


int performCalculations(float**** my_array,int dim_size,int neighbourU,int neighbourD,int neighbourE,int neighbourW,MPI_Comm comm2D,int number_of_processes){

/* The logic behind this function is as follows: As a block I own a partition only of the total array. I check to see if I have neighbours up,down,left and right. With each neighbour I trade my current neighbouring side with them and receive their own with non-block communications Isend,Irecv. While I send and receive with the new sides I update my existing partition using only the cells that are not affected by neighbours. When the neighbouring sides arrive I re-update my array to become the fully updated array */

	int i=0,iz=0,Iam;

	MPI_Request send_r[number_of_processes];
	MPI_Request recv_r[number_of_processes];
	MPI_Status status;

	MPI_Comm_rank(comm2D,&Iam);

	for(i=0; i<STEPS; i++){

		//printf("This is process with ID: %d!\n",Iam);

		if(neighbourU >= 0){
			//printf("ID: %d - I have a NORTH neighbour: %d, sending him my side!\n",Iam,neighbourU);
			sendMyUpSide(*((*my_array)+iz),dim_size,neighbourU,comm2D,&(send_r[neighbourU]));
			//printf("ID: %d - Sent my side successfully!\n",Iam);
			receiveUpSide(*((*my_array)+iz),dim_size,neighbourU,comm2D,&(recv_r[neighbourU]));
			//printf("ID: %d - Received his side successfully!\n",Iam);
		}
		if(neighbourD >= 0){
			//printf("ID: %d - I have a SOUTH neighbour: %d, sending him my side!\n",Iam,neighbourD);
			sendMyDownSide(*((*my_array)+iz),dim_size,neighbourD,comm2D,&(send_r[neighbourD]));
			//printf("ID: %d - Sent my side successfully!\n",Iam);
			receiveDownSide(*((*my_array)+iz),dim_size,neighbourD,comm2D,&(recv_r[neighbourD]));
			//printf("ID: %d - Received his side successfully!\n",Iam);
		}
		if(neighbourE >= 0){
			//printf("ID: %d - I have a EAST neighbour: %d, sending him my side!\n",Iam,neighbourE);
			sendMyRightSide(*((*my_array)+iz),dim_size,neighbourE,comm2D,&(send_r[neighbourE]));
			//printf("ID: %d - Sent my side successfully!\n",Iam);
			receiveRightSide(*((*my_array)+iz),dim_size,neighbourE,comm2D,&(recv_r[neighbourE]));
			//printf("ID: %d - Received his side successfully!\n",Iam);
		}
		if(neighbourW >= 0){
			//printf("ID: %d - I have a WEST neighbour: %d, sending him my side!\n",Iam,neighbourW);
			sendMyLeftSide(*((*my_array)+iz),dim_size,neighbourW,comm2D,&(send_r[neighbourW]));
			//printf("ID: %d - Sent my side successfully!\n",Iam);
			receiveLeftSide(*((*my_array)+iz),dim_size,neighbourW,comm2D,&(recv_r[neighbourW]));
			//printf("ID: %d - Received his side successfully!\n",Iam);
		}

		updateOldArray(*((*my_array)+iz),*((*my_array)+1-iz),dim_size); /*Update all but the neighbouring sides until I receive my neighbours' neighbouring sides*/

		/*if(neighbourU >= 0){ 
			MPI_Wait(&(recv_r[neighbourU]),&status); /*If I haven't received the neighbouring sides*/
		/*}
		if(neighbourD >= 0) MPI_Wait(&recv_r[neighbourD],&status);
		if(neighbourE >= 0) MPI_Wait(&recv_r[neighbourE],&status);
		if(neighbourW >= 0) MPI_Wait(&recv_r[neighbourW],&status);*/

		MPI_Barrier(MPI_COMM_WORLD);
	
		//printf("ID: %d - LET'S GO!\n",Iam);
		
		updateWithNewSides(*((*my_array)+iz),*((*my_array)+1-iz),dim_size); /*Now update the neighbouring sides */
		//if(neighbourU > 0) MPI_Wait(&send_r[neighbourU],&status); /*This is necessary to ensure my neighbours are ready for the next STEP*/
		//if(neighbourD > 0) MPI_Wait(&send_r[neighbourD],&status);
		//if(neighbourE > 0) MPI_Wait(&send_r[neighbourE],&status);
		//if(neighbourW > 0) MPI_Wait(&send_r[neighbourW],&status);

		iz = 1-iz;

		MPI_Barrier(MPI_COMM_WORLD);

	}

}

int updateWithNewSides(float** old_array,float** new_array,int dim_size){ /* Update the non-updated border sides using my neighbours' border sides*/

	struct Parms { 
	  float cx;
	  float cy;
	} parms = {0.1, 0.1};

	int j,i;

	for(j=1; j<dim_size+1; j++){
		new_array[1][j] = new_array[1][j] + parms.cx * (old_array[0][j]);
		new_array[dim_size][j] = new_array[dim_size][j] + parms.cx * (old_array[dim_size+1][j]);
	}

	for(i=1; i<dim_size+1; i++){
		new_array[i][1] = new_array[i][1] + parms.cy * (old_array[i][0]);
		new_array[i][dim_size] = new_array[i][dim_size] + parms.cy * (old_array[i][dim_size+1]);
	}

}

int updateOldArray(float** old_array,float** new_array,int dim_size){ /* Update my array using only the cells that are not affected from neighbours */

	struct Parms { 
	  float cx;
	  float cy;
	} parms = {0.1, 0.1};

	int i,j;

	for(i=1; i<dim_size+1; i++){
		for(j=1; j<dim_size+1; j++){
			new_array[i][j] = old_array[i][j] + parms.cx * ((-1)* old_array[i][j] * 2.0) + parms.cy * ((-1) * old_array[i][j] * 2.0);
			if(i != 1){
				//UseUp
				new_array[i][j] = new_array[i][j] + parms.cx * (old_array[i-1][j]); 
			}
			if(i != dim_size){
				//UseDown
				new_array[i][j] = new_array[i][j] + parms.cx * (old_array[i+1][j]); 
			}
			if(j != 1){
				//UseLeft
				new_array[i][j] = new_array[i][j] + parms.cy * (old_array[i][j-1]); 
			}
			if(j != dim_size){
				//UseRight
				new_array[i][j] = new_array[i][j] + parms.cy * (old_array[i][j+1]); 
			}
		}
	}
}

int sendMyUpSide(float** my_array,int dim_size,int neighbourU,MPI_Comm comm2D,MPI_Request* request){

	MPI_Isend(&(my_array[1][1]), dim_size, MPI_FLOAT, neighbourU, DTAG,comm2D, request);	

}

int sendMyDownSide(float** my_array,int dim_size,int neighbourD,MPI_Comm comm2D,MPI_Request* request){

	MPI_Isend(&(my_array[dim_size][1]), dim_size, MPI_FLOAT, neighbourD, UTAG, comm2D, request);	

}

int sendMyRightSide(float** my_array,int dim_size,int neighbourR,MPI_Comm comm2D,MPI_Request* request){

	int i;
	float* temp_array;

	temp_array = malloc(sizeof(float)*dim_size);
	if(temp_array == NULL){
		printf("sendMyRightSide: Error allocating memory!\n");
		MPI_Finalize();
		exit(1);
	}

	for(i=0; i<dim_size; i++)
		temp_array[i] = my_array[i+1][dim_size];

	MPI_Isend(&(temp_array[0]), dim_size, MPI_FLOAT, neighbourR, LTAG, comm2D, request);	

	if(temp_array != NULL) free(temp_array);

}

int sendMyLeftSide(float** my_array,int dim_size,int neighbourL,MPI_Comm comm2D,MPI_Request* request){

	int i;
	float* temp_array;

	temp_array = malloc(sizeof(float)*dim_size);
	if(temp_array == NULL){
		printf("sendMyLeftSide: Error allocating memory!\n");
		MPI_Finalize();
		exit(1);
	}

	for(i=0; i<dim_size; i++)
		temp_array[i] = my_array[i+1][0];

	MPI_Isend(&(temp_array[0]), dim_size, MPI_FLOAT, neighbourL, RTAG, comm2D, request);	

	if(temp_array != NULL) free(temp_array);

}

int receiveUpSide(float** my_array,int dim_size,int neighbourU,MPI_Comm comm2D,MPI_Request* request){

	MPI_Irecv(&(my_array[0]), dim_size, MPI_FLOAT, neighbourU, DTAG, comm2D, request);	

}

int receiveDownSide(float** my_array,int dim_size,int neighbourD,MPI_Comm comm2D,MPI_Request* request){

	MPI_Irecv(&(my_array[dim_size+1]), dim_size, MPI_FLOAT, neighbourD, UTAG, comm2D, request);	

}

int receiveRightSide(float** my_array,int dim_size,int neighbourR,MPI_Comm comm2D,MPI_Request* request){

	int i;
	float* temp_array;

	temp_array = malloc(sizeof(float)*dim_size);
	if(temp_array == NULL){
		printf("receiveRightSide: Error allocating memory!\n");
		MPI_Finalize();
		exit(1);
	}

	MPI_Irecv(&(temp_array[0]), dim_size, MPI_FLOAT, neighbourR, LTAG, comm2D, request);

	for(i=0; i<dim_size; i++)
		my_array[i+1][dim_size+1] = temp_array[i];

	if(temp_array != NULL) free(temp_array);

}

int receiveLeftSide(float** my_array,int dim_size,int neighbourL,MPI_Comm comm2D,MPI_Request* request){

	int i;
	float* temp_array;

	temp_array = malloc(sizeof(float)*dim_size);
	if(temp_array == NULL){
		printf("receiveLeftSide: Error allocating memory!\n");
		MPI_Finalize();
		exit(1);
	}

	MPI_Irecv(&(temp_array[0]), dim_size, MPI_FLOAT, neighbourL, RTAG, comm2D, request);

	for(i=0; i<dim_size; i++)
		my_array[i+1][0] = temp_array[i];

	if(temp_array != NULL) free(temp_array);

}


/* Credit for Print Code: 
http://stackoverflow.com/questions/9777828/writing-a-matrix-into-a-single-txt-file-with-mpi */

/*void printInTurns(float** my_array,int nx,int ny,int dim_size,MPI_Comm comm2D,char* filename,int coords2D[2]){

    int ierr, rank, size,i, j,source,dest1,k,c;
    MPI_Offset offset;
    MPI_File   file;
    MPI_Status status;
    MPI_Datatype num_as_string;
    MPI_Datatype localarray;
    const int nrows=nx;
    const int ncols=ny;
    float **data;
    char *const fmt="%6.1f ";
    char *const endfmt="%6.1f ";
    int startrow, endrow, startcol, endcol, locnrows, locncols;

    const int charspernum=7;
	char* line_buffer;

	FILE* fp,*fp2;

		if(coords2D[0] == 0){
			if(coords2D[1] == 0){
				fp = fopen(filename,"w");
				if(fp == NULL){
					printf("ID: Failed to open file for print!\n",rank);
					return;
				}
				for(i=1; i<dim_size+1; i++){
					for(j=1; j<dim_size+1; j++){
						if(j == dim_size)
							fprintf(fp,"%6.1f\n",my_array[i][j]);
						else
							fprintf(fp,"%6.1f ",my_array[i][j]);
					}
				}
				fclose(fp);
			}
			MPI_Barrier(MPI_COMM_WORLD);

			if(coords2D[1] == 1){
				line_buffer = malloc(sizeof(char) *(charspernum)*((dim_size*coords2D[1]) + 4) * 2);
				if(line_buffer == NULL){
					printf("ID: Failed to allocate buffer for print!\n",rank);
					return;
				}

				fp = fopen(filename,"r");
				fp2 = fopen(filename,"w");
				for(i=1; i<dim_size+1; i++){
					for(j=1; j<dim_size+1; j++){
						k=0;
						do{
							c = fgetc(fp);
							
							if(c == EOF) break;
							k++;
						} while(c != '\n');
					}
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if(coords2D[1] == 2)
				offsety=2*dim_size;
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		else if(coords2D[0] == 1){
			offsetx=dim_size;
			if(coords2D[1] == 0){
				offsety=0;
			}
			else if(coords2D[1] == 1){
				offsety=dim_size;
			}
			else{
				offsety = 2*dim_size;
			}
		}
		else if(coords2D[0] == 2){
			offsetx=2*dim_size;
			if(coords2D[1] == 0){
				offsety=0;
			}
			else if(coords2D[1] == 1){
				offsety=dim_size;
			}
			else{
				offsety = 2*dim_size;
			}
		}

    ierr|= MPI_Comm_size(comm2D, &size);
    ierr|= MPI_Comm_rank(comm2D, &rank);

    locnrows = dim_size;
	locncols = dim_size;

    startrow = coords2D[0] * locnrows;
    endrow = startrow + locnrows - 1;
	startcol = coords2D[1] * locncols;
	endcol = startcol + locncols - 1;

	int index = 0; /* shift along the 1st index (out of 2)  aka ROWS*/
	/*int displ = 1; /* shift by 1 */
	/*MPI_Cart_shift(comm2D, index, displ, &source, &dest1);

	if(source < 0){
		startrow = 0;
	}
	if(dest1 < 0){
		//printf("This is ID: %d, below me there is nothing!\n",rank);
		endrow = nrows - 1;
		locnrows = endrow - startrow + 1;
	}

	index = 1;
	displ = 1;
	MPI_Cart_shift(comm2D, index, displ, &source, &dest1);

	if(source < 0){
		startcol = 0;
	}
	if(dest1 < 0){
		//printf("This is ID: %d, right from me there is nothing!\n",rank);
		endcol = ncols - 1;
		locncols = endcol - startcol + 1;
	}

    /* allocate local data */
    /*data = alloc2d(locnrows, locncols);

    /* fill local data */
    /*for (i=0; i<locnrows; i++) 
        for (j=0; j<locncols; j++)
            data[i][j] = /*my_array[i+1][j+1]*/ /*rank;*/

    /* each number is represented by charspernum chars */
    /*MPI_Type_contiguous(charspernum, MPI_CHAR, &num_as_string); 
    MPI_Type_commit(&num_as_string); 

    /* convert our data into txt */
   /* char *data_as_txt = malloc(locnrows*locncols*charspernum*sizeof(char));
    int count = 0;
    for (i=0; i<locnrows; i++) {
        for (j=0; j<locncols-1; j++) {
            sprintf(&data_as_txt[count*charspernum], fmt, data[i][j]);
            count++;
        }
        sprintf(&data_as_txt[count*charspernum], endfmt, data[i][locncols-1]);
        count++;
    }

    //printf("%d: %s\n", rank, data_as_txt);

    /* create a type describing our piece of the array */
    /*int globalsizes[2] = {nrows, ncols};
    int localsizes [2] = {locnrows, locncols};
    int starts[2]      = {startrow, startcol};
    int order          = MPI_ORDER_C;

	printf("ID: %d - LR-LC:(%d,%d) - SR-SC:(%d,%d) - ER-EC:(%d,%d)\n",rank,locnrows,locncols,startrow,startcol,endrow,endcol);

    MPI_Type_create_subarray(2, globalsizes, localsizes, starts, order, num_as_string, &localarray);
    MPI_Type_commit(&localarray);

    /* open the file, and set the view */
    /*MPI_File_open(comm2D, filename, 
                  MPI_MODE_CREATE|MPI_MODE_WRONLY,
                  MPI_INFO_NULL, &file);

    MPI_File_set_view(file, 0,  MPI_CHAR, localarray, 
                           "native", MPI_INFO_NULL);*/


	/*for(i=0; i<locnrows; i++){
		fseek (file , coords2D[0]* , SEEK_SET );
		MPI_File_write(file, &(localarray[i]), locncols,  &num_as_strings, &status);
	}*/

    /*MPI_File_write_all(file, data_as_txt, locnrows*locncols, num_as_string, &status);
    MPI_File_close(&file);

    MPI_Type_free(&localarray);
    MPI_Type_free(&num_as_string);

    free(data[0]);
    free(data);

}*/

int destroyArray(float**** my_array,int dim_size,int Iam){
	
	int i,j;

	if(*my_array != NULL){
		for(i=0; i<2; i++){
			if(*((*my_array)+i) != NULL){
				for(j=0; j<dim_size+2; j++){
					if( *(*((*my_array)+i) + j) != NULL)
						free(*(*((*my_array)+i) + j));
				}
				free(*((*my_array)+i));
			}
		}
		free(*my_array);
	}

}

float **alloc2d(int n, int m) {
	int i;

   	float *data = malloc(n*m*sizeof(float));
    float **array = malloc(n*sizeof(float *));
    for (i=0; i<n; i++)
        array[i] = &(data[i*m]);
    return array;
}
