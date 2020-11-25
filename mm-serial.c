/* mm-serial.c
 * Griffin Dube
 *
 * Matrix-Matrix multiplication.
 * mm-serial input_file1 input_file2 output_file
 * Reads 2 input matrix files and computes the product of the 
 * matrix in input_file1 with the matrix in input_file2 (remember 
 * in matrix multiply, order matters), and then writes resulting 
 * product into the output file. The columns of input_file1 and 
 * the rows of input_file2 must be equal, otherwise print an 
 * error and stop.  
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#define DIMS 2
#define R 0
#define C 1

void matrix_multiply(double **M1, double **M2, double ***O, int *dout, int *d1, int *d2);
void print_matrix(int r, int c, double **A);
void read_matrix(char *file_name, int *r, int *c, double ***A);
void write_matrix(char *file_name, int r, int c, double ***A);
void free_matrix(int r, double ***A);

int main(int argc, char * argv[]) {
	char *infile1, *infile2, *outfile;
	int d1[DIMS], d2[DIMS], dout[DIMS];
	int i;
	double **M1, **M2, **O;

	if (argc != 4){
		printf("Usage: ./exec <input_file1> <input_file2> <output_file>");
		exit(0);
	}
	else {
		infile1 = argv[1];
		infile2 = argv[2];
		outfile = argv[3];
	}

	//load input matrices
	read_matrix(infile1, &(d1[R]), &(d1[C]), &M1);
	printf("Matrix one dimensions: %d x %d\n",d1[R],d1[C]);
	read_matrix(infile2, &(d2[R]), &(d2[C]), &M2);
	printf("Matrix two dimensions: %d x %d\n",d2[R],d2[C]);

	//assign output dimensions
	dout[R] = d1[R];
	dout[C] = d2[C];

	//allocate output matrix
	O = (double **) malloc(dout[R] * sizeof(double *));
	if(O == NULL) {
		printf("Output matrix malloc failed exiting...\n");
		exit(0);
	}
	for(i = 0; i < dout[R]; i++) {
		O[i] = (double *)malloc(dout[C] * sizeof(double));
		if(O[i] == NULL) {
			printf("Output failed to allocate matrix columns exiting ..\n");
			exit(0);
		}
	}
  
	//perform MMM
	matrix_multiply(M1,M2,&O,dout,d1,d2);

	write_matrix(outfile, dout[R], dout[C], &O);

	free_matrix(d1[R],&M1);
	free_matrix(d2[R],&M2);
	free_matrix(dout[R],&O);
	return 0;
}

void matrix_multiply(double **M1, double **M2, double ***O, int *dout, int *d1, int *d2) {
int i,j,k;

	for (i = 0; i < dout[R]; i++) {
		for (j = 0; j < dout[C]; j++) {
			for (k = 0; k < d1[C]; k++) {
				(*O)[i][j] += M1[i][k]*M2[k][j];	
			}
		}
	}

	return;
}

void print_matrix(int r, int c, double ** A) {
    int i,j;
    printf("Array is a %d x %d matrix\n\n",r,c);
    for(i = 0; i < r; i++) {
        for(j = 0; j < c; j++) {
					printf("%10.3f ",A[i][j]);
        }
        printf("\n");
    }
	printf("\n");

	return;
}

void read_matrix(char *file_name, int *r, int *c, double ***A) {

	FILE *input;
	double **AStorage;
	int i,j;

	input = fopen(file_name,"r");
	if(input == NULL) {
		printf("encountered error opening file exiting...\n");
		exit(0);
	}	

	if(fread(r,sizeof(int),1,input) != 1) {
		printf("error reading matrix row size exiting...\n");
		exit(0);
	}
	if(fread(c,sizeof(int),1,input) != 1) {
		printf("error reading matrix column size exiting...\n");
		exit(0);
	}

	//read in matrix values
	AStorage = (double **) malloc(*r * sizeof(double *));
	if(AStorage == NULL) {
		printf("matrix malloc failed exiting...\n");
		exit(0);
	}
	for(i = 0; i <*r; i++) {
		AStorage[i] = (double *)malloc(*c * sizeof(double));
		if(AStorage[i] == NULL) {
			printf("failed to allocate matrix columns exiting ..\n");
			exit(0);
		}
	}
	for(i = 0; i < *r; i++) {
		for(j = 0; j < *c; j++) {
			if(fread(&AStorage[i][j],sizeof(double),1,input) != 1) {
				printf("error reading matrix value exiting...\n");
				exit(0);
			}
		}
	}

	*A = AStorage;
	//assign to proper matrix
	fclose(input);
	return;
}

void write_matrix(char *file_name, int r, int c, double ***A) {

	FILE *output = NULL;
	int i,j;
	double temp;
	output = fopen(file_name,"w");
	if(output == NULL) {
		printf("encountered error opening file exiting...\n");
		exit(0);
	}	

	if(fwrite(&r,sizeof(int),1,output) != 1) {
		printf("error writing matrix rows exiting...\n");
		exit(0);
	}
	if(fwrite(&c,sizeof(int),1,output) != 1) {
		printf("error writing matrix columns exiting...\n");
		exit(0);
	}
	//write in matrix values
	for(i = 0; i < r; i++) {
		for(j = 0; j < c; j++) {
			temp = (*A)[i][j];
			if(fwrite(&temp,sizeof(double),1,output) != 1) {
				printf("error reading matrix value exiting...\n");
				exit(0);
			}
		}
	}
	//assign to proper matrix
	fclose(output);
	return;
}

void free_matrix(int r,double ***A){
	for(int i =0; i < r; i++) {
			free((*A)[i]);
		}
		free(*A);
		*A = NULL;
}
