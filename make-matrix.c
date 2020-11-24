/* make-matrix.c
 * Griffin Dube
 *
 * Used to create matrices and store them to a file.
 * This program writes to output_file a binary file with 
 * an n by n (square) matrix representing one of
 * the matrices to be multiplied, where the value of each 
 * element of the matrix is randomly generated based on a
 * uniform random variable  
 * 
 * -r <rows> number of rows in the matrix
 * -c <cols> number of columns in the matrix
 * -l <low> lower bound on the values in the matrix
 * -u <up> upper bound of the values in the matrix
 * -o <fname> output file name
 * Data is 64-bit double precision floating point. For each 
 * value generate a random integer between ‘l’ and ‘u’, convert 
 * to double, and then divide by 1000.0
 *
 * double v = ((double) l + ((double)random() / (u – l))) / 1000.0;
 * 
 * The first two 32-bit integer words of the file should be 
 * ‘r’ and ‘c’ the number of rows and columns of matrix, and the 
 * data should be stored in row major order.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

void write_matrix(char *file_name, int r, int c, double ***A);

int main(int argc, char * argv[]) {

    int opt;
    int i,j,rows,cols,low,up;
		double val;
    char * outputname;
    double **Mstorage;
//    srand(time(NULL));

    if(argc != 11) {
        printf("Usage: ./exec -r <rows> -c <cols> -l <low> -u <up> -o <fname>\n");
        exit(0);
    }

    while((opt = getopt(argc,argv,"r:c:l:u:o:")) != -1) {
        switch(opt) {
            case 'r':
                rows = atoi(optarg);
                break;
            case 'c':
                cols = atoi(optarg);
                break;
            case 'l':
                low = atoi(optarg);
                break;
            case 'u':
                up = atoi(optarg);
                break;
            case 'o':
                outputname = optarg;
                break;
            default:
                printf("Usage: ./exec -r <rows> -c <cols> -l <low> -u <up> -o <fname>\n");
                exit(0);
        }
    }

    Mstorage = (double **) malloc(rows * sizeof(double *));
    if(Mstorage == NULL) {
	printf("failed to allocate matrix exiting ..\n");
	exit(0);
    }
    for(i = 0; i < rows; i++) {
	    Mstorage[i] = (double *) malloc(cols * sizeof(double));
	    if(Mstorage[i] == NULL) {
		    printf("failed to allocate matrix exiting ..\n");
		    exit(0);
	    }
    }

    for(i = 0; i < rows; i++) {
        for(j = 0; j < cols; j++) {
						val = ((double)low + ((double)rand() / (up - low))) / 1000.0;
            Mstorage[i][j] = val;  
        }
    }
    write_matrix(outputname,rows,cols,&Mstorage);	
    return 0;
}

void write_matrix(char *file_name, int r, int c, double ***A) {

	FILE *output;
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
	for(i =0; i < r; i++) {
		free((*A)[i]);
	}
	free(*A);
	*A = NULL;
	return;
}

