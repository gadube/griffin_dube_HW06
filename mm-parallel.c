/* mm-parallel.c
 * Griffin Dube
 *
 * MPI Implementation of Cannons algorithm for MMM using
 * a 2-D cartesian topology.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "matrix_checkerboard_io.h"
#include "MyMPI.h"

#define NDIMS 2

#define LEFT  0
#define RIGHT 1
#define UP    2
#define DOWN  3

#define R     0
#define C     1

void get_args(int argc, char *argv[], char **infile1, char **infile2, char **outfile);
void create_grid_comm(MPI_Comm *comm);
void allocate_matrix(void ***subs, void **storage, MPI_Datatype dtype, int rows, int cols, MPI_Comm grid_comm);
void matrix_multiplication(double **subsA, double *storageA, double **subsB, double *storageB, double **subsC, \
		                       double *storageC, int *dA, int *dB, int *dC, MPI_Comm grid_comm);
void mmm(int lr, int lc, int inner, double **subsA, double **subsB, double **subsC);
void free_array(double ***subs, double **storage);

int main(int argc, char *argv[])
{
	int rank, size, dA[NDIMS], dB[NDIMS], dC[NDIMS];
	double **subsA, **subsB, **subsC, *storageA, *storageB, *storageC;
	char *infile1, *infile2, *outfile;
	MPI_Comm GRID_COMM;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	/* read command line arguments */
	get_args(argc, argv, &infile1, &infile2, &outfile);

	/* create cartesian topology */
	create_grid_comm(&GRID_COMM);
	MPI_Comm_rank(GRID_COMM, &rank);

	/* read input matrices A and B */
	read_checkerboard_matrix(infile1, (void ***)&subsA, (void **)&storageA, MPI_DOUBLE, &(dA[R]), &(dA[C]), GRID_COMM);
	read_checkerboard_matrix(infile2, (void ***)&subsB, (void **)&storageB, MPI_DOUBLE, &(dB[R]), &(dB[C]), GRID_COMM);
	debug("%d: dA = %d x %d; dB = %d x %d;\n",rank, dA[R], dA[C], dB[R], dB[C]);

	/* allocate output matrix C */
	debug("%d: Allocating Output Submatrix...\n", rank);
	dC[R] = dA[R];
	dC[C] = dB[C];
	allocate_matrix((void ***)&subsC, (void **)&storageC, MPI_DOUBLE, dC[R], dC[C], GRID_COMM);
	debug("%d: Completed Allocating Output Submatrix...\n", rank);

	/* perform matrix multiplication */
	double time = MPI_Wtime();
	matrix_multiplication(subsA, storageA, subsB, storageB, subsC, storageC, dA, dB, dC, GRID_COMM);
	time = MPI_Wtime() - time;

	/* find timing values */
	double maxTime;
	MPI_Reduce(&time, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

	if (rank == 0) printf("Num. Processes: %d\nMatrix Size: %dx%d\nCompute Time: %f\n", size, dC[R], dC[C], maxTime);
	
	/* write C output to file */
	write_checkerboard_graph(outfile, (void ***)&subsC, (void **)&storageC, MPI_DOUBLE, dC[R], dC[C], GRID_COMM);
	
	/* clean up memory */
	free_array(&subsA, &storageA);
	free_array(&subsB, &storageB);
	MPI_Comm_free(&GRID_COMM);
	MPI_Finalize();
	return 0;
}

/* frees 2d allocated array */
void free_array(double ***subs, double **storage)
{
	free(*storage);
	*storage = NULL;
	free(*subs);
	*subs = NULL;

	return;
}

/* matrix multiplication kernel */
void mmm(int lr, int lc, int inner, double **subsA, double **subsB, double **subsC)
{
	int i, j, k;
	for (i = 0; i < lr; i++)
	{
		for (j = 0; j < lc; j++)
		{
			for (k = 0; k < inner; k++)
			{
				subsC[i][j] += subsA[i][k] * subsB[k][j];
			}
		}
	}

	return;
}

/* handles communication and performs matrix multiplication */
void matrix_multiplication(double **subsA, double *storageA, double **subsB, double *storageB, double **subsC, \
		                       double *storageC, int *dA, int *dB, int *dC, MPI_Comm grid_comm)
{
	int grid_id, local_rows, local_cols, src, dest, lrA, lrB, lcA, lcB;
	int grid_size[NDIMS], grid_period[NDIMS], grid_coord[NDIMS];
	int nbr[4];
	MPI_Status status;

	MPI_Comm_rank(grid_comm, &grid_id);
	MPI_Cart_get(grid_comm, NDIMS, grid_size, grid_period, grid_coord);

	/* find local row and column sizes for each proc */
	lrA = BLOCK_SIZE(grid_coord[R], grid_size[R], dA[R]);
	lcA = BLOCK_SIZE(grid_coord[C], grid_size[C], dA[C]);
	lrB = BLOCK_SIZE(grid_coord[R], grid_size[R], dB[R]);
	lcB = BLOCK_SIZE(grid_coord[C], grid_size[C], dB[C]);
	local_rows = BLOCK_SIZE(grid_coord[R], grid_size[R], dC[R]);
	local_cols = BLOCK_SIZE(grid_coord[C], grid_size[C], dC[C]);
	debug("%d: local_rows = %d; local_cols = %d;\n",grid_id, local_rows, local_cols);

	/* find neigboring procs */
	MPI_Cart_shift(grid_comm, C, 1, &(nbr[LEFT]), &(nbr[RIGHT]));
	MPI_Cart_shift(grid_comm, R, 1, &(nbr[UP]), &(nbr[DOWN]));
	debug("%d: Neighbors: left = %d right = %d; up = %d; down = %d\n", grid_id,nbr[LEFT],nbr[RIGHT],nbr[UP],nbr[DOWN]);

	/*perfom initial shift */
	MPI_Cart_shift(grid_comm, C, -grid_coord[R], &src, &dest);
	MPI_Sendrecv_replace(storageA, lrA * lcA, MPI_DOUBLE, dest, 1, src, 1, grid_comm, &status);
	debug("%d: Shifting from Proc %d to Proc %d for A\n",grid_id,src,dest);
	MPI_Cart_shift(grid_comm, R, -grid_coord[C], &src, &dest);
	MPI_Sendrecv_replace(storageB, lrB * lcB, MPI_DOUBLE, dest, 1, src, 1, grid_comm, &status);
	debug("%d: Shifting from Proc %d to Proc %d for B\n",grid_id,src,dest);

	/* compute partial multiplication sum and shift */
	for (int i = 0; i < grid_size[R]; i++)
	{
		mmm(local_rows, local_cols, lcA, subsA, subsB, subsC);
		debug("%d: finished iter %d of mmm\n", grid_id, i);
		MPI_Sendrecv_replace(storageA, lrA * lcA, MPI_DOUBLE, nbr[LEFT], 1, nbr[RIGHT], 1, grid_comm, &status);
		MPI_Sendrecv_replace(storageB, lrB * lcB, MPI_DOUBLE, nbr[UP], 1, nbr[DOWN], 1, grid_comm, &status);
		debug("%d: finished comm %d of mmm\n", grid_id, i);
		
	}

	return;
}

/* allocates 2d sub-matrix with storage pointer */
void allocate_matrix(void ***subs, void **storage, MPI_Datatype dtype, int rows, int cols, MPI_Comm grid_comm)
{
	int datum_size, grid_id, local_rows, local_cols;
	void **lptr, *rptr;
	int grid_size[NDIMS], grid_period[NDIMS], grid_coord[NDIMS];

  MPI_Type_size(dtype, &datum_size);	

	MPI_Comm_rank(grid_comm, &grid_id);
	MPI_Cart_get(grid_comm, NDIMS, grid_size, grid_period, grid_coord);

	/*find number of local cols and rows*/
	local_rows = BLOCK_SIZE(grid_coord[R], grid_size[R], rows);
	local_cols = BLOCK_SIZE(grid_coord[C], grid_size[C], cols);
	debug("%d: local_rows = %d; local_cols = %d;\n",grid_id, local_rows, local_cols);
	
	/*alloc and zero memory for subblock*/
	*storage = my_calloc(grid_id, local_rows * local_cols, datum_size);
	*subs = (void **) my_calloc(grid_id, local_rows, PTR_SIZE);
	lptr = (void *) *subs;
	rptr = (void *) *storage;
	for (int i = 0; i < local_rows; i++)
	{
		*(lptr++) = (void *) rptr;
		rptr += local_cols * datum_size;
	}

	return;
}

/* create grid communicator */
void create_grid_comm(MPI_Comm *comm)
{
	int size, rank, dims[NDIMS], period[NDIMS];

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	dims[R] = dims[C] = (int)sqrt(size); //works only with square # procs
	if (dims[R] * dims[C] != size)
	{
		debug("%d: Error in create_grid_comm: invalid nprocs\n", rank);
		MPI_Finalize();
		exit(0);
	}
	period[R] = period[C] = 1;

	MPI_Cart_create(MPI_COMM_WORLD, NDIMS, dims, period, 1, comm);

	return;
}

/* capture command line arguments */
void get_args(int argc, char *argv[], char **infile1, char **infile2, char **outfile)
{
	if (argc != 4)
	{
		printf("Usage: ./exec <infile1> <infile2> <outfile>\n");
		MPI_Finalize();
		exit(0);
	}
	else
	{
		*infile1 = argv[1];
		*infile2 = argv[2];
		*outfile = argv[3];
	}

	return;
}
