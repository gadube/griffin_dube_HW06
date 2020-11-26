/* mm-parallel.c
 * Griffin Dube
 *
 * MPI Implementation of Cannon's Matrix Multiplication.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <mpi.h>
#include "matrix_checkerboard_io.h"

#define NDIMS 2
#define R 0
#define C 1
#define LEFT 0
#define RIGHT 1
#define UP 2
#define DOWN 3

void getArgs(int argc, char *argv[], char **infile1, char **infile2, char **outfile);
void create_communicators(int size, int rank, MPI_Comm *cart_comm);
void matrix_multiply(double **A, double **B, double ***O,int *d1, int *d2, int **dout, MPI_Comm cart_comm);

int main(int argc, char *argv[])
{
  int rank, size, d1[NDIMS], d2[NDIMS], dout[NDIMS];
	double **M1, **M2, **O, *storage1, *storage2, *storageO;
  char *infile1, *infile2, *outfile;
	MPI_Comm GRID_COMM;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
	if (rank == 0) debug("MPI Size: %d\n",size);
	debug("MPI Rank: %d\n",rank);

	/* read in arguments */
 	getArgs(argc, argv, &infile1, &infile2, &outfile);

	/* create communicator for checkerboard pattern */
	create_communicators(size, rank, &GRID_COMM);
	MPI_Barrier(MPI_COMM_WORLD);
	double timeIO = MPI_Wtime();
	/* read in file and distribute among procs */
	read_checkerboard_graph(infile1, (void ***) &M1, (void **) &storage1, MPI_DOUBLE, (int *) d1, GRID_COMM);
	read_checkerboard_graph(infile2, (void ***) &M2, (void **) &storage2, MPI_DOUBLE, (int *) d2, GRID_COMM);
	debug("%d: in1 dims = %dx%d, in2 dims = %dx%d\n",rank,d1[R],d1[C],d2[R],d2[C]);
	dout[R] = d1[R];
	dout[C] = d2[C];

	MPI_Barrier(MPI_COMM_WORLD);
	double timeNOIO = MPI_Wtime();
	/* matrix multiply */
	matrix_multiply(M1, M2, &O, d1, d2, (int **) &dout, GRID_COMM);
	timeNOIO = MPI_Wtime() - timeNOIO;

	/* write matrix to file */
//	write_checkerboard_graph(outfile, (void ***) &O, (void **) &storageO, MPI_DOUBLE, dout, GRID_COMM);
	timeIO = MPI_Wtime() - timeIO;

	double maxTimeIO,maxTimeNOIO;
	MPI_Reduce(&timeIO,&maxTimeIO,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
	MPI_Reduce(&timeNOIO,&maxTimeNOIO,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

	if (rank == 0) printf("Max Time with IO: %4.6f\nMax Compute Time: %4.6f\n",maxTimeIO,maxTimeNOIO);
	
	/* free memory */
	MPI_Comm_free(&GRID_COMM);
  MPI_Finalize();
  return 0;
}

void matrix_multiply(double **A, double **B, double ***O,int *d1, int *d2, int **dout, MPI_Comm cart_comm)
{
	int grid_id, grid_size[NDIMS], grid_coord[NDIMS], grid_period[NDIMS], subsz1[NDIMS], subsz2[NDIMS];
	int nbr[4];
	int ID, ret, src, dest;
	MPI_Status status;

	MPI_Comm_rank(cart_comm, &grid_id); ID = grid_id;

	debug("%d: Entering Matrix Multiplication Function...\n",ID);

	ret = MPI_Cart_get(cart_comm, NDIMS, grid_size, grid_period, grid_coord);
  error_out(ret, ID, NULL);

  debug( "%d: g_size[0] = %d, g_size[1] = %d\n", ID, grid_size[R], grid_size[C] );
  debug( "%d: g_peri[0] = %d, g_peri[1] = %d\n", ID, grid_period[R], grid_period[C] );
  debug( "%d: g_coor[0] = %d, g_coor[1] = %d\n", ID, grid_coord[R], grid_coord[C] );

  subsz1[R] = BLOCK_SIZE(grid_coord[R], grid_size[R], d1[R]);
  subsz1[C] = BLOCK_SIZE(grid_coord[C], grid_size[C], d1[C]);
  debug( "%d: subsz[0] = %d, subsz[1] = %d\n", ID, subsz1[R], subsz1[C] );
  subsz2[R] = BLOCK_SIZE(grid_coord[R], grid_size[R], d2[R]);
  subsz2[C] = BLOCK_SIZE(grid_coord[C], grid_size[C], d2[C]);
  debug( "%d: subsz[0] = %d, subsz[1] = %d\n", ID, subsz2[R], subsz2[C] );


	ret = MPI_Cart_shift(cart_comm, R, -1, &nbr[RIGHT], &nbr[LEFT]); //find neighbors
	error_out(ret, ID, NULL);
	ret = MPI_Cart_shift(cart_comm, C, -1, &nbr[DOWN], &nbr[UP]); //find neighbors
	error_out(ret, ID, NULL);
	debug("%d: Neighbors: left = %d, right = %d,up = %d, down = %d\n",ID,nbr[LEFT],nbr[RIGHT],nbr[UP],nbr[DOWN]);
	
	debug("%d: Perfoming initial alignment...\n",ID);
	ret = MPI_Cart_shift(cart_comm, R, -grid_coord[C], &src, &dest); //initial alignment of A
	error_out(ret, ID, NULL);
	debug("%d: A shift src = %d dest = %d\n",ID,src,dest);
	ret = MPI_Sendrecv_replace(*A, subsz1[R]*subsz1[C], MPI_DOUBLE, dest, 1, src, 1, cart_comm, &status);
	error_out(ret, ID, &status);
	debug("%d: Completed A initial alignment...\n",ID);

	ret = MPI_Cart_shift(cart_comm, C, -grid_coord[R], &src, &dest); //initial alignment of B
	error_out(ret, ID, NULL);
	debug("%d: B shift src = %d dest = %d\n",ID,src,dest);
	ret = MPI_Sendrecv_replace(*B, subsz2[R]*subsz2[C], MPI_DOUBLE, dest, 1, src, 1, cart_comm, &status);
	error_out(ret, ID, &status);
	debug("%d: Completed B initial alignment...\n",ID);




	return;
}

void create_communicators(int size, int rank, MPI_Comm *cart_comm)
{
	debug("%d: Entering Create Communicators Function...\n",rank);
	int dims[NDIMS], periodic[NDIMS], grid_coords[NDIMS];
	int grid_id, ret;

	periodic[0] = periodic[1] = 1;
	dims[0] = dims[1] = 0;
	
	/*find dimensions of grid*/
	MPI_Dims_create(size, NDIMS, dims);
	/* Create Cartesian communicator */
	debug("%d: Communicator dimensions: %d x %d\n",rank,dims[R],dims[C]);
	debug("%d: Communicator periodicity: row=%s,col=%s\n",rank,(periodic[R])?"yes":"no",(periodic[C])?"yes":"no");
	ret = MPI_Cart_create(MPI_COMM_WORLD, NDIMS, dims, periodic, 1, cart_comm);
	error_out(ret,rank,NULL);

	MPI_Comm_rank(*cart_comm, &grid_id);
	MPI_Cart_coords(*cart_comm, grid_id, NDIMS, grid_coords);
	debug("%d: grid_id=%d coord[0]=%d coord[1]=%d\n",rank,grid_id,grid_coords[R],grid_coords[C]);
	//printf("Coords of Rank %d are (%d,%d).\n",grid_id,grid_coords[0],grid_coords[1]);

	return;
}

void getArgs(int argc, char *argv[], char **infile1, char **infile2, char **outfile)
{
  if(argc != 4) {
    printf("Usage: ./exec <infile1> <infile2> <outfile>\n");
    MPI_Finalize();
    exit(0);
  }
	
	debug("Reading %d input arguments...\n", argc);
	*infile1 = argv[1];
	*infile2 = argv[2];
	*outfile = argv[3];
	debug("Input File 1: %s\nInput File 2: %s\nOutput File: %s\n", *infile1,*infile2,*outfile);
}



