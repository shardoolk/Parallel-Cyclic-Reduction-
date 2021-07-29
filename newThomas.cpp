#include<iostream>
#include<mpi.h>
#include<omp.h>
#include<sys/time.h>

/*-------------------------------------------------------------------------------*/
/* IMPORTANT: PLEASE MAKE SURE THAT TOTAL_VOXELS_X IS PERFECTLY DIVISIBLE BY THE */
/* TOTAL NUMBER OF PROCESSES 																										 */
/*-------------------------------------------------------------------------------*/



// #define TOTAL_VOXELS_X 1000									//These will be split up between processes
// #define TOTAL_VOXELS_Y 40											//Y-dimension is not split up (Equivalent to Total Linear Systems)
// #define TOTAL_VOXELS_Z 10	                                    //Z-dimension is not split up

int TOTAL_VOXELS_X,	TOTAL_VOXELS_Y, TOTAL_VOXELS_Z;	

struct timeval begin,end;
struct timeval begin1, end1;
int main(int argc, char *argv[])
{
	

/*----------------------------*/
/* Variable Declaration block */
/*----------------------------*/
	
	int mpi_Provided;														//Requested thread level 
	int mpi_Rank; 															//Rank of a process in MPI_COMM_WORLD
	int mpi_Size; 															//Size of the MPI_COMM_WORLD communicator
	MPI_Comm mpi_Old_comm;											//Initial communicator i.e. MPI_COMM_WORLD of MPI
	MPI_Comm mpi_Cart_comm;											//Handle for new Cartesian communicator
	int mpi_Dims[3];														//Dimensions of Cartesian topology 
	int mpi_Is_periodic[3];											//Periodicity of each dimension
	int mpi_Coords[3];													//Coordinates of each MPI process
	int mpi_Reorder;														//Reorder the processes in the new communicator or not
	
	int x_voxels_pp;														//X voxels per process
	int y_voxels_pp;														//Y voxels per process
	int z_voxels_pp;														//Z voxels per process
	
	
	mpi_Old_comm = MPI_COMM_WORLD; 		//A handle to the communicator
	
/*--------------------*/
/* MPI Initialization */
/*--------------------*/
	
	MPI_Init_thread(NULL, NULL, MPI_THREAD_FUNNELED,&mpi_Provided);
	if(mpi_Provided != MPI_THREAD_FUNNELED)
		MPI_Abort(mpi_Old_comm, -1);
		
	MPI_Comm_size(mpi_Old_comm, &mpi_Size);					//Size of the old communicator i.e. Processes at the beginning 
	MPI_Comm_rank(mpi_Old_comm, &mpi_Rank);					//Rank of each process in the communicator "old_comm"
	
/*----------------------------------------------*/	
/* Specify parameters for creating MPI topology */
/*----------------------------------------------*/
	mpi_Dims[0] = 1;																//Processes in X-direction i.e. Top to down
	mpi_Dims[1] = mpi_Size;                    			//Processes in Y-direction i.e. Left to right
	mpi_Dims[2] = 1;																//Processes in Z-direction i.e. Into the paper


	mpi_Is_periodic[0]=0;														//Aperiodic x-dimension
	mpi_Is_periodic[1]=0;														//Aperiodic y-dimension
	mpi_Is_periodic[2]=0;														//Aperiodic z-dimension
	mpi_Reorder = 0;																//Do not reorder ranks in the new topology
	
/*------------------------------------------------------------*/
/*			Create MPI topology, find coords                  		*/
/*------------------------------------------------------------*/

	MPI_Cart_create(mpi_Old_comm, 3, mpi_Dims, mpi_Is_periodic, mpi_Reorder, &mpi_Cart_comm);
	MPI_Cart_coords(mpi_Cart_comm, mpi_Rank, 3, mpi_Coords);	

/*------------------------------------*/
/* Splitting x-voxels among processes */
/*------------------------------------*/

double sizeinput[3] = {1024000,2048000,4096000};
    int arrSize = sizeof(sizeinput)/sizeof(sizeinput[0]);
    
    for (int k = 0; k < arrSize; k++){
    TOTAL_VOXELS_X = sizeinput[k];									
    TOTAL_VOXELS_Y = sizeinput[k]/200;									
    TOTAL_VOXELS_Z = 10;

   

 x_voxels_pp = TOTAL_VOXELS_X/mpi_Size; 			//Or this can be x_voxels_pp=TOTAL_VOXELS_X/mpi_Dims[1]
 y_voxels_pp = TOTAL_VOXELS_Y;								//Y dimension not split up
 z_voxels_pp = TOTAL_VOXELS_Z; 								//Z dimension not split up

/*------------------------------------------------------------------------------------------------------*/ 
/* These x_voxels_pp are equivalent to the number of unknowns per process per linear system 						*/
/* In the X-direction there are a total of TOTAL_VOXELS_Y linear systems in a single X-Y plane 					*/
/* If we consider 3 dimensions then there are a total of TOTAL_VOXELS_Y * TOTAL_VOXELS_Z systems				*/
/* EACH linear system in the X-direction gives rise to a TOTAL_VOXELS_X * TOTAL_VOXELS_X matrix				  */
/* Each process will have a subdiagonal(a), main diagonal(b) and super diagonal(c) and rhs (d) array 		*/
/* for EVERY linear system (and also for unknowns x) 																										*/
/* Each of these arrays will be of length x_voxels_pp+2 on each process (2 ghost elements)							*/
/* Thus, we have "a[TOTAL_VOXELS_Y][x_voxels_pp+2]" and similarly for b,c,d and x. 											*/
/*------------------------------------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------------------*/
/* Declaration of arrays holding sub-diagonal, main diagonal, super-diagonal and rhs and x */
/* Memory allocation 																																 			 */
/*-----------------------------------------------------------------------------------------*/

 double **a, **b, **c, **d, **x;
 
 a = new double* [TOTAL_VOXELS_Y];
 b = new double* [TOTAL_VOXELS_Y];
 c = new double* [TOTAL_VOXELS_Y];
 d = new double* [TOTAL_VOXELS_Y];
 x = new double* [TOTAL_VOXELS_Y];

 for(int i=0; i<=TOTAL_VOXELS_Y-1; i++)
	 {
		 a[i] = new double[x_voxels_pp+2];
		 b[i] = new double[x_voxels_pp+2];
		 c[i] = new double[x_voxels_pp+2];
		 d[i] = new double[x_voxels_pp+2];
		 x[i] = new double[x_voxels_pp+2];
	 }

/*----------------------------------------------------------------------------------------------------------*/
/* The indexing scheme for the coefficients 1<=i<=x_voxels_pp. Indexes 0 and x_voxels_pp+1 will hold values */	
/* from previous and next neighbour processes (as and when required) 																				*/
/* a[1] on Rank 0 is zero 																																									*/
/* c[x_voxels_pp] on Rank = mpi_Size-1 is 0																																	*/
/* Otherwise a[i]=c[i]=-1 and b[i]=2 and d[i]=1 and x[i]=0																									*/ 
/*----------------------------------------------------------------------------------------------------------*/

 for(int i=0; i<=TOTAL_VOXELS_Y-1; i++)
 {
 	for(int j=1;j<=x_voxels_pp;j++)
 	{
 		a[i][j] = -1.0;
 		c[i][j] = -1.0; 
 		b[i][j] =  2.0;	
 		d[i][j] =  1.0;									//Can later add thread ID + mpi_Rank as well 
 		x[i][j] =  0.0; 								//WHAT IS THE NEED FOR THIS ? REMOVE LATER
 	}
 }

/*--------------------------------------------------------------------------------*/ 
/* On first process i.e. Rank = 0, a[1] should all be zero for all linear systems */
/*--------------------------------------------------------------------------------*/ 
 
 if(mpi_Rank == 0)
 {
 	for(int i=0; i<=TOTAL_VOXELS_Y-1;i++)
 		a[i][1] = 0.0; 													
 }

/*--------------------------------------------------------------------------------------------------*/ 
/* On last process i.e. Rank = mpi_Size-1, c[x_voxels_pp] should all be zero for all linear systems */
/*--------------------------------------------------------------------------------------------------*/ 
 
 if(mpi_Rank == mpi_Size-1)
 {
 	for(int i=0; i<=TOTAL_VOXELS_Y-1;i++)
 		c[i][x_voxels_pp] = 0.0;
 }
 
 

 gettimeofday(&begin,0);

for(int i=0; i<=TOTAL_VOXELS_Y-1;i++)
{
/*---------------------------------------------------------*/ 
/* Step-1(a) Forward pass of the Modified Thomas Algorithm */
/* Elimination of lower diagonal elements 								 */ 
/*---------------------------------------------------------*/

	a[i][1] = a[i][1]/b[i][1];
	c[i][1] = c[i][1]/b[i][1];
	d[i][1] = d[i][1]/b[i][1];

	a[i][2] = a[i][2]/b[i][2];
	c[i][2] = c[i][2]/b[i][2];
	d[i][2] = d[i][2]/b[i][2];
	
	for(int j=3; j<=x_voxels_pp;j++)														//Start from third equation 
		{
			double r = 1.0 / ( b[i][j] - a[i][j] * c[i][j-1] );
			d[i][j]  = r * ( d[i][j] - a[i][j] * d[i][j-1] );
			c[i][j]  = r * c[i][j];
			a[i][j]  = -1.0 * r * a[i][j] * a[i][j-1];  
		}
	
/*----------------------------------------------------------*/ 
/* Step-1(b) Backward pass of the Modified Thomas Algorithm */ 
/* Elimination of Upper Diagonal matrices 									*/
/*----------------------------------------------------------*/		

	for(int j=x_voxels_pp-2; j>=2; j--)												//Start from third last equation
	{
		d[i][j] = d[i][j] - c[i][j] * d[i][j+1];  
		a[i][j] = a[i][j] - c[i][j] * a[i][j+1]; 								//REMEMBER: use c[i][j] BEFORE it is updated ! 
		c[i][j] = -1.0 * c[i][j] * c[i][j+1];
	}
	
/*------------------------------------------------------------------------------------------------------*/	
/* Till now all the b's should have become = 1 (although we have not made the 'b' array explicitly = 1) */	
/* The expressions for 'r' and 'd' are wrong in the paper - consult the Fortran program 								*/
/*------------------------------------------------------------------------------------------------------*/	

	double r = 1.0 / ( 1.0 - a[i][2] * c[i][1] );							//r=1.0/(b[i][2]-a[i][2]*c[i][1] but b[i][2]=1.0 now)
	d[i][1]  = r * (d[i][1] - c[i][1] * d[i][2]);							//This was also incorrect in the paper
	c[i][1]  = -1.0 * r * c[i][1] * c[i][2];
	a[i][1]  = r * a[i][1];  

}


	
/*--------------------------------------------------------------------------------------------------------------------------------*/
/* Now every process will send the 1st and last "local" equation of each linear system that resides on it 												*/
/* i.e. each process collects coefficients of a/b/c/d[0<=i<=TOTAL_VOXELS_Y-1][1] and a/b/c/d[0<=i<=TOTAL_VOXELS_Y-1][x_voxels_pp] */
/* and sends it to "some" chosen process so that all equations can be solved 																											*/	
/*--------------------------------------------------------------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------------------------------------*/
/* Extract the first and last "local" equation coefficients on each process and store in contiguous array */
/*--------------------------------------------------------------------------------------------------------*/

 double *red_abcd;
 red_abcd = new double[TOTAL_VOXELS_Y * 2 * 4]; //Total linear systems = TOTAL_VOXELS_Y, Each linear system contributes 2 equations, coefficients in each eq=4
 int counter=0;

/*-------------------------------------------------------------------------------------------------------------------------------------------*/ 
/* Storage sequence of coefficients : (0th system first) a1 b1 c1 d1 a2 b2 c2 d2 (1st system second) a3 b3 c3 d3 a4 b4 c4 d4 (3rd system) ...*/
/*-------------------------------------------------------------------------------------------------------------------------------------------*/ 


 for(int i=0; i<=TOTAL_VOXELS_Y-1; i++)
 {
 	red_abcd[counter++] = a[i][1]; 
 	red_abcd[counter++] = 1.0; 									//all b's are 1 
 	red_abcd[counter++] = c[i][1];
 	red_abcd[counter++] = d[i][1];
 	
 	red_abcd[counter++] = a[i][x_voxels_pp];
 	red_abcd[counter++] = 1.0; 									//all b's are 1
 	red_abcd[counter++] = c[i][x_voxels_pp];
 	red_abcd[counter++] = d[i][x_voxels_pp];
 }	
	
/*---------------------------------------------------------------------------------------*/
/* Choose a rank at which we will solve these reduced systems of Equations 							 */
/* Then declare an array big enough on this rank to contain these coefficients and 			 */
/* put them in the right format i.e. a, b, c, d <--- easier to apply Thomas solver 			 */
/* Later optimization: Distribute systems on multiple processes so that different systems*/
/* can be solved in parallel - this will be a new thing and we keep multiple cores busy  */
/*---------------------------------------------------------------------------------------*/	

 int chosen_rank = 0;
 double *all_red_abcd;															//Declaration on all processes but memory allocation only on chosen_rank
 double **red_a, **red_b, **red_c, **red_d;					//Used to store COMPLETE reduced systems of equations
 
 
 if(mpi_Rank == chosen_rank)
 {
 	all_red_abcd = new double[mpi_Size * TOTAL_VOXELS_Y * 2 * 4];				//Used in MPI_Gather to receive data
 	
 	red_a = new double* [TOTAL_VOXELS_Y];
 	red_b = new double* [TOTAL_VOXELS_Y];
 	red_c = new double* [TOTAL_VOXELS_Y];
 	red_d = new double* [TOTAL_VOXELS_Y];
 	
 	for(int i=0; i<=TOTAL_VOXELS_Y-1; i++)
 	{
 		red_a[i] = new double[2 * mpi_Size];
 		red_b[i] = new double[2 * mpi_Size];
 		red_c[i] = new double[2 * mpi_Size];
 		red_d[i] = new double[2 * mpi_Size];
 	}	 
 }
 
gettimeofday(&begin1,0);
 MPI_Gather(red_abcd, TOTAL_VOXELS_Y * 2 * 4, MPI_DOUBLE, all_red_abcd, TOTAL_VOXELS_Y * 2 * 4 , MPI_DOUBLE, chosen_rank, mpi_Cart_comm);

 if(mpi_Rank == chosen_rank)
 {
 		for(int i=0; i<=TOTAL_VOXELS_Y-1; i++)
 		{
 			counter = i * 8; 
 			for(int j=0; j<=2 * mpi_Size-1; j=j+2)
 			{
 			red_a[i][j] = all_red_abcd[counter++];
 			red_b[i][j] = all_red_abcd[counter++];
 			red_c[i][j] = all_red_abcd[counter++];
 			red_d[i][j] = all_red_abcd[counter++];
 			
 			red_a[i][j+1] = all_red_abcd[counter++];
 			red_b[i][j+1] = all_red_abcd[counter++];
 			red_c[i][j+1] = all_red_abcd[counter++];
 			red_d[i][j+1] = all_red_abcd[counter++];
 			
 			counter = counter + (TOTAL_VOXELS_Y-1) * 2 * 4; 			 
 			}
 		}
 }
 
 /*------------------------------------------------------------------------------------*/
 /* Serial Thomas algorithm on "chosen_rank" to solve the reduced systems of Equations */
 /* Right now a single thread will be used but later multiple threads 								 */
 /*------------------------------------------------------------------------------------*/

	if(mpi_Rank == chosen_rank)
	{
		 for(int i=0; i<=TOTAL_VOXELS_Y-1;i++)
		 {
		 	red_d[i][0] = red_d[i][0]/red_b[i][0];   //b[i][whatever]=1 for all the equations so later change to 1.0
		 	red_c[i][0] = red_c[i][0]/red_b[i][0]; 	 //Similarly no need for this step 
		 	
		 	for(int j=1; j<=2 * mpi_Size-1; j++)
		 	{
		 		double r = 1.0 / ( red_b[i][j] - red_a[i][j] * red_c[i][j-1] );
		 		red_d[i][j]  = r * ( red_d[i][j] - red_a[i][j] * red_d[i][j-1] );
		 		red_c[i][j]  = r * red_c[i][j]; 
		 	}
		 }
		 
		 for(int i=0; i<=TOTAL_VOXELS_Y-1;i++)
		 {
		    for(int j=2*mpi_Size-2;j>=0;j--)
		    {
		    	red_d[i][j] = red_d[i][j] - red_c[i][j] * red_d[i][j+1]; 
		    }
		 } 
	}

/*-----------------------------------------------------------------------------------------------------------------------------*/  
/* Now the solutions of reduced system of equations is in red_d[][] arrays. 							 																		 */
/* red_d[0<=i<=TOTAL_VOXELS_Y][0] and red_d[0<=i<=TOTAL_VOXELS_Y][1] - should go to rank 0 																		 */
/* red_d[0<=i<=TOTAL_VOXELS_Y][2] and red_d[0<=i<=TOTAL_VOXELS_Y][3] - should go to rank 1 																		 */
/* red_d[0<=i<=TOTAL_VOXELS_Y][4] and red_d[0<=i<=TOTAL_VOXELS_Y][5] - should go to rank 2 																		 */
/* ... 																																										 																		 */
/* red_d[0<=i<=TOTAL_VOXELS_Y][2 * mpi_Size-2] and red_d[0<=i<=TOTAL_VOXELS_Y][2 * mpi_Size-1] - should go to rank (mpi_Size-1)*/	
/*-----------------------------------------------------------------------------------------------------------------------------*/	

  double *contig_red_d;				//Declaration on all but allocation only on chosen_rank 
  
  if(mpi_Rank == chosen_rank)
  {
  	contig_red_d = new double[TOTAL_VOXELS_Y * 2 * mpi_Size];
  	int counter = 0; 
  	for(int j=0; j<=2*mpi_Size-1; j=j+2)
  	{ 
  		for(int i=0; i<=TOTAL_VOXELS_Y-1;i++)
  		{
  			contig_red_d[counter++] = red_d[i][j];
  			contig_red_d[counter++] = red_d[i][j+1]; 		
  		}
  	}
  }
  
/*---------------------------------------------------------------------------------------------------------*/
/* Declare the buffer to contain the incoming solutions on each process 																	 */
/* Length of this buffer is TOTAL_VOXELS_Y * 2 i.e. Each process gets 2 solved variables per linear system */
/*---------------------------------------------------------------------------------------------------------*/

  double *recv_red_d = new double[TOTAL_VOXELS_Y * 2]; 
  
  MPI_Scatter(contig_red_d, TOTAL_VOXELS_Y * 2, MPI_DOUBLE, recv_red_d, TOTAL_VOXELS_Y * 2, MPI_DOUBLE, chosen_rank, mpi_Cart_comm); 
  gettimeofday(&end1,0);

long scs = end1.tv_sec-begin1.tv_sec;
long mcs = end1.tv_usec - begin1.tv_usec;
double commtime = scs + mcs*1e-6;
if(mpi_Rank==0){
	std::cout <<"The time elapsed for MPI scatter and gather is : " << commtime << std::endl;
}
/*----------------------------------------------------------------------------------------------*/
/* Now put the received solutions of the first and last local equation in the local d[][] array */
/*----------------------------------------------------------------------------------------------*/
	
	counter = 0;
	for(int i=0; i<=TOTAL_VOXELS_Y-1;i++)
	{
		d[i][1] 					= recv_red_d[counter++];
		d[i][x_voxels_pp] = recv_red_d[counter++]; 
	}

/*----------------------------------------------------------------------------------------------------------*/	
/* Now the solutions to the first and last local variable on each process is in d[][1] and d[][x_voxels_pp] */
/* Use this on every process to solve the local 2nd to 2nd last variable 																		*/
/* All the solutions then are in 'd' array 																																	*/
/*----------------------------------------------------------------------------------------------------------*/

 for(int i=0; i<=TOTAL_VOXELS_Y-1; i++)
 {
 	for(int j=2; j<=x_voxels_pp-1;j++)
 	{
 		d[i][j] = d[i][j] - a[i][j] * d[i][1] - c[i][j] * d[i][x_voxels_pp]; 
 	}
 }

 
gettimeofday(&end,0);
 long seconds = end.tv_sec - begin.tv_sec;
 long microseconds = end.tv_usec - begin.tv_usec;

 double elapsed = seconds + microseconds*1e-6;

if (mpi_Rank == 0){
    std::cout << "The total time elapsed is : " << elapsed << std::endl;

}

/*----------------------------------------------------------------------------------------------------*/ 
/* To see the final solution for Linear System 'i', print all d[i][1<=j<=x_voxels_pp] on each process */
/*----------------------------------------------------------------------------------------------------*/

 int system_to_be_checked = 0; 
 if (mpi_Rank == system_to_be_checked){
     std::cout<< "The number of processors are " << mpi_Size<< std::endl;
    std::cout <<" The size of the array is  " << TOTAL_VOXELS_X << std::endl;
 }
// for(int rnk_ctr = 0; rnk_ctr <= mpi_Size-1; rnk_ctr++)
// {
// 	if(mpi_Rank == rnk_ctr)

// 	{
 		// for(int i=1; i<=x_voxels_pp; i++)
 		// {
 		// 	std::cout<<d[system_to_be_checked][i]<<" " << std::endl; 
 		// }
 //	}
 //		MPI_Barrier(mpi_Cart_comm); 
 
// }

delete[] a;
delete[] b;
delete[] c;
delete[] d;  
delete[] x;
}


  
	MPI_Finalize();
    
    
}
