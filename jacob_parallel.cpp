# include <math.h>
# include <mpi.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <time.h>

double L = 1.0;			//computaional domain size LxL
int N = 100;			//grid size NxN

double *u, *u_new;		//solution arrays

#define INDEX(i,j) ((N+2)*(i)+(j))

int my_rank;			//rank of this process

int *proc;			    //process index
int *proc_begin, *proc_end;		 //begin and end indices of processes
int *left_proc, *right_proc;	// left and right processes


  //Functions:

int main ( int argc, char *argv[] );
void allocation ( );
void jacobi ( int num_procs, double f[] );
void domain_decomp ( int num_procs );
double *source_term ( );


int main ( int argc, char *argv[] ) 

{
  double change;
  double epsilon = 1.0E-06;
  double *f;
  char file_name[100];
  int i;
  int j;
  double my_change;
  int my_n;
  int n;
  int num_procs;
  int step;
  double *swap;
  double wall_time;

  //MPI initialization.

  MPI_Init ( &argc, &argv );

  MPI_Comm_size ( MPI_COMM_WORLD, &num_procs );

  MPI_Comm_rank ( MPI_COMM_WORLD, &my_rank );


    N = 100;
    epsilon = 1.0E-07;
  

  if ( my_rank == 0 ) 
  {
    printf ( "\n" );
    printf ( "  Number of processes         = %d\n", num_procs );
    printf ( "  Grid size = %d x %d \n", N,N );
    printf ( "  Convergence tolerence = %f\n", epsilon );
    printf ( "\n" );
  }

  allocation ( );
  f = source_term ( );
  domain_decomp ( num_procs );

  step = 0;

  //start time count.

  wall_time = MPI_Wtime ( );

  //Start iteration.

  do 
  {
    jacobi ( num_procs, f );
    ++step;
 
  // error calculation 

    change = 0.0;
    n = 0;

    my_change = 0.0;
    my_n = 0;

    for ( i = proc_begin[my_rank]; i <= proc_end[my_rank]; i++ )
    {
      for ( j = 1; j <= N; j++ )
      {
        if ( u_new[INDEX(i,j)] != 0.0 ) 
        {
          my_change = my_change 
            + fabs ( 1.0 - u[INDEX(i,j)] / u_new[INDEX(i,j)] );

          my_n = my_n + 1;
        }
      }
    }
    MPI_Allreduce ( &my_change, &change, 1, MPI_DOUBLE, MPI_SUM,
      MPI_COMM_WORLD );

    MPI_Allreduce ( &my_n, &n, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );

    if ( n != 0 )
    {
      change = change / n;
    }
    if ( my_rank == 0 && ( step % 10 ) == 0 ) 
    {
      printf ( "  N = %d, n = %d, my_n = %d, Step %4d  Error = %g\n", 
        N, n, my_n, step, change );
    }
 
  //Interchange U and U_NEW.

    swap = u;
    u = u_new;
    u_new = swap;
  } while ( epsilon < change );


   //wallclock time.

  wall_time = MPI_Wtime() - wall_time;
  if ( my_rank == 0 )
  {
    printf ( "\n" );
    printf ( "  Wall clock time = %f secs\n", wall_time );
  }

 //MPI termination 

  MPI_Finalize ( );
  free ( f );
  
  if ( my_rank == 0 )
  {
    printf ( "\n" );
    printf ( "POISSON_MPI:\n" );
    printf ( "  Normal end of execution.\n" );
    printf ( "\n" );
    
  }
 
  return 0;
}


//Function for allocation, initialization and boundary conditions//

void allocation ( ) 

{
  int i;
  int grid_size;

  grid_size = ( N + 2 ) * ( N + 2 );

  u = ( double * ) malloc ( grid_size * sizeof ( double ) );
  for ( i = 0; i < grid_size; i++)
  {
    u[i] = 0.0;
  }

  u_new = ( double * ) malloc ( grid_size * sizeof ( double ) );
  for ( i = 0; i < grid_size; i++ )
  {
    u_new[i] = 0.0;
  }

  return;
}


//Function for jacobi algorithm//

void jacobi ( int num_procs, double f[] ) 



{
  double h;
  int i;
  int j;
  MPI_Request request[4];
  int requests;
  MPI_Status status[4];

  h = L / ( double ) ( N + 1 );

  requests = 0;

  if ( left_proc[my_rank] >= 0 && left_proc[my_rank] < num_procs ) 
  {
    MPI_Irecv ( u + INDEX(proc_begin[my_rank] - 1, 1), N, MPI_DOUBLE,
      left_proc[my_rank], 0, MPI_COMM_WORLD,
      request + requests++ );

    MPI_Isend ( u + INDEX(proc_begin[my_rank], 1), N, MPI_DOUBLE,
      left_proc[my_rank], 1, MPI_COMM_WORLD,
      request + requests++ );
  }

  if ( right_proc[my_rank] >= 0 && right_proc[my_rank] < num_procs ) 
  {
    MPI_Irecv ( u + INDEX(proc_end[my_rank] + 1, 1), N, MPI_DOUBLE,
      right_proc[my_rank], 1, MPI_COMM_WORLD,
      request + requests++ );

    MPI_Isend ( u + INDEX(proc_end[my_rank], 1), N, MPI_DOUBLE,
      right_proc[my_rank], 0, MPI_COMM_WORLD,
      request + requests++ );
  }
 
 // update internal vertices in my domain.

  for ( i = proc_begin[my_rank] + 1; i <= proc_end[my_rank] - 1; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      u_new[INDEX(i,j)] =
        0.25 * ( u[INDEX(i-1,j)] + u[INDEX(i+1,j)] +
                 u[INDEX(i,j-1)] + u[INDEX(i,j+1)] +
                 h * h * f[INDEX(i,j)] );
    }
  }
 
  //Wait for completetion of non-blocking communications 

  MPI_Waitall ( requests, request, status );

 //update boundary vertices in my domain.

  i = proc_begin[my_rank];
  for ( j = 1; j <= N; j++ )
  {
    u_new[INDEX(i,j)] =
      0.25 * ( u[INDEX(i-1,j)] + u[INDEX(i+1,j)] +
               u[INDEX(i,j-1)] + u[INDEX(i,j+1)] +
               h * h * f[INDEX(i,j)] );
  }

  i = proc_end[my_rank];
  if (i != proc_begin[my_rank])
  {
    for (j = 1; j <= N; j++)
    {
      u_new[INDEX(i,j)] =
        0.25 * ( u[INDEX(i-1,j)] + u[INDEX(i+1,j)] +
                 u[INDEX(i,j-1)] + u[INDEX(i,j+1)] +
                 h * h * f[INDEX(i,j)] );
    }
  }

  return;
}



//Function for domain decomposition//

void domain_decomp ( int num_procs ) 


{
  double d;
  double eps;
  int i;
  int p;
  double x_max;
  double x_min;
 
  //Allocate arrays 

  proc = ( int * ) malloc ( ( N + 2 ) * sizeof ( int ) );
  proc_begin = ( int * ) malloc ( num_procs * sizeof ( int ) );
  proc_end = ( int * ) malloc ( num_procs * sizeof ( int ) );
  left_proc = ( int * ) malloc ( num_procs * sizeof ( int ) );
  right_proc = ( int * ) malloc ( num_procs * sizeof ( int ) );

  eps = 0.0001;
  d = ( N - 1.0 + 2.0 * eps ) / ( double ) num_procs;

  for ( p = 0; p < num_procs; p++ )
  {

    x_min = - eps + 1.0 + ( double ) ( p * d );
    x_max = x_min + d;

    for ( i = 1; i <= N; i++ )
    {
      if ( x_min <= i && i < x_max )
      {
        proc[i] = p;
      }
    }
  }

  for ( p = 0; p < num_procs; p++ )
  {
    for ( i = 1; i <= N; i++ )
    {
      if ( proc[i] == p )
      {
        break;
      }
    }
    proc_begin[p] = i;

    for ( i = N; 1 <= i; i-- )
    {
      if ( proc[i] == p )
      {
        break;
      }
    }
    proc_end[p] = i;

 //left and right procs. 

    left_proc[p] = -1;
    right_proc[p] = -1;

    if ( proc[p] != -1 ) 
    {
      if ( 1 < proc_begin[p] && proc_begin[p] <= N )
      {
        left_proc[p] = proc[proc_begin[p] - 1];
      }
      if ( 0 < proc_end[p] && proc_end[p] < N )
      {
        right_proc[p] = proc[proc_end[p] + 1];
      }
    }
  }

  return;
}


//function for source term//


double *source_term ( ) 


{
  double *f;
  int i;
  int j;
  int k;
  double q;

  f = ( double * ) malloc ( ( N + 2 ) * ( N + 2 ) * sizeof ( double ) );

  for ( i = 0; i < ( N + 2 ) * ( N + 2 ); i++ )
  {
    f[i] = 0.0;
  }
 
 // dipole.

  q = 10.0;

  i = 1 + N / 4;
  j = i;
  k = INDEX ( i, j );
  f[k] = q;

  i = 1 + 3 * N / 4;
  j = i;
  k = INDEX ( i, j );
  f[k] = -q;
  

  return f;
}
