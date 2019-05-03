# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>
# include <cmath>

using namespace std;

# define NX 100
# define NY 100

int main ( int argc, char *argv[] );
double rms_norm ( int m, int n, double a[NX][NY] );
void initialize ( int nx, int ny, double f[NX][NY] );
void jacob ( int nx, int ny, double dx, double dy, double f[NX][NY], 
double u[NX][NY], double unew[NX][NY] );
double u_exact ( double x, double y );
double source ( double x, double y );

//main program//

int main ( int argc, char *argv[] )
{
  bool converged;
  double diff, dx, dy, error;
  double f[NX][NY];
  int i, j, it;
  int it_max = 100000;
  int nx = NX;
  int ny = NY;
  double tolerance = 0.000001;
  double u[NX][NY];
  double u_norm;
  double udiff[NX][NY];
  double uexact[NX][NY];
  double unew[NX][NY];
  double unew_norm;
  double x, y;

  dx = 1.0 / ( double ) ( nx - 1 );
  dy = 1.0 / ( double ) ( ny - 1 );
  
  
  //Initialize
  
  initialize(nx,ny,f);
  
  // initial approximation and BC
  
for ( j = 0; j < ny; j++ )
  {
    for ( i = 0; i < nx; i++ )
    {
      if ( i == 0 || i == nx - 1 || j == 0 || j == ny - 1 )
      {
        unew[i][j] = f[i][j];
      }
      else
      {
        unew[i][j] = 0.0;
      }
    }
  }
  
 // RMS norm of vector unew

 unew_norm = rms_norm ( nx, ny, unew );


 //exact solution
 
 for ( j = 0; j < ny; j++ )
  {
    y = ( double ) ( j ) / ( double ) ( ny - 1 );
    for ( i = 0; i < nx; i++ )
    {
      x = ( double ) ( i ) / ( double ) ( nx - 1 );
      uexact[i][j] = u_exact ( x, y );
    }
  }
  u_norm = rms_norm ( nx, ny, uexact );
  
  
  //iteration
  
   converged = false;

  cout << "\n";
  cout << "  Step    ||Unew||     ||Unew-U||     ||Unew-Exact||\n";
  cout << "\n";

  for ( j = 0; j < ny; j++ )
  {
    for ( i = 0; i < nx; i++ )
    {
      udiff[i][j] = unew[i][j] - uexact[i][j];
    }
  }
  error = rms_norm ( nx, ny, udiff );
  cout << "  " << setw(4) << 0
       << "  " << setw(14) << unew_norm
       << "  " << "              "
       << "  " << setw(14) << error << "\n";
	   
//iteration loop

  for ( it = 1; it <= it_max; it++ )
  {
    for ( j = 0; j < ny; j++ )
    {
      for ( i = 0; i < nx; i++ )
      {
        u[i][j] = unew[i][j];
      }
    }
	
	//one jacobi iteration
	
	 jacob ( nx, ny, dx, dy, f, u, unew );
	 
	 //check if solution has converged
	 
	 u_norm = unew_norm;
    unew_norm = rms_norm ( nx, ny, unew );

    for ( j = 0; j < ny; j++ )
    {
      for ( i = 0; i < nx; i++ )
      {
        udiff[i][j] = unew[i][j] - u[i][j];
      }
    }
    diff = rms_norm ( nx, ny, udiff );

    for ( j = 0; j < ny; j++ )
    {
      for ( i = 0; i < nx; i++ )
      {
        udiff[i][j] = unew[i][j] - uexact[i][j];
      }
    }
    error = rms_norm ( nx, ny, udiff );

    cout << "  " << setw(4)  << it
         << "  " << setw(14) << unew_norm
         << "  " << setw(14) << diff
         << "  " << setw(14) << error << "\n";

    if ( diff <= tolerance )
    {
      converged = true;
      break;
    }

  }
  
  //print convergence result

  if ( converged )
  {
    cout << "  solution has converged.\n";
  }
  else
  {
    cout << "  solution has NOT converged.\n";
  }
  
  return 0;
}

//end of main program


// function to calculate the rms norm of a vector//

double rms_norm ( int nx, int ny, double a[NX][NY] )
{
  int i;
  int j;
  double v;

  v = 0.0;

  for ( j = 0; j < ny; j++ )
  {
    for ( i = 0; i < nx; i++ )
    {
      v = v + a[i][j] * a[i][j];
    }
  }
  v = sqrt ( v / ( double ) ( nx * ny )  );

  return v;
}


//function to initialize interior and boundary values//

void initialize ( int nx, int ny, double f[NX][NY] )

{
  double fnorm;
  int i;
  int j;
  double x;
  double y;

  
  for ( j = 0; j < ny; j++ )
  {
    y = ( double ) ( j ) / ( double ) ( ny - 1 );
    for ( i = 0; i < nx; i++ )
    {
      x = ( double ) ( i ) / ( double ) ( nx - 1 );
	 
      if ( i == 0 || i == nx - 1 || j == 0 || j == ny - 1 )
      {
        f[i][j] = u_exact ( x, y );
      }
	  
      else
      {
        f[i][j] = - source( x, y );
      }
    }
  }

  fnorm = rms_norm ( nx, ny, f );

  cout << "  RMS of F = " << fnorm << "\n";

  return;
}


// function to perform jacobi iteration//

void jacob ( int nx, int ny, double dx, double dy, double f[NX][NY], 
  double u[NX][NY], double unew[NX][NY] )

 {
  int i;
  int j;

  for ( j = 0; j < ny; j++ )
  {
    for ( i = 0; i < nx; i++ )
    {
      if ( i == 0 || j == 0 || i == nx - 1 || j == ny - 1 )
      {
        unew[i][j] = f[i][j];
      }
      else
      { 
        unew[i][j] = 0.25 * ( 
          u[i-1][j] + u[i][j+1] + u[i][j-1] + u[i+1][j] + f[i][j] * dx * dy );
      }
    }
  }
  return;
}


//function to calculate exact solution//

double u_exact ( double x, double y )

{
  double pi = 3.141592653589793;
  double value;

  value = sin ( pi * x * y );

  return value;
}

//function to calculate source terms//

double source ( double x, double y )

{
  double pi = 3.141592653589793;
  double value;

  value = - pi * pi * ( x * x + y * y ) * sin ( pi * x * y );

  return value;
}

# undef NX
# undef NY
 
  