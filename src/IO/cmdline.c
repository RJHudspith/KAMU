/**
   @file cmdline.c
   @brief read some arguments from the command line
 */
#include "KAMU.h" // struct integral_args

enum { INPUT_MV = 1 , INPUT_TOL = 2 , INPUT_SYM = 3 ,
       INPUT_YMIN = 4 , INPUT_SHIFT = 5 } ;

// strcmp defaults to 0 if they are equal which is contrary to standard if statements
static int
are_equal( const char *str_1 , const char *str_2 ) {
  return !strcmp( str_1 , str_2 ) ;
}

struct integral_args
read_cmdline( int *err ,
	      const char *argv[] )
{  
  // struct with all the data in
  struct integral_args fdata ;
  *err = 0 ;
  
  fdata.Mv = (double)atof( argv[INPUT_MV] ) ;
  if( fdata.Mv < 0 ) {
    fprintf( stderr , "[KAMU] masses usually cannot be negative %e\n" ,
	     fdata.Mv ) ;
    *err = 1 ;
  }
  fdata.Tol = (double)atof( argv[INPUT_TOL] ) ;
  if( fdata.Tol < 0 || fdata.Tol > 100 ) {
    fprintf( stderr , "[KAMU] Tol is strange %e\n" , fdata.Tol ) ;
    *err = 1 ;
  }
  fdata.xmin[1] = 0.0 ;
  fdata.xmin[2] = (double)atof( argv[INPUT_YMIN] ) ;
  if( fdata.xmin[2] < 0 || fdata.xmin[2] > 6.54 ) {
    fprintf( stderr , "[KAMU] Ymin is strange %e\n" , fdata.xmin[2] ) ;
    *err = 1 ;
  }
  if( are_equal( argv[INPUT_SYM] , "NORMAL" ) ) {
    fdata.Sym = NORMAL ;
  } else if( are_equal( argv[INPUT_SYM] , "MKERN" ) ) {
    fdata.Sym = MKERN ;
  } else if( are_equal( argv[INPUT_SYM] , "SYMXY" ) ) {
    fdata.Sym = SYMXY ;
  } else if( are_equal( argv[INPUT_SYM] , "SYMXY0" ) ) {
    fdata.Sym = SYMXY0 ;
  } else {
    fprintf( stderr , "[KAMU] I do not understand symmetrise option %s\n" ,
	     argv[INPUT_SYM] ) ;
    *err = 1 ;
  }
  fdata.shift = false ;
  if( are_equal( argv[INPUT_SHIFT] , "true" ) ) {
    fdata.shift = true ;
  }
  
  // init the kernel
  initialise( &fdata.t ) ;
  
  fdata.xmin[0] = 0.0 ;
  fdata.xmax[1] = fdata.t.Grid.XX[ fdata.t.Grid.nstpx - 1 ] ;
  fdata.xmax[2] = fdata.t.Grid.YY[ fdata.t.Grid.nstpy - 1 ] ;
  fdata.xmax[0] = M_PI ;
  fdata.K = NULL ;
  
  // summarise the parameters
  fprintf( stdout , "[KAMU] Parameter order is beta, x , y\n" ) ;
  fprintf( stdout , "[KAMU] Xmin being used (%g , %g , %g)\n" ,
	   fdata.xmin[0] , fdata.xmin[1] , fdata.xmin[2] ) ;
  fprintf( stdout , "[KAMU] Xmax being used (%g , %g , %g)\n" ,
	   fdata.xmax[0] , fdata.xmax[1] , fdata.xmax[2] ) ;
  fprintf( stdout , "[KAMU] using Lepton loop mass ratio %g\n" , fdata.Mv ) ;
  fprintf( stdout , "[KAMU] using cubature tolerance %g\n" , fdata.Tol ) ;
  switch( fdata.Sym ) {
  case NORMAL : fprintf( stdout , "[KAMU] using the default kernels\n" ) ;
    fdata.K = malloc( 3*sizeof( struct Kernels ) ) ;
    break ;
  case MKERN : fprintf( stdout , "[KAMU] using the M-kernels\n" ) ;
    fdata.K = malloc( 3*sizeof( struct Kernels ) ) ;
    break ;
  case SYMXY  : fprintf( stdout , "[KAMU] using the SYMXY kernels\n" ) ;
    break ;
  case SYMXY0 : fprintf( stdout , "[KAMU] using the SYMXY0 kernels\n" ) ;
    fdata.K = malloc( sizeof( struct Kernels ) ) ;
    break ;
  }
  if( fdata.shift == true ) {
    fprintf( stdout , "[KAMU] using the shifted kernel without pi-hat\n" ) ;
  } else {
    fprintf( stdout , "[KAMU] using the un-shifted kernel with pi-hat\n" ) ;
  }

  return fdata ;
}

