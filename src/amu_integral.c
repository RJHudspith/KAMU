/**
   Written by En-Hung Chao
   
   Modified by Jamie for c

   uses the hcubature routine from cubature see https://github.com/stevengj/cubature
 */
#include "KAMU.h"

#include "cmdline.h"
#include "integrands.h"
#include "shifted.h"

// main function
int
main( const int argc ,
      const char *argv[] )
{
  // check the command-line arguments
  if( argc != 6 ) {
    fprintf( stderr , "[KAMU] Usage :: ./KAMU Mv Tol Sym Ymin shift\n" ) ;
    fprintf( stderr , "[KAMU] Mv :: ratio of lepton loop to muon mass\n" ) ;
    fprintf( stderr , "[KAMU] Tol :: Tolerance of hcubature\n" ) ;
    fprintf( stderr , "[KAMU] Sym :: NORMAL/SYMXY/SYMXY0\n" ) ;
    fprintf( stderr , "[KAMU] Ymin :: Minimum x or y distance\n" ) ;
    fprintf( stderr , "[KAMU] shift :: true/false use the shifted kernel\n" ) ;
    return 1 ;
  }
  
  int err ;
  struct integral_args fdata = read_cmdline( &err , argv ) ;
  // finally do the integration
  const int NKERNEL = 4 ;
  double valv1[NKERNEL] , errv1[NKERNEL] ;

  int (*integrand)( unsigned int ndim ,
		    const double *x ,
		    void *fdata ,
		    unsigned int fdim ,
		    double *fval ) ;
  if( fdata.shift == true ) {
    integrand = amu_integral_shift ;
  } else {
    integrand = amu_integral ;
  }
  
  if( err ) goto memfree ;
  
  const double prefac =
    2*pow(fdata.Mv,7)*pow(8.0*M_PI*M_PI*AlfQED,3)/3;

  start_timer() ;
  
  // does the integral
  int flag ;
  if( ( flag = hcubature( NKERNEL, integrand, &fdata, 3,
			  fdata.xmin, fdata.xmax, 0, 0, fdata.Tol,
			  ERROR_INDIVIDUAL, valv1, errv1 ) ) != 0 ) {
    fprintf( stderr , "[KAMU] cubature failed %d\n" , flag ) ;
    goto memfree ;
  }

  for(int i = 0 ; i < NKERNEL ; i ++) {
    fprintf( stdout , "[KAMU] amu (L%d) = %e +/- %e\n" ,
	     i , prefac*valv1[i] ,  prefac*errv1[i] ) ;
  }

 memfree :
  print_time() ;

  if( fdata.K != NULL ) {
    free( fdata.K ) ;
  }

  free_QED_temps( &fdata.t ) ;

  return 0 ;
}
