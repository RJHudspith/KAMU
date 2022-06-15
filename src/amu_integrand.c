/**
   Written by Jamie based on the amu_integral code by En-Hung

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
  double y = 0.0 , yprev[4] = {} , nint[4] = {} , nerr[4] = {} , xprev = 0. ;

  for(int i = 0 ; i < NKERNEL ; i ++) {
    fprintf( stdout , "[KAMU] amu_L%d %f = %e +/- %e\n" ,
	     i , y , 0. ,  0. ) ;
  }

  // fine-grained scan
  const double fine_max = 0.5 ;
  for( y = 0.05 ; y < fine_max ; y+=0.02 ) {
    fdata.tmp = y ;
    int flag ;
    if( ( flag = hcubature( NKERNEL, integrand, &fdata, 2,
			    fdata.xmin, fdata.xmax, 0, 0, fdata.Tol,
			    ERROR_L2 , valv1, errv1 ) ) != 0 ) {
      fprintf( stderr , "[KAMU] cubature failed %d\n" , flag ) ;
      goto memfree ;
    }
    const double f = prefac*(y-xprev)/2. ;
    for(int i = 0 ; i < NKERNEL ; i ++) {
      fprintf( stdout , "[KAMU] amu_L%d %f = %e +/- %e\n" ,
	       i , y , prefac*valv1[i] ,  prefac*errv1[i] ) ;
      nint[ i ] += f*(valv1[i] + yprev[i]) ;
      nerr[ i ] = sqrt( nerr[i]*nerr[i] + errv1[i]*errv1[i]) ;
      yprev[i] = valv1[i] ;
      fprintf( stdout , "[KAMU] Nint_L%d %f = %e +/- %e\n" ,
	       i , y , nint[i] , f*nerr[i] ) ;
    }
    xprev = y ;
  }
  // coarse-grained scan
  for( y = fine_max ; y < 5.5 ; y+=0.2 ) {
    fdata.tmp = y ;
    int flag ;
    if( ( flag = hcubature( NKERNEL, integrand, &fdata, 2,
			    fdata.xmin, fdata.xmax, 0, 0, fdata.Tol,
			    ERROR_L2 , valv1, errv1 ) ) != 0 ) {
      fprintf( stderr , "[KAMU] cubature failed %d\n" , flag ) ;
      goto memfree ;
    }
    const double f = prefac*(y-xprev)/2. ;
    for(int i = 0 ; i < NKERNEL ; i ++) {
      fprintf( stdout , "[KAMU] amu_L%d %f = %e +/- %e\n" ,
	       i , y , prefac*valv1[i] ,  prefac*errv1[i] ) ;
      nint[ i ] += ( prefac*(valv1[i] + yprev[i]))*( y-xprev )/2. ;
      nerr[ i ] = sqrt( nerr[i]*nerr[i] + errv1[i]*errv1[i]) ;
      yprev[i] = valv1[i] ;
      fprintf( stdout , "[KAMU] Nint_L%d %f = %e +/- %e\n" ,
	       i , y , nint[i] , f*nerr[i] ) ;
    }
    xprev = y ;
  }
  
 memfree :
  print_time() ;

  if( fdata.K != NULL ) {
    free( fdata.K ) ;
  }

  free_QED_temps( &fdata.t ) ;

  return 0 ;
}
