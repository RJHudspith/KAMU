/**
   @file shifted.c
   @brief computes the contraction with the Method 2 kernel
*/
#define _GNU_SOURCE

#include <math.h>
#include "KAMU.h"

static void
contract_LPI_NORMAL( double fval[4] ,
		     const struct Kernels *K ,
		     const double PI[4][4][4][4][4] ,
		     const double PIr[4][4][4][4] ,
		     const double x[4] )
{
  register double sumL0 = 0.0 , sumL1 = 0.0 , sumL2 = 0.0 , sumL3 = 0.0 ;
  const size_t rhomap[6] = { 0 , 0 , 0 , 1 , 1 , 2 } ;
  const size_t sigmap[6] = { 1 , 2 , 3 , 2 , 3 , 3 } ;

  // point out these structures
  const double *L0 = (const double*)K[2].L0 , *L0yx = (const double*)K[1].L0 ;
  const double *L1 = (const double*)K[2].L1 , *L1yx = (const double*)K[1].L1 ;
  const double *L2 = (const double*)K[2].L2 , *L2yx = (const double*)K[1].L2 ;
  const double *L3 = (const double*)K[2].L3 , *L3yx = (const double*)K[1].L3 ;

  size_t rhosig , mu , nu , lda , rho , sig ;
  for( rhosig = 0 ; rhosig < 6 ; rhosig++ ) {
    rho = rhomap[ rhosig ] ; sig = sigmap[ rhosig ] ;

    for( mu = 0 ; mu < 4 ; mu++ ) {
      for( nu = 0 ; nu < 4 ; nu++ ) {
	for( lda = 0 ; lda < 4 ; lda++ ) {

	  // anti-symmetrise the Pi factors
	  const double f1 = ( PI[rho][mu][nu][lda][sig] -
			      PI[sig][mu][nu][lda][rho] ) ;
	  const double f2 = ( x[rho]*PIr[mu][nu][lda][sig] -
			      x[sig]*PIr[mu][nu][lda][rho] ) ; 

	  sumL0 += *L0*f1 + *L0yx*f2 ; L0++ ; L0yx++ ;
	  sumL1 += *L1*f1 + *L1yx*f2 ; L1++ ; L1yx++ ;
	  sumL2 += *L2*f1 + *L2yx*f2 ; L2++ ; L2yx++ ;
	  sumL3 += *L3*f1 + *L3yx*f2 ; L3++ ; L3yx++ ;
	}
      }
    }
  }

  // set fval with factor of 0.5 (same as ipihatFermLoop_antisym)
  // obviously this factor comes from the antisymmetrisation
  fval[0] = 0.5*sumL0 ; fval[1] = 0.5*sumL1 ;
  fval[2] = 0.5*sumL2 ; fval[3] = 0.5*sumL3 ;

  return ;
}

static void
contract_LPI_SYMXY0( double fval[4] ,
		     const struct Kernels *K ,
		     const double PI[4][4][4][4][4] ,
		     const double PIr[4][4][4][4] ,
		     const double x[4] )
{
  register double sumL0 = 0.0 , sumL1 = 0.0 , sumL2 = 0.0 , sumL3 = 0.0 ;
  const size_t rhomap[6] = { 0 , 0 , 0 , 1 , 1 , 2 } ;
  const size_t sigmap[6] = { 1 , 2 , 3 , 2 , 3 , 3 } ;

  // point out these structures
  const double *L0 = (const double*)K[2].L0 ;
  const double *L1 = (const double*)K[2].L1 ;
  const double *L2 = (const double*)K[2].L2 ;
  const double *L3 = (const double*)K[2].L3 ;

  size_t rhosig , mu , nu , lda , rho , sig ;
  for( rhosig = 0 ; rhosig < 6 ; rhosig++ ) {
    rho = rhomap[ rhosig ] ; sig = sigmap[ rhosig ] ;

    for( mu = 0 ; mu < 4 ; mu++ ) {
      for( nu = 0 ; nu < 4 ; nu++ ) {
	for( lda = 0 ; lda < 4 ; lda++ ) {

	  // anti-symmetrise the Pi factors
	  const double f1 = ( PI[rho][mu][nu][lda][sig] -
			      PI[sig][mu][nu][lda][rho] ) ;
	  const double f2 = ( x[rho]*PIr[mu][nu][lda][sig] -
			      x[sig]*PIr[mu][nu][lda][rho] ) ; 

	  sumL0 += *L0*(3*f1 - f2) ; L0++ ; 
	  sumL1 += *L1*(3*f1 - f2) ; L1++ ; 
	  sumL2 += *L2*(3*f1 - f2) ; L2++ ; 
	  sumL3 += *L3*(3*f1 - f2) ; L3++ ; 
	}
      }
    }
  }

  // set fval with factor of 0.5 (same as ipihatFermLoop_antisym)
  // obviously this factor comes from the antisymmetrisation
  fval[0] = 0.5*sumL0 ; fval[1] = 0.5*sumL1 ;
  fval[2] = 0.5*sumL2 ; fval[3] = 0.5*sumL3 ;

  return ;
}

// this is the integrand
int
amu_integral_shift( unsigned int ndim ,
		    const double *x ,
		    void *fdata ,
		    unsigned int fdim ,
		    double *fval )
{
  struct integral_args *F = (struct integral_args*)fdata ;

    // call to sincos is cheaper than calling sin or cos individually
  double si , ci ;
  sincos( x[0] , &si , &ci ) ;

  const double xn = x[1] ;
  const double yn = ( ndim == 2 ) ? F -> tmp : x[2] ;
  const double xv[4] = { xn*ci , xn*si , 0 , 0 } ;
  const double yv[4] = { yn , 0 , 0 , 0 } ;

  const double xMv[4] = { xv[0]*F->Mv , xv[1]*F->Mv ,
			  xv[2]*F->Mv , xv[3]*F->Mv } ;
  const double yMv[4] = { yv[0]*F->Mv , yv[1]*F->Mv ,
			  yv[2]*F->Mv , yv[3]*F->Mv } ;
  pihat1( xMv , yMv , F->t , F->vpihat1 , F->vpihatr1 ) ;

  const double M[4] = { 0.0 , 0.4 , 0.8 , 1.0 } ;

  switch( F->Sym ) {
  case NORMAL :
    compute_con_kernels_v2( xv , yv , F->t , F->K ) ;
    contract_LPI_NORMAL( fval , F->K , F->vpihat1 , F->vpihatr1 , xMv ) ;
    break ;
  case MKERN :
    compute_con_kernelsM_L2( M , xv , yv , F->t , F->K ) ;
    contract_LPI_NORMAL( fval , F->K , F->vpihat1 , F->vpihatr1 , xMv ) ;
    break ;
  case SYMXY :
    fprintf( stderr , "case SYMXY not implemented" ) ; exit(1) ;
    break ;
  case SYMXY0 :
    compute_all_kernels_SYMXY0_v2( xv , yv , F->t , F->K ) ;
    contract_LPI_SYMXY0( fval , F->K , F->vpihat1 , F->vpihatr1 , xMv ) ;
    break ;
  }

  const double prefac = xn*xn*xn*yn*yn*yn*si*si ;
  fval[0] *= prefac ; fval[1] *= prefac ;
  fval[2] *= prefac ; fval[3] *= prefac ;

  return 0 ; // success
}
