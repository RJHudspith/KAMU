/**
   @file integrands.c
   @brief typical QED kernel integrands
 */
#define _GNU_SOURCE // sincos

#include <math.h>
#include "KAMU.h"

// unshifted integrand
int
amu_integral( unsigned int ndim ,
	      const double *x ,
	      void *fdata ,
	      unsigned int fdim ,
	      double *fval )
{
  const struct integral_args *F = (const struct integral_args*)fdata ;

  // call to sincos is cheaper than calling sin or cos individually
  double si , ci ;
  sincos( x[0] , &si , &ci ) ;

  const double xn = x[1] ;
  const double yn = ( ndim == 2 ) ? F -> tmp : x[2] ;
  const double xv[4] = { xn*ci , xn*si , 0 , 0 } ;
  const double yv[4] = { yn , 0 , 0 , 0 } ;

  struct QED_Kernels K ;
  const double M[4] = { 0.0 , 0.4 , 0.8 , 1.0 } ;

  switch( F -> Sym ) {
  case NORMAL :
    compute_all_kernels( xv , yv , F->t , &K) ;
    break ;
  case MKERN :
    compute_all_Mkernels( M , xv , yv , F->t , &K ) ;
    break ;
  case SYMXY :
    compute_all_kernels_SYMXY( xv , yv , F->t , &K ) ;
    break ;
  case SYMXY0 :
    compute_all_kernels_SYMXY0( xv , yv , F->t , &K ) ;
    break ;
  }

  const double xMv[4] = { xv[0]*F->Mv , xv[1]*F->Mv ,
			  xv[2]*F->Mv , xv[3]*F->Mv } ;
  const double yMv[4] = { yv[0]*F->Mv , yv[1]*F->Mv ,
			  yv[2]*F->Mv , yv[3]*F->Mv } ;
  ipihatFermLoop_antisym( xMv , yMv , F->t , K.L0.yx ) ;

  const double prefac = xn*xn*xn*yn*yn*yn*si*si ;
#if (defined __AVX__)
  const __m256d *K0 = (const __m256d*)K.L0.xy , *K1 = (const __m256d*)K.L1.xy ;
  const __m256d *K2 = (const __m256d*)K.L2.xy , *K3 = (const __m256d*)K.L3.xy ;
  const __m256d *PI = (const __m256d*)K.L0.yx ;
  register __m256d sumL0 = _mm256_setzero_pd() ;
  register __m256d sumL1 = _mm256_setzero_pd() ;
  register __m256d sumL2 = _mm256_setzero_pd() ;
  register __m256d sumL3 = _mm256_setzero_pd() ;
  for( size_t i = 0 ; i < 96 ; i++ ) {
    register const __m256d A = _mm256_mul_pd( *PI , *K0 ) ; K0++ ;
    register const __m256d B = _mm256_mul_pd( *PI , *K1 ) ; K1++ ;
    register const __m256d C = _mm256_mul_pd( *PI , *K2 ) ; K2++ ;
    register const __m256d D = _mm256_mul_pd( *PI , *K3 ) ; K3++ ;
    PI++ ;
    sumL0 = _mm256_add_pd( sumL0 , A ) ;
    sumL1 = _mm256_add_pd( sumL1 , B ) ;
    sumL2 = _mm256_add_pd( sumL2 , C ) ;
    sumL3 = _mm256_add_pd( sumL3 , D ) ;
  }
  const __m256d PF = _mm256_broadcast_sd( &prefac ) ;
  sumL0 = _mm256_mul_pd( sumL0 , PF ) ;
  sumL1 = _mm256_mul_pd( sumL1 , PF ) ;
  sumL2 = _mm256_mul_pd( sumL2 , PF ) ;
  sumL3 = _mm256_mul_pd( sumL3 , PF ) ;

  double t[4] = { 0 } ;
  _mm256_store_pd( t , sumL0 ) ; fval[0] = t[0]+t[1]+t[2]+t[3] ;
  _mm256_store_pd( t , sumL1 ) ; fval[1] = t[0]+t[1]+t[2]+t[3] ;
  _mm256_store_pd( t , sumL2 ) ; fval[2] = t[0]+t[1]+t[2]+t[3] ;
  _mm256_store_pd( t , sumL3 ) ; fval[3] = t[0]+t[1]+t[2]+t[3] ;
#else
  const double *K0 = (const double*)K.L0.xy , *K1 = (const double*)K.L1.xy ;
  const double *K2 = (const double*)K.L2.xy , *K3 = (const double*)K.L3.xy ;
  const double *PI = (const double*)K.L0.yx ;
  register double sumL0 = 0.0 , sumL1 = 0.0 , sumL2 = 0.0 , sumL3 = 0.0 ;
  for( size_t i = 0 ; i < 384 ; i++ ) {
    sumL0 += *K0 * (*PI) ; K0++ ;
    sumL1 += *K1 * (*PI) ; K1++ ;
    sumL2 += *K2 * (*PI) ; K2++ ;
    sumL3 += *K3 * (*PI) ; K3++ ;
    PI++ ;
  }
  fval[0] = prefac*sumL0 ;
  fval[1] = prefac*sumL1 ;
  fval[2] = prefac*sumL2 ;
  fval[3] = prefac*sumL3 ;
#endif

  return 0 ; // success
}
