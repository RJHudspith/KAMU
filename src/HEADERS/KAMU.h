#ifndef KAMU_H
#define KAMU_H

#include "KQED.h"

#include "hcubature.h"

typedef enum { NORMAL = 0 , SYMXY = 1 , SYMXY0 = 2 } symmetrise ; 

struct integral_args{
  double Mv , Tol ;
  double xmin[3] ; 
  double xmax[3] ;
  struct QED_kernel_temps t ;
  struct Kernels *K ;
  double vpihat1[4][4][4][4][4] ;
  double vpihatr1[4][4][4][4] ;
  symmetrise Sym ;
  double tmp ;
  bool shift ;
} ;

#endif
