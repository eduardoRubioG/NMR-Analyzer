//
//  dft.cpp
//  nmr
//
//  Created by Eduardo Rubio on 11/27/19.
//  Copyright © 2019 Eduardo Rubio. All rights reserved.
//
//  This file is where the Dicrete Fourier Transform filtering methods will live 
#include <stdio.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
// #include <math.h> 
#include "include.cpp"
#include "globals.h"

gsl_matrix_complex* invert_matrix_complex( gsl_matrix_complex *matrix, const int& n );
void dft_direct_solver( gsl_matrix_complex* Z, gsl_vector_complex* GC, gsl_vector_complex* Y, const int& n );
void dft_iterative_solver( );
void dft_inverse_solver( gsl_matrix_complex* Z, gsl_vector_complex* C, gsl_vector_complex* Y, const int& n );
void dft_iterative_solver( Data& data, const int& N, gsl_matrix_complex* Z, gsl_vector_complex* GC, gsl_vector_complex* Y, const double& TOL, const int& MAX );
double two_vector_complex_norm( gsl_vector_complex* A, gsl_vector_complex* B, const int& N );

void dft( Data& data ){ 
  /** 
   * Evaluate the following: 
   * c = Z y 
   * 
   * Where: 
   * c: the DFT coefficients 
   * Z: n x n Matrix 
   * y: the y-values that are to be filtered 
   */

  /* Memory allocation for the needed data structures */
  gsl_vector_complex *C  = NULL; 
  gsl_vector_complex *Y  = NULL; 
  gsl_matrix_complex *Z  = NULL; 
  gsl_vector_complex *GC = NULL; 

  C  = gsl_vector_complex_alloc( data.n ); 
  Y  = gsl_vector_complex_alloc( data.n ); 
  GC = gsl_vector_complex_alloc( data.n ); 
  Z  = gsl_matrix_complex_alloc( data.n, data.n ); 

  /* Set the y vector */
  for( int i = 0; i < data.n; i++ ){ 
    gsl_complex val = gsl_complex_rect( data.y[ i ], 0 );
    gsl_vector_complex_set( Y, i, val );
  }

  /* Set the Z matrix  */
  double sqrt_n = sqrt( data.n ); 
  double twoPiDivN = -1 * 2 * M_PI / data.n;
  for( int r = 0; r < data.n; r++ ){ 
    for( int c = 0; c < data.n; c++ ){ 
      double z_val_arg = twoPiDivN * r * c; 
      gsl_complex z_val = gsl_complex_rect( cos( z_val_arg), sin(z_val_arg)); 
      z_val = gsl_complex_div_real( z_val, sqrt_n );
      gsl_matrix_complex_set( Z, r, c, z_val );
    }
  }

  /* Compute c = Z y  and store the DFT coefficients in c */
  gsl_blas_zgemv( CblasNoTrans, GSL_COMPLEX_ONE, Z, Y, GSL_COMPLEX_ZERO, C );

  /* Allocate space for the NMR filter function: c = G c  */
  gsl_matrix_complex *G = NULL; 
  G = gsl_matrix_complex_alloc( data.n, data.n );

  /* Fill out the G matrix */
  double  num = -1 * 4 * log(2); 
  double dnum = pow( data.n, 1.5 );
  for( int r = 0; r < data.n; r++ ){ 
    for( int c = 0; c < data.n; c++ ){ 
      gsl_complex g_val;
      /* Dirac Delta Function */
      if( r == c ){
        double numrc = num * r * c;
        g_val = gsl_complex_rect( exp( numrc / dnum ), 0.0);
      } else { 
        g_val = GSL_COMPLEX_ZERO;
      }
      gsl_matrix_complex_set( G, r, c, g_val );
    }
  }

  /* Compute (gc) = G c */
  gsl_blas_zgemv( CblasNoTrans, GSL_COMPLEX_ONE, G, C, GSL_COMPLEX_ZERO, GC );

  switch (options.filterSize){
    case 0: // inverse recovery 
      dft_inverse_solver( Z, GC, Y, data.n );
      break; 
    case 1: // direct recovery 
      dft_direct_solver( Z, GC, Y, data.n );
      break; 
    case 2: // iterative recovery 
      dft_iterative_solver( data, data.n, Z, GC, Y, options.tol, 30 );
      break; 

    default: 
      printf("Error: DFT filter method was selected but invalid recovery method was requested...");
      exit( 1 );
      break; 
  }

  /* Print out the filtered y values */
  for( int i = 0; i < data.n; i++ ){ 
    data.y[ i ] = GSL_REAL( gsl_vector_complex_get(Y,i) ); 
  }
}

/** 
 * dft_direct_solver: 
 * Uses LU decomposition to solve for the Y values in the form: Z y = (GC)
 * 
 * INPUT:
 *  Z: complex matrix pointer to the Z Fourier matrix 
 *  GC: complex vector pointer to the Fourier coefficients vector
 *  Y: complex vector pointer to the Y values we are solving for 
 *  n: integer of the dimensions needed
 */
void dft_direct_solver( gsl_matrix_complex* Z, gsl_vector_complex* GC, gsl_vector_complex* Y, const int& n ){
  int s; 
  gsl_permutation *p = gsl_permutation_alloc( n ); 
  gsl_linalg_complex_LU_decomp( Z, p, &s ); 
  gsl_linalg_complex_LU_solve( Z, p, GC, Y ); // Solving for y directly in the form Z y = (gc) 
  gsl_permutation_free(p);
}

/** 
 * dft_inverse_solver: 
 * Solve the form y = (GC) Z^-1 
 * 
 * INPUT:
 *  Z: complex matrix pointer to the Z Fourier matrix 
 *  GC: complex vector pointer to the Fourier coefficients vector
 *  Y: complex vector pointer to the Y values we are solving for 
 *  n: integer of the dimensions needed
 */
void dft_inverse_solver( gsl_matrix_complex* Z, gsl_vector_complex* GC, gsl_vector_complex* Y, const int& n ){ 
  gsl_matrix_complex* invZ = gsl_matrix_complex_alloc( n, n );
  gsl_matrix_complex_memcpy( invZ, Z );
  for( int r = 0; r < n; r++ ){
    for( int c = 0; c < n; c++ ){
      gsl_complex CJ = gsl_matrix_complex_get( invZ, r, c );
      CJ = gsl_complex_conjugate( CJ );
      gsl_matrix_complex_set( invZ, r, c, CJ );
    }
  }
  gsl_blas_zgemv( CblasNoTrans, GSL_COMPLEX_ONE, invZ, GC, GSL_COMPLEX_ZERO, Y );
}

double two_vector_complex_norm( gsl_vector_complex* A, gsl_vector_complex* B, const int& N ){ 
  double a_norm, b_norm, norm = 0.0; 
  for( int i = 0; i < N; i++ ){
    a_norm = 0.0; 
    b_norm = 0.0; 
    a_norm = pow(GSL_REAL(gsl_vector_complex_get(A,i)),2) + pow(GSL_IMAG(gsl_vector_complex_get(A,i)),2);
    b_norm = pow(GSL_REAL(gsl_vector_complex_get(B,i)),2) + pow(GSL_IMAG(gsl_vector_complex_get(B,i)),2);
    a_norm = sqrt( a_norm );
    b_norm = sqrt( b_norm );
    norm = std::max( norm, fabs( a_norm - b_norm ));
  }
  return norm; 
}

void dft_iterative_solver( Data& data, const int& N, gsl_matrix_complex* Z, gsl_vector_complex* GC, gsl_vector_complex* Y, const double& TOL, const int& MAX ){ 
  // Step 1 
  bool ok = false; 
  gsl_vector_complex* XO = gsl_vector_complex_calloc( N );
  gsl_vector_complex* x  = gsl_vector_complex_alloc( N );
  gsl_complex sum; 
  int k = 0;
  // Step 2 
  while( k < MAX ){ 
    // Step 3 
    for( int i = 0; i < N; i++ ){ 
      sum = gsl_complex_rect( 0.0, 0.0 );
      for( int j = 0; j < N; j++ )
        if( j != i )
          sum = gsl_complex_add( sum, gsl_complex_mul( gsl_matrix_complex_get(Z,i,j) , gsl_vector_complex_get(XO,j)));
      sum = gsl_complex_mul_real( sum , -1.0 );
      sum = gsl_complex_mul_imag( sum , -1.0 );
      sum = gsl_complex_add( sum, gsl_vector_complex_get(GC,i) );
      sum = gsl_complex_div( sum, gsl_matrix_complex_get(Z,i,i) );
      gsl_vector_complex_set( x, i, sum );
    }
    // Step 4 
    if( two_vector_complex_norm( x, XO, N ) < TOL ){ 
      ok = !ok ;
      gsl_vector_complex_memcpy( x, Y );
      break; 
    }
    // Step 5 
    k++; 
    // Step 6 
    gsl_vector_complex_memcpy(XO,x);
  }
  if( !ok ) {
    printf("Max iterations reached for the iterative method. Defaulting to the inverse method instead...\n");
    dft_inverse_solver( Z, GC, Y, N );
  }
}