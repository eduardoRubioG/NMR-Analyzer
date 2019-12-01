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
#include <math.h> 
#include "include.cpp"
#include "globals.h"

gsl_matrix_complex* invert_matrix_complex( gsl_matrix_complex *matrix, const int& n );
void dft_direct_solver( gsl_matrix_complex* Z, gsl_vector_complex* GC, gsl_vector_complex* Y, const int& n );
void dft_iterative_solver( );
void dft_inverse_solver( gsl_matrix_complex* Z, gsl_vector_complex* C, gsl_vector_complex* Y, const int& n );

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
      printf("BEFORE\n");
      for( int i = 0; i < 10; i++ ) 
        printf("%g + i%g\n",GSL_REAL(gsl_vector_complex_get( Y, i )), GSL_IMAG(gsl_vector_complex_get( Y, i )));
      printf("okok\n");
      dft_inverse_solver( Z, C, Y, data.n );
      printf("AFTER\n");
      for( int i = 0; i < 10; i++ )
        printf("%g + i%g\n",GSL_REAL(gsl_vector_complex_get( Y, i )), GSL_IMAG(gsl_vector_complex_get( Y, i )));
      break; 
    case 1: // direct recovery 
      dft_direct_solver( Z, GC, Y, data.n );
      break; 
    case 2: // iterative recovery 
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

void dft_inverse_solver( gsl_matrix_complex* Z, gsl_vector_complex* C, gsl_vector_complex* Y, const int& n ){ 
  printf( "Doing the inverse\n");
  gsl_matrix_complex* invZ = gsl_matrix_complex_alloc( n, n ); 
  invZ = invert_matrix_complex( Z, n ); // Defined below
  gsl_blas_zgemv( CblasNoTrans, GSL_COMPLEX_ONE, invZ, C, GSL_COMPLEX_ZERO, Y );
}

/**
 * invert_matrix_complex: 
 * INPUT: 
 *   matrix: a complex matrix pointer that you wish to invert 
 *        n: int of the dimensions needed 
 * OUTPUT: 
 *  Returns a complex matrix pointer that is the inverse of the matrix argument 
 */
gsl_matrix_complex* invert_matrix_complex(gsl_matrix_complex* matrix, const int& n ){  

/** 
 * DON"T DO IT LIKE THIS
 * 
 * Find the complex conjugate of all Z(i,j) instead of doing an LU factorization + decomposition
 */


  // Compute the LU decomposition of this matrix
  int s;
  gsl_permutation *p = gsl_permutation_alloc( n );
  gsl_matrix_complex *tmp = gsl_matrix_complex_alloc( n, n ); 
  gsl_matrix_complex_memcpy( tmp, matrix );
  gsl_linalg_complex_LU_decomp(tmp, p, &s);

  // Compute the  inverse using the LU decomposition
  gsl_matrix_complex *inv = gsl_matrix_complex_alloc(n, n);
  gsl_linalg_complex_LU_invert(tmp, p, inv);
  gsl_permutation_free(p);

  return inv;
}
