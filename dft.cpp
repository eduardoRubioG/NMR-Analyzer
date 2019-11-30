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

  FILE* C_FILE; 
  FILE* Z_FILE;
  FILE* G_FILE; 
  FILE* Y_FILE; 
  FILE* SOL; 

  C_FILE = fopen( "c_file.dat", "w"); 
  Z_FILE = fopen( "z_file.dat", "w");
  G_FILE = fopen( "g_file.dat", "w"); 
  Y_FILE = fopen( "y_file.dat", "w"); 
  SOL = fopen( "sol.dat", "w" );

  /* Memory allocation for the needed data structures */
  gsl_vector_complex *c = NULL; 
  gsl_vector_complex *y = NULL; 
  gsl_matrix_complex *Z = NULL; 

  c = gsl_vector_complex_alloc( data.n ); 
  y = gsl_vector_complex_alloc( data.n ); 
  Z = gsl_matrix_complex_alloc( data.n, data.n ); 

  /* Set the y vector */
  for( int i = 0; i < data.n; i++ ){ 
    gsl_complex val = gsl_complex_rect( data.y[ i ], 0 );
    gsl_vector_complex_set( y, i, val );
  }

  /* Set the Z matrix  */
  double sqrt_n = sqrt( data.n ); 
  for( int r = 0; r < data.n; r++ ){ 
    for( int c = 0; c < data.n; c++ ){ 
      double z_val_arg = ( -1 * 2 * M_PI * r * c ) / data.n; 
      gsl_complex z_val = gsl_complex_rect( cos( z_val_arg), sin(z_val_arg)); 
      z_val = gsl_complex_div_real( z_val, sqrt_n );
      // printf("%d,%d: %g + i%g\n", r,c, GSL_REAL(z_val),GSL_IMAG(z_val));
      gsl_matrix_complex_set( Z, r, c, z_val );
    }
  }

  /* Compute c = Z y  and store the DFT coefficients in c */
  gsl_complex alpha = gsl_complex_rect( 1.0, 1.0 );
  gsl_complex  beta = gsl_complex_rect( 0.0, 0.0 );
  gsl_blas_zgemv( CblasNoTrans, alpha, Z, y, beta, c );

  /* Allocate space for the NMR filter function: c = G c  */
  gsl_matrix_complex *G = NULL; 
  G = gsl_matrix_complex_alloc( data.n, data.n );
  for( int r = 0; r < data.n; r++ ){ 
    for( int c = 0; c < data.n; c++ ){ 
      gsl_complex g_val;
      if( r == c ){
        double  num = -1 * 4 * log(2) * r * c; 
        double dnum = pow( data.n, 1.5 );
        g_val = gsl_complex_rect( exp( num / dnum ), 0.0);
      } else { 
        g_val = gsl_complex_rect( 0.0, 0.0 );
      }
      gsl_matrix_complex_set( G, r, c, g_val );
    }
  }

  gsl_blas_zgemv( CblasNoTrans, alpha, G, c, beta, c );

  /* Print out all variables into respective files  */
  gsl_vector_complex_fprintf(C_FILE,c,"%g");
  gsl_vector_complex_fprintf(Y_FILE,y,"%g");
  gsl_matrix_complex_fprintf(G_FILE,G,"%g");
  gsl_matrix_complex_fprintf(Z_FILE,Z,"%g");

  /* Ready to solve  */
  int s; 
  gsl_permutation *p = gsl_permutation_alloc( data.n ); 
  gsl_linalg_complex_LU_decomp( Z, p, &s ); 
  gsl_linalg_complex_LU_solve( Z, p, c, y ); // Solving for y directly in the form Z y = c 

  gsl_vector_complex_fprintf(SOL,y,"%g");
}
