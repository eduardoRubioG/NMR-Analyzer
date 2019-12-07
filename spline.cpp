//
//  spline.cpp
//  nmr
//
//  Created by Eduardo Rubio on 10/17/19.
//  Copyright © 2019 Eduardo Rubio. All rights reserved.
//

#include <stdio.h>

/**
 * INPUT:
 *  int n - Number of points
 *  vector<double> x - all x points
 *  vector<double> y - all y points (referred to as all a points in the book)
 */
void spline( const Data& data, Spline& spline ){
    int n = data.n;
    //double h[ n ];
    //double alpha[ n ];
    //double l[ n ];
    //double mu[ n ];
    //double z[ n ];
    double* h = new double[ n ];
    double* alpha = new double[ n ];
    double* l = new double[ n ];
    double* mu = new double[ n ];
    double* z = new double[ n ];
    spline.resize( n );
    
    //Step one:
    for( int i = 0; i < n; i++ ){
        h[ i ] = data.x[ i + 1 ] - data.x[ i ];
    }
    
    //Step two:
    for( int i = 1; i < n; i++ ){
        alpha[ i ] = (3/h[i])*(data.y[i+1]-data.y[i]) - (3/h[i-1])*(data.y[i]-data.y[i-1]);
    }
    
    // Step 3, 4, 5 and part of 6 solves a tridiagonal linear system
    
    //Step three:
    l[ 0 ] = 1; mu[ 0 ] = 0; z[ 0 ] = 0;
    
    //Step four:
    for( int i = 1; i < n; i++ ){
        l[ i ] = 2*(data.x[i+1]-data.x[i-1]) - h[i-1]*mu[i-1];
        mu[ i ] = h[ i ] / l[ i ];
        z[ i ] = (alpha[ i ] - h[i - 1]*z[ i - 1 ]) / l[ i ];
    }
    
    //Step five:
    l[ n-1 ] = 1; 
    z[ n-1 ] = 0; 
    spline.C[ n-1 ] = 0;
    
    //Step six:
    for( int j = (n-2); j >= 0; j-- ){
        spline.C[ j ] = z[ j ] - mu[ j ]*spline.C[ j + 1 ];
        spline.B[ j ] = (data.y[ j + 1] - data.y[ j ])/h[ j ] - h[ j ]*( spline.C[ j + 1 ] + 2 * spline.C[ j ]) / 3;
        spline.D[ j ] = (spline.C[ j + 1 ] - spline.C[ j ]) / (3 * h[ j ]);
    }

    //Clean up 
    delete[] h;
    delete[] alpha;
    delete[] l;
    delete[] mu;
    delete[] z;
}
