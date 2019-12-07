//
//  roots.cpp
//  nmr
//
//  Created by Eduardo Rubio on 10/21/19.
//  Copyright © 2019 Eduardo Rubio. All rights reserved.
//

#include <stdio.h>
#include "globals.h"
/**
 * Bisection Method:
 *   a:   initial endpoint
 *   b:   ending endpoint
 * tol: tolerance level for finding roots
 *   n:   max number of iterations
 * idx:   integer to determine which spline interpolant to use
 */
double bisection( double a, double b, double tol, int n, const Spline spline, int idx ) {
    int i = 1;
    double p, FP;
    double FA = f( a -_data.x[ idx ],
                  _data.y[ idx ] - options.baseline,
                  spline.B[ idx ],
                  spline.C[ idx ],
                  spline.D[ idx ] );
    
    while( i <= n ){
        p = a + ((b-a)/2);
        FP = f(p - _data.x[ idx ],
               _data.y[ idx ] - options.baseline,
               spline.B[ idx ],
               spline.C[ idx ],
               spline.D[ idx ] );
        
        if( FP == 0 || fabs((b-a)/2) < tol ) return p;
        i++;
        if( FA * FP > 0 ){
            a = p;
            FA = FP;
        } else {
            b = p;
        }
    }
    printf("Error: No root found in index %d\n", idx );
    return 0;
}

/**
 * Runs through the filtered data and find all the spline intersections with the user designated baseline
 */
void findIntersections( const Spline spline ){
    std::ofstream file( "roots.dat" );
    double adjustedA, adjustedB;
    bool newPeak = true;
    Peak peak;
    for( int i = 0; i < _data.n-1; i++ ){
        adjustedA = _data.y[ i ] - options.baseline;
        adjustedB = _data.y[ i + 1 ] - options.baseline;
        if( adjustedA * adjustedB < 0 ){
            if( newPeak ){
                peak.rootA = bisection(_data.x[ i ], _data.x[ i + 1 ], options.tol, 100, spline, i );
                peak.indexA = i;
            } else {
                peak.rootB = bisection(_data.x[ i ], _data.x[ i + 1 ], options.tol, 100, spline, i );
                peak.indexB = i;
                peak.isComplete = true;
                peak.midpoint = (peak.rootB + peak.rootA)/2;
                peaks.push_back( peak );
            }
            newPeak = !newPeak;
        }
    }
}
