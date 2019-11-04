//
//  roots.cpp
//  nmr
//
//  Created by Eduardo Rubio on 10/21/19.
//  Copyright © 2019 Eduardo Rubio. All rights reserved.
//

#include <stdio.h>
/**
 * Bisection Method:
 *   a:   initial endpoint
 *   b:   ending endpoint
 * tol: tolerance level for finding roots
 *   n:   max number of iterations
 * idx:   integer to determine which spline interpolant to use
 */
double bisection( double a, double b, double tol, int n, Spline spline, int idx ) {
    int i = 1;
    double p, FP;
    double FA = f( a-data.x[ idx ],
                    data.y[ idx ],
                   spline.B[ idx ],
                   spline.C[ idx ],
                   spline.D[ idx ] );
    FA = FA -  options.baseline;

    while( i <= n ){
        p = a + (b-a) / 2;
        FP = f(p-data.x[ idx ],
               data.y[ idx ],
               spline.B[ idx ],
               spline.C[ idx ],
               spline.D[ idx ] );
        FP = FP -  options.baseline;
        
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
    return NULL;
}

/**
 * Runs through the filtered data and find all the spline intersections with the user designated baseline
 */
void findIntersections( Spline spline ){
    std::ofstream file( "roots.dat" );
    double adjustedA, adjustedB;
    bool newPeak = true;
    Peak peak;
    for( int i = 0; i < data.n-1; i++ ){
        adjustedA = data.y[ i ] - options.baseline;
        adjustedB = data.y[ i + 1 ] - options.baseline;
        if( adjustedA * adjustedB < 0 ){
            if( newPeak ){
                peak.rootA = bisection(data.x[ i ], data.x[ i + 1 ], options.tol, 100, spline, i );
                peak.indexA = i;
            } else {
                peak.rootB = bisection(data.x[ i ], data.x[ i + 1 ], options.tol, 100, spline, i );
                peak.indexB = i;
                peak.isCompete = true;
                peak.midpoint = (peak.rootB + peak.rootA)/2; 
                peaks.push_back( peak );
            }
            newPeak = !newPeak;
        }
    }
    for( int i = 0; i < peaks.size(); i++ ){
        file << peaks[i].rootA << ' ' << options.baseline << std::endl;
        file << peaks[i].rootB << ' ' << options.baseline << std::endl;
//        printf("PEAK %d\n", i+1);
//        peaks[i].print();
//        printf("\n\n\n");
    }
}
