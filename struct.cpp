//
//  struct.cpp
//  nmr
//
//  Created by Eduardo Rubio on 10/17/19.
//  Copyright © 2019 Eduardo Rubio. All rights reserved.
//

#include <stdio.h>
#include <vector> 

/**
 * Data Struct:
 * Will store the original data given by user
 * x: a vector of doubles that contain x-values
 * y: a vector of doubles that contain y-values
 * n: an int representing the how many x-y value pairs exist in the struct
 */
struct Data {
    std::vector<double> x;
    std::vector<double> y;
    int n = 0; 
};

/**
 * Spline Struct:
 * This struct contains the cubic spline interpolant coefficients as well as a size overhaul
 * for the vectors.
 */
struct Spline {
    std::vector<double> B;
    std::vector<double> C;
    std::vector<double> D;
    
    void resize( int n ){
        B.resize(n);
        C.resize(n);
        D.resize(n);
    }
};

/**
 * Struct options:
 * in_filename: the original data file name
 * out_filename: the file the program will store the results in
 * baseline: double representing the baseline for the NMR interpretation
 * double tol: double representing the tolerance level for calculations
 * filterType: int deciding the filter type
 * filterSize: int representing the filter size
 * filterPasses: int representing the amount of filter passes
 * integrationType: int deciding which integration method will be applied
 */
struct Options {
    std::string in_filename;
    std::string out_filename;
    double baseline;
    double tol;
    int filterType;
    int filterSize;
    int filterPasses;
    int integrationType;
    clock_t time; 
};

/**
 * The Peak struct will serve to store information on the various peaks found in the data
 * indexA: the index of the first spline function
 * indexB: the index of the second spline function
 * rootA: The first root (x-value)
 * rootB: the second root (x-value)
 * midpoint: the midpoint between the two roots
 * manifold: the area of the peak
 * isComplete: whether or not enough degrees of freedom are satisfied for the peak
 */
struct Peak {
    int indexA;
    int indexB;
    double rootA;
    double rootB;
    double midpoint;
    double manifold;
    bool isCompete = false;
    int hydrogens; 
    
    void print( ){
        printf("idxA: %d\nidxB: %d\n rootA: %f\nrootB: %f\n", indexA, indexB, rootA, rootB);
    }
    
};
