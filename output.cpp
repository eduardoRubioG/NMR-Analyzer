//
//  output.cpp
//  nmr
//
//  Created by Eduardo Rubio on 10/17/19.
//  Copyright © 2019 Eduardo Rubio. All rights reserved.
//
// This file will manage all processes that handle outputting results and information on the client side

#include <stdio.h>
#include <fstream>

/**
 * Writes out the data points inside the Data object. Useful after filtering original data.
 */
void dataToFile( const std::string outfile, Data data ){
    std::ofstream file( outfile.c_str() );
    const int n = data.n;
    for( int i = 0; i < n; i++ )
        file << data.x[i] << " " << data.y[i] << std::endl;
}

/**
 * Cubic Spline function with the input being the spline coefficients
 */
double f(double x,double a, double b, double c, double d){
    return (a + b*x + c*pow(x,2) + d*pow(x,3));
}

/**
 * INPUT:
 *  int n - number of (x,y) points given
 *  b, c, d - the spline coefficients calculated with the spline function
 */
void populateSpline( const Data& data, const Spline spline ){
    std::ofstream out( "spline.dat" );
    int splinepoints = 15;
    int n = data.n;
    for (int i = 0; i < n-2; i++){
        double delta = (data.x[i+1] - data.x[i])/splinepoints;
        for ( int j = 1; j <= splinepoints; j++ ){
            out
            << data.x[i] + delta*j
            << " "
            << f(delta*j, data.y[i], spline.B[i], spline.C[i], spline.D[i])
            << std::endl;
        }
    }
    out.close();
}

/**
 * The following three templated function are to assist with the formatting of the output file
 */
template<typename T> void printElement(T t, const int& width, const char& fill, std::ofstream& file){
    file << std::left << std::setw(width) << std::setfill(fill) << t;
}

template<typename T> void printSciNum( T t, const int& width, std::ofstream& file ){
    file << std::scientific << std::setw(width) << std::setfill(' ') << t;
}

/**
 * Output the options and content of the program
 */
void outputAnalysis( ) {
    
    std::ofstream file( options.out_filename.c_str() );
    file << std::setprecision (std::numeric_limits<double>::digits10 + 1);
    file << "============================" << std::endl;
    file << "N M R  I N T E R P R E T E R" << std::endl;
    file << "============================\n" << std::endl;
    file << "O P T I O N S" << std::endl;
    file << "==============" << std::endl;
    printElement("Baseline", 20, '.', file);
    printElement(options.baseline, 25, ' ', file);
    file << std::endl;
    printElement("Tolerance", 20, '.', file);
    printElement(options.tol, 25, ' ', file);
    file << std::endl;
    file << std::endl;
    switch( options.filterType ){
        case 0: file << "No Data Filtering" << std::endl; break;
        case 1: file << "Boxcar Filtering" << std::endl; break;
        case 2: file << "Savitzky-Golay Filtering" << std::endl; break;
    }
    file << "Filter Size (if applicable):\t" << options.filterSize << std::endl;
    file << "Filter Passes (if applicable):\t" << options.filterPasses << std::endl;
    file << std::endl;
    file << "I N T E G R A T I O N  M E T H O D" << std::endl;
    file << "==================================" << std::endl;
    switch ( options.integrationType ) {
        case 0: file << "Apative Quadrature" << std::endl; break;
        case 1: file << "Romrberg's Method" << std::endl; break;
        case 2: file << "Gauss-Laguerre Method" << std::endl; break;
        case 3: file << "Composite Simpson's Method" << std::endl; break;
    }
    file << std::endl;
    file << "P L O T  F I L E  D A T A" << std::endl;
    file << "=========================" << std::endl;
    file << "File: " << options.out_filename << std::endl;
    file << "Plot Shifted " << TMSadjustment << " ppm for TMS calibration" << std::endl;
    file << std::endl;

    printElement("Peak", 8, ' ', file);
    printElement("Begin", 25, ' ', file);
    printElement("End", 25, ' ', file);
    printElement("Location", 25, ' ', file);
    printElement("Area", 25, ' ', file);
    printElement("Hydrogens", 8, ' ', file);
    file << std::endl;
    printElement("^", 120, '^', file);
    file << std::endl;
    for( int i = 0; i < peaks.size(); i++ ){
        Peak tmp = peaks[i];
        printElement(i+1, 8, ' ', file);
        printSciNum(tmp.rootA, 25, file);
        printSciNum(tmp.rootB, 25, file);
        printSciNum(tmp.midpoint, 25, file);
        printSciNum(tmp.manifold, 25, file);
        printElement(tmp.hydrogens, 8, ' ', file);
        file << std::endl;
    }
    file << std::fixed << "\nAnalysis took " << ((float)options.time)/CLOCKS_PER_SEC << " seconds\n" << std::endl;
}
