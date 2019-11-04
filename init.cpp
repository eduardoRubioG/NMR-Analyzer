//  init.cpp
//  nmr
//
//  Created by Eduardo Rubio on 10/17/19.
//  Copyright © 2019 Eduardo Rubio. All rights reserved.
//
//  This file be used to initialize the project. It will be in charge of
//  loading the data and setting all the options 
#include <stdio.h>
#include "include.cpp"
#include <fstream>
#include <sstream>
#include <algorithm>

/**
 * readOptions:
 * filename:The string representation of the options file provided by the user
 */
Options readOptions( const std::string filename ){
    using namespace std;
    Options options;
    std::stringstream ss;
    std::ifstream file( filename.c_str() );
    if( !file.good() ){
        printf("Error opening file \"%s\". Quitting program...\n", filename.c_str());
        exit( 1 );
    }
    std::string buffer;
    
    getline( file, buffer ); options.in_filename = buffer;
    getline( file, buffer ); options.baseline = stod(buffer);
    getline( file, buffer ); options.tol = stod(buffer);
    getline( file, buffer ); options.filterType = stoi(buffer);
    getline( file, buffer ); options.filterSize = stoi(buffer);
    getline( file, buffer ); options.filterPasses = stoi(buffer);
    getline( file, buffer ); options.integrationType = stoi(buffer);
    getline( file, buffer ); options.out_filename = buffer;
    
    return options;
}

/**
 * readData:
 * filename: The string representation of the data file the project will be running
 */
Data readData( const std::string filename ){
    // Initialize variables
    Data data;
    double buffer;
    bool xvalue = true;
    // Open file
    std::ifstream in( filename.c_str(), std::ios::in | std::ios::binary );
    if( !in.good() ){
        printf("Error opening file \"%s\". Quitting program...\n", filename.c_str());
        exit( 1 );
    }
    // Parse the data
    while( in >> buffer ){
        // load the Data struct
        xvalue ? data.x.push_back(buffer) : data.y.push_back(buffer);
        xvalue ? data.n++ : data.n += 0;
        xvalue = !xvalue;
    }
    return data;
}

/**
 * Find the TMS peak and adjust all points accordingly
 */
void peakAdjustment( Data& data ){
    double maxpeak = -__DBL_MAX__;
    for( int i = 0; i < data.n; i++ )
        if( data.y[ i ] > options.baseline && data.x[ i ] > maxpeak )
            maxpeak = data.x[ i ];
    for( int i = 0; i < data.n; i++ )
        data.x[i] = data.x[i] - maxpeak;
    TMSadjustment = maxpeak; 
}
