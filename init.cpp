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

Spline splineCo;
Data _data;
Options options;
std::vector<Peak> peaks;
double TMSadjustment;

/**
 * Essentially a stoi / stod function since Mercer won't update to c++11 :'((
 */
template<typename T> T stringToType( const std::string& str ){
    T toReturn;
    std::stringstream ss( str );
    ss >> toReturn;
    return toReturn;
}

/**
 * readOptions:
 * filename:The string representation of the options file provided by the user
 */
Options readOptions( const std::string filename ){
    using namespace std;
    Options options;
    std::ifstream file( filename.c_str() );
    if( !file.good() ){
        printf("Error opening file \"%s\". Quitting program...\n", filename.c_str());
        exit( 1 );
    }
    
    std::string buffer;
    getline( file, buffer ); options.in_filename = buffer;
    getline( file, buffer ); options.baseline = stringToType<double>(buffer);
    getline( file, buffer ); options.tol = stringToType<double>(buffer);
    getline( file, buffer ); options.filterType = stringToType<int>(buffer);
    getline( file, buffer ); options.filterSize = stringToType<int>(buffer);
    getline( file, buffer ); options.filterPasses = stringToType<int>(buffer);
    getline( file, buffer ); options.integrationType = stringToType<int>(buffer);
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

/**
 * init:
 * Encapsulates all the initializing functions to be called from main
 */
void init( ){
    options = readOptions( "nmr.in" );
    _data = readData( options.in_filename );
}
