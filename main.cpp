//
//  main.cpp
//  nmr
//
//  Created by Eduardo Rubio on 10/17/19.
//  Copyright © 2019 Eduardo Rubio. All rights reserved.
//

#include <iostream>
#include "struct.cpp"
#include "include.cpp"
#include "globals.h"
#include "init.cpp"
#include "output.cpp"
#include "filtering.cpp"
#include "spline.cpp"
#include "roots.cpp"
#include "integration.cpp"
//using namespace std;
int main(int argc, const char * argv[]) {
    
    system( "clear" );
    
    // Read in the data from the command line specified file as well as program options from
    // the nmr.in file
    init();
    
    // Note time in order to calculate analysis time
    options.time = clock( );
    
    // Adjust the data for the TMS peak
    peakAdjustment(_data); //good
    
    // Filter the data in accordance to the user specified options
    nmrFilter(_data, options); //good
    
    //Export filtered data to "filter.txt"
    dataToFile("filter.txt", _data); //good
    
    // Generate the spline coefficients and export spline curve data to "spline.dat"
    spline(_data, splineCo);
    populateSpline(_data, splineCo); // good
    
    // Find all intersections the spline curve has with the baseline and store peaks
    findIntersections(splineCo);// good 
    
    // Integrate the peaks and compute hydrogen counts
    integratePeaks(); // good ?
    findHydrogens(); // good ? 
    
    // Print out program findings
    options.time = clock() - options.time;
    outputAnalysis();
    system("cat analysis.txt");
    
    return 0;
}
