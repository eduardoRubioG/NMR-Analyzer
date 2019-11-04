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
#include "globals.cpp"
#include "init.cpp"
#include "output.cpp"
#include "filtering.cpp"
#include "spline.cpp"
#include "roots.cpp"
#include "integration.cpp"
using namespace std;
int main(int argc, const char * argv[]) {
    
    
    // Read in the data from the command line specified file as well as program options from
    // the nmr.in file
    data = readData( argv[1] );
    options = readOptions("nmr.in");
    
    // Adjust the data for the TMS peak
    peakAdjustment(data);
    
    // Filter the data in accordance to the user specified options
    nmrFilter(data, options);
    
    //Export filtered data to "filter.txt"
    dataToFile("filter.txt", data);
    
    // Generate the spline coefficients and export spline curve data to "spline.dat"
    spline(data, splineCo);
    populateSpline(data, splineCo);
    
    // Find all intersections the spline curve has with the baseline and store peaks
    findIntersections(splineCo);
    
    // Integrate the peaks and compute hydrogen counts
    integratePeaks();
    findHydrogens();
    
    // Print out program findings
    outputAnalysis();
    system("cat analysis.txt");
    
    return 0;
}