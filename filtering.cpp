//
//  filtering.cpp
//  nmr
//
//  Created by Eduardo Rubio on 10/17/19.
//  Copyright © 2019 Eduardo Rubio. All rights reserved.
//
//  This file is where the filtering functions will live
#include <stdio.h>
#include "dft.cpp"
/**
 * Implementation of the Boxcar data filter. Pass in the Data object as well as the filter size
 */
void boxcar( Data& data, int filterSize ){
    if( filterSize > 0 ){
        int half = filterSize/2;
        double sum = 0;
        // Prep the data
        for( int i = 0; i < half; i++ )
            data.y.push_back(data.y[i]);
        for( int i = 0; i < half; i++ )
            data.y.insert(data.y.begin(), data.y[ data.y.size() - half - i - 1]);
        // Begin the boxcar filter
        std::vector<double> tmp = data.y;
        for( int i = half; i < (data.y.size() - half); i++ ){
            sum = 0.0;
            for( int j = i - half; j <= ( i + half ); j++ )
                sum += data.y[j];//tmp[ j ];
            data.y[ i ] = sum / filterSize;
        }
        // Adjust vector sizes
        data.y.erase(data.y.begin(),data.y.begin()+half);
        data.y.erase(data.y.end()-half,data.y.end());
    }
}

/**
 * Savitzky-Golay Helper
 */
void savitzkyGolay( Data& data, int SG[], int C, int size ){
    int half = size/2;
    for( int i = 0; i < half; i++ )
        data.y.push_back(data.y[i]);
    for( int i = 0; i < half; i++ )
        data.y.insert(data.y.begin(),
                      data.y[data.y[data.y.size() - half - i - 1]]);
    // Begin the filtering process
    double sum = 0;
    std::vector<double> tmp = data.y;
    for( int i = half; i < (data.n - half ); i++ ){
        sum = 0;
        for( int j = i - half; j <= (i + half); j++ )
            sum += tmp[ j ]*SG[ j - i + half ];
        data.y[ i ] = sum / C;
    }
    data.y.erase(data.y.begin(), data.y.begin()+half);
    data.y.erase(data.y.end()-half, data.y.end());
}

/**
 * Handler for filter options in the NMR analyzer
 * Dependent on the options file the user inputs
 */
void nmrFilter( Data& data, const Options& options ){
    int size = options.filterSize;
    switch (options.filterType) {
        case 3: // DFT 
            dft( data );
            break;
        case 2: //Savitzky-Golay
            switch (size) {
                case 5  :
                    for( int i = 0; i < options.filterPasses; i++ )
                        savitzkyGolay(data, SG5, 35, options.filterSize);
                    break;
                case 11 :
                    for( int i = 0; i < options.filterPasses; i++ )
                        savitzkyGolay(data, SG11, 429, options.filterSize);
                    break;
                case 17 :
                    for( int i = 0; i < options.filterPasses; i++ )
                        savitzkyGolay(data, SG17, 323, options.filterSize);
                    break;
                default :
                    printf("Error: Invalid filter size for Savitzky-Golay method\n");
                    break;
            }
            break;
        case 1: //Boxcar
            if( size == 0 ) printf("Boxcar method with filter size 0 requested; Filtering is, in essence, off");
            for( int i = 0; i < options.filterPasses; i++ )
                boxcar(data, size);
            break;
        case 0: break; // No filter 
        default : // Invalid option
            printf("Error: Invalid filtering option requested\n");
            break;
    }
}


