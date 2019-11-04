//
//  globals.cpp
//  nmr
//
//  Created by Eduardo Rubio on 10/19/19.
//  Copyright © 2019 Eduardo Rubio. All rights reserved.
//

/**
 * Global Variables needed throughout the program 
 */

#include <stdio.h>

int SG5[5] = {-3,12,17,12,-3};
int SG11[11] = {-36,9,44,69,84,89,84,69,44,9,-36};
int SG17[17] = {-21,-6,7,18,27,34,39,42,43,42,39,34,27,18,7,-6,-21};

Spline splineCo;
Data data;
Options options;
std::vector<Peak> peaks;
double TMSadjustment; 
