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

// #include <stdio.h>
#ifndef globals_h
#define globals_h
#include <vector>
#include "struct.cpp"

int SG5[5] = {-3,12,17,12,-3};
int SG11[11] = {-36,9,44,69,84,89,84,69,44,9,-36};
int SG17[17] = {-21,-6,7,18,27,34,39,42,43,42,39,34,27,18,7,-6,-21};

extern Spline splineCo;
extern Data _data;
extern Options options;
extern std::vector<Peak> peaks;
extern double TMSadjustment;

#endif 
