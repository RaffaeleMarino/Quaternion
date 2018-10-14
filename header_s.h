/*
 *  header_s.h
 *  Porj_EffectiveDiffusionNEW
 *
 *  Created by Marino Raffaele on 19/10/15.
 *  Copyright 2015 Raffaele Marino. All rights reserved.
 *
 */

#ifndef HEADER_S_H
#define HEADER_S_H
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iterator>
#include <algorithm>
#include <valarray>
#include <numeric>
#include <complex>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <random>


#define MAX_SIZE 50000000



using namespace std;
enum{CONF=1000}; //MAX1024*1000000 with -mcmodel=medium
enum{PASSI=50000};
const double PI=3.1415926535897933;
const double numberlc=5;
//const size_t  MAX_SIZE_ARRAY=CONF*PASSI;
const double k_B=+1.3806E-2; //femtoNewton*Micrometr*Kelvin{-1}
const double dt=0.001;
const double deltat=1;
const double the_time=50000000.;
const double myzero=1E-15;
const double tau=10000000000;
#endif

