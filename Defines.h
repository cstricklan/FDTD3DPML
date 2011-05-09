#ifndef DEFINES_H_
#define DEFINES_H_

#define TIXML_USE_TICPP 1

#include<stdlib.h>
//#include<ticpp.h>
#include<map>
#include<string>
#include<iostream>
#include<vector>
/***********************************
 *
 * Constants
 *
 ***********************************
*/


const double e0=8.8541878176203898505365630317107502606083701665994498081024171524053950954599821142852891607182008932e-12;  // F/m
const double u0=1.2566370614359172953850573533118011536788677597500423283899778369231265625144835994512139301368468271e-6;   // H/m
const double c0=299792458; // m/s
const double pi=3.14159;

const double meters = 1;
const double decimeters = 1e-1 * meters;
const double centimeters = 1e-2 * meters;
const double millimeters = 1e-3 * meters;
const double inches = 2.54 * centimeters;
const double feet = 12 * inches;
const double seconds = 1;
const double hertz = 1/seconds;
const double kilohertz = 1e3 * hertz;
const double megahertz = 1e6 * hertz;
const double gigahertz = 1e9 * hertz;

#endif /* DEFINES_H_ */
