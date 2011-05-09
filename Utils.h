/*
 * Utils.h
 *
 *  Created on: May 1, 2011
 *      Author: linus
 */

#ifndef UTILS_H_
#define UTILS_H_

#include <complex>

double* Zeros(int Nx, int Ny, int Nz);
double* Zeros(int Size);
std::complex<double>* ZerosC(int Nx, int Ny, int Nz);
std::complex<double>* ZerosC(int Size);


double* Ones(int Nx, int Ny, int Nz);
double* linspace(long Low, long High, int Num);
double* FillArray(int Low, int High, int Diff=1);
double* MultArray(double* Array, int Size, const int value);
double* MultArray(double* Array, int Size, double value);

double* MeshGridx(double* m, int X);
double* MeshGridy(double* n, int Y);

std::complex<double>* sft(std::complex<double>* Array, int Length);
std::complex<double>* stf(std::complex<double>* Array, int Nx, int Ny);
void CopyArray(std::complex<double> *a, std::complex<double> *b, int Length);
#endif /* UTILS_H_ */
