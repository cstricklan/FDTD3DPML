/*
 * Utils.cpp
 *
 *  Created on: May 1, 2011
 *      Author: linus
 */

#include <complex>

double* Zeros(int Nx, int Ny, int Nz)
{
	double* Matrix = new double[Nx*Ny*Nz];

#pragma omp for
	for(int i=0; i<Nx; ++i)
		for(int j=0; j<Ny; ++j)
			for(int k=0; k<Nz; ++k)
				Matrix[i+j*Nx+k*Nx*Ny] = 0;

	return Matrix;
}

double* Zeros(int Size)
{
	double* Matrix = new double[Size];

#pragma omp for
	for(int i=0; i<Size; ++i)
		Matrix[i] = 0;

	return Matrix;
}


std::complex<double>* ZerosC(int Nx, int Ny, int Nz)
{
	std::complex<double>* Matrix = new std::complex<double>[Nx*Ny*Nz];

#pragma omp for
	for(int i=0; i<Nx; ++i)
		for(int j=0; j<Ny; ++j)
			for(int k=0; k<Nz; ++k)
				Matrix[i+j*Nx+k*Nx*Ny] = std::complex<double>(0,0);

	return Matrix;
}

std::complex<double>* ZerosC(int Size)
{
	std::complex<double>* Matrix = new std::complex<double>[Size];

#pragma omp for
	for(int i=0; i<Size; ++i)
		Matrix[i] = std::complex<double>(0,0);

	return Matrix;
}


double* Ones(int Nx, int Ny, int Nz)
{
	double* Matrix = new double[Nx*Ny*Nz];

#pragma omp for
	for(int i=0; i<Nx; ++i)
		for(int j=0; j<Ny; ++j)
			for(int k=0; k<Nz; ++k)
				Matrix[i+j*Nx+k*Nx*Ny] = 1;

	return Matrix;
}



double* linspace(long Low, long High, int Num)
{
	double* Array = new double[Num];

	double diff = High-Low;
	diff = diff/(double)Num;

	Array[0] = Low;
	Array[Num-1] = High;

	for(int i=1;i<Num-1;++i)
		Array[i] = Array[i-1]+diff;

	return Array;
}

double* FillArray(int Low, int High, int Diff=1)
{
	int Size = (High-Low)*Diff;
	double* Array = new double[Size];

	for(int i=0; i<Size; ++i)
		Array[i] = Low + i*Diff;

	return Array;
}



double* MultArray(double* Array, int Size, const int value)
{
	double* Result = new double[Size];

#pragma omp for
	for(int i=0;i<Size;++i)
		Result[i] = Array[i] * value;

	return Result;
}

double* MultArray(double* Array, int Size, double value)
{
	double* Result = new double[Size];

#pragma omp for
	for(int i=0;i<Size;++i)
		Result[i] = Array[i] * value;

	return Result;
}


double* MeshGridx(double* m, int X)
{
	double *mt = Zeros(X*X);

	for(int i=0; i<X; ++i)
	{
		for(int j=0; j<X; ++j)
		{
			mt[i+j*X] = m[i];
		}
	}

	return mt;
}

double* MeshGridy(double* n, int Y)
{
	double *nt = Zeros(Y*Y);

	for(int i=0; i<Y; ++i)
	{
		for(int j=0; j<Y; ++j)
		{
			nt[i*Y + j] = n[i];
		}
	}

	return nt;
}


std::complex<double>* sft(std::complex<double>* Array, int Length)
{
	int N = Length;
	std::complex<double> i = std::complex<double>(0,1);
	std::complex<double> seed = std::complex<double>(2*M_PI/N);

	std::complex<double> *F = ZerosC(N);

	for (int k=-1;k<N-1;++k)
	{
		for (int n=-1;n<N-1;++n)
	  {
			F[k+1] = F[k+1] + Array[n+1]*exp(-i*seed*std::complex<double>(k*n));
	  }
	}

	return F;
}
//
//
//std::complex<double>* sft(std::complex<double>* Array, int Nx, int Ny)
//{
//	F = f;
//	for nx = 1 : Nx
//	    F(nx,:) = sft(F(nx,:));
//	end
//	for ny = 1 : Ny
//	    F(:,ny) = sft(F(:,ny));
//	end
//
//}


void CopyArray(std::complex<double> *a, std::complex<double> *b, int Length)
{
	for(int i=0; i<Length; ++i)
	{
		b[i] = a[i];
	}
}


