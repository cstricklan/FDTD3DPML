/*
 * FDTD.h
 *
 *  Created on: May 1, 2011
 *      Author: linus
 */

#include <complex>
#include "OutputSlice.h"
#include "OutputFile.h"

#ifndef FDTD_H_
#define FDTD_H_

class FDTD
{
		/***************** Device Params ************************/
		  double Lx;
			double Ly;
		  double r;
		  double d;
		  double dsub;
		  double dm;
		  double edge;
		  double er1;
		  double erm;
		  double erd;
		/***************** Device Params ************************/

		int NFREQ;
		long fmax;
		long fmin;
		double* FREQ;

		double nmax;

		int Nx, Ny, Nz;
		int NPML;
		int NRES;
		double BUFZ;

		double dx,dy,dz;
		double dt;
		double* ta;

		// Physical Size of Grid
		double Sx;
		double Sy;
		double Sz;

		int STEPS;

		// Record Planes
		int nz_ref; //position of reflection record plane
		int nz_trn; //position of transmission record plane


		// Grid Axis
		double *xa, *ya, *za;

		// PMLS
		double *sigz;

		// Materials
		double *Materials;
		double *URxx, *URyy, *URzz;
		double *ERxx, *ERyy, *ERzz;

		// Source
		int nz_src; 		// Source Plane
		double nref;
		double ntrn;
		double *g;
		double *Ex_src;
		double *Ey_src;
		double *Hx_src;
		double *Hy_src;

		// Update Coefficients
		double *mHx0, *mHx1, *mHx2;
		double *mHy0, *mHy1, *mHy2;
		double *mHz2, *mHz3;

		double *mDx0, *mDx1, *mDx2;
		double *mDy0, *mDy1, *mDy2;
		double *mDz2, *mDz3;

		double *mEx1, *mEy1, *mEz1;




		// Fields
		double *Hx, *Hy, *Hz;
		double *Dx, *Dy, *Dz;
		double *Ex, *Ey, *Ez;

		//Curl Terms
		double *CEx, *CEy, *CEz;
		double *CHx, *CHy, *CHz;

		//Integration Terms
		double* ICHz;
		double* ICEz;

		// REF/TRN Terms
		std::complex<double> *K;
		std::complex<double> *Exref, *Eyref, *Extrn, *Eytrn;
		std::complex<double> *SRC;


		OutputSlice** Out;
		OutputFile** Outf;

		void ComputeCurlEx();
		void ComputeCurlEy();
		void ComputeCurlEz();

		void ComputeCurlHx();
		void ComputeCurlHy();
		void ComputeCurlHz();

		void UpdateIntegration(double*, double*);
		inline int Index(int, int, int);
		void ComputeSource();
		void AddDevice();
		void ComputePML();
		void ComputeUpdateCoeff();
		void InitializeFFT();
		void UpdateFFT(int t);
		void OutputFFT(int t);
		double* GetSig(int, int, int);

	public:
		int ocount;
		void init();
		void exec();

		FDTD();
		virtual	~FDTD();
};

#endif /* FDTD_H_ */
