/*
 * FDTD.cpp
 *
 *  Created on: May 1, 2011
 *      Author: linus
 */

#include "FDTD.h"
#include "Defines.h"
#include <iostream>
#include "Utils.h"
#include <algorithm>
#include<math.h>
#include"OutputPower.h"
#include<complex>
#include<fftw3.h>


using namespace std;

void FDTD::init()
{
	// Values;
	FREQ = linspace(fmin,fmax,NFREQ);

	double lam_min  = c0/fmax;
	dx = lam_min/nmax/NRES;
	dy = lam_min/nmax/NRES;
	dz = lam_min/nmax/NRES;

	//SNAP GRID TO CRITICAL DIMENSIONS
	Nx = 2*ceil(Lx/dx/2) + 1; //make Nx odd
	dx = Lx/Nx;

	Ny = 2*ceil(Ly/dy/2) + 1; //make Ny odd
	dy = Ly/Ny;

	Nz = ceil(d/dz); //Nz does not have to be odd
	dz = d/Nz;

	// COMPUTE SIZE OF GRID
	Sx = Lx;
	Sy = Ly;

	Sz = BUFZ + d + BUFZ;
	Nz = NPML + 3 + ceil(Sz/dz) + 2 + NPML;
	Sz = Nz*dz;


	// COMPUTE STABLE TIME STEP
	double ds[] = {dx, dy, dz};
	double dmin = *std::min_element(ds, ds+3);
	dt   = 0.5*dmin/c0;

	cout << "Steps: " << STEPS << endl;

	// Time Axis
	ta = MultArray(FillArray(0,STEPS-1), STEPS,dt);

	// Grid Axis
	xa = MultArray(FillArray(0,Nx-1),Nx,dx);
	ya = MultArray(FillArray(0,Ny-1),Ny,dy);
	za = MultArray(FillArray(0,Nz-1),Nz,dz);


	// Material
	URxx = Ones(Nx,Ny,Nz);
	URyy = Ones(Nx,Ny,Nz);
	URzz = Ones(Nx,Ny,Nz);
	ERxx = Ones(Nx,Ny,Nz);
	ERyy = Ones(Nx,Ny,Nz);
	ERzz = Ones(Nx,Ny,Nz);
	AddDevice();

	// Setup Record Planes
	nz_ref = NPML + 1;
	nz_trn = Nz - NPML;

	ComputeSource();
	InitializeFFT();
	ComputePML();
	ComputeUpdateCoeff();

	ocount = 2;
	Outf = new OutputFile*[ocount];
	Outf[0] = new OutputFile("/home/linus/Source/FDTD/FDTDPML/Debug/Results/file/", nz_ref);
	Outf[1] = new OutputFile("/home/linus/Source/FDTD/FDTDPML/Debug/Results/file/", nz_trn);

//
//	for(int i=0; i<ocount;++i)
//	{
//		Outf[i] = new OutputFile("/home/linus/Source/FDTD/FDTDPML/Debug/Results/File/", i+nz_src);
//	}
}

void FDTD::exec()
{
#pragma omp parallel
	{
	// Fields
	Hx = Zeros(Nx,Ny,Nz);	Hy = Zeros(Nx,Ny,Nz); Hz = Zeros(Nx,Ny,Nz);
	Dx = Zeros(Nx,Ny,Nz); Dy = Zeros(Nx,Ny,Nz);	Dz = Zeros(Nx,Ny,Nz);
	Ex = Zeros(Nx,Ny,Nz);	Ey = Zeros(Nx,Ny,Nz);	Ez = Zeros(Nx,Ny,Nz);

	//Curl Terms
	CEx = Zeros(Nx,Ny,Nz); CEy = Zeros(Nx,Ny,Nz);	CEz = Zeros(Nx,Ny,Nz);
	CHx = Zeros(Nx,Ny,Nz); CHy = Zeros(Nx,Ny,Nz);	CHz = Zeros(Nx,Ny,Nz);

	//Integration Terms
	ICHz = Zeros(Nx,Ny,Nz); ICEz = Zeros(Nx,Ny,Nz);
	}
	for(int T=0;T<STEPS;++T)
	{
		#pragma omp parallel
		{
		// Compute Curls
		ComputeCurlEx();
		ComputeCurlEy();
		ComputeCurlEz();


		#pragma omp barrier
		// TF/SF

		#pragma omp for
		for(int i=0;i<Nx;++i)
		{
			for(int j=0;j<Ny;++j)
			{
				CEx[Index(i,j,nz_src-1)] = CEx[Index(i,j,nz_src-1)] + Ey_src[T]/dz;
				CEy[Index(i,j,nz_src-1)] = CEy[Index(i,j,nz_src-1)] - Ex_src[T]/dz;
			}
		}
		// TF/SF

		#pragma omp barrier
		// Update H Integrations
		UpdateIntegration(ICEz, CEz);

		// UPdate H Field

		#pragma omp for
		for(int i=0;i<Nx;++i)
		{
			for(int j=0;j<Ny;++j)
			{
				for(int k=0;k<Nz;++k)
				{
					Hx[Index(i,j,k)] = mHx1[Index(i,j,k)]*Hx[Index(i,j,k)] + mHx2[Index(i,j,k)]*CEx[Index(i,j,k)];
					Hy[Index(i,j,k)] = mHy1[Index(i,j,k)]*Hy[Index(i,j,k)] + mHy2[Index(i,j,k)]*CEy[Index(i,j,k)];
					Hz[Index(i,j,k)] = Hz[Index(i,j,k)] + mHz2[Index(i,j,k)]*CEz[Index(i,j,k)] + mHz3[Index(i,j,k)]*ICEz[Index(i,j,k)];
				}
			}
		}

		#pragma omp barrier
		// Compute Curls
		ComputeCurlHx();
		ComputeCurlHy();
		ComputeCurlHz();

		#pragma omp barrier
		// TF/SF
		#pragma omp for
		for(int i=0;i<Nx;++i)
		{
			for(int j=0;j<Ny;++j)
			{
				CHx[Index(i,j,nz_src-1)] = CHx[Index(i,j,nz_src-1)] + Hy_src[T]/dz;
				CHy[Index(i,j,nz_src-1)] = CHy[Index(i,j,nz_src-1)] - Hx_src[T]/dz;
			}
		}
		// TF/SF

		#pragma omp barrier
	  // Update D Integrations
	  UpdateIntegration(ICHz,CHz);

		#pragma omp barrier

	  // Update Dz
		#pragma omp for
		for(int i=0;i<Nx;++i)
		{
			for(int j=0;j<Ny;++j)
			{
				for(int k=0;k<Nz;++k)
				{
					Dx[Index(i,j,k)] = mDx1[Index(i,j,k)]*Dx[Index(i,j,k)] + mDx2[Index(i,j,k)]*CHx[Index(i,j,k)];
					Dy[Index(i,j,k)] = mDy1[Index(i,j,k)]*Dy[Index(i,j,k)] + mDy2[Index(i,j,k)]*CHy[Index(i,j,k)];
					Dz[Index(i,j,k)] = Dz[Index(i,j,k)] + mDz2[Index(i,j,k)]*CHz[Index(i,j,k)] + mDz3[Index(i,j,k)]*ICHz[Index(i,j,k)];
				}
			}
		}

		#pragma omp barrier

	  // Update Ez
		#pragma omp for
		for(int i=0;i<Nx;++i)
		{
			for(int j=0;j<Ny;++j)
			{
				for(int k=0;k<Nz;++k)
				{
					Ex[Index(i,j,k)] = mEx1[Index(i,j,k)]*Dx[Index(i,j,k)];
					Ey[Index(i,j,k)] = mEy1[Index(i,j,k)]*Dy[Index(i,j,k)];
					Ez[Index(i,j,k)] = mEz1[Index(i,j,k)]*Dz[Index(i,j,k)];
				}
			}
		}

		#pragma omp barrier
		UpdateFFT(T);
		}

		if(T%1000==0)
	  {
	  	cout << "Step: " << T << endl;
  		Out[0]->Output(T, Ex, Nx, Ny, Nz);
  		Out[1]->Output(T, Ey, Nx, Ny, Nz);
  		Out[2]->Output(T, Ez, Nx, Ny, Nz);
	  }
		if(T%1000==0)
		{
			OutputFFT(T);
		}
		if(T%1000==0)
		{
			for(int i=0;i<ocount;++i)
			{
				Outf[i]->Output(T, Ex, "Ex", Nx,Ny);
				Outf[i]->Output(T, Ey, "Ey", Nx,Ny);
				Outf[i]->Output(T, Ez, "Ez", Nx,Ny);
			}
		}
	}

	OutputFFT(STEPS+1);

}

/************************************** CURLS **************************************/
void FDTD::ComputeCurlEx()
{
	#pragma omp for
	for (int nx=0 ; nx<Nx ; nx++)
	{
		// ny=Ny, nz=Nz boundary
		CEx[Index(nx,(Ny-1),(Nz-1))] = (Ez[Index(nx,0,(Nz-1))] - Ez[Index(nx,(Ny-1),(Nz-1))])/dy - (Ey[Index(nx,(Ny-1),0)] - Ey[Index(nx,(Ny-1),(Nz-1))])/dz;

		// nz=Nz boundary
		for (int ny=0 ; ny<Ny-1 ; ny++)
		{
			CEx[Index(nx,ny,(Nz-1))] = (Ez[Index(nx,(ny+1),(Nz-1))] - Ez[Index(nx,ny,(Nz-1))])/dy - (Ey[Index(nx,ny,0)] - Ey[Index(nx,ny,(Nz-1))])/dz;
		}

		// ny=Ny boundary
		for (int nz=0 ; nz<Nz-1 ; nz++)
		{
			CEx[Index(nx,(Ny-1),nz)] = (Ez[Index(nx,0,nz)] - Ez[Index(nx,(Ny-1),nz)])/dy - (Ey[Index(nx,(Ny-1),(nz+1))] - Ey[Index(nx,(Ny-1),nz)])/dz;
		}

		// rest of array
		for (int ny=0 ; ny<Ny-1 ; ny++)
		{
			for (int nz=0 ; nz<Nz-1 ; nz++)
			{
				CEx[Index(nx,ny,nz)] = (Ez[Index(nx,(ny+1),nz)] - Ez[Index(nx,ny,nz)])/dy - (Ey[Index(nx,ny,(nz+1))] - Ey[Index(nx,ny,nz)])/dz;
			}
		}
	}
}
void FDTD::ComputeCurlEy()
{
	#pragma omp for
	for (int ny=0 ; ny<Ny ; ny++)
	{
		// nx=Nx, nz=Nz boundary
		CEy[Index(Nx-1,ny,Nz-1)] = (Ex[Index(Nx-1,ny,0)] - Ex[Index(Nx-1,ny,Nz-1)])/dz - (Ez[Index(0,ny,Nz-1)] - Ez[Index(Nx-1,ny,Nz-1)])/dx;

		// nz=Nz boundary
		for (int nx=0 ; nx<Nx-1 ; nx++)
		{
			CEy[Index(nx,ny,Nz-1)] = (Ex[Index(nx,ny,0)] - Ex[Index(nx,ny,Nz-1)])/dz - (Ez[Index(nx+1,ny,Nz-1)] - Ez[Index(nx,ny,Nz-1)])/dx;
		}

		// nx=Nx boundary
		for (int nz=0 ; nz<Nz-1 ; nz++)
		{
			CEy[Index(Nx-1,ny,nz)] = (Ex[Index(Nx-1,ny,nz+1)] - Ex[Index(Nx-1,ny,nz)])/dz - (Ez[Index(0,ny,nz)] - Ez[Index(Nx-1,ny,nz)])/dx;
		}

		// rest of array
		for (int nx=0 ; nx<Nx-1 ; nx++)
		{
			for (int nz=0 ; nz<Nz-1 ; nz++)
			{
				CEy[Index(nx,ny,nz)] = (Ex[Index(nx,ny,nz+1)] - Ex[Index(nx,ny,nz)])/dz - (Ez[Index(nx+1,ny,nz)] - Ez[Index(nx,ny,nz)])/dx;
			}
		}
	}
}
void FDTD::ComputeCurlEz()
{
#pragma omp for

	for (int nz=0 ; nz<Nz ; nz++)
	{
		// nx=Nx, ny=Ny boundary
		CEz[Index(Nx-1,Ny-1,nz)] = (Ey[Index(0,Ny-1,nz)] - Ey[Index(Nx-1,Ny-1,nz)])/dx - (Ex[Index(Nx-1,0,nz)] - Ex[Index(Nx-1,Ny-1,nz)])/dy;

		// ny=Ny boundary
		for (int nx=0 ; nx<Nx-1 ; nx++)
		{
			CEz[Index(nx,Ny-1,nz)] = (Ey[Index(nx+1,Ny-1,nz)] - Ey[Index(nx,Ny-1,nz)])/dx - (Ex[Index(nx,0,nz)] - Ex[Index(nx,Ny-1,nz)])/dy;
		}

		// nx=Nx boundary
		for (int ny=0 ; ny<Ny-1 ; ny++)
		{
			CEz[Index(Nx-1,ny,nz)] = (Ey[Index(0,ny,nz)] - Ey[Index(Nx-1,ny,nz)])/dx - (Ex[Index(Nx-1,ny+1,nz)] - Ex[Index(Nx-1,ny,nz)])/dy;
		}

		// rest of array
		for (int nx=0 ; nx<Nx-1 ; nx++)
		{
			for (int ny=0 ; ny<Ny-1 ; ny++)
			{
				CEz[Index(nx,ny,nz)] = (Ey[Index(nx+1,ny,nz)] - Ey[Index(nx,ny,nz)])/dx - (Ex[Index(nx,ny+1,nz)] - Ex[Index(nx,ny,nz)])/dy;
			}
		}
	}

}

void FDTD::ComputeCurlHx()
{
#pragma omp for

	for (int nx=0 ; nx<Nx ; nx++)
		{
			// ny=nz=0 boundary
			CHx[Index(nx,0,0)] = (Hz[Index(nx,0,0)] - Hz[Index(nx,Ny-1,0)])/dy - (Hy[Index(nx,0,0)] - Hy[Index(nx,0,Nz-1)])/dz;

			// nz=0 boundary
			for (int ny=1 ; ny<Ny ; ny++)
			{
				CHx[Index(nx,ny,0)] = (Hz[Index(nx,ny,0)] - Hz[Index(nx,ny-1,0)])/dy - (Hy[Index(nx,ny,0)] - Hy[Index(nx,ny,Nz-1)])/dz;
			}

			// ny=0 boundary
			for (int nz=1 ; nz<Nz ; nz++)
			{
				CHx[Index(nx,0,nz)] = (Hz[Index(nx,0,nz)] - Hz[Index(nx,Ny-1,nz)])/dy - (Hy[Index(nx,0,nz)] - Hy[Index(nx,0,nz-1)])/dz;
			}

			// rest of array
			for (int ny=1 ; ny<Ny ; ny++)
			{
				for (int nz=1 ; nz<Nz ; nz++)
				{
					CHx[Index(nx,ny,nz)] = (Hz[Index(nx,ny,nz)] - Hz[Index(nx,ny-1,nz)])/dy - (Hy[Index(nx,ny,nz)] - Hy[Index(nx,ny,nz-1)])/dz;
				}
			}
		}
}
void FDTD::ComputeCurlHy()
{
#pragma omp for

	for (int ny=0 ; ny<Ny ; ny++)
		{
			// nx=nz=0 boundary
			CHy[Index(0,ny,0)] = (Hx[Index(0,ny,0)] - Hx[Index(0,ny,Nz-1)])/dz - (Hz[Index(0,ny,0)] - Hz[Index(Nx-1,ny,0)])/dx;

			// nz=0 boundary
			for (int nx=1 ; nx<Nx ; nx++)
			{
				CHy[Index(nx,ny,0)] = (Hx[Index(nx,ny,0)] - Hx[Index(nx,ny,Nz-1)])/dz - (Hz[Index(nx,ny,0)] - Hz[Index(nx-1,ny,0)])/dx;
			}

			// nx=0 boundary
			for (int nz=1 ; nz<Nz ; nz++)
			{
				CHy[Index(0,ny,nz)] = (Hx[Index(0,ny,nz)] - Hx[Index(0,ny,nz-1)])/dz - (Hz[Index(0,ny,nz)] - Hz[Index(Nx-1,ny,nz)])/dx;
			}

			// rest of array
			for (int nx=1 ; nx<Nx ; nx++)
			{
				for (int nz=1 ; nz<Nz ; nz++)
				{
					CHy[Index(nx,ny,nz)] = (Hx[Index(nx,ny,nz)] - Hx[Index(nx,ny,nz-1)])/dz - (Hz[Index(nx,ny,nz)] - Hz[Index(nx-1,ny,nz)])/dx;
				}
			}
		}
}
void FDTD::ComputeCurlHz()
{
#pragma omp for

	for (int nz=0 ; nz<Nz ; nz++)
		{
			// nx=ny=0 boundary
			CHz[Index(0,0,nz)] = (Hy[Index(0,0,nz)] - Hy[Index(Nx-1,0,nz)])/dx - (Hx[Index(0,0,nz)] - Hx[Index(0,Ny-1,nz)])/dy;

			// ny=0 boundary
			for (int nx=1 ; nx<Nx ; nx++)
			{
				CHz[Index(nx,0,nz)] = (Hy[Index(nx,0,nz)] - Hy[Index(nx-1,0,nz)])/dx - (Hx[Index(nx,0,nz)] - Hx[Index(nx,Ny-1,nz)])/dy;
			}

			// nx=0 boundary
			for (int ny=1 ; ny<Ny ; ny++)
			{
				CHz[Index(0,ny,nz)] = (Hy[Index(0,ny,nz)] - Hy[Index(Nx-1,ny,nz)])/dx - (Hx[Index(0,ny,nz)] - Hx[Index(0,ny-1,nz)])/dy;
			}

			// rest of array
			for (int nx=1 ; nx<Nx ; nx++)
			{
				for (int ny=1 ; ny<Ny ; ny++)
				{
					CHz[Index(nx,ny,nz)] = (Hy[Index(nx,ny,nz)] - Hy[Index(nx-1,ny,nz)])/dx - (Hx[Index(nx,ny,nz)] - Hx[Index(nx,ny-1,nz)])/dy;
				}
			}
		}
}
/************************************** CURLS **************************************/


void FDTD::UpdateIntegration(double* Int, double* C)
{
#pragma omp for
	for(int i=0;i<Nx;++i)
		for(int j=0;j<Ny;++j)
			for(int k=0;k<Nz;++k)
				Int[i+j*Nx+k*Nx*Ny] += C[i+j*Nx+k*Nx*Ny];
}


inline int FDTD::Index(int x, int y, int z)
{
	return x+y*Nx+z*Nx*Ny;
}


void FDTD::ComputeSource()
{
	// Source Position
	nz_src = NPML + 2;

	int P[]={1,0};

	// COMPUTE REFRACTIVE INDICES IN UTILITY PLANES
	nref = sqrt(ERxx[Index(1,1,nz_ref)]);        //refractive index above grating
	ntrn = sqrt(ERxx[Index(1,1,nz_trn)]);        //refractive index below grating

	double tau  = 0.5/fmax;
	double t0   = 6*tau;

	double A    = -sqrt(ERxx[Index(1,1,nz_src)]/URxx[Index(1,1,nz_src)]);
	double delt = 0.5*dz/c0 + dt/2;

	g = new double[STEPS];

	for(int i=0;i<STEPS;++i)
		g[i] = exp(-1*(pow((ta[i]-t0)/tau,2)));

	Ex_src = MultArray(g,STEPS,P[0]);
	Ey_src = MultArray(g,STEPS,P[1]);

	double s = -nref*dz/(2*c0) + dt/2;
	double hx = -P[1]*sqrt(ERxx[Index(1,1,nz_ref)]*URxx[Index(1,1,nz_ref)]);
	double hy = P[0]*sqrt(ERxx[Index(1,1,nz_ref)]*URxx[Index(1,1,nz_ref)]);

	Hx_src = Zeros(STEPS);
	Hy_src = Zeros(STEPS);

	for(int i=0;i<STEPS;++i)
	{
		Hx_src[i] = hx*exp(-1*pow((ta[i]-t0+s)/tau,2));
		Hy_src[i] = hy*exp(-1*pow((ta[i]-t0+s)/tau,2));
	}

}



void FDTD::AddDevice()
{
//	COMPUTE POSITION INDICES
	int nz  = round(dsub/dz);
	int nza = NPML + round(BUFZ/dz) + 1;
	int nzb = nza + nz - 1;

	//ADD DIELECTRIC
	for(int i=0;i<Nx;++i)
	{
		for(int j=0;j<Ny;++j)
		{
			for(int k=nza;k<=nzb;++k)
			{
				ERxx[i+j*Nx+k*Nx*Ny] = erd;
				ERyy[i+j*Nx+k*Nx*Ny] = erd;
				ERzz[i+j*Nx+k*Nx*Ny] = erd;
			}
		}
	}

	int nzbm = nza-1;
	nz = round(dm/dz);
	int nzam = nzbm-nz;

	int rx = round(r/dx);
	int srx = round((Lx-r)/dx);

	int ex = round(edge/dx);

	// Top Edges
	for(int i=0;i<ex;++i)
	{
		for(int j=0;j<Ny;++j)
		{
			for( int k=nzam;k<=nzbm;++k)
			{
				ERxx[i+j*Nx+k*Nx*Ny] = erm;
				ERyy[i+j*Nx+k*Nx*Ny] = erm;
				ERzz[i+j*Nx+k*Nx*Ny] = erm;
			}
		}
	}

	// Bottom Edges
	for(int i=Nx-ex-1;i<Nx;++i)
	{
		for(int j=0;j<Ny;++j)
		{
			for( int k=nzam;k<=nzbm;++k)
			{
				ERxx[i+j*Nx+k*Nx*Ny] = erm;
				ERyy[i+j*Nx+k*Nx*Ny] = erm;
				ERzz[i+j*Nx+k*Nx*Ny] = erm;
			}
		}
	}

	// Left Edges
	for(int i=0;i<Nx;++i)
	{
		for(int j=0;j<ex;++j)
		{
			for( int k=nzam;k<=nzbm;++k)
			{
				ERxx[i+j*Nx+k*Nx*Ny] = erm;
				ERyy[i+j*Nx+k*Nx*Ny] = erm;
				ERzz[i+j*Nx+k*Nx*Ny] = erm;
			}
		}
	}

	// right Edges
	for(int i=0;i<Nx;++i)
	{
		for(int j=Ny-ex-1;j<Ny;++j)
		{
			for( int k=nzam;k<=nzbm;++k)
			{
				ERxx[i+j*Nx+k*Nx*Ny] = erm;
				ERyy[i+j*Nx+k*Nx*Ny] = erm;
				ERzz[i+j*Nx+k*Nx*Ny] = erm;
			}
		}
	}

	// Middle
	for(int i=srx; i<srx+rx; ++i)
	{
		for(int j=srx; j<srx+rx; ++j)
		{
			for( int k=nzam;k<=nzbm;++k)
			{
				ERxx[i+j*Nx+k*Nx*Ny] = erm;
				ERyy[i+j*Nx+k*Nx*Ny] = erm;
				ERzz[i+j*Nx+k*Nx*Ny] = erm;
			}
		}
	}


}



void FDTD::ComputePML()
{
	//	COMPUTE PML PARAMETERS
	//	COMPUTE SIZE OF 2X GRID
	int Nx2 = 2*Nx;
	int Ny2 = 2*Ny;
	int Nz2 = 2*Nz;

	// COMPUTE sigz
	sigz = Zeros(Nx2,Ny2,Nz2);

#pragma omp for
	for (int i=0;i<Nx2;++i)
	{
		for(int j=0;j<Ny2;++j)
		{
			for(int k=0;k<2*NPML;++k)
			{
				int nz = 2*NPML - k + 1;
				sigz[i+j*Nx2+nz*Nx2*Ny2] = (0.5*e0/dt)*pow(((double)k/(double)2/(double)NPML),3);
			}

			for(int k=0;k<2*NPML;++k)
			{
				int nz = Nz2 - 2*NPML + k;
			  sigz[i+j*Nx2+nz*Nx2*Ny2] = (0.5*e0/dt)*pow(((double)k/(double)2/(double)NPML),3);
			}
		}
	}
}

double* FDTD::GetSig(int xstart, int ystart, int zstart)
{
	double *sig = Zeros(Nx,Ny,Nz);

	int Nx2 = 2*Nx;
	int Ny2 = 2*Ny;
	int Nz2 = 2*Nz;

	for(int i=0;i<Nx;++i)
		for(int j=0;j<Ny;++j)
			for(int k=0;k<Nz;++k)
				sig[Index(i,j,k)] = sigz[((i*2)+xstart) + (((j*2)+ystart)*Nx2) + (((k*2)+zstart)*Nx2*Ny2)];

	return sig;


}

void FDTD::ComputeUpdateCoeff()
{
	//COMPUTE UPDATE COEFFICIENTS
	mHx0 = Zeros(Nx,Ny,Nz); mHx1 = Zeros(Nx,Ny,Nz); mHx2 = Zeros(Nx,Ny,Nz);
	mHy0 = Zeros(Nx,Ny,Nz); mHy1 = Zeros(Nx,Ny,Nz); mHy2 = Zeros(Nx,Ny,Nz);
	mHz2 = Zeros(Nx,Ny,Nz); mHz3 = Zeros(Nx,Ny,Nz);

	mDx0 = Zeros(Nx,Ny,Nz); mDx1 = Zeros(Nx,Ny,Nz); mDx2 = Zeros(Nx,Ny,Nz);
	mDy0 = Zeros(Nx,Ny,Nz); mDy1 = Zeros(Nx,Ny,Nz); mDy2 = Zeros(Nx,Ny,Nz);
	mDz2 = Zeros(Nx,Ny,Nz); mDz3 = Zeros(Nx,Ny,Nz);

	mEx1 = Zeros(Nx,Ny,Nz); mEy1 = Zeros(Nx,Ny,Nz); mEz1 = Zeros(Nx,Ny,Nz);

	double *sig;
	// Hx
	sig  = GetSig(0,1,1);
	for(int i=0;i<Nx;++i)
	{
		for(int j=0;j<Ny;++j)
		{
			for(int k=0;k<Nz;++k)
			{
				mHx0[Index(i,j,k)] = (1/dt) + sig[Index(i,j,k)]/(2*e0);
				mHx1[Index(i,j,k)] = ((1/dt) - sig[Index(i,j,k)]/(2*e0))/mHx0[Index(i,j,k)];
				mHx2[Index(i,j,k)] = -c0/URxx[Index(i,j,k)]/mHx0[Index(i,j,k)];
			}
		}
	}
	delete sig;

	// Hy
	sig  = GetSig(1,0,1);
	for(int i=0;i<Nx;++i)
	{
		for(int j=0;j<Ny;++j)
		{
			for(int k=0;k<Nz;++k)
			{
				mHy0[Index(i,j,k)] = (1/dt) + sig[Index(i,j,k)]/(2*e0);
				mHy1[Index(i,j,k)] = ((1/dt) - sig[Index(i,j,k)]/(2*e0))/mHy0[Index(i,j,k)];
				mHy2[Index(i,j,k)] = - c0/URyy[Index(i,j,k)]/mHy0[Index(i,j,k)];
			}
		}
	}
	delete sig;

	// Hz
	sig  = GetSig(1,1,0);
	for(int i=0;i<Nx;++i)
	{
		for(int j=0;j<Ny;++j)
		{
			for(int k=0;k<Nz;++k)
			{
				mHz2[Index(i,j,k)] = - c0*dt/URzz[Index(i,j,k)];
				mHz3[Index(i,j,k)] = -(c0*pow(dt,2)/e0)*sig[Index(i,j,k)]/URzz[Index(i,j,k)];
			}
		}
	}
	delete sig;

	// Dx
	sig  = GetSig(1,0,0);
	for(int i=0;i<Nx;++i)
	{
		for(int j=0;j<Ny;++j)
		{
			for(int k=0;k<Nz;++k)
			{
				mDx0[Index(i,j,k)] = (1/dt) + sig[Index(i,j,k)]/(2*e0);
				mDx1[Index(i,j,k)] = ((1/dt) - sig[Index(i,j,k)]/(2*e0))/mDx0[Index(i,j,k)];
				mDx2[Index(i,j,k)] = c0/mDx0[Index(i,j,k)];
			}
		}
	}
	delete sig;

	// Dy
	sig  = GetSig(0,1,0);
	for(int i=0;i<Nx;++i)
	{
		for(int j=0;j<Ny;++j)
		{
			for(int k=0;k<Nz;++k)
			{
				mDy0[Index(i,j,k)] = (1/dt) + sig[Index(i,j,k)]/(2*e0);
				mDy1[Index(i,j,k)] = ((1/dt) - sig[Index(i,j,k)]/(2*e0))/mDx0[Index(i,j,k)];
				mDy2[Index(i,j,k)] = c0/mDx0[Index(i,j,k)];
			}
		}
	}
	delete sig;

	// Dz
	sig  = GetSig(0,0,1);
	for(int i=0;i<Nx;++i)
	{
		for(int j=0;j<Ny;++j)
		{
			for(int k=0;k<Nz;++k)
			{
				mDz2[Index(i,j,k)] = c0*dt;
				mDz3[Index(i,j,k)] = (c0*pow(dt,2)/e0)*sig[Index(i,j,k)];
			}
		}
	}
	delete sig;

	// Ex, Ey, and Ez
	for(int i=0;i<Nx;++i)
	{
		for(int j=0;j<Ny;++j)
		{
			for(int k=0;k<Nz;++k)
			{
				mEx1[Index(i,j,k)]  = 1/ERxx[Index(i,j,k)];
				mEy1[Index(i,j,k)]  = 1/ERyy[Index(i,j,k)];
				mEz1[Index(i,j,k)]  = 1/ERzz[Index(i,j,k)];
			}
		}
	}
}

void FDTD::InitializeFFT()
{
	//	INITIALIZE FOURIER TRANSFORMS
	//	KERNELS FOR REFLECTANCE AND TRANSMITTANCE

	complex<double> i = complex<double>(0,1);
	complex<double> seed = complex<double>(2*M_PI*dt);

	K = new complex<double>[NFREQ];

	for(int f=0;f<NFREQ;++f)
		K[f] = exp(-i*seed*complex<double>(FREQ[f]));  //kernels for sweep

	Exref = ZerosC(Nx,Ny,NFREQ);     //steady-state reflected field for Ex
	Eyref = ZerosC(Nx,Ny,NFREQ);     //steady-state reflected field for Ey
	Extrn = ZerosC(Nx,Ny,NFREQ);     //steady-state transmitted field for Ex
	Eytrn = ZerosC(Nx,Ny,NFREQ);     //steady-state transmitted field for Ey
	SRC   = ZerosC(NFREQ);         //source transform
}


void FDTD::UpdateFFT(int t)
{
	// Update Fourier Transforms

	for(int i=0;i<Nx;++i)
	{
		for(int j=0;j<Ny;++j)
		{
			for(int f=0;f<NFREQ;++f)
			{
				Exref[i+j*Nx+f*Nx*Ny] = Exref[i+j*Nx+f*Nx*Ny] + pow(K[f],t)*Ex[Index(i,j,nz_ref)]*dt;
				Eyref[i+j*Nx+f*Nx*Ny] = Eyref[i+j*Nx+f*Nx*Ny] + pow(K[f],t)*Ey[Index(i,j,nz_ref)]*dt;
				Extrn[i+j*Nx+f*Nx*Ny] = Extrn[i+j*Nx+f*Nx*Ny] + pow(K[f],t)*Ex[Index(i,j,nz_trn)]*dt;
				Eytrn[i+j*Nx+f*Nx*Ny] = Eytrn[i+j*Nx+f*Nx*Ny] + pow(K[f],t)*Ey[Index(i,j,nz_trn)]*dt;
			}
		}
	}

	for(int f=0;f<NFREQ;++f)
	{
		SRC[f] = SRC[f] + (pow(K[f],t))*g[t]*dt;
	}
}



void FDTD::OutputFFT(int t)
{
	//	COMPcalUTE REFLECTANCE AND TRANSMITTANCE
	//	INITIALIZE REFLECTANCE AND TRANSMITTANCE


//	complex<double> *REF  = ZerosC(NFREQ);
//	complex<double> *TRN  = ZerosC(NFREQ);
//	complex<double> *CON  = ZerosC(NFREQ);
//
//
//	//trans ex/ey
//  // ref ey/hx
//
//	/********************************* SRC *******************************************/
//	complex<double> SRCSUM = 0;
//	for(int nfreq=0; nfreq<NFREQ;++nfreq)
//	{
//		SRCSUM += SRC[nfreq];
//	}
//	/********************************* SRC *******************************************/
//
//
//	/********************************* TRANS *******************************************/
//	complex<double> *Pz = ZerosC(Nx*Ny*NFREQ);
//	for(int nfreq=0;nfreq<NFREQ;++nfreq)
//	{
//		for(int i=0;i<Nx;++i)
//		{
//			for(int j=0;j<Ny;++j)
//			{
//				Pz[i+j*Nx+nfreq*Nx*Ny] = Ex[i+j*Nx+nfreq*Nx*Ny] * Hy[i+j*Nx+nfreq*Nx*Ny];
//			}
//		}
//	}
//
//	for(int nfreq=0;nfreq<NFREQ;++nfreq)
//	{
//		for(int i=0;i<Nx;++i)
//		{
//			for(int j=0;j<Ny;++j)
//			{
//				TRN[nfreq] += Pz[i+j*Nx+nfreq*Nx*Ny];
//			}
//		}
//	}
//
//	for(int nfreq=0; nfreq<NFREQ;++nfreq)
//	{
//		TRN[nfreq] /= SRCSUM;
//	}
//
//	/********************************* TRANS *******************************************/
//
//
//
//	/********************************* REF *******************************************/
//	complex<double> *Pz1 = ZerosC(Nx*Ny*NFREQ);
//	for(int nfreq=0;nfreq<NFREQ;++nfreq)
//	{
//		for(int i=0;i<Nx;++i)
//		{
//			for(int j=0;j<Ny;++j)
//			{
//				Pz1[i+j*Nx+nfreq*Nx*Ny] = Ey[i+j*Nx+nfreq*Nx*Ny] * Hx[i+j*Nx+nfreq*Nx*Ny];
//			}
//		}
//	}
//
//	for(int nfreq=0;nfreq<NFREQ;++nfreq)
//	{
//		for(int i=0;i<Nx;++i)
//		{
//			for(int j=0;j<Ny;++j)
//			{
//				REF[nfreq] += Pz1[i+j*Nx+nfreq*Nx*Ny];
//			}
//		}
//	}
//
//	for(int nfreq=0; nfreq<NFREQ;++nfreq)
//	{
//		REF[nfreq] /= SRCSUM;
//	}
//
//	/********************************* REF *******************************************/
//
//	for(int nfreq=0; nfreq<NFREQ;++nfreq)
//	{
//		CON[nfreq] = REF[nfreq]+TRN[nfreq];
//	}
//
//
//	double* REFd = new double[NFREQ];
//	double* TRNd = new double[NFREQ];
//	double* CONd = new double[NFREQ];
//
//	for(int nfreq=0; nfreq<NFREQ;++nfreq)
//	{
//		REFd[nfreq] = REF[nfreq].real();
//		TRNd[nfreq] = TRN[nfreq].real();
//		CONd[nfreq] = CON[nfreq].real();
//	}
//
//	OutputPower* Output = new OutputPower("/home/linus/Source/FDTD/FDTDPML/Debug/Results/", 1280, 1024);
//	Output->Output(NFREQ, REFd, TRNd, CONd);

	double* REF = new double[NFREQ];
	double* TRN = new double[NFREQ];
	double* CON = new double[NFREQ];



	// TRANSVERSE WAVE VECTOR EXPANSION
	int *m = new int[Nx];
	int *n = new int[Ny];
	for(int i=0; i<Nx; ++i)
		m[i] = i - floor(Nx/2);
	for(int i=0; i<Ny; ++i)
		n[i] = i - floor(Ny/2);

	double *kx1 = new double[Nx];
	for(int i=0; i<Nx; ++i)
		kx1[i] = -2*M_PI*m[i]/Sx;
	double *ky1 = new double[Ny];
	for(int i=0; i<Ny; ++i)
		ky1[i] = -2*M_PI*n[i]/Sy;

	double *kx = MeshGridx(kx1,Nx);
	double *ky = MeshGridy(ky1, Ny);

	delete kx1;
	delete ky1;

	// LOOP OVER FREQUENCY
	for(int nfreq=0;nfreq<NFREQ;++nfreq)
	{
		// Compute Longitudinal Wave Vector Components
		double lam = c0/FREQ[nfreq];
		double k0 = 2*pi/lam;
		double kzinc = k0*nref;

	  complex<double> *kzR = new complex<double>[Nx*Ny];
	  complex<double> *kzT = new complex<double>[Nx*Ny];

		for(int i=0; i<Nx*Ny; ++i)
		{
			kzR[i] = sqrt(complex<double>(pow((k0*nref),2) - pow(kx[i],2) - pow(ky[i],2)));
			kzT[i] = sqrt(complex<double>(pow((k0*ntrn),2) - pow(kx[i],2) - pow(ky[i],2)));
		}


	  // Interpolate Fields at Origin
    complex<double> *Exr = ZerosC(Nx*Ny);
    for(int k=0; k<Ny; ++k)
    	Exr[0+k*Nx] = (Exref[Nx-1+k*Nx+nfreq*Nx*Ny] + Exref[0+k*Nx+nfreq*Nx*Ny])/complex<double>(2);
    for(int i=1; i<Nx; ++i)
    	for(int k=0; k<Ny; ++k)
    		Exr[i+k*Nx] = (Exref[(i-1)+k*Nx+nfreq*Nx*Ny] + Exref[i+k*Nx+nfreq*Nx*Ny])/complex<double>(2);

		complex<double> *Eyr = ZerosC(Nx*Ny);
    for(int i=0; i<Ny; ++i)
    	Eyr[i+0*Nx] = (Eyref[i+(Ny-1)+nfreq*Nx*Ny] + Eyref[i+0*Nx+nfreq*Nx*Ny])/complex<double>(2);
    for(int i=0; i<Nx; ++i)
    	for(int k=1; k<Ny; ++k)
    		Eyr[i+k*Nx] = (Eyref[i+(k-1)*Nx+nfreq*Nx*Ny] + Eyref[i+k*Nx+nfreq*Nx*Ny])/complex<double>(2);

		complex<double> *Ext = ZerosC(Nx*Ny);
    for(int k=0; k<Ny; ++k)
    	Ext[0+k*Nx] = (Extrn[Nx-1+k*Nx+nfreq*Nx*Ny] + Extrn[1+k*Nx+nfreq*Nx*Ny])/complex<double>(2);
    for(int i=1; i<Nx; ++i)
    	for(int k=0; k<Ny; ++k)
    		Ext[i+5*Nx] = (Extrn[(i-1)+k*Nx+nfreq*Nx*Ny] + Extrn[i+k*Nx+nfreq*Nx*Ny])/complex<double>(2);

		complex<double> *Eyt = ZerosC(Nx*Ny);
    for(int i=0; i<Ny; ++i)
    	Eyt[i+0*Nx] = (Eytrn[i+(Ny-1)+nfreq*Nx*Ny] + Eytrn[i+0*Nx+nfreq*Nx*Ny])/complex<double>(2);
    for(int i=0; i<Nx; ++i)
    	for(int k=1; k<Ny; ++k)
    		Eyt[i+k*Nx] = (Eytrn[i+(k-1)*Nx+nfreq*Nx*Ny] + Eytrn[i+k*Nx+nfreq*Nx*Ny])/complex<double>(2);

    // Normalize to Source
    for(int i=0; i<Nx*Ny; ++i)
    {
    	Exr[i] = Exr[i] / SRC[nfreq];
			Eyr[i] = Eyr[i] / SRC[nfreq];
			Ext[i] = Ext[i] / SRC[nfreq];
			Eyt[i] = Eyt[i] / SRC[nfreq];
    }

    	complex<double>* fftin = ZerosC(Nx*Ny);
      complex<double>* fftout = ZerosC(Nx*Ny);

  		fftw_plan fftplan;
  		fftplan = fftw_plan_dft_2d(Nx, Ny, reinterpret_cast<fftw_complex*>(fftin), reinterpret_cast<fftw_complex*>(fftout), FFTW_FORWARD, FFTW_MEASURE);

  		// Calculate Amplitude Components of Spatial Harmonics
  		int Length = Nx*Ny;
  		CopyArray(Exr, fftin, Length);
	    fftw_execute(fftplan);
	    CopyArray(fftout, Exr, Length);
	    for(int i=0;i<Length;++i)
	    	Exr[i] /= complex<double>(Nx*Ny);

  		CopyArray(Eyr, fftin, Length);
	    fftw_execute(fftplan);
	    CopyArray(fftout, Eyr, Length);
	    for(int i=0;i<Length;++i)
	    	Eyr[i] /= complex<double>(Nx*Ny);

  		CopyArray(Ext, fftin, Length);
	    fftw_execute(fftplan);
	    CopyArray(fftout, Ext, Length);
	    for(int i=0;i<Length;++i)
	    	Ext[i] /= complex<double>(Nx*Ny);

  		CopyArray(Eyt, fftin, Length);
	    fftw_execute(fftplan);
	    CopyArray(fftout, Eyt, Length);
	    for(int i=0;i<Length;++i)
	    	Eyt[i] /= complex<double>(Nx*Ny);

	    fftw_destroy_plan(fftplan);

	    complex<double>* Ezr = ZerosC(Nx*Ny);
	    complex<double>* Ezt = ZerosC(Nx*Ny);
	    // Calculate Longitudinal Components
	    for(int i=0;i<Nx;++i)
	    {
	    	for(int j=0;j<Ny;++j)
	    	{
					Ezr[i+j*Nx] = -(kx[i+j*Nx]*Exr[i+j*Nx] + ky[i+j*Nx]*Eyr[i+j*Nx])/kzR[i+j*Nx];
					Ezt[i+j*Nx] = -(kx[i+j*Nx]*Ext[i+j*Nx] + ky[i+j*Nx]*Eyt[i+j*Nx])/kzT[i+j*Nx];
	    	}
	    }

	    complex<double>* ref = ZerosC(Nx*Ny);
	    complex<double>* trn = ZerosC(Nx*Ny);
	    // Calculate Amplitude of Spatial Harmonics
	    for(int i=0;i<Nx;++i)
	    {
	    	for(int j=0;j<Ny;++j)
	    	{
	    		ref[i+j*Nx] = pow(abs(Exr[i+j*Nx]),2) + pow(abs(Eyr[i+j*Nx]),2) + pow(abs(Ezr[i+j*Nx]),2);
	    		trn[i+j*Nx] = pow(abs(Ext[i+j*Nx]),2) + pow(abs(Eyt[i+j*Nx]),2) + pow(abs(Ezt[i+j*Nx]),2);
	    	}
	    }

	    int x = 0;
    	// Calculate Diffraction Efficiencies
	    for(int i=0;i<Nx;++i)
	    {
	    	for(int j=0;j<Ny;++j)
	    	{
	    		ref[i+j*Nx] = kzR[i+j*Nx].real()/kzinc * ref[i+j*Nx].real();
	    		trn[i+j*Nx] = kzT[i+j*Nx].real()/kzinc * trn[i+j*Nx].real();
	    	}
	    }
x++;
			// Compute Reflectance and Transmittance
			for(int i=0;i<Nx;++i)
			{
				for(int j=0;j<Ny;++j)
				{
					REF[nfreq] += abs(ref[i+j*Nx].real());
					TRN[nfreq] += abs(trn[i+j*Nx].real());
				}
			}
			x++;
	}


	// CALCULATE CONSERVATION OF ENERGY

	for(int i=0;i<NFREQ;++i)
	{
		CON[i] = REF[i] + TRN[i];
	}

	OutputPower* Output = new OutputPower("/home/linus/Source/FDTD/FDTDPML/Debug/Results/", 1280, 1024);
	Output->Output(t, NFREQ, REF, TRN, CON);

}

FDTD::FDTD()
{
	STEPS = 100000;
	NFREQ = 250;
	fmax = 20 * gigahertz;
	fmin = 10 * gigahertz;

  // DEVICE PARAMETERS
  Lx  = 8 * millimeters;
	Ly  = 8 * millimeters;
  r   = 6 * millimeters;
  edge = 1 * millimeters;

//	Lx  = 1.5 * centimeters;
//	Ly  = 1.5 * centimeters;
//	r   = 0.5 * centimeters;
//	d   = 0.5 * centimeters;
//	er1 = 1.0;
//	erm = 10000;
//	erd = 6.0;

  dsub   = 0.0508 * millimeters;
  dm = 0.035 * millimeters;
  d = dsub + dm;
  er1 = 1.0;
  erm = 10000;
  erd = 2.9;

	nmax = sqrt(erd); //maximum refractive index
	NRES = 10;      //wavelength resolution
	NPML = 20; //PML size [zlo zhi]
	BUFZ = (0.05*centimeters);

	Out = new OutputSlice*[3];
	Out[0] = new OutputSlice("/home/linus/Source/FDTD/FDTDPML/Debug/Results/", "Ex", 1280, 1024, 60, 40, 0, 250);
	Out[1] = new OutputSlice("/home/linus/Source/FDTD/FDTDPML/Debug/Results/", "Ey", 1280, 1024, 60, 40, 0, 250);
	Out[2] = new OutputSlice("/home/linus/Source/FDTD/FDTDPML/Debug/Results/", "Ez", 1280, 1024, 60, 40, 0, 250);



}

FDTD::~FDTD()
{
	delete xa;
	delete ya;
	delete za;

	delete sigz;

	delete Materials;
	delete URxx;
	delete URyy;
	delete URzz;

	delete ERxx;
	delete ERyy;
	delete ERzz;

	delete g;
	delete Ex_src;
	delete Ey_src;
	delete Hx_src;
	delete Hy_src;

	delete mHx0;
	delete mHx1;
	delete mHx2;
	delete mHy0;
	delete mHy1;
	delete mHy2;
	delete mHz2;
	delete mHz3;

	delete mDx0;
	delete mDx1;
	delete mDx2;
	delete mDy0;
	delete mDy1;
	delete mDy2;
	delete mDz2;
	delete mDz3;

	delete mEx1;
	delete mEy1;
	delete mEz1;

	delete Hx;
	delete Hy;
	delete Hz;
	delete Dx;
	delete Dy;
	delete Dz;
	delete Ex;
	delete Ey;
	delete Ez;

	delete CEx;
	delete CEy;
	delete CEz;
	delete CHx;
	delete CHy;
	delete CHz;

	delete ICHz;
	delete ICEz;

	delete K;
	delete Exref;
	delete Eyref;
	delete Extrn;
	delete Eytrn;
	delete SRC;

	delete Out;
	delete Outf;
}
