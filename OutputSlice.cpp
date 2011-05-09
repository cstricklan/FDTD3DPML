/*
 * OutputSlice.cpp
 *
 *  Created on: May 2, 2011
 *      Author: linus
 */

#include "OutputSlice.h"
#include <sstream>

using namespace std;

void OutputSlice::Output(int Step, double* Points, int Nx, int Ny, int Nz)
{
//	stringstream filename;
//	filename << this->_Location << this->_Field <<  "-" << Nx <<  "-" << Ny <<  "-" << Nz << "-" << Step+1 << ".png";

	mglData DataSet(Nx, Ny, Nz);

	#pragma omp for
	for(int nx=0; nx<Nx; ++nx)
		for(int ny=0; ny<Ny; ++ny)
			for(int nz=0; nz<Nz; ++nz)
				DataSet.a[nx+ny*Nx+nz*Nx*Ny] = Points[nx+ny*Nx+nz*Nx*Ny];

	std::stringstream Title;
	Title << "Step: " << Step;


	graph->NewFrame();
	graph->Text(mglPoint(0,1.3,1),Title.str().c_str());
//	mglGraphZB graph(this->_Width, this->_Height);
	graph->Cmin =  DataSet.Maximal();
	graph->Cmax =  DataSet.Minimal();
	graph->Org = mglPoint(0,0,0);
	graph->Alpha(true);
	graph->Light(true);
	graph->Light(0,mglPoint(1,0,-1));
	graph->Rotate(this->_RotX, this->_RotY, this->_RotZ);
	//graph.Box();
	graph->DensA(DataSet);
	//graph.Axis();
	//graph->Colorbar();
	graph->EndFrame();

//	graph.WritePNG(filename.str().c_str());
}


OutputSlice::OutputSlice(std::string Location, std::string Field, int Width, int Height, float RotX, float RotY, float RotZ, int Delay)
{
	this->_Location = Location;
	this->_Field = Field;
	this->_Width = Width;
	this->_Height = Height;
	this->_RotX = RotX;
	this->_RotY = RotY;
	this->_RotZ = RotZ;
	this->_Delay = Delay;

	stringstream filename;
	filename << this->_Location << this->_Field << ".gif";


	this->graph = new mglGraphZB(this->_Width, this->_Height);
	this->graph->StartGIF(filename.str().c_str(), this->_Delay);
}


OutputSlice::~OutputSlice()
{
	this->graph->CloseGIF();
	delete this->graph;
}
