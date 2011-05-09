/*
 * OutputFile.cpp
 *
 *  Created on: May 3, 2011
 *      Author: linus
 */

#include "OutputFile.h"
#include<sstream>
#include<fstream>
#include<iomanip>

using namespace std;
void OutputFile::Output(int t, double* Field, string FieldName, int Nx, int Ny)
{
	stringstream filename;
	filename << this->_Location << FieldName << "-" << this->_Slice << "-" << t;

	int nz = this->_Slice;

	ofstream crossfile;
	crossfile.open (filename.str().c_str());

	crossfile.setf(ios_base::left, ios_base::dec);

	for(int nx=0; nx<Nx;++nx)
	{
		for(int ny=0; ny<Ny;++ny)
		{
			crossfile.setf(ios::fixed, ios::floatfield);
			crossfile << setprecision(10) << setw(15) << Field[nx+ny*Nx+nz*Nx*Ny] << " ";
		}

		crossfile << endl;
	}

	crossfile.close();
}

OutputFile::OutputFile(string Location, int Slice)
{
	this->_Location = Location;
	this->_Slice = Slice;

}

OutputFile::~OutputFile()
{
	// TODO Auto-generated destructor stub
}
