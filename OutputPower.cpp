/*
 * OutputPower.cpp
 *
 *  Created on: May 2, 2011
 *      Author: linus
 */

#include "OutputPower.h"
#include <sstream>
#include<mgl/mgl_zb.h>

void OutputPower::Output(int t, int NFREQ, double* REF, double* TRN,double* CON)
{
	std::stringstream filename;
	filename << this->_Location << "Power-" << t << ".jpg";

	mglData dTRN(NFREQ);
	mglData dREF(NFREQ);
	mglData dCON(NFREQ);

	for(int f=0;f<NFREQ;++f)
	{
		dTRN.a[f] = REF[f];
		dREF.a[f] = TRN[f];
		dCON.a[f] = CON[f];
	}

		mglGraphZB graph(this->_Width, this->_Height);
//		graph.Cmin =  dCON.Maximal();
//		graph.Cmax =  dCON.Minimal();
		graph.AutoOrg = true;

		graph.Box();
		graph.Plot(dCON);
		graph.Plot(dREF);
		graph.Plot(dTRN);
		//graph.Axis();

		graph.WriteJPEG(filename.str().c_str());
}


OutputPower::OutputPower(std::string Location, int Width, int Height)
{
	this->_Location = Location;
	this->_Width = Width;
	this->_Height = Height;
}

OutputPower::~OutputPower()
{
	// TODO Auto-generated destructor stub
}
