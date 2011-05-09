/*
 * OutputSlice.h
 *
 *  Created on: May 2, 2011
 *      Author: linus
 */

#ifndef OUTPUTSLICE_H_
#define OUTPUTSLICE_H_

#include"Defines.h"
#include<string.h>
#include<mgl/mgl_zb.h>


using namespace std;

class OutputSlice
{
	private:
		int _Width;
		int _Height;
		int _Delay;

		float _RotX;
		float _RotY;
		float _RotZ;

		string _Location;
		string _Field;

		mglGraphZB* graph;

	protected:

	public:
		void Output(int, double*, int,int,int);
		OutputSlice(std::string, std::string, int, int, float, float, float, int);
		virtual ~OutputSlice();
};

#endif /* OUTPUTSLICE_H_ */
