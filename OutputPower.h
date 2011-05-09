/*
 * OutputPower.h
 *
 *  Created on: May 2, 2011
 *      Author: linus
 */

#ifndef OUTPUTPOWER_H_
#define OUTPUTPOWER_H_

#include"Defines.h"
#include<string.h>

class OutputPower
{
		int _Width;
		int _Height;

		std::string _Location;

	public:

		void Output(int, int, double*, double*,double*);

		OutputPower(std::string, int, int);
		virtual
		~OutputPower();
};

#endif /* OUTPUTPOWER_H_ */
