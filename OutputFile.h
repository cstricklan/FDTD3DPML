/*
 * OutputFile.h
 *
 *  Created on: May 3, 2011
 *      Author: linus
 */

#ifndef OUTPUTFILE_H_
#define OUTPUTFILE_H_

#include<string>

class OutputFile
{
	std::string _Location;
	int _Slice;


	public:
	void Output(int t, double* Field, std::string FieldName, int Nx, int Ny);
		OutputFile(std::string, int);
		virtual
		~OutputFile();
};

#endif /* OUTPUTFILE_H_ */
