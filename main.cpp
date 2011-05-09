#include"OpenMP.h"
#include"Defines.h"
#include"FDTD.h"
#include<string.h>
#include<iostream>
using namespace std;

int main(int argc, char* argv[])
{
	int Threads = -1;
  int ocount = 0;
	/****************************************************
	 *
	 * Get Command Line Parameters
	 *
	 ****************************************************/
	if (argc == 3)
	{
		Threads = atoi(argv[1]);
		ocount = atoi(argv[2]);
	}
	else
	{
		Threads = 1;
	}

	/****************************************************
	 *
	 * Initialize Parallel
	 *
	 ****************************************************/
	cout << "  Initialize Parallel " << endl;
	omp_set_num_threads(Threads);

	#pragma omp parallel default(none) shared(Threads)
	{
		Threads = omp_get_num_threads();
	}

	cout << "    Running Threads: " << Threads << endl;
	cout << "  End Initialize Parallel" << endl << endl;

  FDTD* engine = new FDTD();
  engine->ocount = ocount;
  engine->init();
  engine->exec();
  delete engine;




}
