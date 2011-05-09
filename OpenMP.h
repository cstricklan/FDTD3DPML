/*
 * OpenMP.h
 *
 *  Created on: Feb 21, 2010
 *      Author: linus
 */

#ifndef OPENMP_H_
#define OPENMP_H_

#ifdef _OPENMP
	#include <omp.h>
#else
	#define omp_get_thread_num() 0
	#define omp_get_num_threads() 0
  #define omp_set_num_threads(int) 0
	#define omp_get_num_procs() 0
#endif



#endif /* OPENMP_H_ */
