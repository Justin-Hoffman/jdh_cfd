/*
 * util.h
 *
 *  Created on: Feb 20, 2014
 *      Author: justin
 */

#ifndef UTIL_H_
#define UTIL_H_
void heapsort(double* R, int* Ni, int* Nj, int n );
int G_from_flux(int i, int j, double** G, int** Markers,double* Close, int* Closei, int* Closej, int nClose, int nghost, double dx, double dy);

#endif /* UTIL_H_ */
