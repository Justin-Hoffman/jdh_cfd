/*
 * util.c
 *
 *  Created on: Feb 20, 2014
 *      Author: justin
 */
#define M_PI 3.14159265358979323846
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>
#include <unistd.h>
void heapsort(double* R, int* Ni, int* Nj, int n ){
	R = R-1;
	Ni = Ni-1;
	Nj = Nj-1;
	if (n < 2) {
		return;
	}
	int l = n/2+1;
	int r = n;
	int i,j;
	double RR = 0;
	int Nii, Njj;
	while(1){
		if (l > 1) {
			l = l-1;
			RR = R[l];
			Nii = Ni[l];
			Njj = Nj[l];
		} else {
			RR = R[r];
			Nii = Ni[r];
			Njj = Nj[r];
			R[r] = R[1];
			Ni[r] = Ni[1];
			Nj[r] = Nj[1];
			r -= 1;
			if (r == 1) {
				R[1] = RR;
				Ni[1] = Nii;
				Nj[1] = Njj;
				break;
			}
		}
		j = l;
		while(1){
			i = j;
			j *= 2;
			if (j>r){
				R[i] = RR;
				Ni[i] = Nii;
				Nj[i] = Njj;
				break;
			}
			if (j<r){
				if(R[j]<R[j+1]){
					j=j+1;
				}
			}
			if(RR>=R[j]){
				R[i]=RR;
				Ni[i] = Nii;
				Nj[i] = Njj;
				break;
			}
			R[i]=R[j];
			Ni[i] = Ni[j];
			Nj[i] = Nj[j];
		}
	}
}

int G_from_flux(int i, int j, double** G, int** Markers,double* Close, int* Closei, int* Closej, int nClose, int nghost, double dx, double dy){
	printf("Attempting %i,%i:",i,j);
	double alpha,beta,a,b,A,B,C,Gtemp,ip,im,jp,jm;
	ip = (abs(Markers[i+nghost+1][j+nghost]) > 1) ? 100000.0 : 1.0;
	im = (abs(Markers[i+nghost-1][j+nghost]) > 1) ? 100000.0 : 1.0;
	jp = (abs(Markers[i+nghost][j+nghost+1]) > 1) ? 100000.0 : 1.0;
	jm = (abs(Markers[i+nghost][j+nghost-1]) > 1) ? 100000.0 : 1.0;
	a = fmin(fabs(G[i+nghost+1][j+nghost]*ip),fabs(G[i+nghost-1][j+nghost]*im));
	b = fmin(fabs(G[i+nghost][j+nghost+1]*jp),fabs(G[i+nghost][j+nghost-1]*jm));

	if(ip == 1.0 || im  == 1.0 ){
		alpha = 1.0;
	} else {
		alpha = 0.0;
	}
	if(jp == 1.0 || jm  == 1.0 ){
		beta = 1.0;
	} else {
		beta = 0.0;
	}
	A = alpha+beta;
	B = -2*(alpha*a+beta*b);
	C = alpha*a*a+beta*b*b-dx*dy;
	Gtemp = (-B+sqrt(B*B-4*A*C))/(2*A);\
	printf("\t Found Gtemp\n");
	if (Gtemp < 0) {
		printf("\n******* NEGATIVE GTEMP *******\n");
		printf("A = %f, B = %f, C = %f \n\n", A, B, C);
	}
	if (isnan(Gtemp)) {
		printf("\n******* NAN GTEMP *******\n");
		printf("Alpha = %f, Beta = %f,\n\n", alpha,beta);
		printf("A = %f, B = %f, C = %f \n\n", A, B, C);
		Gtemp = fabs(Gtemp);
	}


	if(abs(Markers[i+nghost][j+nghost]) > 2){ // Point was far
		Close[nClose] = Gtemp;
		G[i+nghost][j+nghost] = (G[i+nghost][j+nghost] >= 0) ? Gtemp : -Gtemp;
		Markers[i+nghost][j+nghost] = (Markers[i+nghost][j+nghost] > 0) ? 2 : -2;
		Closei[nClose] = i;
		Closej[nClose] = j;
		printf("Updated G=%f at %i,%i \n",Gtemp,i,j);
		return nClose+1;
	} else { //Point was close
		int ii = 0;
		while(Closei[ii]!=i && Closej[ii] != j){
			ii++;
		}
		G[i+nghost][j+nghost] = (G[i+nghost][j+nghost] >= 0) ? Gtemp : -Gtemp;
		Close[ii] = Gtemp;
		printf("Updated Previous G=%f at %i,%i \n",Gtemp,i,j);
		return nClose;
	}
}



