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
#include "util.h"
#include "lvlset.h"


void heapsort(double* R, int* Ni, int* Nj, int** MarkPoint, int n,int nghost){
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
	int Nii, Njj, MM;
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
			MarkPoint[Ni[r]+nghost][Nj[r]+nghost] = r;

			r -= 1;
			if (r == 1) {
				R[1] = RR;
				Ni[1] = Nii;
				Nj[1] = Njj;
				MarkPoint[Nii+nghost][Njj+nghost]=1;
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
				MarkPoint[Ni[i]+nghost][Nj[i]+nghost] = i;
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
				MarkPoint[Ni[i]+nghost][Nj[i]+nghost] = i;
				break;
			}
			R[i]=R[j];
			Ni[i] = Ni[j];
			Nj[i] = Nj[j];
			MarkPoint[Ni[i]+nghost][Nj[i]+nghost] = j;
		}
	}
}

int G_from_flux(int i, int j, double** G, int** Markers, int** MarkPoint, double* Close, int* Closei, int* Closej, int nClose, int nghost, double dx, double dy){
#ifdef DBGGFF
	printf("Attempting %i,%i: with G = %f , Marker = %i",i,j,G[i+nghost][j+nghost],Markers[i+nghost][j+nghost]);
#endif
	double alpha,beta,a,b,A,B,C,Gtemp,ip,im,jp,jm;
	ip = (abs(Markers[i+nghost+1][j+nghost]) > 1) ? 100000.0 : 1.0;
	im = (abs(Markers[i+nghost-1][j+nghost]) > 1) ? 100000.0 : 1.0;
	jp = (abs(Markers[i+nghost][j+nghost+1]) > 1) ? 100000.0 : 1.0;
	jm = (abs(Markers[i+nghost][j+nghost-1]) > 1) ? 100000.0 : 1.0;
	a = fabs(fmax(0.0,fmin(-fabs(G[i+nghost+1][j+nghost]*ip),-fabs(G[i+nghost-1][j+nghost]*im))));
	b = fabs(fmax(0.0,fmin(-fabs(G[i+nghost][j+nghost+1]*jp),-fabs(G[i+nghost][j+nghost-1]*jm))));
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
	Gtemp = (-B+sqrt(B*B-4*A*C))/(2*A);
#ifdef DBGGFF
	printf("\t Found Gtemp = %f \n",Gtemp);
#endif
	if (Gtemp < 0) {
		printf("\n******* NEGATIVE GTEMP *******\n");
	}
	if (isnan(Gtemp)) {
		printf("\n******* NAN GTEMP *******\n");
		Gtemp = fabs(Gtemp);
	}
#ifdef DBGGFF
	printf("Alpha = %f, Beta = %f,\n", alpha,beta);
	printf("A = %f, B = %f, C = %f \n", A, B, C);
	printf("a = %f, b = %f, ip = %f, im = %f, jp = %f, jm = %f \n", a,b,ip,im,jp,jm);
#endif

	if(abs(Markers[i+nghost][j+nghost]) > 2){ // Point was far
		Close[nClose] = fabs(Gtemp);
		G[i+nghost][j+nghost] = (G[i+nghost][j+nghost] >= 0) ? Gtemp : -Gtemp;
		Markers[i+nghost][j+nghost] = (Markers[i+nghost][j+nghost] > 0) ? 2 : -2;
		MarkPoint[i+nghost][j+nghost] = nClose-1;
		Closei[nClose] = i;
		Closej[nClose] = j;
#ifdef DBGGFF
		printf("Updated G=%f at %i,%i \n\n",Gtemp,i,j);
#endif
		return nClose+1;
	} else { //Point was close
		G[i+nghost][j+nghost] = (G[i+nghost][j+nghost] >= 0) ? Gtemp : -Gtemp;
		Close[MarkPoint[i+nghost][j+nghost]] = fabs(Gtemp);
#ifdef DBGGFF
		printf("Updated Previous G=%f at %i,%i \n\n",Gtemp,i,j);
#endif
		return nClose;
	}
}

void get_alpha(double** G, double** alpha,  int nghost, double dx, double dy,int nx, int ny){
	double a = (dx+dy)/3.0;

#pragma omp parallel for
	for(int i = 0; i<nx; i++){
		for(int j = 0; j<ny; j++){
			alpha[i+nghost][j+nghost] = alpha_H(G[i+nghost][j+nghost], a);
		}
	}
}
double get_alpha_sum(double** alpha,  int nghost, double dx, double dy,int nx, int ny){
	double sum = 0.0;
	for(int i = 0; i<nx; i++){
		for(int j = 0; j<ny; j++){
			sum+= alpha[i+nghost][j+nghost]*dx*dy;
		}
	}
	return sum;
}
double get_alpha_diff(double** alpha, double** alpha2,  int nghost, double dx, double dy,int nx, int ny){
	double sum = 0.0;
	double a = (dx+dy)/2.0;
	for(int i = 0; i<nx; i++){
		for(int j = 0; j<ny; j++){
			sum+= fabs(alpha[i+nghost][j+nghost]-alpha2[i+nghost][j+nghost])*dx*dy;
		}
	}
	return sum;
}

V2D vof_alpha_grad(double** alpha, int i, int j, int nghost, double dx, double dy){
	V2D grad;
	double mag;
	double eps = 0.00000000001;
	grad.x = 0.125*(alpha[i+nghost+1][j+nghost+1]-alpha[i+nghost-1][j+nghost+1])/dx*2.0+
		     0.250*(alpha[i+nghost+1][j+nghost  ]-alpha[i+nghost-1][j+nghost  ])  /dx*2.0+
		     0.125*(alpha[i+nghost+1][j+nghost-1]-alpha[i+nghost-1][j+nghost-1])/dx*2.0;

	grad.y = 0.125*(alpha[i+nghost+1][j+nghost+1]-alpha[i+nghost+1][j+nghost-1])/dy*2.0+
		     0.250*(alpha[i+nghost  ][j+nghost+1]-alpha[i+nghost  ][j+nghost-1])/dy*2.0+
		     0.125*(alpha[i+nghost-1][j+nghost+1]-alpha[i+nghost-1][j+nghost-1])/dy*2.0;
	mag = fmax(sqrt(grad.x*grad.x+grad.y*grad.y),eps);
	grad.x = -grad.x/mag;
	grad.y = -grad.y/mag;
	if(mag == eps){
		grad.y = 1.0;
		grad.x = 0.0;
	}
	return grad;
}

V2D rot90(V2D v){
	double temp;
	temp = v.x;
	v.x = -v.y;
	v.y = temp;
	return v;
}

V2D rotn90(V2D v){
	double temp;
	temp = v.x;
	v.x = v.y;
	v.y = -temp;
	return v;
}

Wall cycle_wall(Wall w){
	double tmp;
	tmp = w.l;
	w.l = w.d;
	w.d = w.r;
	w.r = w.u;
	w.u = tmp;
	return w;
}
void vof_get_surfaces(double** alpha,double** darr, V2D** marr,lseg** segs, int nghost, double dx, double dy,int nx, int ny){
	for(int i = 0; i<nx; i++){
		for(int j = 0; j<ny; j++){
			vof_get_surface(alpha,darr,marr,segs,i,j,nghost,dx,dy);
		}
	}
}

void vof_construct_ds(double** alpha,double** darr, V2D** marr, int nghost, double dx, double dy,int nx, int ny){
	for(int i = 0; i<nx; i++){
		for(int j = 0; j<ny; j++){
			double a = alpha[i+nghost][j+nghost];
			V2D m = vof_alpha_grad(alpha,i,j,nghost,dx,dy);
			marr[i+nghost][j+nghost] = m;
			darr[i+nghost][j+nghost] = vof_d_from_alpha(a,m,dx,dy);
		}
	}
}
double H(double x){
	if (x<0) {
		return 0.0;
	} else {
		return 1.0;
	}
}
double vof_alpha_from_surf(double d, V2D m, double hx, double hy){
	double alpha;
	m.y = fabs(m.y);
	m.x = fabs(m.x);
	if (d < 0.0){
		printf("Negative D!*******\n");
	}
	if(d == 0.0){
		alpha = 0.0;
	}else if ((d >= m.x*hx+m.y*hy) && (m.x > 0.0 || m.y > 0.0)){
		return 1.0;
	} else if (m.x < 0.00000001 || m.y < 0.00000001){
		if (m.y < 0.00000001) {
			alpha = d/hx;
		} else if(m.x < 0.00000001){
			alpha = d/hy;
		}
		return alpha;
	}else {
		alpha = d*d/(2.0*m.x*m.y*hx*hy)*(1.0-H(d-m.x*hx)*((d-m.x*hx)/d)*((d-m.x*hx)/d)
				-H(d-m.y*hy)*((d-m.y*hy)/d)*((d-m.y*hy)/d));
	}
	return alpha;
}
double vof_d_from_alpha(double a, V2D m, double hx, double hy){
	double a1, a2, d,d1,d2,A,B,C;
	double at1,at2;
	m.y = fabs(m.y);
	m.x = fabs(m.x);
	d1 = hy*m.y;
	d2 = hx*m.x;
	a1 = vof_alpha_from_surf(d1,m,hx,hy);
	a2 = vof_alpha_from_surf(d2,m,hx,hy);
	/*
	printf("alpha = %f, a1,a2 = %f,%f\n",a,a1,a2);
	printf("MxMy = %f,%f\n",m.x,m.y);
	printf("d1d2 = %f,%f\n",d1,d2);
	printf("hxhy = %f,%f\n",hx,hy);
	*/
	if(a1 < 0 || a2 < 0){
		printf("BAD ALPHA12 %f,%f\n",a1,a2);
	}
	if (a <= 0 || (hx == 0 || hy == 0)){
		d = 0.0;
	} else if (a >= (1.0-0.00000001)){
		d = hx*m.x+hy*m.y;
	} else if ((m.x < 0.00000001 || m.y < 0.00000001) ) {
		//printf("rectangle \n");
		if (m.x < 0.00000001){
			d  = a*hx;
		} else if (m.y < 0.00000001){
			d  = a*hy;
		}
	} else if (a <= fmin(a1,a2)){
		//printf("simple \n");
		d = sqrt(2.0*m.x*m.y*hx*hy*a);
	} else if (a < a1) {
		//printf("a < a1 \n");
		B = +2.0*m.x*hx/(2.0*m.x*m.y);
		C = -m.x*m.x*hx*hx/(2.0*m.x*m.y)-a*hx*hy;
		d = -C/B;

	} else if (a < a2) {
		//printf("a < a1 \n");
		B = +2.0*m.y*hy/(2.0*m.x*m.y);
		C = -m.y*m.y*hy*hy/(2.0*m.x*m.y)-a*hx*hy;
		d = -C/B;

	} else {
		//printf("quadratic \n");
		A = 1.0/(2.0*m.x*m.y)-1.0/(2.0*m.x*m.y)-1.0/(2.0*m.x*m.y);
		B = +2.0*m.x*hx/(2.0*m.x*m.y)+2.0*m.y*hy/(2.0*m.x*m.y);
		C = -m.x*m.x*hx*hx/(2.0*m.x*m.y)-m.y*m.y*hy*hy/(2.0*m.x*m.y)-a*hx*hy;
		if((B*B-4.0*A*C) > 0){
			double r1,r2;
			r1 = sqrt(B*B-4.0*A*C);
			r2 = (-B-r1)/(2.0*A);
			r1 = (-B+r1)/(2.0*A);

			at1 = r1*r1/(2.0*m.x*m.y*hx*hy)*(1-((r1-m.x*hx)/r1)*((r1-m.x*hx)/r1)-((r1-m.y*hy)/r1)*((r1-m.y*hy)/r1));
			at2 = r2*r2/(2.0*m.x*m.y*hx*hy)*(1-((r2-m.x*hx)/r2)*((r2-m.x*hx)/r2)-((r2-m.y*hy)/r2)*((r2-m.y*hy)/r2));
			if ((r1 >= fmax(d1,d2) && r1 <= (m.x*hx+m.y*hy)) && (r2 >= fmax(d1,d2) && r2 <= (m.x*hx+m.y*hy))){
				printf("BOTH ROOTS FIT!!!!\n");
				d = r1;
			} else if ((r2 >= fmax(d1,d2) && r2 <= (m.x*hx+m.y*hy))){
				//printf("R2!!!!\n");
				d = r2;
			} else if ((r1 >= fmax(d1,d2) && r1 <= (m.x*hx+m.y*hy))){
				//printf("R1!!!!\n");
				d = r1;
			} else {
				printf("NO ROOTS FIT (1)!!!!******\n");
				d = r1;
			}
		} else {
			printf("COMPLEX ROOTS IN D (1) %f,%f,%f!!!!******\n",a,a1,a2);
			printf("alpha = %f, a1,a2 = %f,%f\n",a,a1,a2);
			printf("MxMy = %f,%f\n",m.x,m.y);
			printf("d1d2 = %f,%f\n",d1,d2);
			printf("hxhy = %f,%f\n",hx,hy);
			printf("A,B,C = %f, %f, %f\n",A,B,C);
			d = a*hx;
		}
	}
	double acheck;
	acheck = vof_alpha_from_surf(d,m,hx,hy);

	/*if (fabs(a-acheck) > 0.01) {
		printf("BAD ALPHA a - a1,a2 - acheck %f- %f,%f - %f for d = %f\n",a,a1,a2,acheck,d);
		printf("alpha = %f, a1,a2 = %f,%f\n",a,a1,a2);
		printf("MxMy = %f,%f\n",m.x,m.y);
		printf("d1d2 = %f,%f\n",d1,d2);
		printf("hxhy = %f,%f\n",hx,hy);
		printf("A,B,C = %f, %f, %f\n",A,B,C);
	}
	if (d > m.x*hx+m.y*hy && (a > 0.0 && a < 1.0)){
		printf("bad d=%f, a,a1,a2= %f,%f,%f!!!!******\n",d,a,a1,a2);
	}*/


	return d;
}
void vof_get_surface(double** alpha, double** darr, V2D** marr, lseg** segs, int i, int j, int nghost, double dx, double dy){
	V2D m =	vof_alpha_grad(alpha,i,j,nghost,dx,dy);
	marr[i+nghost][j+nghost] = m;
	Wall w;
	lseg l;
	double a, a1, a2, d,d1,d2;
	double at1,at2;
	a = alpha[i+nghost][j+nghost];
	int mark;
	int rot = 0;
	//d = vof_d_from_alpha(a,m,dx,dy);
	d = darr[i+nghost][j+nghost];
	while(m.x < 0.0 || m.y < 0.0){
		rot += 1;
		m = rot90(m);
	}
	if (a < 1.0 && a > 0.0){
		// found d, get fill
		if (m.y == 0 || m.x == 0.0){
			if (m.y == 0) {
				w.l = d/dy;
				w.r = d/dy;
				w.u = 0.0;
				w.d = 1.0;
			} else {
				w.l = 1.0;
				w.r = 0.0;
				w.u = d/dx;
			    w.d = d/dx;
			}
		} else {
			if (d/m.y < dy){//doesn't fill left wall
				w.l = d/m.y/dy;
				w.u = 0.0;
			} else { //over left wall
				w.l = 1.0;
				w.u = (d-m.y*dy)/(m.x*dx);
			}
			if (d/m.x < dx){ //doesn't fill bottom wall
				w.d = d/m.x/dx;
				w.r = 0.0;
			} else { // fills bottom wall
				w.d = 1.0;
				w.r = (d-m.x*dx)/(m.y*dx);
			}
		}
		l.start.x = w.u*dx-dx/2.0;
		l.start.y = w.l*dy-dy/2.0;
		l.end.x = w.d*dx-dx/2.0;
		l.end.y = w.r*dy-dy/2.0;
		while (rot > 0) { //rotate wall back to normal
			w = cycle_wall(w);
			l.start = rotn90(l.start);
			l.end = rotn90(l.end);
			rot -= 1;
		}
		l.start.x = l.start.x+dx/2.0;
		l.start.y = l.start.y+dy/2.0;
		l.end.x = l.end.x+dx/2.0;
		l.end.y = l.end.y+dy/2.0;
	} else if (a >= 1){
		w.l = 1.0;
		w.u = 1.0;
		w.r = 1.0;
		w.d = 1.0;
		l.start.x = dx;l.start.y = dx;
		l.end.x=dy;l.end.y=dy;
	} else {
		w.l = 0.0;
		w.u = 0.0;
		w.r = 0.0;
		w.d = 0.0;
		l.start.x = 0.0;l.start.y = 0.0;
		l.end.x=0.0;l.end.y=0.0;
	}
	l.start.x = l.start.x+(i*dx);l.start.y = l.start.y+(j*dy);
	l.end.x = l.end.x+(i*dx);l.end.y = l.end.y+(j*dy);
	segs[i+nghost][j+nghost] = l;
	darr[i+nghost][j+nghost] = d;

}
void vof_get_fluxes(double** u, double** v, double** alpha,double** darr, V2D** marr, double dt, int nghost, double dx, double dy, int nx, int ny){
	double dalpha,d,dp,uf,vf,a;
	double hx,hy;
	V2D m;
	// flux u
	hy = dy;
	vof_construct_ds(alpha,darr,marr,nghost, dx, dy, nx, ny);
	for(int i = 0; i<nx; i++){
		for(int j = 0; j<ny; j++){
			uf = u[i+nghost][j+nghost];
			hx = fabs(uf*dt);
			if(uf > 0){
				d = darr[i+nghost-1][j+nghost];
				m = marr[i+nghost-1][j+nghost];
			} else {
				d = darr[i+nghost][j+nghost];
				m = marr[i+nghost][j+nghost];
			}
			if (uf*m.x > 0){
				dp = fmax(0.0,d-fabs(m.x)*(dx-hx));
			} else {
				dp = d;
			}
			dalpha = vof_alpha_from_surf(dp,m,hx,hy);
			//printf("dalpha = %f,d = %f, dp = %f, uf = %f, dt = %f, m = %f,%f hx = %f, hy = %f     -------   %f\n",dalpha,d,dp,uf,dt,m.x,m.y,hx,hy,dalpha*hx*hy/(dx*dy));
			dalpha = dalpha*hx*hy/(dx*dy);
			if(uf > 0) {
				//dalpha = fmin(dalpha*hx*hy/(dx*dy),alpha[i+nghost-1][j+nghost]);
				alpha[i+nghost-1][j+nghost] = alpha[i+nghost-1][j+nghost]-dalpha;
				alpha[i+nghost][j+nghost] = alpha[i+nghost][j+nghost]+dalpha;
			} else if (uf < 0){
				//dalpha = fmin(dalpha*hx*hy/(dx*dy),alpha[i+nghost][j+nghost]);
				alpha[i+nghost-1][j+nghost] = alpha[i+nghost-1][j+nghost]+dalpha;
				alpha[i+nghost][j+nghost] = alpha[i+nghost][j+nghost]-dalpha;
			}
		}
	}
	vof_construct_ds(alpha,darr, marr,nghost, dx, dy, nx, ny);

	// flux v
	hx = dx;
	for(int i = 0; i<nx; i++){
		for(int j = 0; j<ny; j++){
			vf = v[i+nghost][j+nghost];
			hy = fabs(vf*dt);
			if(vf > 0){
				d = darr[i+nghost][j+nghost-1];
				m = marr[i+nghost][j+nghost-1];
			} else {
				d = darr[i+nghost][j+nghost];
				m = marr[i+nghost][j+nghost];
			}
			if (vf*m.y > 0){
				dp = fmax(0.0,d-fabs(m.y)*(dy-hy));
			} else {
				dp = d;
			}
			dalpha = vof_alpha_from_surf(fmax(dp,0.0),m,hx,hy);
			dalpha = dalpha*hx*hy/(dx*dy);
			if(vf > 0) {
				alpha[i+nghost][j+nghost-1] = alpha[i+nghost][j+nghost-1]-dalpha;
				alpha[i+nghost][j+nghost] = alpha[i+nghost][j+nghost]+dalpha;
			} else if(vf < 0) {
				alpha[i+nghost][j+nghost-1] = alpha[i+nghost][j+nghost-1]+dalpha;
				alpha[i+nghost][j+nghost] = alpha[i+nghost][j+nghost]-dalpha;
			}
		}
	}
	vof_construct_ds(alpha,darr, marr,nghost, dx, dy, nx, ny);
}





