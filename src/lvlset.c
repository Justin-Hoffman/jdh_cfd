/*
 * lvlset.c
 *
 *  Created on: Jan 30, 2014
 *      Author: justin
 */
#define M_PI 3.14159265358979323846
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>
#include <unistd.h>
#include "lvlset.h"
#include "slv.h"

double dist_from_arc(double x, double y, double xc, double yc, double R, double t0, double t1){
	double tc = atan2((y-yc),(x-xc)); // Angle of current pt
	if (tc < 0.0) {
		tc = 2.0*M_PI+tc;
	}
	if ((tc > t0) && (tc < t1)){ //Point within blanked arc
		double xp0 = xc + R*cos(t0);
		double yp0 = yc + R*sin(t0);
		double xp1 = xc + R*cos(t1);
		double yp1 = yc + R*sin(t1);
		return -fmin(dist_from_pt(x,y,xp0,yp0), dist_from_pt(x,y,xp1,yp1));
	} else { //Point outside of blanked arc
		return dist_from_circ(x,y,xc,yc,R);
	}
}
double dist_from_pt(double x, double y, double xc, double yc){
	return sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc));
}
double dist_from_circ(double x, double y, double xc, double yc, double R){
	return(-dist_from_pt(x,y,xc,yc) + R);
}
double dist_from_line(double x, double y, double x0, double y0, double x1, double y1){
	double xo0 = x0-(y1-y0); //Orthogonal vector
	double yo0 = y0+(x1-x0); //Orthogonal vector
	double xo1 = x1-(y1-y0); //
	double yo1 = y1+(x1-x0); //

	double a0 = vec_angle((x-x0),(y-y0), (xo0-x0),(yo0-y0)); // x -> xo0 with x0 as vertex
	double a1 = vec_angle((x-x0),(y-y0), (x1-x0),(y1-y0)); // x -> x1 with x0 as vertex
	double a2 = vec_angle((x-x1),(y-y1), (xo1-x1),(yo1-y1)); // x -> x0 with x1 ax vertex;

	if ((a0>=0.0 && (a1 <= 0.0 && a2 <= 0.0) ) || (a0 >= 0.0 && a1 >= 0.0 && a2 < 0.0) ){
		return dist_from_vec((x-x0),(y-y0),(x1-x0),(y1-y0));
	} else {
		return fmin(dist_from_pt(x,y,x0,y0),dist_from_pt(x,y,x1,y1));
	}
}
double dist_from_vec(double x, double y, double x0, double y0){
	if ((fabs(x0) > 0.0 || fabs(y0) > 0.0)){ //If I won't divide by zero
		double sc = vec_dot(x,y,x0,y0)/(vec_mag(x0,y0)*vec_mag(x0,y0)); //vector scale
		return (dist_from_pt(x,y,x0*sc,y0*sc));
	} else {
		if ((fabs(x0) > 0.0 || fabs(y0) > 0.0)) {
			return vec_mag(x0,y0);
		} else {
			return vec_mag(x,y);
		}
	}
}
double vec_dot(double x0, double y0, double x1, double y1){
	return (x0*x1+y0*y1);
}
double vec_cross(double x0, double y0, double x1, double y1){
	return(x0*y1-x1*y0);
}
double vec_angle(double x0, double y0, double x1, double y1){
	double mag1 = vec_mag(x0,y0);
	double mag2 = vec_mag(x1,y1);
	if (mag1 != 0 && mag2 != 0){
		return asin(vec_cross(x0,y0,x1,y1)/(mag1*mag2));
	} else {
		return 0.0;
	}
}
double vec_mag(double x0, double y0){
	return sqrt(x0*x0+y0*y0);
}

double dist_from_box(double x, double y, double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3){
	double angle0 = vec_angle((x-x0),(y-y0), (x1-x0),(y1-y0));
	double angle1 = vec_angle((x-x1),(y-y1), (x2-x1),(y2-y1));
	double angle2 = vec_angle((x-x2),(y-y2), (x3-x2),(y3-y2));
	double angle3 = vec_angle((x-x3),(y-y3), (x0-x3),(y0-y3));
	double C = M_PI/2.0;
	double d0 = dist_from_line(x,y,x0,y0,x1,y1);
	double d1 = dist_from_line(x,y,x1,y1,x2,y2);
	double d2 = dist_from_line(x,y,x2,y2,x3,y3);
	double d3 = dist_from_line(x,y,x3,y3,x0,y0);
	double m0 = fmin(d0,d1);
	double m1 = fmin(d2,d3);
	if (angle0 >= 0.0 && angle1 >= 0.0 && angle2 >= 0.0 && angle3 >= 0.0){ // If all cross product angles are greater than 0.0 then point HAS TO BE IN BOXX
		return -fmin(m0,m1);
	} else { // Must be OUTSIDE BOX
		return fmin(m0,m1);
	}
}
double dist_from_notch(double x, double y, double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3){ //Open sided box
	/*(--++)---(+-++)---(+--+)-
	 *-------------------------
	 *-------1========2--------
	 *(-+++)-|--------|--------
	 *-------|-(++++)-|-(++-+)-
	 *-------|--------|-------
	 *-------|--------|--------
	 *-------0--------3--------
	 *-------------------------
	 *(-++-)--(+++-)----(++--)-
	 *-------------------------
	 */
	// = or | notch boundary
	// * location of pt


	double angle0 = vec_angle((x-x0),(y-y0), (x1-x0),(y1-y0));
	double angle1 = vec_angle((x-x1),(y-y1), (x2-x1),(y2-y1));
	double angle2 = vec_angle((x-x2),(y-y2), (x3-x2),(y3-y2));
	double angle3 = vec_angle((x-x3),(y-y3), (x0-x3),(y0-y3));
	double C = M_PI/2.0;
	double d0 = dist_from_line(x,y,x0,y0,x1,y1);
	double d1 = dist_from_line(x,y,x1,y1,x2,y2);
	double d2 = dist_from_line(x,y,x2,y2,x3,y3);
	double d3 = dist_from_line(x,y,x3,y3,x0,y0);
	double m0 = fmin(d0,d1);

	/*
	printf("X%i = %f, Y%i = %f\n",0,x0,0,y0);
	printf("X%i = %f, Y%i = %f\n",1,x1,1,y1);
	printf("X%i = %f, Y%i = %f\n",2,x2,2,y2);
	printf("X%i = %f, Y%i = %f\n",3,x3,3,y3);

	printf("Angles: %f\t %f\t %f\t %f\t\n",angle0, angle1, angle2, angle3);
	printf("Distances: %f\t %f\t %f\t %f\t\n",d0, d1, d2, d3);
	*/

	if (angle0 >= 0.0 && angle1 >= 0.0 && angle2 >= 0.0 && angle3 >= 0.0){ // If all cross product angles are greater than 0.0 then point HAS TO BE IN BOX
		return -fmin(m0,d2);
	} else if (angle0 <= 0.0 && angle1 >= 0.0 && angle2 >= 0.0 && angle3 <= 0.0){ // Lower Left
		return -d0;
	} else if (angle0 <= 0.0 && angle1 >= 0.0 && angle2 >= 0.0 && angle3 >= 0.0){ // Left
		return d0;
	} else if (angle0 <= 0.0 && angle1 <= 0.0 && angle2 >= 0.0 && angle3 >= 0.0){ // Top Left
		return d0;
	} else if (angle0 >= 0.0 && angle1 <= 0.0 && angle2 >= 0.0 && angle3 >= 0.0){ // Top
		return d1;
	} else if( angle0 >= 0.0 && angle1 <= 0.0 && angle2 <= 0.0 && angle3 >= 0.0 ) { // Top Right
		return d1;
	} else if (angle0 >= 0.0 && angle1 >= 0.0 && angle2 <= 0.0 && angle3 >= 0.0){ // Right
		return d2;
	} else if (angle0 >= 0.0 && angle1 >= 0.0 && angle2 <= 0.0 && angle3 <= 0.0){ // Lower Right
		return -d3;
	} else { //Lower
		return -fmin(d0, d3);
	}
}

double dist_from_zalesak(double x,double y,double xc,double yc, double R, double cwidth, double cdepth, double cangle){
	double halfcut = atan2(cwidth/2.0,R); //Trig to identify half angle of notch
	double t0 = fmod((cangle-halfcut),(2.0*M_PI)); //Arc angle at start of cut
	double t1 = fmod((cangle+halfcut),(2.0*M_PI)); //Arc angle at end of cut
	double theta = fmod(atan2((x-xc),(y-yc)),(2.0*M_PI)); // Angle of current pt

	if (dist_from_pt(x,y,xc,yc) < R){ // Inside Disk
		double xn0 = xc+R*cos(t0); double yn0 = yc+R*sin(t0);
		double xn1 = xn0-cos(cangle)*cdepth; double yn1 = yn0-sin(cangle)*cdepth;
		double xn3 = xc+R*cos(t1); double yn3 = yc+R*sin(t1);
		double xn2 = xn3-cos(cangle)*cdepth; double yn2 = yn3-sin(cangle)*cdepth;

		double dist2notch = dist_from_notch(x, y, xn0, yn0, xn1, yn1, xn2, yn2, xn3, yn3);
		double dist2arc = dist_from_arc(x,y,xc,yc,R,t0,t1);
		if (dist2notch < 0.0){
			return dist2notch;
		} else {
			return fmin(dist2notch,dist2arc);
		}
	} else {//Outside disk
		return (dist_from_arc(x,y,xc,yc,R,t0,t1));
	}
}

void init_zalesak(double** G, int nx, int ny, int nghost, double dx, double dy){
	double xc = 0.5;
	double yc = 0.75;
	double R = 0.15;
	double cwidth = 0.05;
	double cdepth = 0.25;
	double cangle = M_PI*3.0/2.0;
	#pragma omp parallel for
	for (int i = -nghost; i<nx+nghost-1; i++){
		for (int j = -nghost; j<ny+nghost-1; j++){
			G[i+nghost][j+nghost] = dist_from_zalesak((i+nghost)*dx+dx/2.0f-nghost*dx,(j+nghost)*dy+dy/2.0f-nghost*dy,xc,yc,R, cwidth, cdepth, cangle);
		}
	}
	//printf("\n\n Zalesak Test: D = %f \n\n\n", dist_from_zalesak(0.0, -0.25, 0.0, 0.0, 1.0, .5, 1.5, 3.0/2.0*M_PI) );
}
void init_circle(double** G, int nx, int ny, int nghost, double dx, double dy){
	double xc = 0.5;
	double yc = 0.75;
	double R = 0.15;
	#pragma omp parallel for
	for (int i = -nghost; i<nx+nghost-1; i++){
		for (int j = -nghost; j<ny+nghost-1; j++){
			G[i+nghost][j+nghost] = dist_from_circ((i+nghost)*dx+dx/2.0f-nghost*dx,(j+nghost)*dy+dy/2.0f-nghost*dy,xc,yc,R);
		}
	}
}
void init_uv_test(double** u, double** v, int nx, int ny, int nghost, double dx, double dy){
	for (int i = -nghost; i<nx+nghost+1; i++){
		for (int j = -nghost; j<ny+nghost+1; j++){
			if(i != nx+nghost){
				u[i+nghost][j+nghost] = (0.5-(j*dy+dy/2));
				//printf("u(%i, %i) = %7.7f\n",i+nghost,j+nghost,u[i+nghost][j+nghost]);
			}
			if( j!= ny+nghost){
				v[i+nghost][j+nghost] = -(0.5-(i*dx+dx/2));
			}
		}
	}
}

void init_uv_test_static(double** u, double** v, int nx, int ny, int nghost, double dx, double dy){
	for (int i = -nghost; i<nx+nghost+1; i++){
		for (int j = -nghost; j<ny+nghost+1; j++){
			if(i != nx+nghost){
				u[i+nghost][j+nghost] = 0.0;
				//printf("u(%i, %i) = %7.7f\n",i+nghost,j+nghost,u[i+nghost][j+nghost]);
			}
			if( j!= ny+nghost){
				v[i+nghost][j+nghost] = 0.0;
			}
		}
	}
}
void init_uv_test_t3(double** u, double** v, int nx, int ny, int nghost, double dx, double dy, double t){
	double x;
	double y;
	#pragma omp parallel for private(x,y)
	for (int i = -nghost; i<nx+nghost+1; i++){
		x = (i*dx+dx/2);
		for (int j = -nghost; j<ny+nghost+1; j++){
			y = (j*dy+dy/2)-0.5;
			if(i != nx+nghost){
				u[i+nghost][j+nghost] = -cos(M_PI*t/8.0)*sin(M_PI*x)*sin(M_PI*x)*sin(2.0*M_PI*y);
			}
			if( j!= ny+nghost){
				v[i+nghost][j+nghost] = -cos(M_PI*t/8.0)*sin(2.0*M_PI*x)*cos(M_PI*y)*cos(M_PI*y);
			}
		}
	}
}



void levelset_diff(double** restrict G,double** restrict dGdt,double** restrict dGdxp,double** restrict dGdxm,double** restrict dGdyp,double** restrict dGdym, double** restrict u, double** restrict v, double dx, double dy, int nx, int ny, int nghost){
#ifdef DBGLVLST
	printf("\t Calculating dGdx\n");
#endif
#pragma omp parallel for
	for(int i = -nghost; i<nx+nghost-2; i++){
		for(int j = -nghost; j<ny+nghost-1; j++){
			dGdxp[i+nghost][j+nghost] = (G[i+nghost+1][j+nghost]-G[i+nghost][j+nghost]);
			dGdxm[i+nghost+1][j+nghost] = dGdxp[i+nghost][j+nghost];
		}
	}
#ifdef DBGLVLST
	printf("\t Calculating dGdy\n");
#endif
	#pragma omp parallel for
	for(int i = -nghost; i<nx+nghost-1; i++){
		for(int j = -nghost; j<ny+nghost-2; j++){
			dGdyp[i+nghost][j+nghost] = (G[i+nghost][j+nghost+1]-G[i+nghost][j+nghost]);
			dGdym[i+nghost][j+nghost+1] = dGdxp[i+nghost][j+nghost];
		}
	}
#ifdef DBGLVLST
	printf("\t Calculating dGdt\n");
#endif
	#pragma omp parallel for
	for(int i = 0; i<nx-1; i++){
		for(int j = 0; j<ny-1; j++){

			double uh = (u[i+nghost+1][j+nghost]+u[i+nghost][j+nghost])/2.0;
			double vh = (v[i+nghost][j+nghost+1]+v[i+nghost][j+nghost])/2.0;

			dGdt[i+nghost][j+nghost] = uh/(12.0*dx)*(-dGdxp[i+nghost-2][j+nghost]+7.0*dGdxp[i+nghost-1][j+nghost]+7.0*dGdxp[i+nghost][j+nghost]-dGdxp[i+nghost+1][j+nghost]) +
										vh/(12.0*dy)*(-dGdyp[i+nghost][j+nghost-2]+7.0*dGdyp[i+nghost][j+nghost-1]+7.0*dGdyp[i+nghost][j+nghost]-dGdyp[i+nghost][j+nghost+1]);
			if (uh >= 0){
				dGdt[i+nghost][j+nghost] -= uh*slv_psi_weno(
						(dGdxp[i+nghost-2][j+nghost]-dGdxp[i+nghost-3][j+nghost])/dx,
						(dGdxp[i+nghost-1][j+nghost]-dGdxp[i+nghost-2][j+nghost])/dx,
						(dGdxp[i+nghost-0][j+nghost]-dGdxp[i+nghost-1][j+nghost])/dx,
						(dGdxp[i+nghost+1][j+nghost]-dGdxp[i+nghost-0][j+nghost])/dx);
			} else {
				dGdt[i+nghost][j+nghost] += uh*slv_psi_weno(
						(dGdxp[i+nghost+2][j+nghost]-dGdxp[i+nghost+1][j+nghost])/dx,
						(dGdxp[i+nghost+1][j+nghost]-dGdxp[i+nghost-0][j+nghost])/dx,
						(dGdxp[i+nghost+0][j+nghost]-dGdxp[i+nghost-1][j+nghost])/dx,
						(dGdxp[i+nghost-1][j+nghost]-dGdxp[i+nghost-2][j+nghost])/dx);
			}
			if (vh >= 0.0){
				dGdt[i+nghost][j+nghost] -= vh*slv_psi_weno(
						(dGdyp[i+nghost][j+nghost-2]-dGdyp[i+nghost][j+nghost-3])/dy,
						(dGdyp[i+nghost][j+nghost-1]-dGdyp[i+nghost][j+nghost-2])/dy,
						(dGdyp[i+nghost][j+nghost-0]-dGdyp[i+nghost][j+nghost-1])/dy,
						(dGdyp[i+nghost][j+nghost+1]-dGdyp[i+nghost][j+nghost-0])/dy);
			} else {
				dGdt[i+nghost][j+nghost] += vh*slv_psi_weno(
						(dGdyp[i+nghost][j+nghost+2]-dGdyp[i+nghost][j+nghost+1])/dy,
						(dGdyp[i+nghost][j+nghost+1]-dGdyp[i+nghost][j+nghost+0])/dy,
						(dGdyp[i+nghost][j+nghost+0]-dGdyp[i+nghost][j+nghost-1])/dy,
						(dGdyp[i+nghost][j+nghost-1]-dGdyp[i+nghost][j+nghost-2])/dy);
			}
		}
	}
#ifdef DBGLVLST
	printf("\t Done with dGdt \n");
#endif
}


void reinit_diff(double** G, double** G0, double** restrict dGdt,double** restrict dGdxp,double** restrict dGdxm,double** restrict dGdyp,double** restrict dGdym, double dx, double dy, int nx, int ny, int nghost){
	double eps = (dx+dy)/2.0;
	double txp, txm, typ, tym, S;
#ifdef DBGLVLST
	printf("\t Calculating dGdx\n");
#endif
#pragma omp parallel for
	for(int i = -nghost; i<nx+nghost-2; i++){
		for(int j = -nghost; j<ny+nghost-1; j++){
			dGdxp[i+nghost][j+nghost] = (G[i+nghost+1][j+nghost]-G[i+nghost][j+nghost]);
			dGdxm[i+nghost+1][j+nghost] = dGdxp[i+nghost][j+nghost];
		}
	}
#ifdef DBGLVLST
	printf("\t Calculating dGdy\n");
#endif
	#pragma omp parallel for
	for(int i = -nghost; i<nx+nghost-1; i++){
		for(int j = -nghost; j<ny+nghost-2; j++){
			dGdyp[i+nghost][j+nghost] = (G[i+nghost][j+nghost+1]-G[i+nghost][j+nghost]);
			dGdym[i+nghost][j+nghost+1] = dGdxp[i+nghost][j+nghost];
		}
	}
#ifdef DBGLVLST
	printf("\t Calculating dGdt\n");
#endif
	#pragma omp parallel for private(txp,txm,typ,tym,S)
	for(int i = 0; i<nx-1; i++){
		for(int j = 0; j<ny-1; j++){

			S = G0[i+nghost][j+nghost]/sqrt((G0[i+nghost][j+nghost]*G0[i+nghost][j+nghost]+eps*eps));
			txp = 1.0/(12.0*dx)*(-dGdxp[i+nghost-2][j+nghost]+7.0*dGdxp[i+nghost-1][j+nghost]+7.0*dGdxp[i+nghost][j+nghost]-dGdxp[i+nghost+1][j+nghost]);
			txm = txp;
			typ = 1.0/(12.0*dy)*(-dGdyp[i+nghost][j+nghost-2]+7.0*dGdyp[i+nghost][j+nghost-1]+7.0*dGdyp[i+nghost][j+nghost]-dGdyp[i+nghost][j+nghost+1]);
			tym = typ;

			txm -=slv_psi_weno(
					(dGdxp[i+nghost-2][j+nghost]-dGdxp[i+nghost-3][j+nghost])/dx,
					(dGdxp[i+nghost-1][j+nghost]-dGdxp[i+nghost-2][j+nghost])/dx,
					(dGdxp[i+nghost-0][j+nghost]-dGdxp[i+nghost-1][j+nghost])/dx,
					(dGdxp[i+nghost+1][j+nghost]-dGdxp[i+nghost-0][j+nghost])/dx);
			txp += slv_psi_weno(
					(dGdxp[i+nghost+2][j+nghost]-dGdxp[i+nghost+1][j+nghost])/dx,
					(dGdxp[i+nghost+1][j+nghost]-dGdxp[i+nghost-0][j+nghost])/dx,
					(dGdxp[i+nghost+0][j+nghost]-dGdxp[i+nghost-1][j+nghost])/dx,
					(dGdxp[i+nghost-1][j+nghost]-dGdxp[i+nghost-2][j+nghost])/dx);
			tym -= slv_psi_weno(
					(dGdyp[i+nghost][j+nghost-2]-dGdyp[i+nghost][j+nghost-3])/dy,
					(dGdyp[i+nghost][j+nghost-1]-dGdyp[i+nghost][j+nghost-2])/dy,
					(dGdyp[i+nghost][j+nghost-0]-dGdyp[i+nghost][j+nghost-1])/dy,
					(dGdyp[i+nghost][j+nghost+1]-dGdyp[i+nghost][j+nghost-0])/dy);
			typ += slv_psi_weno(
					(dGdyp[i+nghost][j+nghost+2]-dGdyp[i+nghost][j+nghost+1])/dy,
					(dGdyp[i+nghost][j+nghost+1]-dGdyp[i+nghost][j+nghost+0])/dy,
					(dGdyp[i+nghost][j+nghost+0]-dGdyp[i+nghost][j+nghost-1])/dy,
					(dGdyp[i+nghost][j+nghost-1]-dGdyp[i+nghost][j+nghost-2])/dy);

			if (G[i+nghost][j+nghost] > 0.0){
				dGdt[i+nghost][j+nghost] = S*(sqrt(fmax( fmax(txm,0.0)*fmax(txm,0.0) , fmin(txp,0.0)*fmin(txp,0.0) ) + fmax( fmax(tym,0.0)*fmax(tym,0.0),fmin(typ,0.0)*fmin(typ,0.0)))-1.0);
			} else if (G[i+nghost][j+nghost] < 0.0) {
				dGdt[i+nghost][j+nghost] = S*(sqrt(fmax( fmin(txm,0.0)*fmin(txm,0.0) , fmax(txp,0.0)*fmax(txp,0.0) ) + fmax( fmin(tym,0.0)*fmin(tym,0.0),fmax(typ,0.0)*fmax(typ,0.0)))-1.0);
			} else {
				dGdt[i+nghost][j+nghost] = 0.0;
			}
		}
	}
#ifdef DBGLVLST
	printf("\t Done with dGdt \n");
#endif
}


void levelset_advect_euler(double** restrict G,double** restrict dGdt,double** restrict dGdxp,double** restrict dGdxm,double** restrict dGdyp,double** restrict dGdym, double** restrict u, double** restrict v, double dx, double dy, double dt, int nx, int ny, int nghost){
	levelset_diff(G,dGdt,dGdxp,dGdxm,dGdyp,dGdym,u,v,dx,dy,nx, ny,nghost);
	#pragma omp parallel for
	for(int i = 0; i<nx; i++){
		for(int j = 0; j<ny; j++){
			G[i+nghost][j+nghost] -= dt*dGdt[i+nghost][j+nghost];
		}
	}
}
void levelset_advect_TVDRK3(double** restrict G,double** restrict G1,double** restrict G2,double** restrict dGdt,double** restrict dGdt1,double** restrict dGdt2,double** restrict dGdxp,double** restrict dGdxm,double** restrict dGdyp,double** restrict dGdym, double** restrict u, double** restrict v, double dx, double dy, double dt, double time, int nx, int ny, int nghost){
	//Step1
	init_uv_test_t3(u, v, nx, ny, nghost, dx, dy,time);
	levelset_diff(G,dGdt,dGdxp,dGdxm,dGdyp,dGdym,u,v,dx,dy,nx, ny,nghost);
	#pragma omp parallel for
	for(int i = 0; i<nx; i++){
		for(int j = 0; j<ny; j++){
			G1[i+nghost][j+nghost] = G[i+nghost][j+nghost]-dt*dGdt[i+nghost][j+nghost];
		}
	}
	set_all_bcs_neumann(G1,dx,dy,nx,ny,nghost,nghost);
	//Step2
	init_uv_test_t3(u, v, nx, ny, nghost, dx, dy,time+1.0/4.0*dt);
	levelset_diff(G1,dGdt1,dGdxp,dGdxm,dGdyp,dGdym,u,v,dx,dy,nx, ny,nghost);
	#pragma omp parallel for
	for(int i = 0; i<nx; i++){
		for(int j = 0; j<ny; j++){
			G2[i+nghost][j+nghost] = G1[i+nghost][j+nghost]+3.0/4.0*dt*dGdt[i+nghost][j+nghost]
			                                               -1.0/4.0*dt*dGdt1[i+nghost][j+nghost];
		}
	}
	set_all_bcs_neumann(G2,dx,dy,nx,ny,nghost,nghost);
	//Step3
	init_uv_test_t3(u, v, nx, ny, nghost, dx, dy,time+2.0/3.0*dt);
	levelset_diff(G2,dGdt2,dGdxp,dGdxm,dGdyp,dGdym,u,v,dx,dy,nx,ny,nghost);
	#pragma omp parallel for
	for(int i = 0; i<nx; i++){
		for(int j = 0; j<ny; j++){
			G[i+nghost][j+nghost] = G2[i+nghost][j+nghost]+1.0/12.0*dt*dGdt[i+nghost][j+nghost]
			                                              +1.0/12.0*dt*dGdt1[i+nghost][j+nghost]
			                                              -2.0/3.0 *dt*dGdt2[i+nghost][j+nghost];
		}
	}
}

double reinit_advect_TVDRK3(double** restrict G,double** restrict G0, double** restrict G1,double** restrict G2,double** restrict dGdt,double** restrict dGdt1,double** restrict dGdt2,double** restrict dGdxp,double** restrict dGdxm,double** restrict dGdyp,double** restrict dGdym, double dx, double dy, double dt, int nx, int ny, int nghost){
	//Step1
	reinit_diff(G,G0,dGdt,dGdxp, dGdxm,dGdyp,dGdym, dx, dy, nx, ny, nghost);
	#pragma omp parallel for
	for(int i = 0; i<nx; i++){
		for(int j = 0; j<ny; j++){
			G1[i+nghost][j+nghost] = G[i+nghost][j+nghost]-dt*dGdt[i+nghost][j+nghost];
		}
	}
	set_all_bcs_neumann(G1,dx,dy,nx,ny,nghost,nghost);
	//Step2
	//levelset_diff(G1,dGdt1,dGdxp,dGdxm,dGdyp,dGdym,u,v,dx,dy,nx, ny,nghost);
	reinit_diff(G1,G0,dGdt1,dGdxp,dGdxm,dGdyp,dGdym, dx, dy, nx, ny, nghost);
	#pragma omp parallel for
	for(int i = 0; i<nx; i++){
		for(int j = 0; j<ny; j++){
			G2[i+nghost][j+nghost] = G1[i+nghost][j+nghost]+3.0/4.0*dt*dGdt[i+nghost][j+nghost]
			                                               -1.0/4.0*dt*dGdt1[i+nghost][j+nghost];
		}
	}
	set_all_bcs_neumann(G2,dx,dy,nx,ny,nghost,nghost);
	//Step3
	//levelset_diff(G2,dGdt2,dGdxp,dGdxm,dGdyp,dGdym,u,v,dx,dy,nx,ny,nghost);
	reinit_diff(G2,G0,dGdt2,dGdxp, dGdxm,dGdyp,dGdym, dx, dy, nx, ny, nghost);
	double err;
	double alpha = 6.0*dx;
	int count = 0;
	//#pragma omp parallel for
	for(int i = 0; i<nx; i++){
		for(int j = 0; j<ny; j++){
			if (fabs(G[i+nghost][j+nghost]) < alpha){
			err += fabs((G[i+nghost][j+nghost])-(G2[i+nghost][j+nghost]+1.0/12.0*dt*dGdt[i+nghost][j+nghost]
			                   																	  +1.0/12.0*dt*dGdt1[i+nghost][j+nghost]
			                   																	  -2.0/3.0 *dt*dGdt2[i+nghost][j+nghost]));
			count ++;
			}
			G[i+nghost][j+nghost] = G2[i+nghost][j+nghost]+1.0/12.0*dt*dGdt[i+nghost][j+nghost]
																	  +1.0/12.0*dt*dGdt1[i+nghost][j+nghost]
																	  -2.0/3.0 *dt*dGdt2[i+nghost][j+nghost];
		}
	}
	return (err/((double)count));
}

void reinit_FMM(double** restrict G,double** restrict G0, int** restrict Markers, double* restrict FMM, int* restrict FMMi, int* restrict FMMj, double dx, double dy, double dt, int nx, int ny, int nghost){
	int nAccept = 0;
	int nClose = 0;
	int nFar = 0;
	int nTot = (nx+2*nghost-1)*(ny+2*nghost-1);
	printf("\n nAccept = %i, \t nClose = %i, \t nTot = %i, \t\n ", nAccept, nClose,nTot);
	double* Close;
	int* Closei;
	int* Closej;
	double* Far;
	int* Fari;
	int* Farj;
	double Large = 1000000.0;//Million
	printf("Marker 1 \n");
	//Find Accepted
	for (int i = 0; i < nx; i++){
		for(int j = 0; j < ny; j++){

			if(G[i+nghost][j+nghost] >= 0.0){//Positive
				//If sign changed
				if((G[i+nghost+1][j+nghost] < 0.0) || (G[i+nghost-1][j+nghost] < 0.0) ||
						(G[i+nghost][j+nghost+1] < 0.0) || (G[i+nghost][j+nghost-1] < 0.0)){
					Markers[i+nghost][j+nghost]  = +1;
					FMM[nAccept] = G[i+nghost][j+nghost];
					FMMi[nAccept] = i;
					FMMj[nAccept] = j;
					nAccept++;
				}
			} else { //Negative
				//If sign changed
				if((G[i+nghost+1][j+nghost] > 0.0) || (G[i+nghost-1][j+nghost] > 0.0) ||
						(G[i+nghost][j+nghost+1] > 0.0) || (G[i+nghost][j+nghost-1] > 0.0)){
					Markers[i+nghost][j+nghost]  = -1;
					FMM[nAccept] = -G[i+nghost][j+nghost];
					FMMi[nAccept] = i;
					FMMj[nAccept] = j;
					nAccept++;
				}
			}
		}
	}
	printf("Marker 2 \n");
	//Find Close and Far
	for (int i = 0; i < nx; i++){
		for(int j = 0; j < ny; j++){
			if(G[i+nghost][j+nghost] >= 0.0 && Markers[i+nghost][j+nghost] != 1 && Markers[i+nghost][j+nghost] != -1){//Positive
				//If Close
				if(((Markers[i+nghost+1][j+nghost] == 1) || (Markers[i+nghost-1][j+nghost] == 1) ||
						(Markers[i+nghost][j+nghost+1] == 1) || (Markers[i+nghost][j+nghost-1] == 1))){
					Markers[i+nghost][j+nghost]  = +2;
					FMM[nAccept+nClose] = G[i+nghost][j+nghost];
					FMMi[nAccept+nClose] = i;
					FMMj[nAccept+nClose] = j;
					nClose++;
				} else {//point is far
					Markers[i+nghost][j+nghost]  = +3;
					FMM[nTot-1-nFar] = Large;
					FMMi[nTot-1-nFar] = i;
					FMMj[nTot-1-nFar] = j;
					nFar++;
				}
			} else if(Markers[i+nghost][j+nghost] != 1 && Markers[i+nghost][j+nghost] != -1) { //Negative
				//If close
				if(((Markers[i+nghost+1][j+nghost] == -1) || (Markers[i+nghost-1][j+nghost] == -1) ||
						(Markers[i+nghost][j+nghost+1] == -1) || (Markers[i+nghost][j+nghost-1] == -1))){
					Markers[i+nghost][j+nghost]  = -2;
					FMM[nAccept+nClose] = -G[i+nghost][j+nghost];
					FMMi[nAccept+nClose] = i;
					FMMj[nAccept+nClose] = j;
					nClose++;
				} else {//point is far
					Markers[i+nghost][j+nghost]  = -3;
					FMM[nTot-1-nFar] = Large;
					FMMi[nTot-1-nFar] = i;
					FMMj[nTot-1-nFar] = j;
					nFar++;
				}
			}
		}
	}
	printf("Marker 3 \n");
	Close = &FMM[nAccept-1];
	Closei = &FMMi[nAccept-1];
	Closej = &FMMj[nAccept-1];

	Far = &FMM[nTot-1-nFar];
	Fari = &FMMi[nTot-1-nFar];
	Farj = &FMMj[nTot-1-nFar];

	printf("Marker 4 \n");
	printf("\n nAccept = %i, \t nClose = %i, \t nFar = %i, \t nTot = %i, \t Close-FMM %i \t\n ", nAccept, nClose,nFar,nTot,Close-FMM);
	print_array(FMM, nTot);
	print_array(Close, nClose);
	printf("Marker 5 \n");
	heapsort(Close,Closei,Closej,nClose);
	printf("Marker 6 \n");
	print_array(Close, nClose);
	//print_array_int(FMMi, nTot);
	//print_array_int(FMMj, nTot);

	int i,j;
	double alpha,beta,a,b,A,B,C,Gtemp;
	int count = 0;
	write_matrix_2d_int(Markers,nx+2*nghost-1,ny+2*nghost-1,"Markers.0");
/*
	while (nClose > 0){
		i = Closei[0];
		j = Closej[0];

		a = fmin(fabs(G[i+nghost+1][j+nghost]),fabs(G[i+nghost-1][j+nghost]));
		b = fmin(fabs(G[i+nghost][j+nghost+1]),fabs(G[i+nghost][j+nghost+1]));

		if(Markers[i+nghost+1][j+nghost] == 1 || Markers[i+nghost-1][j+nghost] == 1 || Markers[i+nghost+1][j+nghost] == -1 || Markers[i+nghost-1][j+nghost] == -1){
			alpha = 1.0;
		} else {
			alpha = 0.0;
		}
		if(Markers[i+nghost][j+nghost+1] == 1 || Markers[i+nghost][j+nghost-1] == 1 || Markers[i+nghost][j+nghost+1] == -1 || Markers[i+nghost][j+nghost-1] == -1){
					beta = 1.0;
		} else {
			beta = 0.0;
		}
		A = alpha+beta;
		B = -2*(alpha*a+beta*b);
		C = alpha*a*a+beta*b*b-dx*dy;
		Gtemp = (-B+sqrt(B*B-4*A*C))/(2*A);
		if (Gtemp < 0) {
			printf("\n******* NEGATIVE GTEMP *******\n");
			printf("A = %f, B = %f, C = %f \n\n", A, B, C);
		}


		if(G[i+nghost][j+nghost]  > 0.0){
			Close[0] = Gtemp;
			G[i+nghost][j+nghost] = Gtemp;
		} else {
			Close[0] = -Gtemp;
			G[i+nghost][j+nghost] = -Gtemp;
		}
		nAccept++;
		Close = &Close[1];
		Closei = &Closei[1];
		Closej = &Closej[1];
		nClose -= 1;

		//Fix +-i
		if(i>0 && i<(nx-1)){
			if(Markers[i+nghost+1][j+nghost] == 3){
				Markers[i+nghost+1][j+nghost] = 2;
				Close[nClose] = G[i+nghost+1][j+nghost];
				Closei[nClose] = i+1;
				Closej[nClose] = j;
				nClose++;
			} else if(Markers[i+nghost+1][j] == -3){
				Markers[i+nghost+1][j+nghost] = -2;
				Close[nClose] = -G[i+nghost+1][j+nghost];
				Closei[nClose] = i+1;
				Closej[nClose] = j;
				nClose++;
			}
			if(Markers[i+nghost-1][j+nghost] == 3){
				Markers[i+nghost+1][j+nghost] = 2;
				Close[nClose] = G[i+nghost+1][j+nghost];
				Closei[nClose] = i-1;
				Closej[nClose] = j;
				nClose++;
			} else if(Markers[i+nghost-1][j+nghost] == -3){
				Markers[i+nghost-1][j+nghost] = -2;
				Close[nClose] = -G[i+nghost-1][j+nghost];
				Closei[nClose] = i-1;
				Closej[nClose] = j;
				nClose++;
			}
		}
		//Fix +-j
		if(j > 0 && j < (ny-1)){
			if(Markers[i+nghost][j+nghost+1] == 3){
				Markers[i+nghost][j+nghost+1] = 2;
				Close[nClose] = G[i+nghost][j+nghost+1];
				Closei[nClose] = i;
				Closej[nClose] = j+1;
				nClose++;
			} else if(Markers[i+nghost+1][j+nghost+1] == -3){
				Markers[i+nghost][j+nghost+1] = -2;
				Close[nClose] = -G[i+nghost][j+nghost+1];
				Closei[nClose] = i;
				Closej[nClose] = j+1;
				nClose++;
			}
			if(Markers[i+nghost][j+nghost-1] == 3){
				Markers[i+nghost][j+nghost-1] = 2;
				Close[nClose] = G[i+nghost][j+nghost-1];
				Closei[nClose] = i;
				Closej[nClose] = j-1;
				nClose++;
			} else if(Markers[i+nghost][j+nghost-1] == -3){
				Markers[i+nghost][j+nghost-1] = -2;
				Close[nClose] = -G[i+nghost][j+nghost-1];
				Closei[nClose] = i;
				Closej[nClose] = j-1;
				nClose++;
			}
		}
		heapsort(Close,Closei,Closej,nClose);//Resort
		count++;
		if (count == 1){
			write_matrix_2d_int(Markers,nx+2*nghost-1,ny+2*nghost-1,"Markers.1");
		}
		if (count == 2){
			write_matrix_2d_int(Markers,nx+2*nghost-1,ny+2*nghost-1,"Markers.2");
		}
		if (count == 100){
			write_matrix_2d_int(Markers,nx+2*nghost-1,ny+2*nghost-1,"Markers.100");
		}
		printf("\tCount = %i, nAccept = %i, nClose = %i \n",count,nAccept, nClose);
		printf("i,j = %i,%i\n",i,j);
	} */

}



double lvl_H (double fx, double a){
	if(fx >= a){
		return 1.0;
	} else if (fx <= -a) {
		return 0.0;
	} else {
		return 1.0/2.0*(fx/a+1/M_PI*sin(M_PI*fx/a))+0.5;
	}
}

double get_shape_err(double** restrict G, double** restrict Gt, double volGt, double a, int nx, int ny, int nghost,double dx, double dy){

	double err;
	for(int i = 0; i<nx; i++){
		for(int j = 0; j<ny; j++){
			err += fabs(lvl_H(G[i+nghost][j+nghost], a) - lvl_H(Gt[i+nghost][j+nghost],a))*dx*dy;
		}
	}
	return err/volGt;
}
double get_vol(double** restrict G, double a, int nx, int ny, int nghost,double dx, double dy){

	double vol;
	for(int i = 0; i<nx; i++){
		for(int j = 0; j<ny; j++){
			vol += lvl_H(G[i+nghost][j+nghost], a)*dx*dy;
		}
	}
	return vol;
}

