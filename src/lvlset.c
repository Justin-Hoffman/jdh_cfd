/*
 * lvlset.c
 *
 *  Created on: Jan 30, 2014
 *      Author: justin
 */
#define M_PI 3.14159265358979323846L
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>
#include <unistd.h>
#include "lvlset.h"

double dist_from_arc(double x, double y, double xc, double yc, double R, double t0, double t1){
	double tc = atan2((y-yc),(x-xc)); // Angle of current pt
	if (tc < 0.0) {
		tc = 2.0L*M_PI+tc;
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

	if ((a0>=0.0 && a1 <= 0.0 & a2 <= 0.0) || (a0 >= 0.0 && a1 >= 0.0 && a2 < 0.0) ){
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
		return 0.0l;
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
	double C = M_PI/2.0l;
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
	double C = M_PI/2.0l;
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
	double halfcut = atan2(cwidth/2.0l,R); //Trig to identify half angle of notch
	double t0 = fmod((cangle-halfcut),(2.0l*M_PI)); //Arc angle at start of cut
	double t1 = fmod((cangle+halfcut),(2.0l*M_PI)); //Arc angle at end of cut
	double theta = fmod(atan2((x-xc),(y-yc)),(2.0l*M_PI)); // Angle of current pt

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
	double yc = 0.5;
	double R = 0.15;
	double cwidth = 0.05;
	double cdepth = 0.25;
	double cangle = M_PI*2.5l/2.0l;
	for (int i = -nghost; i<nx+nghost-1; i++){
		for (int j = -nghost; j<ny+nghost; j++){
			G[i+nghost][j+nghost] = dist_from_zalesak((i+nghost)*dx+dx/2.0f-nghost*dx,(j+nghost)*dy+dy/2.0f-nghost*dy,xc,yc,R, cwidth, cdepth, cangle);
		}
	}

	//printf("\n\n Zalesak Test: D = %f \n\n\n", dist_from_zalesak(0.0, -0.25, 0.0, 0.0, 1.0, .5, 1.5, 3.0L/2.0L*M_PI) );

}
