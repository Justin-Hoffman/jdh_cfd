/*
 * util.h
 *
 *  Created on: Feb 20, 2014
 *      Author: justin
 */

#ifndef UTIL_H_
#define UTIL_H_
typedef struct {
	double x;
	double y;
} V2D;
typedef struct {
	V2D start;
	V2D end;
} lseg;
typedef struct {
	double l;
	double u;
	double r;
	double d;
} Wall;
void heapsort(double* R, int* Ni, int* Nj, int** MarkPoint, int n, int nghost );
int G_from_flux(int i, int j, double** G, int** Markers,int** MarkPoint, double* Close, int* Closei, int* Closej, int nClose, int nghost, double dx, double dy);
void get_alpha(double** G, double** alpha,  int nghost, double dx, double dy,int nx, int ny);
double get_alpha_sum(double** alpha,  int nghost, double dx, double dy,int nx, int ny);
double get_alpha_diff(double** alpha, double** alpha2,  int nghost, double dx, double dy,int nx, int ny);
double get_tke(double** alpha, double** u, double** v, int nghost, double dx, double dy,int nx, int ny, struct slv_settings st);
double get_minp(double** alpha, double** phi, int nghost, double dx, double dy,int nx, int ny, struct slv_settings st);
double get_maxp(double** alpha, double** phi, int nghost, double dx, double dy,int nx, int ny, struct slv_settings st);
void get_curv(double** K, double** G, int nx, int ny, int nghost, double dx, double dy, struct slv_settings st);
V2D vof_alpha_grad(double** alpha, int i, int j, int nghost, double dx, double dy);
double vof_d_from_alpha(double a, V2D m, double hx, double hy);
double vof_alpha_from_surf(double d, V2D m, double hx, double hy);
double H(double x);
void vof_construct_ds(double** alpha,double** darr, V2D** marr, int nghost, double dx, double dy,int nx, int ny);
void vof_get_surfaces(double** alpha,double** darr, V2D** marr,lseg** segs, int nghost, double dx, double dy,int nx, int ny);
void vof_get_surface(double** alpha, double** darr, V2D** marr, lseg** segs, int i, int j, int nghost, double dx, double dy);
void vof_get_fluxes(double** u, double** v, double** alpha,double** darr, V2D** marr, double dt, int nghost, double dx, double dy, int nx, int ny,int t);

#endif /* UTIL_H_ */
