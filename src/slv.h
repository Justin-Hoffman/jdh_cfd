/*
 * slv.h
 *
 *  Created on: Jan 21, 2014
 *      Author: justin
 */

#ifndef SLV_H_
#define SLV_H_

void solve_matrix(int n, double* a, double* b, double* c, double* v, double* x);

/* Set Boundary Conditions */
void set_bcs(double** restrict u, double** restrict v, double dx, double dy, int  nx, int ny, double Re);

/* Get u and V */
void get_uv(double** restrict u, double** restrict v, double** restrict usv,double** restrict vsv,int nx, int ny);

/* Solve for vorticity */
void get_vorticity(double** u, double** v,double** omega, double dx, double dy, int nx, int ny);

/* Solve for Stream Function */
void get_stream(double** restrict strm, double** restrict strmnext, double** restrict omega, double dx, double dy, int nx, int ny);

/* Viscous Burgers Equation Solver */
/* Uses AB2 + ADI */
void slv_vbe(double** restrict u, double** restrict us, double** restrict dus, double** restrict duss, double** restrict v, double** restrict vs, double** restrict dvs, double** restrict dvss, double** restrict hu, double** restrict huold,double** restrict hv, double** restrict hvold, double dx, double dy, int nx, int ny, double Re, double dt);

/* ftcs solver */
void slv_vbeftcs(double** restrict u, double** restrict us, double** restrict dus, double** restrict duss, double** restrict v, double** restrict vs, double** restrict dvs, double** restrict dvss, double** restrict hu, double** restrict huold,double** restrict hv, double** restrict hvold, double dx, double dy, int nx, int ny, double Re, double dt);

/* Poisson Solver */
void slv_pssn(double** restrict phi,double** restrict phinext, double** restrict us, double** restrict vs, double dx, double dy, int nx, int ny,double dt,double min);

/* Apply Solution Projection into Solenoidal */
void apply_projection(double** restrict phi, double** restrict u, double** restrict us, double** restrict v, double** restrict vs, double dx, double dy, int nx, int ny, double dt);

/* Write 2D Matrix to File */
void write_matrix_2d(double** mat, int nx, int ny, char* filename);


#endif /* SLV_H_ */
