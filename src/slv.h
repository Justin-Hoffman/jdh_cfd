/*
 * slv.h
 *
 *  Created on: Jan 21, 2014
 *      Author: justin
 */

#ifndef SLV_H_
#define SLV_H_

typedef enum {WALL,PERIODIC} jdh_BC;

	struct slv_settings {
		int nx;
		int ny;
		int nt;
		double Re;
		double dt;
		jdh_BC XBC;
		jdh_BC YBC;
		double XUV;
		double XLV;
		double YLV;
		double YRV;
		double XConst;
		double YConst;
	};

void solve_matrix(int n, double* a, double* b, double* c, double* v, double* x);

void copy_2D(double** target, double** source, int nx, int ny);

/* Set Boundary Conditions */
void set_bcs(double** restrict u, double** restrict v, double dx, double dy, int  nx, int ny, double Re,struct slv_settings st);

/* Set Neumann Boundary Conditions */
void set_all_bcs_neumann(double** restrict x, double dx, double dy, int  nx, int ny, int nghostx, int nghosty);

/* Get u and V */
void get_uv(double** restrict u, double** restrict v, double** restrict usv,double** restrict vsv,int nx, int ny);

/* Solve for vorticity */
void get_vorticity(double** u, double** v,double** omega, double dx, double dy, int nx, int ny);

/* Solve for Stream Function */
void get_stream(double** restrict strm, double** restrict strmnext, double** restrict omega, double dx, double dy, int nx, int ny);

/* Viscous Burgers Equation Solver */
/* Uses AB2 + ADI */
void slv_vbe(double** restrict u, double** restrict us, double** restrict dus, double** restrict duss, double** restrict v, double** restrict vs, double** restrict dvs, double** restrict dvss, double** restrict hu, double** restrict huold,double** restrict hv, double** restrict hvold, double dx, double dy, int nx, int ny, double Re, double dt,struct slv_settings st);

/* ftcs solver */
void slv_vbeftcs(double** restrict u, double** restrict us, double** restrict dus, double** restrict duss, double** restrict v, double** restrict vs, double** restrict dvs, double** restrict dvss, double** restrict hu, double** restrict huold,double** restrict hv, double** restrict hvold, double dx, double dy, int nx, int ny, double Re, double dt,struct slv_settings st);

/* Poisson Solver */
void slv_pssn(double** restrict phi,double** restrict phinext, double** restrict us, double** restrict vs, double dx, double dy, int nx, int ny,double dt,double min);

/* Poisson Solver */
void slv_pssn_gmres(double** restrict phi,double** restrict phinext, double** restrict us, double** restrict vs, double dx, double dy, int nx, int ny,double dt,double min, struct slv_settings st);

/* find value in sparse triple pair array */
double find_in_A (double* A, int* iA, int* jA, int npt, int i, int j);

/* Apply Solution Projection into Solenoidal */
void apply_projection(double** restrict phi, double** restrict u, double** restrict us, double** restrict v, double** restrict vs, double dx, double dy, int nx, int ny, double dt);


/* Write 2D Matrix to File */
void write_matrix_2d(double** mat, int nx, int ny, char* filename);

/* Write 2D Matrix to File */
void write_matrix_2d_int(int** mat, int nx, int ny, char* filename);

/*Init Settings*/
struct slv_settings init_settings();

/*Get weno5 weight*/
double slv_psi_weno(double a, double b, double c, double d);

void print_array(double* X, int n);

void print_array_int(int* X, int n);


#endif /* SLV_H_ */
