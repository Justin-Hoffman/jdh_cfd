#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "slv.h"
#include "mgmres.h"
#include "lvlset.h"
#include "util.h"


/* Tri-Diagonal Solver */
void solve_matrix(int n, double* a, double* b, double* c, double* v, double* x)
{ 
	for(int i = 1; i<n; i++){
		double m = a[i]/b[i-1];
		b[i] = b[i]-m*c[i-1];
		v[i] = v[i]-m*v[i-1];
	}
	x[n-1] = v[n-1]/b[n-1];
	for(int i = n-2;i>=0;i--){
		x[i] = (v[i]-c[i]*x[i+1])/b[i];
	}
}

//2D Array Copy
void copy_2D(double** target, double** source, int nx, int ny){
	size_t data_size = sizeof(double)*ny;
	for(int i = 0; i<nx; i++){
		memcpy(target[i], source[i], data_size);
	}
}

int main(int argc, char** argv)
{
	/* Initialize some variables */
	int nx = 40;
	int ny = 40;	
	int nt = 1000;
	int nsave = nt-1;
	int saves = 0;
	int memcntr = 0;
	int ninit = 1;
	double Re = 100.0;
	double dt = 0.001;
	double time = 0.0;
	struct slv_settings st;
	st = init_settings();
	st.YBC = WALL;
	st.XBC = WALL;
	st.XUV = 0.0;
	st.XLV = 0.0;
	st.XConst = 0.0;
	st.YConst = 0.0;

	/*Heapsort Test
	int Ni[6] = {0,1,2,3,4,5};
	int Nj[6] = {0,1,2,3,4,5};
	int ntest = 5;
	double R[6] = {0.0,2.0,3.0,1.0,1.1,2.5,-2.2};
	print_array(R, ntest+1);
	heapsort(R, Ni, Nj, ntest);
	print_array(R, ntest+1);
	*/

	/* Parse Input Arguments */
	for(int i = 0; i<argc-1; i++){
		if(strcmp(argv[i],"-nx")==0){
			nx = atoi(argv[i+1]);
			i++;
		} else if(strcmp(argv[i],"-ny")==0){
			ny = atoi(argv[i+1]);
			i++;
		} else if (strcmp(argv[i],"-nt")==0){
			nt = atoi(argv[i+1]);
			i++;
		} else if (strcmp(argv[i],"-nsave")==0){
			nsave = atoi(argv[i+1]);
			i++;
		} else if (strcmp(argv[i],"-Re")==0){
			Re = atof(argv[i+1]);
			i++;
		} else if (strcmp(argv[i],"-dt")==0){
			dt = atof(argv[i+1]);
			i++;
		} else if (strcmp(argv[i],"-XBC")==0){
			dt = atof(argv[i+1]);
			i++;
		} else if (strcmp(argv[i],"-ninit")==0){
			ninit = atoi(argv[i+1]);
			i++;
		}
	}

	/* Calcs for Input */
	double dx = 1.0/((double)(nx-1));
	double dy = 1.0/((double)(ny-1));
	double reinitl2 = 0.0;

	/* Echo to User */
	printf("***\t2D Incompressible Navier Stokes Lid Driven Cavity Solver\t***\n");
	printf("\n\tGrid size in X: %i\n",nx);
	printf("\tGrid size in Y: %i\n",ny);
	printf("\tReynolds Number: %f\n",Re);
	printf("\tTime Step: %f\n",dt);
	printf("\tTime Iterations: %i\n",nt);


	/*Allocate Memory for Spacial Variables */
#ifdef DEBUG
	printf("\t Begin Memory Allocation \n");
#endif

	int nghost = 3;
	double** G = malloc((nx+2*nghost-1)*sizeof(double*));memcntr += (nx+2*nghost-1)*sizeof(double*);
	double** G0 = malloc((nx+2*nghost-1)*sizeof(double*));memcntr += (nx+2*nghost-1)*sizeof(double*);
	double** G1 = malloc((nx+2*nghost-1)*sizeof(double*));memcntr += (nx+2*nghost-1)*sizeof(double*);
	double** G2 = malloc((nx+2*nghost-1)*sizeof(double*));memcntr += (nx+2*nghost-1)*sizeof(double*);
	double** dGdt = malloc((nx+2*nghost-1)*sizeof(double*));memcntr += (nx+2*nghost-1)*sizeof(double*);
	double** dGdt1 = malloc((nx+2*nghost-1)*sizeof(double*));memcntr += (nx+2*nghost-1)*sizeof(double*);
	double** dGdt2 = malloc((nx+2*nghost-1)*sizeof(double*));memcntr += (nx+2*nghost-1)*sizeof(double*);
	double** dGdxp = malloc((nx+2*nghost-1)*sizeof(double*));memcntr += (nx+2*nghost-1)*sizeof(double*);
	double** dGdxm = malloc((nx+2*nghost-1)*sizeof(double*));memcntr += (nx+2*nghost-1)*sizeof(double*);
	double** dGdyp = malloc((nx+2*nghost-1)*sizeof(double*));memcntr += (nx+2*nghost-1)*sizeof(double*);
	double** dGdym = malloc((nx+2*nghost-1)*sizeof(double*));memcntr += (nx+2*nghost-1)*sizeof(double*);

	double* FMM = malloc((nx+2*nghost-1)*(ny+2*nghost-1)*sizeof(double)); memcntr += (nx+2*nghost-1)*(ny+2*nghost-1)*sizeof(double);
	int** Markers = malloc((nx+2*nghost-1)*sizeof(int*));memcntr += (nx+2*nghost-1)*sizeof(int*);
	int* FMMi = malloc((nx+2*nghost-1)*(ny+2*nghost-1)*sizeof(int)); memcntr += (nx+2*nghost-1)*(ny+2*nghost-1)*sizeof(int);
	int* FMMj = malloc((nx+2*nghost-1)*(ny+2*nghost-1)*sizeof(int)); memcntr += (nx+2*nghost-1)*(ny+2*nghost-1)*sizeof(int);


	double** omega = malloc(nx*sizeof(double*));memcntr += nx*sizeof(double*);

	double** phi = malloc((nx+1)*sizeof(double*));memcntr += (nx+1)*sizeof(double*);
	double** phinext = malloc((nx+1)*sizeof(double*));memcntr += (nx+1)*sizeof(double*);
	double** strm = malloc((nx+3)*sizeof(double*));memcntr += (nx+1)*sizeof(double*);
	double** strmnext = malloc((nx+3)*sizeof(double*));memcntr += (nx+3)*sizeof(double*);
	
	double** u = malloc((nx+2*nghost)*sizeof(double*));memcntr += (nx+2*nghost)*sizeof(double*);
	double** us = malloc((nx+2*nghost)*sizeof(double*));memcntr += (nx+2*nghost)*sizeof(double*);
	double** dus = malloc((nx+2*nghost)*sizeof(double*));memcntr += (nx+2*nghost)*sizeof(double*);
	double** duss = malloc((nx+2*nghost)*sizeof(double*));memcntr += (nx+2*nghost)*sizeof(double*);

	double** v = malloc((nx+1+2*nghost)*sizeof(double*));memcntr += (nx+1+2*nghost)*sizeof(double*);
	double** vs = malloc((nx+1+2*nghost)*sizeof(double*));memcntr += (nx+1+2*nghost)*sizeof(double*);
	double** dvs = malloc((nx+1+2*nghost)*sizeof(double*));memcntr += (nx+1+2*nghost)*sizeof(double*);
	double** dvss = malloc((nx+1+2*nghost)*sizeof(double*));memcntr += (nx+1+2*nghost)*sizeof(double*);
		
	double** usv = malloc((nx+2*nghost)*sizeof(double*));memcntr += (nx+2*nghost)*sizeof(double*);
	double** vsv = malloc((nx+2*nghost)*sizeof(double*));memcntr += (nx+2*nghost)*sizeof(double*);

	double** hu = malloc((nx+2*nghost)*sizeof(double*));memcntr += (nx+2*nghost)*sizeof(double*);
	double** huold = malloc((nx+2*nghost)*sizeof(double*));memcntr += (nx+2*nghost)*sizeof(double*);

	double** hv = malloc((nx+1+2*nghost)*sizeof(double*));memcntr += (nx+1+2*nghost)*sizeof(double*);
	double** hvold = malloc((nx+1+2*nghost)*sizeof(double*));memcntr += (nx+1+2*nghost)*sizeof(double*);
	
#ifdef DBGMEM
	printf("\t Malloc Levelset Arrays\n");
#endif
	//Allocate Levelset Array
	for(int i = -nghost; i<nx+nghost-1; i++){
		G[i+nghost] = malloc((ny+2*nghost-1)*sizeof(double));memcntr += (nx+2*nghost-1)*sizeof(double);
		G0[i+nghost] = malloc((ny+2*nghost-1)*sizeof(double));memcntr += (nx+2*nghost-1)*sizeof(double);
		G1[i+nghost] = malloc((ny+2*nghost-1)*sizeof(double));memcntr += (nx+2*nghost-1)*sizeof(double);
		G2[i+nghost] = malloc((ny+2*nghost-1)*sizeof(double));memcntr += (nx+2*nghost-1)*sizeof(double);
		dGdt[i+nghost] = malloc((ny+2*nghost-1)*sizeof(double));memcntr += (nx+2*nghost-1)*sizeof(double);
		dGdt1[i+nghost] = malloc((ny+2*nghost-1)*sizeof(double));memcntr += (nx+2*nghost-1)*sizeof(double);
		dGdt2[i+nghost] = malloc((ny+2*nghost-1)*sizeof(double));memcntr += (nx+2*nghost-1)*sizeof(double);
		dGdxp[i+nghost] = malloc((ny+2*nghost-1)*sizeof(double));memcntr += (nx+2*nghost-1)*sizeof(double);
		dGdxm[i+nghost] = malloc((ny+2*nghost-1)*sizeof(double));memcntr += (nx+2*nghost-1)*sizeof(double);
		dGdyp[i+nghost] = malloc((ny+2*nghost-1)*sizeof(double));memcntr += (nx+2*nghost-1)*sizeof(double);
		dGdym[i+nghost] = malloc((ny+2*nghost-1)*sizeof(double));memcntr += (nx+2*nghost-1)*sizeof(double);
		Markers[i+nghost] = malloc((ny+2*nghost-1)*sizeof(int));memcntr += (ny+2*nghost-1)*sizeof(int);
	}
#ifdef DBGMEM
	printf("\t Malloc Normal Arrays\n");
#endif
	for(int i = -nghost; i<nx+nghost+1; i++){

		if(i < (nx+nghost+1)){
			u[i+nghost] = malloc((ny+1+2*nghost)*sizeof(double));memcntr += (ny+1+2*nghost)*sizeof(double);
			us[i+nghost] = malloc((ny+1+2*nghost)*sizeof(double));memcntr += (ny+1+2*nghost)*sizeof(double);
			dus[i+nghost] = malloc((ny+1+2*nghost)*sizeof(double));memcntr += (ny+1+2*nghost)*sizeof(double);
			duss[i+nghost] = malloc((ny+1+2*nghost)*sizeof(double));memcntr += (ny+1+2*nghost)*sizeof(double);
					
			hu[i+nghost] = malloc((ny+1+2*nghost)*sizeof(double));memcntr += (ny+1+2*nghost)*sizeof(double);
			huold[i+nghost] = malloc((ny+1+2*nghost)*sizeof(double));memcntr += (ny+1+2*nghost)*sizeof(double);

			usv[i+nghost] = malloc((ny+2*nghost)*sizeof(double));memcntr += (ny+2*nghost)*sizeof(double);
			vsv[i+nghost] = malloc((ny+2*nghost)*sizeof(double));memcntr += (ny+2*nghost)*sizeof(double);
			omega[i+nghost] = malloc((ny+2*nghost)*sizeof(double));memcntr += (ny+2*nghost)*sizeof(double);
		}
		if(i < (nx+nghost+2)){
			phi[i+nghost] = malloc((ny+1+2*nghost)*sizeof(double));memcntr += (ny+1+2*nghost)*sizeof(double);
			phinext[i+nghost] = malloc((ny+1+2*nghost)*sizeof(double));memcntr += (ny+1+2*nghost)*sizeof(double);
	
			v[i+nghost] = malloc((ny+2*nghost)*sizeof(double));memcntr += (ny+2*nghost)*sizeof(double);
			vs[i+nghost] = malloc((ny+2*nghost)*sizeof(double));memcntr += (ny+2*nghost)*sizeof(double);
			dvs[i+nghost] = malloc((ny+2*nghost)*sizeof(double));memcntr += (ny+2*nghost)*sizeof(double);
			dvss[i+nghost] = malloc((ny+2*nghost)*sizeof(double));memcntr += (ny+2*nghost)*sizeof(double);
	
			hv[i+nghost] = malloc((ny+2*nghost)*sizeof(double));memcntr += (ny+2*nghost)*sizeof(double);
			hvold[i+nghost] = malloc((ny+2*nghost)*sizeof(double));memcntr += (ny+2*nghost)*sizeof(double);
		}
		
		strm[i+nghost] = malloc((ny+2+2*nghost)*sizeof(double));memcntr += (ny+2+2*nghost)*sizeof(double);
		strmnext[i+nghost] = malloc((ny+2+2*nghost)*sizeof(double));memcntr += (ny+2+2*nghost)*sizeof(double);
		
		for(int j=-nghost;j<ny+nghost+1;j++){
			//printf("\t (%i, %i) \n",i,j);
			if (i < nx+1){
				u[i+nghost][j+nghost] = 0.0;
				us[i+nghost][j+nghost] = 0.0;
				dus[i+nghost][j+nghost] = 0.0;
				duss[i+nghost][j+nghost] = 0.0;
				hu[i+nghost][j+nghost] = 0.0;
				huold[i+nghost][j+nghost] = 0.0;
				if(j != ny && j != (ny+1)){
					usv[i+nghost][j+nghost]=0.0;
					vsv[i+nghost][j+nghost]=0.0;
					omega[i+nghost][j+nghost] = 0.0;
				}
			}
			if (j < (ny+1) && i < (nx+1)){
				v[i+nghost][j+nghost] = 0.0;
				vs[i+nghost][j+nghost] = 0.0;
				dvs[i+nghost][j+nghost] = 0.0;
				dvss[i+nghost][j+nghost] = 0.0;
				hv[i+nghost][j+nghost] = 0.0;
				hvold[i+nghost][j+nghost] = 0.0;
			}

			if (i < (nx+1)){
				if(j <(ny+1)){
					phi[i+nghost][j+nghost] = 0.0;
					phinext[i+nghost][j+nghost] = 0.0;
				}
			}
			strm[i+nghost][j+nghost] = 0.0;
			strmnext[i+nghost][j+nghost] = 0.0;
		}
	}
	printf("\n\tEstimated Memory Usage: %i (words)\n", memcntr);
#ifdef DEBUG
	printf("\t\t MEMORY ALLOCATION DONE \n\n");

#endif
	printf("\n\t Initializing Zalesak's Disk into G \n");
	//init_zalesak(G, nx, ny, nghost, dx, dy);
	init_circle(G, nx, ny, nghost, dx, dy);
	printf("\n\t Done Initializing!\n\n");
#ifdef DEBUG
	printf("\n\t Setting Neumann BC for G\n\n");
#endif DEBUG
	set_all_bcs_neumann(G,dx,dy,nx,ny,nghost,nghost);
	//init_uv_test(u, v, nx, ny, nghost, dx, dy);

	//init_uv_test_t3(u, v, nx, ny, nghost, dx, dy,0.0);

	double min = 0.00001;
	/* Time Iteration Loop */
	for(int t = 0;t<nt;t++){
		time = t*dt;
		if(t > nt*9.0/10.0){
			min = 0.00001;
		}
		if(t%100==0){
			printf("\t %i / %i \n",t,nt);
		}
#ifdef DEBUG
		printf("\t Simulation Loop. Time: %f / %f \n",(double)t*dt,dt*nt);
#endif

		/* Set BCS */
#ifdef DEBUG
		printf("\t Setting BCS \n");
#endif
		//set_bcs(u, v, dx, dy, nx, ny, Re, st);

		/* Solv Viscous Bergers Equation */
#ifdef DEBUG
		printf("\t Solving Viscous Bergers \n");
#endif
		//slv_vbe(u,us,dus,duss,v,vs,dvs,dvss,hu,huold,hv,hvold,dx,dy,nx,ny,Re,dt,st);

		/* Calculate phi */	
#ifdef DEBUG
		printf("\t Solving Poisson Equation \n");
#endif
		//slv_pssn(phi,phinext,us,vs,dx,dy,nx,ny,dt,min);
		//slv_pssn_gmres(phi,phinext,us,vs,dx,dy,nx,ny,dt,min,st);

		/* Apply Projection */
#ifdef DEBUG
		printf("\t Applying Projection \n");
#endif
		//apply_projection(phi,u,us,v,vs,dx,dy,nx,ny,dt);

#ifdef DEBUG
		printf("\t Advecting G \n");
#endif
		reinit_FMM(G,G0,Markers,FMM,FMMi,FMMj, dx, dy, dt, nx, ny,nghost);
		/*
		levelset_advect_TVDRK3(G,G1,G2,dGdt,dGdt1,dGdt2,dGdxp,dGdxm,dGdyp,dGdym,u,v,dx,dy,dt,time,nx,ny,nghost);
		set_all_bcs_neumann(G,dx,dy,nx,ny,nghost,nghost);
		if((t+1)%ninit == 0){
			copy_2D(G0, G, nx+2*nghost-1, ny+2*nghost-1);
			reinitl2 = reinit_advect_TVDRK3(G,G0,G1,G2,dGdt,dGdt1,dGdt2,dGdxp,dGdxm,dGdyp,dGdym,dx,dy,dx/4.0,nx,ny,nghost);
			int ind = 1;
			while(reinitl2 > dt*dx && ind < 150){
				reinitl2 = reinit_advect_TVDRK3(G,G0,G1,G2,dGdt,dGdt1,dGdt2,dGdxp,dGdxm,dGdyp,dGdym,dx,dy,dx/4.0,nx,ny,nghost);
				ind++;
			}
			//printf("Reinit converged at step %i:\t L2 err = %7.7E\n",ind, reinitl2);
			set_all_bcs_neumann(G,dx,dy,nx,ny,nghost,nghost);
		}

		if( ((t+1)%nsave) == 0) {
			char fname[50];
			//strcopy(fname,"u.");
			sprintf(&fname,"G.%i",saves);
			write_matrix_2d(G,nx+2*nghost-1,ny+2*nghost-1,fname);
			saves++;
		}

		/* Repeat */
		if (t==nt-1){
			get_uv(u,v,usv,vsv,nx+nghost,ny+nghost);
			printf("\t Getting Vorticity\n");
			//get_vorticity(u,v,omega,dx,dy,nx,ny);
			printf("\t Getting Stream Function\n");
			//get_stream(strm,strmnext,omega,dx,dy,nx,ny);
			printf("\t Writing Data to Files\n");

			//Zalesak Disk Volumes
			double a = 0.005;
			init_circle(G1, nx, ny, nghost, dx, dy);
			double vol =  get_vol(G, a, nx, ny, nghost,dx,dy);
			double vol_act =  get_vol(G1, a, nx, ny, nghost,dx,dy);
			double shp_err = get_shape_err(G, G1, vol_act, a, nx, ny, nghost, dx, dy);
			printf("\t volume 1 rot = %7.7E \n", vol);
			printf("\t volume true = %7.7E \n", vol_act);
			printf("\t volume error = %7.7E \n", fabs(vol-vol_act));
			printf("\t shape error = %7.7E \n", shp_err);

			write_matrix_2d(usv, nx+2*nghost-1, ny+2*nghost-1, "u.dat");
			write_matrix_2d(vsv, nx+2*nghost-1, ny+2*nghost-1, "v.dat");
			write_matrix_2d(u,nx+2*nghost,ny+2*nghost+1,"uraw.dat");
			write_matrix_2d(v,nx+2*nghost+1,ny+2*nghost,"vraw.dat");
			write_matrix_2d(us,nx,ny+1,"us.dat");
			write_matrix_2d(vs,nx+1,ny,"vs.dat");
			write_matrix_2d(omega, nx, ny, "omega.dat");
			write_matrix_2d(phi,nx+1,ny+1,"phi.dat");
			write_matrix_2d(strm,nx+2,ny+2,"stream.dat");
			write_matrix_2d(G,nx+2*nghost-1,ny+2*nghost-1,"G.dat");
			write_matrix_2d(G1,nx+2*nghost-1,ny+2*nghost-1,"G1.dat");
			write_matrix_2d(dGdt,nx+2*nghost-1,ny+2*nghost-1,"dGdt.dat");
			write_matrix_2d_int(Markers,nx+2*nghost-1,ny+2*nghost-1,"Markers.dat");
		}
	} 
	
	printf("***\tDONE!\t***\n");
	return EXIT_SUCCESS;
}
	
/* Set Boundary Conditions */
void set_bcs(double** restrict u, double** restrict v, double dx, double dy, int  nx, int ny, double Re, struct slv_settings st)
{
#ifdef DBGBCS
		printf("\t Setting BCS \n");
		if(st.XBC == PERIODIC){
			printf("\t BCS is Periodic\n");
		}
#endif
	int i = 0;
	int j = 0;
#ifdef DBGBCS
		printf("\t Setting U Lower BC \n");
#endif
	#pragma omp parallel for
	for(i = 0; i<nx; i++){ //Lower Wall
		if(st.YBC == PERIODIC){
			u[i][j]=u[i][ny-1];
		} else {
			u[i][j]=u[i][j+1]-2*(u[i][j+1]-st.XLV);
		}
	}
#ifdef DBGBCS
		printf("\t Setting U Upper BC \n");
#endif
	j = ny;
	#pragma omp parallel for
	for(i = 0; i<nx; i++){ //Upper Wall
		if(st.YBC == PERIODIC){
			u[i][j]=u[i][1];
		} else {
			u[i][j]=u[i][j-1]-2*(u[i][j-1]-st.XUV);
		}
	}
#ifdef DBGBCS
		printf("\t Setting U Left BC \n");
#endif
	j = 0;
	i = 0;
	#pragma omp parallel for
	for(j = 0; j<ny+1; j++){ //Left Wall
		if(st.XBC == PERIODIC){
			u[i][j]=(u[1][j]+u[nx-1][j])/2;
		} else {
			u[i][j]=0.0;
		}
	}
#ifdef DBGBCS
		printf("\t Setting U Right BC \n");
#endif
	i = nx-1;
	#pragma omp parallel for
	for(j = 0; j<ny; j++){ //Right Wall
		if(st.XBC == PERIODIC){
			u[i][j]=(u[1][j]+u[nx-1][j])/2;
		} else {
			u[i][j]=0.0;
		}
	}

#ifdef DBGBCS
		printf("\t Setting V Lower BC \n");
#endif
	//V Velocity
	i = 0;
	j = 0;
	#pragma omp parallel for
	for(i = 0; i<nx+1; i++){ //Lower Wall
		if(st.YBC == PERIODIC){
			v[i][j]=(v[i][1]+v[i][ny])/2;
		} else {
			v[i][j]=0.0;
		}
	}
	j = ny-1;
#ifdef DBGBCS
		printf("\t Setting V Upper BC \n");
#endif
	#pragma omp parallel for
	for(i = 0; i<nx+1; i++){ //Upper Wall
		if(st.YBC == PERIODIC){
			v[i][j]=(v[i][1]+v[i][ny])/2;
		} else {
			v[i][j]=0.0;
		}
	}
#ifdef DBGBCS
		printf("\t Setting V Left BC \n");
#endif
	i = 0;
	#pragma omp parallel for
	for(j = 0; j<ny; j++){ //Left Wall
		if(st.XBC == PERIODIC){
			v[i][j] = v[nx][j];
		} else {
			v[i][j]=v[i+1][j]-2*(v[i+1][j]-st.YLV);
		}
	}
#ifdef DBGBCS
		printf("\t Setting V Right BC \n");
#endif
	i = nx;
	#pragma omp parallel for
	for(j = 0; j<ny; j++){ //Right Wall
		if(st.XBC == PERIODIC){
			v[i][j] = v[1][j];
		} else {
			v[i][j]=v[i-1][j]-2*(v[i-1][j]-st.YRV);
		}
	}
#ifdef DBGBCS
		printf("\t BCS Done \n");
#endif

	/*
	#pragma omp parallel for
	for(int i=0;i<nx+1;i++){
		for(int j=0;j<ny+1;j++){
			if (j==0){
				if(i!=nx){
					u[i][j]=-u[i][j+1];
				}
				v[i][j]=0.0;
			}
			if (j==ny-1){
				v[i][j]=0.0;
			}
			if(i==0){
				if(j!=ny){
					v[i][j]=-v[1][j];
				}
				u[i][j]=0.0;
			}
			if(i==nx){
				if(j!=ny){
					v[i][j]=-v[i-1][j];
				}
			}
			if(i==nx-1){
				u[i][j]=0.0;
			}	
			if (j==ny){
				if(i!=nx){
					u[i][j]=2-u[i][j-1];				
				}	
			}
		}
	}
	*/
}
/* Set Boundary Conditions */
void set_all_bcs_neumann(double** restrict x, double dx, double dy, int  nx, int ny, int nghostx, int nghosty){
#ifdef DBGBCS
		printf("\t Neumann BC Lower/Upper \n");
#endif
	#pragma omp parallel for
	for (int i = -nghostx; i<nx+nghostx-1; i++){ //Lower and Upper Boundary
		for (int j = 0;j < nghosty; j++){
			x[i+nghostx][j] = x[i+nghostx][2*nghosty-j]; //lower
			x[i+nghostx][ny+2*nghosty-2-j] = x[i+nghostx][ny-2+j]; //upper
		}
	}
#ifdef DBGBCS
		printf("\t Neumann BC Left/Right \n");
#endif
	#pragma omp parallel for
	for (int i = 0; i<nghostx; i++){ //Left and Right Boundary
		//printf("\t i = %i \n",i);
		for (int j = -nghosty; j<ny+nghosty-1; j++){
			//printf("\t x(%i, %i) = x(%i, %i) \n",nx+2*nghostx-2-i,j+nghosty,nx-1+i,j+nghosty);
			x[i][j+nghosty] = x[2*nghostx-i][j+nghosty];//left
			x[nx+2*nghostx-2-i][j+nghosty] = x[nx-2+i][j+nghosty];//right
		}
	}
#ifdef DBGBCS
		printf("\t Neumann BC Done! \n");
#endif
}

/* Get u and V */
void get_uv(double** restrict u, double** restrict v, double** restrict usv,double** restrict vsv,int nx, int ny)
{
	#pragma omp parallel for
	for(int i = 0; i<nx; i++){
		for(int j = 0;j<ny; j++){
			usv[i][j] = (u[i+1][j]+u[i][j])/2.0;
			vsv[i][j] = (v[i][j+1]+v[i][j])/2.0;
		}
	}	
}

/* Solve for vorticity */
void get_vorticity(double** u, double** v,double** omega, double dx, double dy, int nx, int ny)
{
	#pragma omp parallel for
	for(int i = 0;i<nx;i++){
		for(int j=0;j<ny;j++){
			omega[i][j] = (v[i+1][j]-v[i][j])/dx-(u[i][j+1]-u[i][j])/dy;
		}
	}
}

/* Solve for Stream Function */
void get_stream(double** restrict strm, double** restrict strmnext, double** restrict omega, double dx, double dy, int nx, int ny)
{	
	double dy2 = pow(dy,2.0);
	double dx2 = pow(dx,2.0);
	double resid = 100.0;
	for(int cnt=0;cnt<100000;cnt++){
		#pragma omp parallel for
		for(int i=1;i<nx;i++){
			for(int j=1;j<ny;j++){
				strmnext[i][j] = 1.0/(2.0/dx2+2.0/dy2)*((strm[i+1][j]+strm[i-1][j])/(dx2)+(strm[i][j+1]+strm[i][j-1])/(dy2))+1.0/(2.0/dx2+2.0/dy2)*omega[i-1][j-1];
			}
		}
		#pragma omp parallel for
		for(int i=1;i<nx;i++){
			for(int j=1;j<ny;j++){
				strm[i][j] = 1.0/(2.0/dx2+2.0/dy2)*((strmnext[i+1][j]+strmnext[i-1][j])/(dx2)+(strmnext[i][j+1]+strmnext[i][j-1])/(dy2))+1.0/(2.0/dx2+2.0/dy2)*omega[i-1][j-1];
			}
		}
		resid = 0.0;
		for(int i=1;i<nx;i++){
			for(int j=1;j<ny;j++){
				resid  = resid+fabs(strm[i][j]-strmnext[i][j]);
			}
		}
	}	
	printf("residual = %f \n",resid);
}


/* Viscous Burgers Equation Solver */
/* Uses AB2 + ADI */
void slv_vbe(double** restrict u, double** restrict us, double** restrict dus, double** restrict duss, double** restrict v, double** restrict vs, double** restrict dvs, double** restrict dvss, double** restrict hu, double** restrict huold,double** restrict hv, double** restrict hvold, double dx, double dy, int nx, int ny, double Re, double dt, struct slv_settings st)
{
	double dx2 = pow(dx,2.0);
	double dy2 = pow(dy,2.0);
	double ce = dt/(2.0*Re);
	#pragma omp parallel
	{	
		double* bum = malloc((nx-2)*sizeof(double));
		double* lum = malloc((nx-2)*sizeof(double));
		double* mum = malloc((nx-2)*sizeof(double));
		double* uum = malloc((nx-2)*sizeof(double));
		double* xum = malloc((nx-2)*sizeof(double));


		double* bvm = malloc((nx-1)*sizeof(double));
		double* lvm = malloc((nx-1)*sizeof(double));
		double* mvm = malloc((nx-1)*sizeof(double));
		double* uvm = malloc((nx-1)*sizeof(double));
		double* xvm = malloc((nx-1)*sizeof(double));
		
		#pragma omp parallel for
		for(int i=1;i<nx;i++){
			for(int j=1;j<ny;j++){	
				/* Do this for all u */					
				if(i!=nx-1){
					huold[i][j] = hu[i][j];
					hu[i][j] = -((u[i+1][j]+u[i][j])*(u[i+1][j]+u[i][j])-(u[i-1][j]+u[i][j])*(u[i-1][j]+u[i][j]))/(4.0*dx)-((u[i][j+1]+u[i][j])*(v[i][j]+v[i+1][j])-(u[i][j-1]+u[i][j])*(v[i][j-1]+v[i+1][j-1]))/(4.0*dy) + st.XConst;
				}		
				/* Do this for all v */
				if(j!=ny-1){
					hvold[i][j] = hv[i][j];
					hv[i][j] = -((v[i][j+1]+v[i][j])*(v[i][j+1]+v[i][j])-(v[i][j-1]+v[i][j])*(v[i][j-1]+v[i][j]))/(4.0*dy)-((v[i+1][j]+v[i][j])*(u[i][j]+u[i][j+1])-(v[i-1][j]+v[i][j])*(u[i-1][j]+u[i-1][j+1]))/(4.0*dx) + st.YConst;
				}
			}
		}

		
		#pragma omp parallel for
		for(int j=1;j<ny;j++){
			for(int i = 1; i<nx;i++){
			
				/* limit range for u */
				if (i == 1){
					bum[i-1] = dt/2.0*(3.0*hu[i][j]-huold[i][j])+2.0*ce*((u[i+1][j]-2.0*u[i][j]+u[i-1][j])/dx2+(u[i][j+1]-2.0*u[i][j]+u[i][j-1])/dy2);
					lum[i-1] = -ce/dx2;
					mum[i-1] = (1+2*ce/dx2);
					
					if(j!=ny-1){
						bvm[i-1] = dt/2.0*(3.0*hv[i][j]-hvold[i][j])+2.0*ce*((v[i+1][j]-2.0*v[i][j]+v[i-1][j])/dx2+(v[i][j+1]-2.0*v[i][j]+v[i][j-1])/dy2);
						lvm[i-1] = -ce/dx2;
						mvm[i-1] = (1.0+2.0*ce/dx2);
					}
				} else if (i == nx-2){
					bum[i-1] = dt/1.0*(2.0*hu[i][j]-huold[i][j])+2.0*ce*((u[i+1][j]-2.0*u[i][j]+u[i-1][j])/dx2+(u[i][j+1]-2.0*u[i][j]+u[i][j-1])/dy2);
					mum[i-1] = (1+2*ce/dx2);
					uum[i-1] = -ce/dx2;
					
					if(j!=ny-1){
						bvm[i-1] = dt/2.0*(3.0*hv[i][j]-hvold[i][j])+2.0*ce*((v[i+1][j]-2.0*v[i][j]+v[i-1][j])/dx2+(v[i][j+1]-2.0*v[i][j]+v[i][j-1])/dy2);
						lvm[i-1] = -ce/dx2;
						mvm[i-1] = (1.0+2.0*ce/dx2);
						uvm[i-1] = -ce/dx2;
					}
				} else if (i == nx-1){
					if(j!=ny-1){
						bvm[i-1] = dt/2.0*(3.0*hv[i][j]-hvold[i][j])+2.0*ce*((v[i+1][j]-2.0*v[i][j]+v[i-1][j])/dx2+(v[i][j+1]-2.0*v[i][j]+v[i][j-1])/dy2);
						mvm[i-1] = (1.0+2.0*ce/dx2);
						uvm[i-1] = -ce/dx2;
					}
				} else {
					
					bum[i-1] = dt/2.0*(3.0*hu[i][j]-huold[i][j])+2.0*ce*((u[i+1][j]-2.0*u[i][j]+u[i-1][j])/dx2+(u[i][j+1]-2.0*u[i][j]+u[i][j-1])/dy2);
					lum[i-1] = -ce/dx2;
					mum[i-1] = (1.0+2.0*ce/dx2);
					uum[i-1] = -ce/dx2;
	
					if(j!=ny-1){
						bvm[i-1] = dt/2.0*(3.0*hv[i][j]-hvold[i][j])+2.0*ce*((v[i+1][j]-2.0*v[i][j]+v[i-1][j])/dx2+(v[i][j+1]-2.0*v[i][j]+v[i][j-1])/dy2);
						lvm[i-1] = -ce/dx2;
						mvm[i-1] = (1.0+2.0*ce/dx2);
						uvm[i-1] = -ce/dx2;
					}
				}
			}
			
			solve_matrix(nx-2,lum,mum,uum,bum,xum);
			if(j!=ny-1){
				solve_matrix(nx-1,lvm,mvm,uvm,bvm,xvm);
			}			

			for(int i = 1;i<nx;i++){
				if(i!=nx-1){
					duss[i][j] = xum[i-1];
				}
				if(j!=ny-1){
					dvss[i][j] = xvm[i-1];
				}
			}	
		}
	}
	#pragma omp parallel for
	for(int i = 1; i<nx; i++){
		for(int j=1; j<ny; j++){
			if(j!=ny-1){
				vs[i][j] = v[i][j]+dvss[i][j];
			}
			if(i!=nx-1){
				us[i][j] = u[i][j]+duss[i][j];
			}
		}
	}
	
	set_bcs(us, vs, dx, dy, nx, ny, Re, st);
	
	#pragma omp parallel
	{
				double* bumy = malloc((nx-1)*sizeof(double));
                double* lumy = malloc((nx-1)*sizeof(double));
                double* mumy = malloc((nx-1)*sizeof(double));
                double* uumy = malloc((nx-1)*sizeof(double));
                double* xumy = malloc((nx-1)*sizeof(double));

                double* bvmy = malloc((nx-2)*sizeof(double));
                double* lvmy = malloc((nx-2)*sizeof(double));
                double* mvmy = malloc((nx-2)*sizeof(double));
                double* uvmy = malloc((nx-2)*sizeof(double));
                double* xvmy = malloc((nx-2)*sizeof(double));		

		for(int i = 1; i<nx; i++){
			for(int j = 1; j<ny; j++){
				/* limit range for v */
				if (j == 1){
					if(i!=nx-1){	
						bumy[j-1] = duss[i][j];	
						bumy[j-1] = dt/2.0*(3.0*hu[i][j]-huold[i][j])+2.0*ce*((us[i+1][j]-2.0*us[i][j]+us[i-1][j])/dx2+(us[i][j+1]-2.0*us[i][j]+us[i][j-1])/dy2);
						lumy[j-1] = -ce/dy2;
						mumy[j-1] = (1.0+2.0*ce/dy2);
					}
					bvmy[j-1] = dvss[i][j];		
					bvmy[j-1] = dt/2.0*(3.0*hv[i][j]-hvold[i][j])+2.0*ce*((vs[i+1][j]-2.0*vs[i][j]+vs[i-1][j])/dx2+(vs[i][j+1]-2.0*vs[i][j]+vs[i][j-1])/dy2);
					lvmy[j-1] = -ce/dy2;
					mvmy[j-1] = (1.0+2.0*ce/dy2);
				} else if (j == ny-2){
					if(i!=nx-1){
						bumy[j-1] = duss[i][j];	
						bumy[j-1] = dt/2.0*(3.0*hu[i][j]-huold[i][j])+2.0*ce*((us[i+1][j]-2.0*us[i][j]+us[i-1][j])/dx2+(us[i][j+1]-2.0*us[i][j]+us[i][j-1])/dy2);
						lumy[j-1] = -ce/dy2;
						mumy[j-1] = (1.0+2.0*ce/dy2);
						uumy[j-1] = -ce/dy2;
					}
					bvmy[j-1] = dvss[i][j];		
					bvmy[j-1] = dt/2.0*(3.0*hv[i][j]-hvold[i][j])+2.0*ce*((vs[i+1][j]-2.0*vs[i][j]+vs[i-1][j])/dx2+(vs[i][j+1]-2.0*vs[i][j]+vs[i][j-1])/dy2);
					mvmy[j-1] = (1.0+2.0*ce/dy2);
					uvmy[j-1] = -ce/dy2;
				} else if (j == ny-1){
					if(i!=nx-1){
						bumy[j-1] = duss[i][j];	
						bumy[j-1] = dt/2.0*(3.0*hu[i][j]-huold[i][j])+2.0*ce*((us[i+1][j]-2.0*us[i][j]+us[i-1][j])/dx2+(us[i][j+1]-2.0*us[i][j]+us[i][j-1])/dy2);
						mumy[j-1] = (1.0+2.0*ce/dy2);
						uumy[j-1] = -ce/dy2;
					}
				} else {
					if(i!=nx-1){	
						bumy[j-1] = duss[i][j];	
						bumy[j-1] = dt/2.0*(3.0*hu[i][j]-huold[i][j])+2.0*ce*((us[i+1][j]-2.0*us[i][j]+us[i-1][j])/dx2+(us[i][j+1]-2.0*us[i][j]+us[i][j-1])/dy2);
						lumy[j-1] = -ce/dy2;
						mumy[j-1] = (1.0+2.0*ce/dy2);
						uumy[j-1] = -ce/dy2;
					}
					bvmy[j-1] = dvss[i][j];	
					bvmy[j-1] = dt/2.0*(3.0*hv[i][j]-hvold[i][j])+2.0*ce*((vs[i+1][j]-2.0*vs[i][j]+vs[i-1][j])/dx2+(vs[i][j+1]-2.0*vs[i][j]+vs[i][j-1])/dy2);
					lvmy[j-1] = -ce/dy2;
					mvmy[j-1] = (1.0+2.0*ce/dy2);
					uvmy[j-1] = -ce/dy2;
				}
			}
			if(i!=nx-1){
				solve_matrix(ny-1,lumy,mumy,uumy,bumy,xumy);
			}
			solve_matrix(ny-2,lvmy,mvmy,uvmy,bvmy,xvmy);		
			for(int j = 1;j<ny;j++){
				if(j!=ny-1){
					dvs[i][j] = xvmy[j-1];
				}
				if(i!=nx-1){
					dus[i][j] = xumy[j-1];
				}
			}	
		}
	}
	#pragma omp parallel for
	for(int i = 1; i<nx; i++){
		for(int j=1; j<ny; j++){
			if(j!=ny-1){
				vs[i][j] = vs[i][j]+dvs[i][j];
			}
			if(i!=nx-1){
				us[i][j] = us[i][j]+dus[i][j];
			}
		}
	}
	/*
	#pragma omp parallel for
	for(int i = 0;i<nx+1;i++){
		vs[i][0] = v[i][0];	
		vs[i][ny-1] = v[i][ny-1];
		if(i!=nx){
			us[i][0] = u[i][0];
			us[i][ny] = u[i][ny];
		}
	}
	for(int j=0;j<ny+1;j++){
		us[0][j] = u[0][j];
		us[nx-1][j] = u[nx-1][j];
		if(j!=ny){
			vs[0][j] = v[0][j];
			vs[nx][j] = v[nx][j];
		}
	}
	*/
}

/* ftcs solver */
void slv_vbeftcs(double** restrict u, double** restrict us, double** restrict dus, double** restrict duss, double** restrict v, double** restrict vs, double** restrict dvs, double** restrict dvss, double** restrict hu, double** restrict huold,double** restrict hv, double** restrict hvold, double dx, double dy, int nx, int ny, double Re, double dt, struct slv_settings st)
{
	
	double dx2 = pow(dx,2.0);
	double dy2 = pow(dy,2.0);
	double ce = dt/(Re);
	#pragma omp parallel for
	for(int i=1;i<nx;i++){
		for(int j=1;j<ny;j++){
			/* Do this for all u */					
			if(i!=nx-1){
				huold[i][j] = hu[i][j];
				hu[i][j] = -((u[i+1][j]+u[i][j])*(u[i+1][j]+u[i][j])-(u[i-1][j]+u[i][j])*(u[i-1][j]+u[i][j]))/(4.0*dx)-((u[i][j+1]+u[i][j])*(v[i][j]+v[i+1][j])-(u[i][j-1]+u[i][j])*(v[i][j-1]+v[i+1][j-1]))/(4.0*dy) + st.XConst;
			}		
			/* Do this for all v */
			if(j!=ny-1){
				hvold[i][j] = hv[i][j];
				hv[i][j] = -((v[i][j+1]+v[i][j])*(v[i][j+1]+v[i][j])-(v[i][j-1]+v[i][j])*(v[i][j-1]+v[i][j]))/(4.0*dy)-((v[i+1][j]+v[i][j])*(u[i][j]+u[i][j+1])-(v[i-1][j]+v[i][j])*(u[i-1][j]+u[i-1][j+1]))/(4.0*dx) + st.YConst;
			}
		}
	}
	#pragma omp parallel for
	for(int i=1;i<nx;i++){
		for(int j = 1; j<nx;j++){
			if(i!=nx-1){	
				us[i][j] =u[i][j]+dt*(hu[i][j])+ce*((u[i+1][j]-2.0*u[i][j]+u[i-1][j])/dx2+(u[i][j+1]-2.0*u[i][j]+u[i][j-1])/dy2);
			}
			if(j!=ny-1){
				vs[i][j] = v[i][j]+dt*(hv[i][j])+ce*((v[i+1][j]-2.0*v[i][j]+v[i-1][j])/dx2+(v[i][j+1]-2.0*v[i][j]+v[i][j-1])/dy2);
			}
		}	
	}
	#pragma omp parallel for
	for(int i = 0;i<nx+1;i++){
		vs[i][0] = v[i][0];	
		vs[i][ny-1] = v[i][ny-1];
		if(i!=nx){
			us[i][0] = u[i][0];
			us[i][ny] = u[i][ny];
		}
	}
	#pragma omp parallel for
	for(int j=0;j<ny+1;j++){
		us[0][j] = u[0][j];
		us[nx-1][j] = u[nx-1][j];
		if(j!=ny){
			vs[0][j] = v[0][j];
			vs[nx][j] = v[nx][j];
		}
	}
}


/* Poisson Solver */
void slv_pssn(double** restrict phi,double** restrict phinext, double** restrict us, double** restrict vs, double dx, double dy, int nx, int ny,double dt,double min)
{
	double dx2 = pow(dx,2.0);
	double dy2 = pow(dy,2.0);
	double save = 0.00;
	double resid = 100.0;
	int cnt = 0;
	while(resid>min){
		cnt = cnt+1;
		#pragma omp parallel for
		for(int i = 0;i<nx+1;i++){
			phi[i][0] = phi[i][1];
			phi[i][ny] = phi[i][ny-1];
		}
		#pragma omp parallel for 
		for(int j = 0;j<ny+1;j++){
			phi[0][j] = phi[1][j];
			phi[nx][j] = phi[nx-1][j];
		}
		#pragma omp parallel for
		for(int i = 1; i<nx;i++){
			for(int j = 1; j<ny; j++){
				phinext[i][j] = 1.0/(2.0/dx2+2.0/dy2)*((phi[i+1][j]+phi[i-1][j])/dx2+(phi[i][j+1]+phi[i][j-1])/dy2)-1.0/(dt*(2.0/dx2+2.0/dy2))*((us[i][j]-us[i-1][j])/dx+(vs[i][j]-vs[i][j-1])/dy);
			}
		}
		#pragma omp parallel for
		for(int i = 0;i<nx+1;i++){
			phinext[i][0] = phinext[i][1];
			phinext[i][ny] = phinext[i][ny-1];
		}
		#pragma omp parallel for 
		for(int j = 0;j<ny+1;j++){
			phinext[0][j] = phinext[1][j];
			phinext[nx][j] = phinext[nx-1][j];
		}
		#pragma omp parallel for //Solve system
		for(int i = 1; i<nx;i++){
				for(int j = 1; j<ny; j++){
					phi[i][j] = 1.0/(2.0/dx2+2.0/dy2)*((phinext[i+1][j]+phinext[i-1][j])/dx2+(phinext[i][j+1]+phinext[i][j-1])/dy2)-1.0/(dt*(2.0/dx2+2.0/dy2))*((us[i][j]-us[i-1][j])/dx+(vs[i][j]-vs[i][j-1])/dy);
				}
		}
		resid = 0.0;
		for(int i=0; i<nx; i++){
			for(int j=0;j<ny;j++){
				save = fabs(phi[i][j]-phinext[i][j]); //Infinity Norm
				if (save>resid){
					resid = save;
				} 	
			}
		}
	}
#ifdef DEBUG
	printf("\t\t Poisson Iterations: %i\n",cnt);
#endif
}

/* Poisson Solver */
void slv_pssn_gmres(double** restrict phi,double** restrict phinext, double** restrict us, double** restrict vs, double dx, double dy, int nx, int ny,double dt,double min, struct slv_settings st)
{
#ifdef DEBUG
	printf("\n Mgmres poisson call \n");
#endif
	double dx2 = pow(dx,2.0);
	double dy2 = pow(dy,2.0);
	int ic = 0;
	int jc = 0;
	int n = (nx+1)*(ny+1); //Rank of Linear System
	int nz_num = (2*nx+2*ny) + nx*ny*5; //the number of nonzero matrix values.
	double* A = malloc(nz_num*sizeof(double)); //Matrix Values
	int* iA = malloc(nz_num*sizeof(int)); //Matrix Row Index
	int* jA = malloc(nz_num*sizeof(int)); //Matrix Column Index
	double* x = malloc(n*sizeof(double)); //An approximation to the solution.  On output, an improved approximation.
    double* rhs = malloc(n*sizeof(double)); //The right hand side of the linear system.
    int itr_max = 5; //The maximum number of (outer) iterations to take.
    int mr = 5; //the maximum number of (inner) iterations to take.
    double tol_abs = 0.001; //An absolute tolerance applied to the current residual.
    double tol_rel = 0.001; //A relative tolerance comparing the current residual to the initial residual.
	int cnt = 0;
	int cntrhs =0;
	//Ground lower corner
	A[cnt] = 1.0l; iA[cnt] = 0; jA[cnt] = 0; cnt++;
	rhs[cntrhs] = 0.0;cntrhs++;
	int i = 0;
	int j = 0;
	for (j = 0; j < ny+1 ; j++){ //Bottom Wall
		if (st.YBC == PERIODIC){
			A[cnt] =  1.0l; iA[cnt] = i*(ny+1)+j; jA[cnt] = i*(ny+1)+j;    cnt++;//This value...
			A[cnt] = -1.0l; iA[cnt] = i*(ny+1)+j; jA[cnt] = (nx)*(ny+1)+j; cnt++;//and its periodic neighbor...
			rhs[cntrhs] = 0.0; cntrhs++;                                         //are equal
		} else {
			A[cnt] = 1.0; iA[cnt] = i*(ny+1)+j; jA[cnt] = i*(ny+1)+j; cnt++; //This value...
			A[cnt] = 1.0; iA[cnt] = i*(ny+1)+j; jA[cnt] = (i+1)*(ny+1)+j; cnt++;   //and its inside neighbor...
			rhs[cntrhs] = 0.00; cntrhs++;                                    //are equal
		}
	}
	i = 0;
	j = 0;
	for (i = 0; i < nx+1 ; i++){//Left Wall
		if (st.XBC == PERIODIC){
			A[cnt] =  1.0l; iA[cnt] = i*(ny+1)+j; jA[cnt] = i*(ny+1)+j;  cnt++;//This value...
			A[cnt] = -1.0l; iA[cnt] = i*(ny+1)+j; jA[cnt] = i*(ny+1)+ny; cnt++;//and its periodic neighbor...
			rhs[cntrhs] = 0.0; cntrhs++;                                       //are equal
		} else {
			A[cnt] = 1.0; iA[cnt] = i*(ny+1)+j; jA[cnt] = i*(ny+1)+j;   cnt++; //This value...
			A[cnt] = 1.0; iA[cnt] = i*(ny+1)+j; jA[cnt] = i*(ny+1)+j+1; cnt++; //and its inside neighbor...
			rhs[cntrhs] = 0.00; cntrhs++;                                      //are equal
		}
	}
#ifdef DEBUG
	printf("\n Equations from BCs - %i", cntrhs);
#endif
	double cnst = -1.0/(2.0/dx2+2.0/dy2);
	for (i = 1; i<nx; i++){//Inside Points
		for (j =1; j<ny; j++){
			A[cnt] =  1.0l; iA[cnt] = (i)*(ny+1)+j; jA[cnt] = i*(ny+1)+(j);  cnt++;  //This value...
			A[cnt] = cnst/dy2; iA[cnt] = (i)*(ny+1)+j; jA[cnt] = i*(ny+1)+(j+1); cnt++; //Point +j...
			A[cnt] = cnst/dy2; iA[cnt] = (i)*(ny+1)+j; jA[cnt] = i*(ny+1)+(j-1); cnt++; //Point -j...
			A[cnt] = cnst/dx2; iA[cnt] = (i)*(ny+1)+j; jA[cnt] = (i+1)*(ny+1)+(j); cnt++; //Point +i...
			A[cnt] = cnst/dx2; iA[cnt] = (i)*(ny+1)+j; jA[cnt] = (i-1)*(ny+1)+(j-1); cnt++; //Point -i...
			rhs[cntrhs] = -cnst/dt*(((us[i][j]-us[i-1][j])/dx+(vs[i][j]-vs[i][j-1])/dy)); cntrhs++; //Equal to ustar and vstar forcing
		}
	}
	//mgmres_st (int n, int nz_num, int ia[], int ja[], double a[], double x[], double rhs[], int itr_max, int mr, double tol_abs, double tol_rel )
#ifdef DEBUG
	printf("\n Calling mgmres \n");
	printf("\n Predicted %i rank - counted %i", n, cntrhs);
	printf("\n Predicted %i A entries - counted %i \n", nz_num, cnt);
#endif
#ifdef DBGMGM
	for (int i = 0; i<n; i++){
		printf("\t|");
		for (int j = 0; j < n; j++){
			printf("%+2.1f, ",find_in_A(A,iA,jA,cnt,i,j));
		}
		printf("| = |%+5.3f| \n", rhs[i]);
	}
#endif
	mgmres_st (n, nz_num, iA, jA, A, x, rhs, itr_max, mr, tol_abs,tol_rel);
#ifdef DEBUG
	printf("\n Returned from mgmres call \n");
#endif

	free(A);
	free(iA);
	free(jA);
	free(x);
	free(rhs);
}
double find_in_A (double* A, int* iA, int* jA, int npt, int i, int j){
	for(int a = 0; a < npt; a++){
		if(iA[a] == i && jA[a] == j){
			return A[a];
		}
	}
	return 0.0l;
}
void apply_projection(double** restrict phi, double** restrict u, double** restrict us, double** restrict v, double** restrict vs, double dx, double dy, int nx, int ny, double dt)
{
	#pragma omp parallel for
	for(int i=1;i<nx;i++){
		for(int j=1;j<ny;j++){
			if(i!=nx-1){
				u[i][j] = us[i][j]-dt*(phi[i+1][j]-phi[i][j])/dx;
			}
			if(j!=ny-1){
				v[i][j] = vs[i][j]-dt*(phi[i][j+1]-phi[i][j])/dy;
			}
		}
	}	
}

/* Get weight for WENO5 term */
double slv_psi_weno(double a, double b, double c, double d){
#ifdef DBGWENO
	printf("\t Calculating Psi WENO\n");
#endif
	double is0 = 13.0*(a-b)*(a-b)+3.0*(a-3.0*b)*(a-3.0*b);
	double is1 = 13.0*(b-c)*(b-c)+3.0*(b+c)*(b+c);
	double is2 = 13.0*(c-d)*(c-d)+3.0*(3.0*c-d)*(3.0*c-d);

	double eps = 0.0000000001; //1^-10

	double a0 = 1.0/((eps+is0)*(eps+is0));
	double a1 = 6.0/((eps+is1)*(eps+is1));
	double a2 = 3.0/((eps+is2)*(eps+is2));

	double w0 = a0/(a0+a1+a2);
	double w2 = a2/(a0+a1+a2);

	return (1.0/3.0)*w0*(a-2.0*b+c)+1.0/6.0*(w2-1.0/2.0)*(b-2.0*c+d);
}
/* Init Settings */
struct slv_settings init_settings(){
	struct slv_settings st;
	st.nx = 40;
	st.ny = 40;
	st.nt = 1000;
	st.Re = 100.0;
	st.dt = 0.001;
	st.XBC = WALL;
	st.YBC = WALL;
	st.XUV = 1.0l;
	st.XLV = 0.0l;
	st.YLV = 0.0l;
	st.YRV = 0.0l;
	return st;
}
	
/* Write 2D Matrix to File */
void write_matrix_2d(double** mat, int nx, int ny, char* filename)
{
	FILE *file = fopen(filename,"w");
	size_t count;
	for (int i = 0; i<nx; i++){
		count = fwrite(mat[i],sizeof(double),ny,file);
	}
}

/* Write 2D Matrix to File */
void write_matrix_2d_int(int** mat, int nx, int ny, char* filename)
{
	FILE *file = fopen(filename,"w");
	size_t count;
	for (int i = 0; i<nx; i++){
		count = fwrite(mat[i],sizeof(int),ny,file);
	}
}

void print_array(double* X, int n){
	printf("\n X = {");
	for(int i = 0; i < n; i++){
		printf("%f,",X[i]);
	}
	printf("}\n");
}
void print_array_int(int* X, int n){
	printf("\n X = {");
	for(int i = 0; i < n; i++){
		printf("%2.2i,",X[i]);
	}
	printf("}\n");
}
