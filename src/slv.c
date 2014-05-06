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
	int ninit = 100;
	double Re = 100.0;
	double dt = 0.001;
	double time = 0.0;
	double tend = 0.0001;
	double maxv = 0.0;
	struct slv_settings st;
	st = init_settings();
	st.Re = 100.0;
	st.YBC = SLIPWALL;
	st.XBC = PERIODIC;
	st.XUV =0.0;
	st.XLV = 0.0;
	st.YLV = 0.0;
	st.YRV = 0.0;
	st.XConst = 0.0;
	st.YConst = 0.0;
	//st.rhol = 1.225;
	st.rhol = 1.0;
	//st.rhog = 0.1694;
	st.rhog = 0.1;
	//st.mul = 0.00313;
	st.mul = 1.0;
	//st.mug = 0.00313;
	st.mug = .00181;
	st.Reg = st.rhog/st.mug;
	st.Rel = st.rhol/st.mul;

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
			st.Re = atof(argv[i+1]);
			i++;
		} else if (strcmp(argv[i],"-tf")==0){
			tend = atof(argv[i+1]);
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
	printf("\tUnit Velocity Re Fluid: %f\n",st.rhol/st.mul);
	printf("\tUnit Velocity Re Gas: %f\n",st.rhog/st.mug);
	printf("\tTime Final: %f\n",tend);
	printf("\tTime Iterations: %f\n", tend/dt);


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

	double** alpha = malloc((nx+2*nghost-1)*sizeof(double*));memcntr += (nx+2*nghost-1)*sizeof(double*);
	double** d = malloc((nx+2*nghost-1)*sizeof(double*));memcntr += (nx+2*nghost-1)*sizeof(double*);
	V2D**  marr   = malloc((nx+2*nghost-1)*sizeof(V2D*));memcntr += (nx+2*nghost-1)*sizeof(V2D*);
	lseg**  segs   = malloc((nx+2*nghost-1)*sizeof(lseg*));memcntr += (nx+2*nghost-1)*sizeof(lseg*);

	double* FMM = malloc((nx+2*nghost-1)*(ny+2*nghost-1)*sizeof(double)); memcntr += (nx+2*nghost-1)*(ny+2*nghost-1)*sizeof(double);
	int** Markers = malloc((nx+2*nghost-1)*sizeof(int*));memcntr += (nx+2*nghost-1)*sizeof(int*);
	int** MarkPoint = malloc((nx+2*nghost-1)*sizeof(int*));memcntr += (nx+2*nghost-1)*sizeof(int*);
	int* FMMi = malloc((nx+2*nghost-1)*(ny+2*nghost-1)*sizeof(int)); memcntr += (nx+2*nghost-1)*(ny+2*nghost-1)*sizeof(int);
	int* FMMj = malloc((nx+2*nghost-1)*(ny+2*nghost-1)*sizeof(int)); memcntr += (nx+2*nghost-1)*(ny+2*nghost-1)*sizeof(int);


	double** omega = malloc(nx*sizeof(double*));memcntr += nx*sizeof(double*);

	double* sysx = malloc((nx-1)*(ny-1)*sizeof(double));memcntr += (nx-1)*(ny-1)*sizeof(double*);

	double** phi = malloc((nx+2*nghost-1)*sizeof(double*));memcntr += (nx+2*nghost-1)*sizeof(double*);
	double** phinext = malloc((nx+2*nghost-1)*sizeof(double*));memcntr += (nx+2*nghost-1)*sizeof(double*);
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
		MarkPoint[i+nghost] = malloc((ny+2*nghost-1)*sizeof(int));memcntr += (ny+2*nghost-1)*sizeof(int);

		alpha[i+nghost] = malloc((nx+2*nghost-1)*sizeof(double));memcntr += (nx+2*nghost-1)*sizeof(double);
		d[i+nghost] = malloc((nx+2*nghost-1)*sizeof(double));memcntr += (nx+2*nghost-1)*sizeof(double);
		marr[i+nghost] = malloc((nx+2*nghost-1)*sizeof(V2D));memcntr += (nx+2*nghost-1)*sizeof(V2D);
		segs[i+nghost] = malloc((nx+2*nghost-1)*sizeof(lseg));memcntr += (nx+2*nghost-1)*sizeof(lseg);
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
			phi[i+nghost] = malloc((ny+2*nghost-1)*sizeof(double));memcntr += (ny+2*nghost-1)*sizeof(double);
			phinext[i+nghost] = malloc((ny+2*nghost-1)*sizeof(double));memcntr += (ny+2*nghost-1)*sizeof(double);
	
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
	//init_circle(G, nx, ny, nghost, dx, dy);
	//init_RT(G,nx,ny,nghost,dx,dy);
	init_drop(G,nx,ny,nghost,dx,dy);
	//init_allgas(G,nx,ny,nghost,dx,dy);
	printf("\n\t Done Initializing!\n\n");
#ifdef DEBUG
	printf("\n\t Setting Neumann BC for G\n\n");
#endif DEBUG
	set_all_bcs_neumann(G,dx,dy,nx,ny,nghost,nghost);
	//set_all_BCS(G,dx,dy,nx,ny,nghost,st);
	//init_uv_test(u, v, nx, ny, nghost, dx, dy);
	get_alpha(G, alpha,  nghost, dx, dy, nx,  ny);
	//set_all_bcs_neumann(alpha,dx,dy,nx,ny,nghost,nghost);
	//init_uv_test_t3(u, v, nx, ny, nghost, dx, dy,0.0);
	//init_uv_test_down(u, v, nx, ny, nghost, dx, dy);
	init_uv_test_static(u, v, nx, ny, nghost, dx, dy);


	double min = 0.0000000001;
	int t = 0;
	/* Time Iteration Loop */
	//for(int t = 0;t<nt;t++){
	while (tend >= time && fabs(tend-time) > 0.0000000001 ){//|| t > 1000){
		if(t > 1){
			min = 0.0000000001;
		}
		if (t == 0){
			dt = 0.00001;
		}else {
			dt = fmin(0.01,0.5/(2.0)*fmin(st.Reg,st.Rel)*fmin(dx*dx,dy*dy)); // Stability limit for FTCS
			dt = fmin(dt, 0.25/(2.0*maxv)*fmin(dx,dy)); //Stability limit for Flux
			//dt = fmin(0.01, 0.5/(2.0*maxv)*fmin(dx,dy)); //Stability limit for Flux
		}
		if((tend - time) < dt){
			dt = tend-time;
		}
		printf("\t T = %9.9f, dt =  %9.9f \n",time, dt);
		if(t%10==0){
			printf("\t %i / %i \n",t,nt);
			printf("\t Max Velocity is %f, CFL(Flux) =  %f, CFL(viscous) =  %f \n",maxv,2.0*maxv*dt/fmin(dx,dy),2.0*dt/(fmin(st.Reg,st.Rel)*fmin(dx*dx,dy*dy)));

		}
#ifdef DEBUG
		printf("\t Simulation Loop. Time: %9.9f / %9.9f \n",(double)t*dt,dt*nt);
#endif

		/* Set BCS */
#ifdef DEBUG
		printf("\t Setting BCS \n");
#endif
		set_bcs(u, v, dx, dy, nx, ny,nghost, st);

		/* Solv Viscous Bergers Equation */
#ifdef DEBUG
		printf("\t Solving Viscous Bergers \n");
#endif
		//slv_vbe(alpha,u,us,dus,duss,v,vs,dvs,dvss,hu,huold,hv,hvold,dx,dy,nx,ny,nghost,dt,st);
		//slv_vbe(alpha,u,us,dus,duss,v,vs,dvs,dvss,hu,huold,hv,hvold,dx,dy,nx,ny,nghost,dt,st);
		slv_vbeftcs(alpha,u,us,dus,duss,v,vs,dvs,dvss,hu,huold,hv,hvold,dx,dy,nx,ny,nghost,dt,st);
		/* Calculate phi */	
#ifdef DEBUG
		printf("\t Solving Poisson Equation \n");
#endif
		//slv_pssn(phi,phinext,alpha,us,vs,dx,dy,nx,ny,nghost,dt,min,st);
		slv_pssn_gmres(phi,sysx,alpha,us,vs,dx,dy,nx,ny,nghost,dt,min,st);
		//write_matrix_2d(phi,nx+2*nghost-1,ny+2*nghost-1,"phi.dat");
		//write_matrix_2d(us,nx,ny+1,"us.dat");
		//write_matrix_2d(vs,nx+1,ny,"vs.dat");
		//slv_pssn_gmres(phi,phinext,us,vs,dx,dy,nx,ny,dt,min,st);

		/* Apply Projection */
#ifdef DEBUG
		printf("\t Applying Projection \n");
#endif
		maxv = apply_projection(phi,alpha,u,us,v,vs,dx,dy,nx,ny,nghost,dt,st);
		set_bcs(u, v, dx, dy, nx, ny,nghost, st);
		//apply_projection(phi,alpha,u,us,v,vs,dx,dy,nx,ny,nghost,dt,st);
#ifdef DEBUG
		printf("\t Advecting G \n");
#endif
		//init_uv_test_t3(u, v, nx, ny, nghost, dx, dy,time);
		//vof_get_fluxes(u,v,alpha,d,marr,dt,nghost, dx, dy, nx, ny,t);
		levelset_advect_TVDRK3(G,G1,G2,dGdt,dGdt1,dGdt2,dGdxp,dGdxm,dGdyp,dGdym,u,v,dx,dy,dt,time,nx,ny,nghost,st);
		set_all_bcs_neumann(G,dx,dy,nx,ny,nghost,nghost);
		//set_all_BCS(G,dx,dy,nx,ny,nghost,st);
		//set_all_bcs_neumann(G,dx,dy,nx,ny,nghost,nghost);
		get_alpha(G, alpha,  nghost, dx, dy, nx,  ny);
		//set_all_bcs_neumann(alpha,dx,dy,nx,ny,nghost,nghost);
		//set_all_BCS(alpha,dx,dy,nx,ny,nghost,st);

		if((t+1)%ninit == 0){
			copy_2D(G0, G, nx+2*nghost-1, ny+2*nghost-1);

			reinitl2 = reinit_advect_TVDRK3(G,G0,G1,G2,dGdt,dGdt1,dGdt2,dGdxp,dGdxm,dGdyp,dGdym,dx,dy,fmin(dx,dy)/8.0,nx,ny,nghost,st);
			int ind = 1;
			while(reinitl2 > dt*fmin(dx,dy)*fmin(dx,dy) && ind < 5){
				reinitl2 = reinit_advect_TVDRK3(G,G0,G1,G2,dGdt,dGdt1,dGdt2,dGdxp,dGdxm,dGdyp,dGdym,dx,dy,fmin(dx,dy)/8.0,nx,ny,nghost,st);
				set_all_bcs_neumann(G,dx,dy,nx,ny,nghost,nghost);
				ind++;
			}
			printf("Reinit converged at step %i:\t L2 err = %7.7E\n",ind, reinitl2);
			//reinit_FMM(G,G0,Markers, MarkPoint,FMM,FMMi,FMMj, dx, dy, dt, nx, ny,nghost);
			set_all_bcs_neumann(G,dx,dy,nx,ny,nghost,nghost);
			//printf("Ran FMM\n");
		}

		//get_alpha(G, alpha,  nghost, dx, dy, nx, ny);
		//vof_get_surfaces(alpha,d,marr,segs,nghost,dx,dy,nx,ny);

		if( ((t+1)%nsave) == 0) {
			char fname[50];
			//strcopy(fname,"u.");
			sprintf(&fname,"G.%i",saves);
			write_matrix_2d(G,nx+2*nghost-1,ny+2*nghost-1,fname);
			saves++;
		}

		/* Repeat */
		if (time+dt == tend){
			get_uv(u,v,usv,vsv,nx+2*nghost-1,ny+2*nghost-1);
			printf("\t Getting Vorticity\n");
			//get_vorticity(u,v,omega,dx,dy,nx,ny);
			printf("\t Getting Stream Function\n");
			//get_stream(strm,strmnext,omega,dx,dy,nx,ny);
			printf("\t Writing Data to Files\n");

			//Zalesak Disk Volumes
			double a = 0.005;
			//init_zalesak(G1, nx, ny, nghost, dx, dy);
			init_circle(G1, nx, ny, nghost, dx, dy);
			double vol =  get_vol(G, a, nx, ny, nghost,dx,dy);
			vol = get_alpha_sum(alpha,nghost,dx,dy,nx,ny);


			double vol_act =  get_vol(G1, a, nx, ny, nghost,dx,dy);
			get_alpha(G1, G2,  nghost, dx, dy, nx,  ny);
			vol_act =  get_alpha_sum(G2,nghost,dx,dy,nx,ny);

			double shp_err = get_shape_err(G, G1, vol_act, a, nx, ny, nghost, dx, dy);
			shp_err = get_alpha_diff(alpha, G2,nghost, dx,dy, nx, ny);

			printf("\t volume 1 rot = %16.16E \n", vol);
			printf("\t volume true = %16.16E \n", vol_act);
			printf("\t volume error = %16.16E \n", fabs(vol-vol_act));
			printf("\t shape error = %16.16E \n", shp_err);

			write_matrix_2d(usv, nx+2*nghost-1, ny+2*nghost-1, "u.dat");
			write_matrix_2d(vsv, nx+2*nghost-1, ny+2*nghost-1, "v.dat");
			write_matrix_2d(u,nx+2*nghost,ny+2*nghost+1,"uraw.dat");
			write_matrix_2d(v,nx+2*nghost+1,ny+2*nghost,"vraw.dat");
			write_matrix_2d(us,nx+2*nghost,ny+2*nghost+1,"us.dat");
			write_matrix_2d(vs,nx+2*nghost+1,ny+2*nghost,"vs.dat");

			write_matrix_2d(omega, nx, ny, "omega.dat");
			//write_matrix_2d(phi,nx+1,ny+1,"phi.dat");
			write_matrix_2d(phi,nx+2*nghost-1,ny+2*nghost-1,"phi.dat");
			write_matrix_2d(phinext,nx+2*nghost-1,ny+2*nghost-1,"phinext.dat");
			write_matrix_2d(strm,nx+2,ny+2,"stream.dat");

			write_matrix_2d(G,nx+2*nghost-1,ny+2*nghost-1,"G.dat");
			write_matrix_2d(G1,nx+2*nghost-1,ny+2*nghost-1,"G1.dat");
			write_matrix_2d(dGdt,nx+2*nghost-1,ny+2*nghost-1,"dGdt.dat");

			write_matrix_2d(alpha,nx+2*nghost-1,ny+2*nghost-1,"alpha.dat");
			write_matrix_2d(d,nx+2*nghost-1,ny+2*nghost-1,"d.dat");
			write_matrix_2d_lseg(segs, nx+2*nghost-1,ny+2*nghost-1, "lsegs.dat");

			write_matrix_2d_int(Markers,nx+2*nghost-1,ny+2*nghost-1,"Markers.dat");

		}
		time = time+dt;
	    t = t+1;
	} 
	
	printf("***\tDONE!\t***\n");
	return EXIT_SUCCESS;
}
	
/* Set Boundary Conditions */
void set_bcs(double** restrict u, double** restrict v, double dx, double dy, int  nx, int ny, int nghost, struct slv_settings st)
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
			u[i+nghost][j+nghost-1]=u[i+nghost][ny-2+nghost]; // lower out of bounds is upper in bounds
		} else if (st.YBC == SLIPWALL){
			u[i+nghost][j+nghost-1]=u[i+nghost][j+nghost]; // lower out of bounds is lower in bounds
		}else {
			//u[i+nghost][j+nghost]=u[i+nghost][j+1+nghost]-2*(u[i+nghost][j+1+nghost]-st.XLV);
			u[i+nghost][j+nghost-1]=-u[i+nghost][j+nghost]+2.0*st.XLV; //average of lowers is XLV
		}
	}
#ifdef DBGBCS
		printf("\t Setting U Upper BC \n");
#endif
	j = ny-1;
	#pragma omp parallel for
	for(i = 0; i<nx; i++){ //Upper Wall
		if(st.YBC == PERIODIC){
			u[i+nghost][j+nghost]=u[i+nghost][nghost];  //upper out of bounds is lower in bounds
		} else if (st.YBC == SLIPWALL){
			u[i+nghost][j+nghost]=u[i+nghost][j+nghost-1];  //upper out of bounds is upper in bounds
		} else {
			//u[i+nghost][j+nghost]=u[i+nghost][j-1+nghost]-2*(u[i+nghost][j-1+nghost]-st.XUV);
			u[i+nghost][j+nghost]=-u[i+nghost][j+nghost-1]+2.0*st.XUV; //average of uppers is XUV
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
			u[i+nghost][j+nghost]=(u[1+nghost][j+nghost]+u[nx-1+nghost][j+nghost])/2;
			u[i+nghost][j+nghost]=0.0;
		} else if (st.XBC == SLIPWALL){
			u[i+nghost][j+nghost]=0.0;
		} else {
			u[i+nghost][j+nghost]=0.0;
		}
	}
#ifdef DBGBCS
		printf("\t Setting U Right BC \n");
#endif
	i = nx-1;
	#pragma omp parallel for
	for(j = 0; j<ny; j++){ //Right Wall
		if(st.XBC == PERIODIC){
			u[i+nghost][j+nghost]=(u[1+nghost][j+nghost]+u[nx-1+nghost][j+nghost])/2;
			u[i+nghost][j+nghost]=0.0;
		} else if (st.XBC == SLIPWALL){
			u[i+nghost][j+nghost]=0.0;
		} else {
			u[i+nghost][j+nghost]=0.0;
		}
	}

#ifdef DBGBCS
		printf("\t Setting V Lower BC \n");
#endif
	//V Velocity
	i = 0;
	j = 0;
	#pragma omp parallel for
	for(i = -1; i<nx+1; i++){ //Lower Wall
		if(st.YBC == PERIODIC){
			v[i+nghost][j+nghost]=(v[i+nghost][1+nghost]+v[i+nghost][ny+nghost-2])/2;
		} else if (st.YBC == SLIPWALL){
			v[i+nghost][j+nghost]=0.0;
		} else {
			v[i+nghost][j+nghost]=0.0;
		}
	}
	j = ny-1;
#ifdef DBGBCS
		printf("\t Setting V Upper BC \n");
#endif
	#pragma omp parallel for
	for(i = -1; i<nx+1; i++){ //Upper Wall
		if(st.YBC == PERIODIC){
			v[i+nghost][j+nghost]=(v[i+nghost][1+nghost]+v[i+nghost][ny+nghost-2])/2;
		} else if (st.YBC == SLIPWALL){
			v[i+nghost][j+nghost]=0.0;
		} else {
			v[i+nghost][j+nghost]=0.0;
		}
	}
#ifdef DBGBCS
		printf("\t Setting V Left BC \n");
#endif
	i = 0;
	#pragma omp parallel for
	for(j = 0; j<ny; j++){ //Left Wall
		if(st.XBC == PERIODIC){
			v[i+nghost-1][j+nghost] = v[i+nghost][j+nghost]; //left oob is left ib %%%%%%%%%%%%%%%%%%%%
			//v[i+nghost-1][j+nghost] = v[i+nghost+nx-2][j+nghost]; //left oob is right ib %%%%%%%%%%%%%%%%%%%%%%%%%%%%5
			//v[i+nghost][j+nghost] = v[i+nghost+nx-2][j+nghost];
			for (int ii = 0; ii < nx/2; ii++){
				v[i+nghost+ii][j+nghost] = v[i+nghost+nx-2-ii][j+nghost];

			}
			for (int ii = 0; ii<ny/2; ii++){
				u[i+nghost+ii][j+nghost] = -u[i+nghost+nx-2-ii][j+nghost];
			}
			//v[i+nghost+1][j+nghost] = v[i+nghost+nx-3][j+nghost];
			//v[i+nghost+2][j+nghost] = v[i+nghost+nx-4][j+nghost];
		} else if (st.XBC == SLIPWALL){
			v[i+nghost-1][j+nghost] = v[i+nghost][j+nghost]; //left oob is left ib
		} else {
			//v[i+nghost][j+nghost]=v[i+1+nghost][j+nghost]-2*(v[i+1+nghost][j+nghost]-st.YLV);
			v[i+nghost-1][j+nghost]=-v[i+nghost][j+nghost]+2.0*st.YLV; //left averages

		}
	}
#ifdef DBGBCS
		printf("\t Setting V Right BC \n");
#endif
	i = nx-1;
	#pragma omp parallel for
	for(j = 0; j<ny; j++){ //Right Wall
		if(st.XBC == PERIODIC){
			v[i+nghost][j+nghost] = v[nghost][j+nghost]; //right oob is left ib
		} else if (st.XBC == SLIPWALL){
			v[i+nghost][j+nghost] = v[i+nghost-1][j+nghost]; //right oob is right ib
		} else {
			//v[i+nghost][j+nghost]=v[i-1+nghost][j+nghost]-2*(v[i-1+nghost][j+nghost]-st.YRV);
			v[i+nghost-1][j+nghost]=-v[i+nghost-1][j+nghost]+2.0*st.YRV; //right averages
			//v[i+nghost][j+nghost]=v[0+nghost-1][j+nghost];
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
void set_all_BCS(double** restrict x, double dx, double dy, int  nx, int ny, int nghost,struct slv_settings st){
	int i,j;
	i = 0;
	for(int j = -nghost; j<ny+nghost-1; j++){//left wall right wall
		if(st.XBC == PERIODIC){
			x[i+nghost-1][j+nghost] = x[i+nghost  ][j+nghost]; //left oob is right ib
			x[i+nghost-2][j+nghost] = x[i+nghost+1][j+nghost]; //left oob is right ib
			x[i+nghost-3][j+nghost] = x[i+nghost+2][j+nghost]; //left oob is right ib
			x[i+nghost+nx-1][j+nghost] = x[i+nghost+nx-2][j+nghost]; //right oob is left ib
			x[i+nghost+nx  ][j+nghost] = x[i+nghost+nx-3][j+nghost]; //right oob is left ib
			x[i+nghost+nx+1][j+nghost] = x[i+nghost+nx-4][j+nghost]; //right oob is left ib
		} /*else if (st.XBC == SLIPWALL) {
			x[i+nghost-1][j+nghost] = x[i+nghost+nx-2][j+nghost]; //left oob is right ib
			x[i+nghost-2][j+nghost] = x[i+nghost+nx-3][j+nghost]; //left oob is right ib
			x[i+nghost-3][j+nghost] = x[i+nghost+nx-3][j+nghost]; //left oob is right ib
			x[i+nghost+nx-1][j+nghost] = x[i+nghost  ][j+nghost]; //right oob is left ib
			x[i+nghost+nx  ][j+nghost] = x[i+nghost+1][j+nghost]; //right oob is left ib
			x[i+nghost+nx+1][j+nghost] = x[i+nghost+2][j+nghost]; //right oob is left ib
		} else {
			x[i+nghost-1][j+nghost] = x[i+nghost  ][j+nghost]; //left oob is right ib
			x[i+nghost-2][j+nghost] = x[i+nghost+1][j+nghost]; //left oob is right ib
			x[i+nghost-3][j+nghost] = x[i+nghost+2][j+nghost]; //left oob is right ib
			x[i+nghost+nx-1][j+nghost] = x[i+nghost+nx-2][j+nghost]; //right oob is left ib
			x[i+nghost+nx  ][j+nghost] = x[i+nghost+nx-3][j+nghost]; //right oob is left ib
			x[i+nghost+nx+1][j+nghost] = x[i+nghost+nx-4][j+nghost]; //right oob is left ib
		}*/
	}
	j = 0;
	for(int i =-nghost; i<nx+nghost-1; i++){
		if(st.YBC == PERIODIC){
			x[i+nghost][j+nghost-1] = x[i+nghost][j+nghost+ny-1]; //bottom oob is top ib
			x[i+nghost][j+nghost-2] = x[i+nghost][j+nghost+ny-2]; //bottom oob is top ib
			x[i+nghost][j+nghost-3] = x[i+nghost][j+nghost+ny-3]; //bottom oob is top ib
			x[i+nghost][j+nghost+ny-1] = x[i+nghost][j+nghost  ]; //top oob is bottom ib
			x[i+nghost][j+nghost+ny  ] = x[i+nghost][j+nghost+1]; //top oob is bottom ib
			x[i+nghost][j+nghost+ny+1] = x[i+nghost][j+nghost+2]; //top oob is bottom ib
		} /*else if (st.YBC == SLIPWALL) {
			x[i+nghost][j+nghost-1] = x[i+nghost][j+nghost+ny-2]; //bottom oob is top ib
			x[i+nghost][j+nghost-2] = x[i+nghost][j+nghost+ny-3]; //bottom oob is top ib
			x[i+nghost][j+nghost-3] = x[i+nghost][j+nghost+ny-4]; //bottom oob is top ib
			x[i+nghost][j+nghost+ny-1] = x[i+nghost][j+nghost  ]; //top oob is bottom ib
			x[i+nghost][j+nghost+ny  ] = x[i+nghost][j+nghost+1]; //top oob is bottom ib
			x[i+nghost][j+nghost+ny+1] = x[i+nghost][j+nghost+2]; //top oob is bottom ib
		} */else {
			x[i+nghost][j+nghost-1] = x[i+nghost][j+nghost]; //bottom oob is top ib
			x[i+nghost][j+nghost-2] = x[i+nghost][j+nghost+1]; //bottom oob is top ib
			x[i+nghost][j+nghost-3] = x[i+nghost][j+nghost+2]; //bottom oob is top ib
			x[i+nghost][j+nghost+ny-2] = x[i+nghost][j+nghost+ny-1]; //top oob is bottom ib
			x[i+nghost][j+nghost+ny-3] = x[i+nghost][j+nghost+ny  ]; //top oob is bottom ib
			x[i+nghost][j+nghost+ny-4] = x[i+nghost][j+nghost+ny+1]; //top oob is bottom ib
		}

	}

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
void slv_vbe(double** restrict a, double** restrict u, double** restrict us, double** restrict dus, double** restrict duss, double** restrict v, double** restrict vs, double** restrict dvs, double** restrict dvss, double** restrict hu, double** restrict huold,double** restrict hv, double** restrict hvold, double dx, double dy, int nx, int ny, int nghost, double dt, struct slv_settings st)
{
	double dx2 = pow(dx,2.0);
	double dy2 = pow(dy,2.0);
	double ce = dt/(2.0);
#ifdef DBGVBE
	printf("\n slv_vbe call \n");
#endif
	//#pragma omp parallel
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
#ifdef DBGVBE
	printf("\n slv_vbe malloc 1 done \n");
#endif
		
		#pragma omp parallel for
		for(int i=0;i<nx;i++){
			for(int j=0;j<ny;j++){

				/* Do this for all u */					
				if(i!=0){
					huold[i][j] = hu[i][j];
					hu[i][j] = -((u[i+1+nghost][j+nghost]+u[i+nghost][j+nghost])*(u[i+1+nghost][j+nghost]+u[i+nghost][j+nghost])-
							     (u[i-1+nghost][j+nghost]+u[i+nghost][j+nghost])*(u[i-1+nghost][j+nghost]+u[i+nghost][j+nghost]))/(4.0*dx)-
						        ((u[i+nghost][j+1+nghost]+u[i+nghost][j+nghost])*(v[i+nghost][j+nghost]+v[i+1+nghost][j+nghost])-
						         (u[i+nghost][j-1+nghost]+u[i+nghost][j+nghost])*(v[i+nghost][j-1+nghost]+v[i+1+nghost][j-1+nghost]))/(4.0*dy) + st.XConst ;
				}		
				/* Do this for all v */
				if(j!=0){
					hvold[i][j] = hv[i][j];
					hv[i][j] = -((v[i+nghost][j+1+nghost]+v[i+nghost][j+nghost])*(v[i+nghost][j+1+nghost]+v[i+nghost][j+nghost])-
							     (v[i+nghost][j-1+nghost]+v[i+nghost][j+nghost])*(v[i+nghost][j-1+nghost]+v[i+nghost][j+nghost]))/(4.0*dy)-
							    ((v[i+1+nghost][j+nghost]+v[i+nghost][j+nghost])*(u[i+nghost][j+nghost]+u[i+nghost][j+1+nghost])-
							     (v[i-1+nghost][j+nghost]+v[i+nghost][j+nghost])*(u[i-1+nghost][j+nghost]+u[i-1+nghost][j+1+nghost]))/(4.0*dx) + st.YConst +st.g;
				}
			}

		}
#ifdef DBGVBE
	printf("\n slv_vbe h 1 done \n");
#endif

		
		//#pragma omp parallel for
		for(int j=0;j<ny-1;j++){
			for(int i = 0; i<nx-1;i++){
				double rho = a[i+nghost][j+nghost]*st.rhol+(1-a[i+nghost][j+nghost])*st.rhog;
				double mu = a[i+nghost][j+nghost]*st.mul+(1-a[i+nghost][j+nghost])*st.mug;
				/* limit range for u */
				if (i == 0){
					bum[i] = dt/2.0*(3.0*hu[i][j]-huold[i][j])+2.0*ce*mu/rho*((u[i+1+nghost][j+nghost]-2.0*u[i+nghost][j+nghost]+u[i-1+nghost][j+nghost])/dx2
							+(u[i+nghost][j+1+nghost]-2.0*u[i+nghost][j+nghost]+u[i+nghost][j-1+nghost])/dy2);
					lum[i] = -ce/dx2*mu/rho;
					mum[i] = (1+2*ce*mu/rho/dx2);
					
					if(j!=ny){
						bvm[i] = dt/2.0*(3.0*hv[i][j]-hvold[i][j])+2.0*ce*mu/rho*((v[i+1+nghost][j+nghost]-2.0*v[i+nghost][j+nghost]+v[i-1+nghost][j+nghost])/dx2
								+(v[i+nghost][j+1+nghost]-2.0*v[i+nghost][j+nghost]+v[i+nghost][j-1+nghost])/dy2);
						lvm[i] = -ce/dx2*mu/rho;
						mvm[i] = (1.0+2.0*ce/dx2*mu/rho);
					}
				} else if (i == nx-3){
					bum[i] = dt/1.0*(2.0*hu[i][j]-huold[i][j])+2.0*ce*mu/rho*((u[i+1+nghost][j+nghost]-2.0*u[i+nghost][j+nghost]+u[i-1+nghost][j+nghost])/dx2
							+(u[i+nghost][j+1+nghost]-2.0*u[i+nghost][j+nghost]+u[i+nghost][j-1+nghost])/dy2);
					mum[i] = (1+2*ce/dx2*mu/rho);
					uum[i] = -ce/dx2*mu/rho;
					
					if(j!=0){
						bvm[i] = dt/2.0*(3.0*hv[i][j]-hvold[i][j])+2.0*ce*mu/rho*((v[i+1+nghost][j+nghost]-2.0*v[i+nghost][j+nghost]+v[i-1+nghost][j+nghost])/dx2+
								(v[i+nghost][j+1+nghost]-2.0*v[i+nghost][j+nghost]+v[i+nghost][j-1+nghost])/dy2);
						lvm[i] = -ce/dx2*mu/rho;
						mvm[i] = (1.0+2.0*ce/dx2*mu/rho);
						uvm[i] = -ce/dx2*mu/rho;
					}
				} else if (i == nx-2){
					if(j!=0){
						bvm[i] = dt/2.0*(3.0*hv[i][j]-hvold[i][j])+2.0*ce*mu/rho*((v[i+1+nghost][j+nghost]-2.0*v[i+nghost][j+nghost]+v[i-1+nghost][j+nghost])/dx2+
								(v[i+nghost][j+1+nghost]-2.0*v[i+nghost][j+nghost]+v[i+nghost][j-1+nghost])/dy2);
						mvm[i] = (1.0+2.0*ce/dx2*mu/rho);
						uvm[i] = -ce/dx2*mu/rho;
					}
				} else {
					
					bum[i] = dt/2.0*(3.0*hu[i][j]-huold[i][j])+2.0*ce*mu/rho*((u[i+1+nghost][j+nghost]-2.0*u[i+nghost][j+nghost]+u[i-1+nghost][j+nghost])/dx2+
							(u[i+nghost][j+1+nghost]-2.0*u[i+nghost][j+nghost]+u[i+nghost][j-1+nghost])/dy2);
					lum[i] = -ce/dx2*mu/rho;
					mum[i] = (1.0+2.0*ce/dx2*mu/rho);
					uum[i] = -ce/dx2*mu/rho;
	
					if(j!=0){
						bvm[i] = dt/2.0*(3.0*hv[i][j]-hvold[i][j])+2.0*ce*mu/rho*((v[i+1+nghost][j+nghost]-2.0*v[i+nghost][j+nghost]+v[i-1][j])/dx2+
								(v[i+nghost][j+1+nghost]-2.0*v[i+nghost][j+nghost]+v[i+nghost][j-1+nghost])/dy2);
						lvm[i] = -ce/dx2*mu/rho;
						mvm[i] = (1.0+2.0*ce/dx2*mu/rho);
						uvm[i] = -ce/dx2*mu/rho;
					}
				}
			}
/*#ifdef DBGVBE
	printf("\n slv_vbe post setup done \n");
#endif*/
			
			solve_matrix(nx-2,lum,mum,uum,bum,xum);
			if(j!=ny-1){
				solve_matrix(nx-1,lvm,mvm,uvm,bvm,xvm);
			}			
/*#ifdef DBGVBE
	printf("\n slv_vbe post solve 1\n");
#endif*/
			for(int i = 0;i<nx-1;i++){
				if(i!=0){
					duss[i][j] = xum[i];
				}
				if(j!=0){
					dvss[i][j] = xvm[i];
				}
			}	
		}/*
		free(bum);
		free(lum);
		free(mum);
		free(uum);
		free(xum);
		free(bvm);
		free(lvm);
		free(mvm);
		free(uvm);
		free(xvm);*/
	}
//#pragma omp end parallel
#pragma omp barrier
	#pragma omp parallel for
	for(int i = 0; i<nx-1; i++){
		for(int j=0; j<ny-1; j++){
			if(j!=0){
				vs[i+nghost][j+nghost] = v[i+nghost][j+nghost]+dvss[i][j];
			}
			if(i!=0){
				us[i+nghost][j+nghost] = u[i+nghost][j+nghost]+duss[i][j];
			}
		}
	}
#pragma omp barrier
#ifdef DBGVBE
	printf("\n slv_vbe post set solution 1\n");
#endif
#ifdef DBGVBE
	printf("\n slv_vbe set intermediate BCS \n");
#endif
	set_bcs(us, vs, dx, dy, nx, ny,nghost, st);
#ifdef DBGVBE
	printf("\n slv_vbe post set intermediate BCS \n");
#endif
	//#pragma omp parallel
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
#ifdef DBGVBE
	printf("\n slv_vbe post malloc 2 \n");
#endif
		//#pragma omp parallel for
		for(int i = 0; i<nx-1; i++){
			for(int j = 0; j<ny-1; j++){
				double rho = a[i+nghost][j+nghost]*st.rhol+(1.0-a[i+nghost][j+nghost])*st.rhog;
				double mu = a[i+nghost][j+nghost]*st.mul+(1.0-a[i+nghost][j+nghost])*st.mug;
				//printf("%i,%i ",i,j);
				/* limit range for v */
				if (j == 0){
					if(i!=0){
						bumy[j] = duss[i][j];
						bumy[j] = dt/2.0*(3.0*hu[i][j]-huold[i][j])+2.0*ce*mu/rho*((us[i+1+nghost][j+nghost]-2.0*us[i+nghost][j+nghost]+us[i-1+nghost][j+nghost])/dx2+
								(us[i+nghost][j+1+nghost]-2.0*us[i+nghost][j+nghost]+us[i+nghost][j-1+nghost])/dy2);
						lumy[j] = -ce/dy2*mu/rho;
						mumy[j] = (1.0+2.0*ce/dy2*mu/rho);
					}
					bvmy[j] = dvss[i][j];
					bvmy[j] = dt/2.0*(3.0*hv[i][j]-hvold[i][j])+2.0*ce*mu/rho*((vs[i+1+nghost][j+nghost]-2.0*vs[i+nghost][j+nghost]+vs[i-1+nghost][j+nghost])/dx2+
							(vs[i+nghost][j+1+nghost]-2.0*vs[i+nghost][j+nghost]+vs[i+nghost][j-1+nghost])/dy2);
					lvmy[j] = -ce/dy2*mu/rho;
					mvmy[j] = (1.0+2.0*ce/dy2*mu/rho);
				} else if (j == ny-3){
					if(i!=0){
						bumy[j] = duss[i][j];
						bumy[j] = dt/2.0*(3.0*hu[i][j]-huold[i][j])+2.0*ce*mu/rho*((us[i+1+nghost][j+nghost]-2.0*us[i+nghost][j+nghost]+us[i-1+nghost][j+nghost])/dx2+
								(us[i+nghost][j+1+nghost]-2.0*us[i+nghost][j+nghost]+us[i+nghost][j-1+nghost])/dy2);
						lumy[j] = -ce/dy2*mu/rho;
						mumy[j] = (1.0+2.0*ce/dy2*mu/rho);
						uumy[j] = -ce/dy2*mu/rho;
					}
					bvmy[j] = dvss[i][j];
					bvmy[j] = dt/2.0*(3.0*hv[i][j]-hvold[i][j])+2.0*ce*mu/rho*((vs[i+1+nghost][j+nghost]-2.0*vs[i+nghost][j+nghost]+vs[i-1+nghost][j+nghost])/dx2+
							(vs[i+nghost][j+1+nghost]-2.0*vs[i+nghost][j+nghost]+vs[i+nghost][j-1+nghost])/dy2);
					mvmy[j] = (1.0+2.0*ce/dy2*mu/rho);
					uvmy[j] = -ce/dy2*mu/rho;
				} else if (j == ny-2){
					if(i!=0){
						bumy[j] = duss[i][j];
						bumy[j] = dt/2.0*(3.0*hu[i][j]-huold[i][j])+2.0*ce*mu/rho*((us[i+1+nghost][j+nghost]-2.0*us[i+nghost][j+nghost]+us[i-1+nghost][j+nghost])/dx2+
								(us[i+nghost][j+1+nghost]-2.0*us[i+nghost][j+nghost]+us[i+nghost][j-1+nghost])/dy2);
						mumy[j] = (1.0+2.0*ce/dy2*mu/rho);
						uumy[j] = -ce/dy2*mu/rho;
					}
				} else {
					if(i!=0){
						bumy[j] = duss[i][j];
						bumy[j] = dt/2.0*(3.0*hu[i][j]-huold[i][j])+2.0*ce*mu/rho*((us[i+1+nghost][j+nghost]-2.0*us[i+nghost][j+nghost]+us[i-1+nghost][j+nghost])/dx2+
								(us[i+nghost][j+1+nghost]-2.0*us[i+nghost][j+nghost]+us[i+nghost][j-1+nghost])/dy2);
						lumy[j] = -ce/dy2*mu/rho;
						mumy[j] = (1.0+2.0*ce/dy2*mu/rho);
						uumy[j] = -ce/dy2*mu/rho;
					}
					bvmy[j] = dvss[i][j];
					bvmy[j] = dt/2.0*(3.0*hv[i][j]-hvold[i][j])+2.0*ce*mu/rho*((vs[i+1+nghost][j+nghost]-2.0*vs[i+nghost][j+nghost]+vs[i-1+nghost][j+nghost])/dx2+
							(vs[i+nghost][j+1+nghost]-2.0*vs[i+nghost][j+nghost]+vs[i+nghost][j-1+nghost])/dy2);
					lvmy[j] = -ce/dy2*mu/rho;
					mvmy[j] = (1.0+2.0*ce/dy2*mu/rho);
					uvmy[j] = -ce/dy2*mu/rho;
				}
			}
/*#ifdef DBGVBE
	printf("\n slv_vbe post loop 2 \n");
#endif*/
			if(i!=0){
				solve_matrix(ny-1,lumy,mumy,uumy,bumy,xumy);
			}
			solve_matrix(ny-2,lvmy,mvmy,uvmy,bvmy,xvmy);		
			for(int j = 0;j<ny-2;j++){
				if(j!=0){
					dvs[i][j] = xvmy[j];
				}
				if(i!=0){
					dus[i][j] = xumy[j];
				}
			}	
		}
#ifdef DBGVBE
	printf("\n slv_vbe post set solution \n");
#endif
		/*free(bumy);
		free(lumy);
		free(mumy);
		free(uumy);
		free(xumy);
		free(bvmy);
		free(lvmy);
		free(mvmy);
		free(uvmy);
		free(xvmy);*/
#ifdef DBGVBE
	printf("\n slv_vbe post free 2 \n");
#endif
	}
#pragma omp barrier
	#pragma omp parallel for
	for(int i = 0; i<nx-1; i++){
		for(int j=0; j<ny-1; j++){
			if(j!=0){
				vs[i+nghost][j+nghost] = vs[i+nghost][j+nghost]+dvs[i][j];
			}
			if(i!=0){
				us[i+nghost][j+nghost] = us[i+nghost][j+nghost]+dus[i][j];
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
	set_bcs(us, vs, dx, dy, nx, ny,nghost, st);
}

/* ftcs solver */
void slv_vbeftcs(double** restrict a, double** restrict u, double** restrict us, double** restrict dus, double** restrict duss, double** restrict v, double** restrict vs, double** restrict dvs, double** restrict dvss, double** restrict hu, double** restrict huold,double** restrict hv, double** restrict hvold, double dx, double dy, int nx, int ny, int nghost, double dt, struct slv_settings st)
{
	
	double dx2 = pow(dx,2.0);
	double dy2 = pow(dy,2.0);
	double ce = dt;
	int xp = 0;
	int yp = 0;
	if (st.XBC == PERIODIC){
		xp = 1;
	}
	if (st.YBC == PERIODIC){
		yp = 1;
	}

	#pragma omp parallel for
	for(int i=0;i<nx-1;i++){
		for(int j=0;j<ny-1;j++){
			/* Do this for all u */					
			if(i!=0){
				huold[i][j] = hu[i][j];
				hu[i][j] = -((u[i+1+nghost][j+nghost]+u[i+nghost][j+nghost])*(u[i+1+nghost][j+nghost]+u[i+nghost][j+nghost])
				        -(u[i-1+nghost][j+nghost]+u[i+nghost][j+nghost])*(u[i-1+nghost][j+nghost]+u[i+nghost][j+nghost]))/(4.0*dx)
						-((u[i+nghost][j+1+nghost]+u[i+nghost][j+nghost])*(v[i+nghost][j+nghost]+v[i+1+nghost][j+nghost])
								-(u[i+nghost][j-1+nghost]+u[i+nghost][j+nghost])*(v[i+nghost][j-1+nghost]+v[i+1+nghost][j-1+nghost]))/(4.0*dy) + st.XConst;
			}		
			/* Do this for all v */
			if(j!=0){
				hvold[i][j] = hv[i][j];
				hv[i][j] = -((v[i+nghost][j+1+nghost]+v[i+nghost][j+nghost])*(v[i+nghost][j+1+nghost]+v[i+nghost][j+nghost])
						-(v[i+nghost][j-1+nghost]+v[i+nghost][j+nghost])*(v[i+nghost][j-1+nghost]+v[i+nghost][j+nghost]))/(4.0*dy)
						-((v[i+1+nghost][j+nghost]+v[i+nghost][j+nghost])*(u[i+nghost][j+nghost]+u[i+nghost][j+1+nghost])
								-(v[i-1+nghost][j+nghost]+v[i+nghost][j+nghost])*(u[i-1+nghost][j+nghost]+u[i-1+nghost][j+1+nghost]))/(4.0*dx) + st.YConst + st.g ;
			}
		}
	}
	#pragma omp parallel for
	for(int i=0;i<nx-1;i++){
		for(int j = 0; j<ny-1;j++){
			double rpj, rmj, rpi, rmi, rmmi, rho;
			rho = a[i+nghost][j+nghost]*st.rhol+(1.0-a[i+nghost][j+nghost])*st.rhog;
			rpj = a[i+nghost][j+nghost+1]*st.rhol+(1.0-a[i+nghost][j+nghost+1])*st.rhog;
			rmj = a[i+nghost][j+nghost-1]*st.rhol+(1.0-a[i+nghost][j+nghost-1])*st.rhog;
			rpi = a[i+nghost+1][j+nghost]*st.rhol+(1.0-a[i+nghost+1][j+nghost])*st.rhog;
			rmi = a[i+nghost-1][j+nghost]*st.rhol+(1.0-a[i+nghost-1][j+nghost])*st.rhog;
			double ru = (rho+rmi)/2.0;
			double rv = (rho+rmj)/2.0;
			double mpj, mmj, mpi,mmimj,mpimj, mmi,mmipj, mhu;
			mhu = a[i+nghost][j+nghost]*st.mul+(1.0-a[i+nghost][j+nghost])*st.mug;
			mpj = a[i+nghost][j+nghost+1]*st.mul+(1.0-a[i+nghost][j+nghost+1])*st.mug;
			mmimj = a[i+nghost-1][j+nghost-1]*st.mul+(1.0-a[i+nghost-1][j+nghost-1])*st.mug;
			mmipj = a[i+nghost-1][j+nghost+1]*st.mul+(1.0-a[i+nghost-1][j+nghost+1])*st.mug;
			mmj = a[i+nghost][j+nghost-1]*st.mul+(1.0-a[i+nghost][j+nghost-1])*st.mug;
			mpimj = a[i+nghost+1][j+nghost-1]*st.mul+(1.0-a[i+nghost+1][j+nghost-1])*st.mug;
			mpi = a[i+nghost+1][j+nghost]*st.mul+(1.0-a[i+nghost+1][j+nghost])*st.mug;
			mmi = a[i+nghost-1][j+nghost]*st.mul+(1.0-a[i+nghost-1][j+nghost])*st.mug;
			if(i!=0){
				us[i+nghost][j+nghost] = u[i+nghost][j+nghost]+dt*(1.5*hu[i][j]-0.5*huold[i][j])
						+ce/ru*(1.0*(mhu)*(u[i+1+nghost][j+nghost]-u[i+nghost][j+nghost])/dx2 -                   1.0*(mmi)*(u[i+nghost][j+nghost]-u[i-1+nghost][j+nghost])/dx2
                 +(mhu+mpj+mmipj+mmi)*.25*(u[i+nghost][j+1+nghost]-u[i+nghost][j+nghost])/dy2 -     (mhu+mmj+mmimj+mmi)*.25*(u[i+nghost][j+nghost]-u[i+nghost][j-1+nghost])/dy2);
						           //+(mpi+mmipj+mhu+mmi)/4.0*(v[i+nghost+1][j+nghost]-v[i+nghost][j+nghost])/(dx*dy)
						           //-(mmi+mmimj+mhu+mmi)/4.0*(v[i+nghost][j+nghost+1]-v[i+nghost-1][j+nghost+1])/(dy*dx) );
			}
			if(j!=0 ){
				vs[i+nghost][j+nghost] = v[i+nghost][j+nghost]+dt*(1.5*hv[i][j]-0.5*hvold[i][j])
						+ce/rv*((mhu+mpi+mmj+mpimj)*.25*(v[i+1+nghost][j+nghost]-v[i+nghost][j+nghost])/dx2 - (mhu+mmi+mmj+mmimj)*.25*(v[i+nghost][j+nghost]-v[i-1+nghost][j+nghost])/dx2
						                     +1.0*(mhu)*(v[i+nghost][j+1+nghost]-v[i+nghost][j+nghost])/dy2 -           1.0*(mmj)*(v[i+nghost][j+nghost]-v[i+nghost][j-1+nghost])/dy2);
						           //+(mhu+mpi+mpimj+mmj)/4.0*(u[i+nghost+1][j+nghost]-u[i+nghost+1][j+nghost-1])/(dx*dy)
							       //-(mhu+mmi+mmimj+mmj)/4.0*(u[i+nghost][j+nghost]-u[i+nghost][j+nghost-1] )/(dy*dx) );
			}
		}	
	}
	/*
	#pragma omp parallel for
	for(int i = 0;i<nx+1;i++){//lower and upper
		vs[i+nghost][0+nghost] = v[i+nghost][0+nghost];
		vs[i+nghost][ny-1+nghost] = v[i+nghost][ny-1+nghost];
		if(i!=nx){
			us[i+nghost][nghost-1] = u[i+nghost][nghost-1];
			us[i+nghost][ny+nghost-1] = u[i+nghost][ny+nghost-1];
		}
	}
	#pragma omp parallel for
	for(int j=0;j<ny+1;j++){//left and right
		us[0+nghost][j+nghost] = u[0+nghost][j+nghost];
		us[nx-1+nghost][j+nghost] = u[nx-1+nghost][j+nghost];
		if(j!=ny){
			vs[0+nghost-1][j+nghost] = v[0+nghost-1][j+nghost];
			vs[nx+nghost-1][j+nghost] = v[nx+nghost-1][j+nghost];
		}
	}*/
	set_bcs(us, vs, dx, dy, nx, ny,nghost, st);
}


/* Poisson Solver */
void slv_pssn(double** restrict phi,double** restrict phinext, double** restrict a, double** restrict us, double** restrict vs, double dx, double dy, int nx, int ny,int nghost, double dt,double min, struct slv_settings st)
{
	double dx2 = pow(dx,2.0);
	double dy2 = pow(dy,2.0);
	double save = 0.00;
	double resid = 100.0;
	int cnt = 0;
	while( resid > min && cnt < nx*ny*10){
		cnt = cnt+1;
		//printf("count %i \n",cnt);
		#pragma omp parallel for
		for(int i = 0;i<nx;i++){
			phi[i+nghost][nghost-1] = phi[i+nghost][nghost]; //lower wall
			phi[i+nghost][ny+nghost-1] = phi[i+nghost][ny-2+nghost]; //upper wall
		}
		#pragma omp parallel for 
		for(int j = 0;j<ny;j++){
			phi[nghost-1][j+nghost] = phi[nghost][j+nghost];
			phi[nx+nghost-1][j+nghost] = phi[nx-2+nghost][j+nghost];
		}
		phinext[nghost][nghost-1] = 0.0;
		phi[nghost][nghost-1] = 0.0;

		//set_all_bcs_neumann(phi,dx,dy,nx,ny,nghost,nghost);
		#pragma omp parallel for
		for(int i = 0; i<nx-1;i++){
			for(int j = 0; j<ny-1; j++){
				//if(i == 0 && j == 0){ j = j+1;}
				//else {
					double rpj, rmj, rpi, rmi, rho;
					rho = a[i+nghost][j+nghost]*st.rhol+(1.0-a[i+nghost][j+nghost])*st.rhog;
					rpj = a[i+nghost][j+nghost+1]*st.rhol+(1.0-a[i+nghost][j+nghost+1])*st.rhog;
					rmj = a[i+nghost][j+nghost-1]*st.rhol+(1.0-a[i+nghost][j+nghost-1])*st.rhog;
					rpi = a[i+nghost+1][j+nghost]*st.rhol+(1.0-a[i+nghost+1][j+nghost])*st.rhog;
					rmi = a[i+nghost-1][j+nghost]*st.rhol+(1.0-a[i+nghost-1][j+nghost])*st.rhog;
					double rpj2,rmj2,rpi2,rmi2;
					rpj2 = (rpj+rho)/2.0;
					rmj2 = (rmj+rho)/2.0;
					rpi2 = (rpi+rho)/2.0;
					rmi2 = (rmi+rho)/2.0;
					//phinext[i][j] = 1.0/rho*1.0/(2.0/dx2+2.0/dy2)*
					//		((rpi*phi[i+1][j]+rmi*phi[i-1][j])/dx2+(rpj*phi[i][j+1]+rmj*phi[i][j-1])/dy2)
					//		-1.0/(dt*(2.0/dx2+2.0/dy2))*((rho*us[i][j]-rmi*us[i-1][j])/dx+(rho*vs[i][j]-rmj*vs[i][j-1])/dy);
					phinext[i+nghost][j+nghost] = 1.0/((-1.0/(rpj2*dy2)-1.0/(rmj2*dy2)-1.0/(rpi2*dx2)-1.0/(rmi2*dx2)))*
							( (us[i+nghost+1][j+nghost]-us[i+nghost][j+nghost])/(dt*dx) + (vs[i+nghost][j+nghost+1]-vs[i+nghost][j+nghost])/(dt*dy)
									-(1.0/(rpi2*dx2)*phi[i+1+nghost][j+nghost] + 1.0/(rmi2*dx2)*phinext[i-1+nghost][j+nghost]
									 +1.0/(rpj2*dy2)*phi[i+nghost][j+1+nghost] + 1.0/(rmj2*dy2)*phinext[i+nghost][j-1+nghost]) );
				//}
			}
		}

		#pragma omp parallel for
		for(int i = 0;i<nx;i++){
			phinext[i+nghost][nghost-1] = phinext[i+nghost][nghost];
			phinext[i+nghost][ny+nghost-1] = phinext[i+nghost][ny-2+nghost];
		}
		#pragma omp parallel for 
		for(int j = 0;j<ny;j++){
			phinext[nghost-1][j+nghost] = phinext[nghost][j+nghost];
			phinext[nx+nghost-1][j+nghost] = phinext[nx-2+nghost][j+nghost];
		}
		phinext[nghost-1][nghost] = 0.0;
		phi[nghost-1][nghost] = 0.0;

		//set_all_bcs_neumann(phinext,dx,dy,nx-1,ny-1,nghost,nghost);
		#pragma omp parallel for //Solve system
		for(int i = 0; i<nx-1;i++){
					for(int j = 0; j<ny-1; j++){
					//if(i == 0 && j == 0){ j = j+1;}
					//else {
						double rpj, rmj, rpi, rmi, rho;
						rho = a[i+nghost][j+nghost]*st.rhol+(1.0-a[i+nghost][j+nghost])*st.rhog;
						rpj = a[i+nghost][j+nghost+1]*st.rhol+(1.0-a[i+nghost][j+nghost+1])*st.rhog;
						rmj = a[i+nghost][j+nghost-1]*st.rhol+(1.0-a[i+nghost][j+nghost-1])*st.rhog;
						rpi = a[i+nghost+1][j+nghost]*st.rhol+(1.0-a[i+nghost+1][j+nghost])*st.rhog;
						rmi = a[i+nghost-1][j+nghost]*st.rhol+(1.0-a[i+nghost-1][j+nghost])*st.rhog;
						double rpj2,rmj2,rpi2,rmi2;
						rpj2 = (rpj+rho)/2.0;
						rmj2 = (rmj+rho)/2.0;
						rpi2 = (rpi+rho)/2.0;
						rmi2 = (rmi+rho)/2.0;
						/*if(i == 40-nghost && j == 32-nghost){
							double e1 = a[i+nghost][j+nghost-1];
							double e2 = st.rhol;
							double e3 = (1.0-a[i+nghost][j+nghost]);
							double e4 = st.rhog;
							double dsxd = 2.0;

						}*/
						//phi[i][j] = 1.0/(rho)*1.0/(2.0/dx2+2.0/dy2)*((rpi*phinext[i+1][j]+rmi*phinext[i-1][j])/dx2+(rpj*phinext[i][j+1]+rmj*phinext[i][j-1])/dy2)
						//		-1.0/(dt*(2.0/dx2+2.0/dy2))*((rho*us[i][j]-rmi*us[i-1][j])/dx+(rho*vs[i][j]-rmj*vs[i][j-1])/dy);
						phi[i+nghost][j+nghost] = 1.0/((-1.0/(rpj2*dy2)-1.0/(rmj2*dy2)-1.0/(rpi2*dx2)-1.0/(rmi2*dx2)))*
												( (us[i+nghost+1][j+nghost]-us[i+nghost][j+nghost])/(dt*dx) + (vs[i+nghost][j+nghost+1]-vs[i+nghost][j+nghost])/(dt*dy)
														-(1.0/(rpi2*dx2)*phinext[i+1+nghost][j+nghost] + 1.0/(rmi2*dx2)*phi[i-1+nghost][j+nghost]
														 +1.0/(rpj2*dy2)*phinext[i+nghost][j+1+nghost] + 1.0/(rmj2*dy2)*phi[i+nghost][j-1+nghost]) );
						//printf("%i,%i\n",i,j);
				//}
			}
		}

		resid = 0.0;
		for(int i=0; i<nx-1; i++){
			for(int j=0;j<ny-1;j++){
				double rpj, rmj, rpi, rmi, rho;
				rho = a[i+nghost][j+nghost]*st.rhol+(1.0-a[i+nghost][j+nghost])*st.rhog;
				rpj = a[i+nghost][j+nghost+1]*st.rhol+(1.0-a[i+nghost][j+nghost+1])*st.rhog;
				rmj = a[i+nghost][j+nghost-1]*st.rhol+(1.0-a[i+nghost][j+nghost-1])*st.rhog;
				rpi = a[i+nghost+1][j+nghost]*st.rhol+(1.0-a[i+nghost+1][j+nghost])*st.rhog;
				rmi = a[i+nghost-1][j+nghost]*st.rhol+(1.0-a[i+nghost-1][j+nghost])*st.rhog;
				double rpj2,rmj2,rpi2,rmi2;
				rpj2 = (rpj+rho)/2.0;
				rmj2 = (rmj+rho)/2.0;
				rpi2 = (rpi+rho)/2.0;
				rmi2 = (rmi+rho)/2.0;
				save = fabs(phi[i+nghost][j+nghost]-1.0/((-1.0/(rpj2*dy2)-1.0/(rmj2*dy2)-1.0/(rpi2*dx2)-1.0/(rmi2*dx2)))*
						( (us[i+nghost+1][j+nghost]-us[i+nghost][j+nghost])/(dt*dx) + (vs[i+nghost][j+nghost+1]-vs[i+nghost][j+nghost])/(dt*dy)
								-(1.0/(rpi2*dx2)*phi[i+1+nghost][j+nghost] + 1.0/(rmi2*dx2)*phi[i-1+nghost][j+nghost]
								 +1.0/(rpj2*dy2)*phi[i+nghost][j+1+nghost] + 1.0/(rmj2*dy2)*phi[i+nghost][j-1+nghost]) )); //Infinity Norm
				if (save>resid){
					resid = save;
				} 	
			}
		}
	}
	printf("count %i \n",cnt);
#ifdef DEBUG
	printf("\t\t Poisson Iterations: %i\n",cnt);
#endif
}
int get_pssn_ind(int i, int j, int nx, int ny){
	int ind;
	ind = j*(nx-1)+i;
	return ind;
}
int get_pssn_i(int ind, int nx, int ny){
	int i = ind%(nx-1);
	return i;
}
int get_pssn_j(int ind, int nx, int ny){
	int i = get_pssn_i(ind, nx, ny);
	int j = (ind-i)/(nx-1);
	return j;
}
/* Poisson Solver */
void slv_pssn_gmres(double** restrict phi,double* restrict x, double** restrict a, double** restrict us, double** restrict vs, double dx, double dy, int nx, int ny, int nghost, double dt,double min, struct slv_settings st)
{
#ifdef DEBUG
	printf("\n Mgmres poisson call \n");
#endif
	double dx2 = pow(dx,2.0);
	double dy2 = pow(dy,2.0);
	int ic = 0;
	int jc = 0;
	int n = (nx-1)*(ny-1); //Rank of Linear System
	int nz_num = (5*n)-2*(ny-1)-2*(nx-1); //the number of nonzero matrix values. 5 per entry minus out of bounds left/right top/bottom

	double* A = malloc(nz_num*sizeof(double)); //Matrix Values
	int* iA = malloc(nz_num*sizeof(int)); //Matrix Row Index
	int* iAcr = malloc((n+1)*sizeof(int)); // Matrix Row Index - Compressed Row;
	int* jA = malloc(nz_num*sizeof(int)); //Matrix Column Index

	//double* x = malloc(n*sizeof(double)); //An approximation to the solution.  On output, an improved approximation.
    double* rhs = malloc(n*sizeof(double)); //The right hand side of the linear system.
    int itr_max = 10; //The maximum number of (outer) iterations to take.

    int mr = n/5; //the maximum number of (inner) iterations to take.
    double tol_abs = 0.0000000000001; //An absolute tolerance applied to the current residual.

    double tol_rel = 1.0; //A relative tolerance comparing the current residual to the initial residual.
#ifdef DBGMGM
	printf("\n Finished MGMRES Malloc");
#endif
    int count = 0;
    for (int j = 0; j<ny-1; j++){
    	for(int i = 0; i<nx-1; i++){
    		double rpj, rmj, rpi, rmi, rho;
			rho = a[i+nghost][j+nghost]*st.rhol+(1.0-a[i+nghost][j+nghost])*st.rhog;
			rpj = a[i+nghost][j+nghost+1]*st.rhol+(1.0-a[i+nghost][j+nghost+1])*st.rhog;
			rmj = a[i+nghost][j+nghost-1]*st.rhol+(1.0-a[i+nghost][j+nghost-1])*st.rhog;
			rpi = a[i+nghost+1][j+nghost]*st.rhol+(1.0-a[i+nghost+1][j+nghost])*st.rhog;
			rmi = a[i+nghost-1][j+nghost]*st.rhol+(1.0-a[i+nghost-1][j+nghost])*st.rhog;
			double rpj2,rmj2,rpi2,rmi2;
			rpj2 = (rpj+rho)/2.0;
			rmj2 = (rmj+rho)/2.0;
			rpi2 = (rpi+rho)/2.0;
			rmi2 = (rmi+rho)/2.0;
			int indpi, indmi, indpj, indmj, ind;
			int cpi, cmi, cpj, cmj,ch;
			ind = get_pssn_ind(i,j,nx,ny);
			//printf("\n Count(%i,%i) =  %i, cmax = %i",i,j,count, nz_num);
			ch = count; count++;
			if(i > 0){
				indmi = get_pssn_ind(i-1,j,nx,ny);
				cmi = count; count++;
			} else {
				indmi = get_pssn_ind(i+1,j,nx,ny);
				cmi = count; rmi2 = rpi2;
			}
			if(i < nx-2){
				indpi = get_pssn_ind(i+1,j,nx,ny);
				cpi = count;count++;
			} else {
				indpi = get_pssn_ind(i-1,j,nx,ny);
				count--;
				cpi = count; count++; rpi2 = rmi2;
			}
			if(j > 0){
				indmj = get_pssn_ind(i,j-1,nx,ny);
				cmj = count; count++;
			} else {
				indmj = get_pssn_ind(i,j+1,nx,ny);
				cmj = count; rmj2 = rpj2;
			}
			if(j < ny-2){
				indpj = get_pssn_ind(i,j+1,nx,ny);
				cpj = count; count++;
			} else {
				indpj = get_pssn_ind(i,j-1,nx,ny);
				count--;
				cpj = count; count++; rpj2 = rmj2;
			}
			A[ch] = 0.0; A[cmi] = 0.0; A[cpi] = 0.0; A[cmj] = 0.0; A[cpj] = 0.0;
			iAcr[ind] = ch; iAcr[ind+1] = ch+3;
			if(ch == 0){
				A[ch] = 1.0; iA[ch] = ind; jA[ch] = ind;
				rhs[ind] = 0.0;
			} else {
				A[ch] = ( -1.0/(rpj2*dy2) -1.0/(rmj2*dy2) -1.0/(rpi2*dx2) -1.0/(rmi2*dx2) ); iA[ch] = ind; jA[ch] = ind;
				rhs[ind] = ( (us[i+nghost+1][j+nghost]-us[i+nghost][j+nghost]) / (dt*dx) + (vs[i+nghost][j+nghost+1]-vs[i+nghost][j+nghost])/(dt*dy) );
			}

				A[cpi] += (1.0/(rpi2*dx2));iA[cpi] = ind; jA[cpi] = indpi;
				A[cmi] += (1.0/(rmi2*dx2));iA[cmi] = ind; jA[cmi] = indmi;
				A[cpj] += (1.0/(rpj2*dy2));iA[cpj] = ind; jA[cpj] = indpj;
				A[cmj] += (1.0/(rmj2*dy2));iA[cmj] = ind; jA[cmj] = indmj;

			/*phi[i+nghost][j+nghost] = 1.0/((-1.0/(rpj2*dy2)-1.0/(rmj2*dy2)-1.0/(rpi2*dx2)-1.0/(rmi2*dx2)))*
															( (us[i+nghost+1][j+nghost]-us[i+nghost][j+nghost])/(dt*dx) + (vs[i+nghost][j+nghost+1]-vs[i+nghost][j+nghost])/(dt*dy)
																	-(1.0/(rpi2*dx2)*phinext[i+1+nghost][j+nghost] + 1.0/(rmi2*dx2)*phi[i-1+nghost][j+nghost]
																	 +1.0/(rpj2*dy2)*phinext[i+nghost][j+1+nghost] + 1.0/(rmj2*dy2)*phi[i+nghost][j-1+nghost]) );*/
    	}
    }

#ifdef DBGMGM
	printf("\n Equations from BCs - %i", count);
#endif
	//mgmres_st (int n, int nz_num, int ia[], int ja[], double a[], double x[], double rhs[], int itr_max, int mr, double tol_abs, double tol_rel )
#ifdef DBGMGM
	printf("\n Calling mgmres \n");

#endif

	//mgmres_st (n, nz_num, iA, jA, A, x, rhs, itr_max, mr, tol_abs,tol_rel);
	pmgmres_ilu_cr (n,nz_num,iAcr,jA,A,x,rhs,itr_max,mr,tol_abs,tol_rel);
#ifdef DBGMGM
	printf("\n Returned from mgmres call \n");
#endif
	//#pragma omp parallel for
	 for (int i = 0; i<nx-1; i++){
	    	for(int j = 0; j<ny-1; j++){
	    		int ind = get_pssn_ind(i,j,nx,ny);
	    		phi[i+nghost][j+nghost] = x[ind];
	    	}
	 }

		for(int i = 0;i<nx;i++){
			phi[i+nghost][nghost-1] = phi[i+nghost][nghost+1]; //lower wall
			phi[i+nghost][ny+nghost-1] = phi[i+nghost][ny-3+nghost]; //upper wall
		}

		for(int j = 0;j<ny;j++){
			phi[nghost-1][j+nghost] = phi[nghost+1][j+nghost];
			phi[nx+nghost-1][j+nghost] = phi[nx-3+nghost][j+nghost];
		}
#ifdef DBGMGM
	printf("\n Values assigned to phi \n A=\n[");
	for(int i = 0; i<n; i++){
		printf("\n");
		for(int j = 0; j<n; j++){
			printf("%+05.3f, ",find_in_A (A,iA,jA,nz_num,i,j));
		}
		printf(";");
	}
	printf("]\n[");
	for(int i = 0; i<n; i++){
		printf("%+05.3f; ",rhs[i]);
	}
	printf("]\n[");
	for(int i = 0; i<n; i++){
		printf("%+05.3f; ",x[i]);
	}
#endif
	double save;
	double resid = 0.0;
	for(int i=1; i<nx-2; i++){
		for(int j=1;j<ny-2;j++){
			double rpj, rmj, rpi, rmi, rho;
			rho = a[i+nghost][j+nghost]*st.rhol+(1.0-a[i+nghost][j+nghost])*st.rhog;
			rpj = a[i+nghost][j+nghost+1]*st.rhol+(1.0-a[i+nghost][j+nghost+1])*st.rhog;
			rmj = a[i+nghost][j+nghost-1]*st.rhol+(1.0-a[i+nghost][j+nghost-1])*st.rhog;
			rpi = a[i+nghost+1][j+nghost]*st.rhol+(1.0-a[i+nghost+1][j+nghost])*st.rhog;
			rmi = a[i+nghost-1][j+nghost]*st.rhol+(1.0-a[i+nghost-1][j+nghost])*st.rhog;
			double rpj2,rmj2,rpi2,rmi2;
			rpj2 = (rpj+rho)/2.0;
			rmj2 = (rmj+rho)/2.0;
			rpi2 = (rpi+rho)/2.0;
			rmi2 = (rmi+rho)/2.0;
			save = fabs( phi[i+nghost][j+nghost]*((-1.0/(rpj2*dy2)-1.0/(rmj2*dy2)-1.0/(rpi2*dx2)-1.0/(rmi2*dx2)))-
					( (us[i+nghost+1][j+nghost] - us[i+nghost][j+nghost] )/(dt*dx) + (vs[i+nghost][j+nghost+1]-vs[i+nghost][j+nghost])/(dt*dy)
							-(1.0/(rpi2*dx2)*phi[i+1+nghost][j+nghost] + 1.0/(rmi2*dx2)*phi[i-1+nghost][j+nghost]
							 +1.0/(rpj2*dy2)*phi[i+nghost][j+1+nghost] + 1.0/(rmj2*dy2)*phi[i+nghost][j-1+nghost]) ) ); //Infinity Norm
			if (save>resid){
				resid = save;
			}
		}
	}
	printf("resid = %E \n",resid);

	free(A);
	free(iA);
	free(iAcr);
	free(jA);
	//free(x);
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
double apply_projection(double** restrict phi, double** restrict a, double** restrict u, double** restrict us, double** restrict v, double** restrict vs, double dx, double dy, int nx, int ny, int nghost, double dt,struct slv_settings st)
{
	double maxv = 0.0;
	for(int i=0;i<nx-1;i++){
		for(int j=0;j<ny-1;j++){
			double rhoi = a[i+nghost-1][j+nghost]*st.rhol+(1.0-a[i+nghost-1][j+nghost])*st.rhog;
			double rpi = a[i+nghost][j+nghost]*st.rhol+(1.0-a[i+nghost][j+nghost])*st.rhog;
			double rhoj = a[i+nghost][j+nghost-1]*st.rhol+(1.0-a[i+nghost][j+nghost-1])*st.rhog;
			double rpj = a[i+nghost][j+nghost]*st.rhol+(1.0-a[i+nghost][j+nghost])*st.rhog;
			if(i!=0){
				u[i+nghost][j+nghost] = us[i+nghost][j+nghost]-dt*(2.0)/(rhoi+rpi)*(phi[i+nghost][j+nghost]-phi[i+nghost-1][j+nghost])/dx;
			}
			if(j!=0){
				v[i+nghost][j+nghost] = vs[i+nghost][j+nghost]-dt*(2.0)/(rhoj+rpj)*(phi[i+nghost][j+nghost]-phi[i+nghost][j+nghost-1])/dy;

			}
			if( fmax(fabs(u[i+nghost][j+nghost]),fabs(v[i+nghost][j+nghost])) > maxv){
				maxv = fmax(fabs(u[i+nghost][j+nghost]),fabs(v[i+nghost][j+nghost]));

			}
		}
	}	
	return maxv;
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
	st.g = -9.81;
	st.g = -1000.0;
	st.dt = 0.001;
	st.XBC = WALL;
	st.YBC = WALL;
	st.XUV = 1.0l;
	st.XLV = 0.0l;
	st.YLV = 0.0l;
	st.YRV = 0.0l;
	st.rhog = 1.0;
	st.rhol = 1.0;
	st.mug  = 1.0/100.0;
	st.mul  = 1.0/100.0;
	st.Reg = st.rhog/st.mug;
	st.Rel = st.rhol/st.mul;
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

/* Write 2D Matrix to File */
void write_matrix_2d_lseg(lseg** segs, int nx, int ny, char* filename)
{
	FILE *file = fopen(filename,"w");
	size_t count;
	for (int i = 0; i<nx; i++){
		count = fwrite(segs[i],sizeof(lseg),ny,file);
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
