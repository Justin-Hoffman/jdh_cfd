/*
 * lvlset.h
 *
 *  Created on: Jan 30, 2014
 *      Author: justin
 */

#ifndef LVLSET_H_
#define LVLSET_H_

	double dist_from_arc(double x, double y, double xc, double yc, double R, double t0, double t1);

	double dist_from_pt(double x, double y, double xc, double yc);

	double dist_from_circ(double x, double y, double xc, double yc, double R);

	double dist_from_line(double x, double y, double x0, double y0, double x1, double y1);

	double dist_from_vec(double x, double y, double x0, double y0);

	double vec_dot(double x0, double y0, double x1, double y1);

	double vec_cross(double x0, double y0, double x1, double y1);

	double vec_angle(double x0, double y0, double x1, double y1);

	double vec_mag(double x0, double y0);

	double dist_from_box(double x, double y, double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3);

	double dist_from_notch(double x, double y, double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3);

	double dist_from_zalesak(double x,double y,double xc,double yc, double R, double cwidth, double cdepth, double cangle);

	void init_zalesak(double** G, int nx, int ny, int nghost, double dx, double dy);

	void init_circle(double** G, int nx, int ny, int nghost, double dx, double dy);

	void init_uv_test(double** u, double** v, int nx, int ny, int nghost, double dx, double dy);


	void init_uv_test_static(double** u, double** v, int nx, int ny, int nghost, double dx, double dy);

	void levelset_diff(double** restrict G,double** restrict dGdt,double** restrict dGdxp,double** restrict dGdxm, double** restrict dGdyp, double** restrict dGdym, double** restrict u, double** restrict v, double dx, double dy, int nx, int ny, int nghost);

	void levelset_advect_euler(double** restrict G,double** restrict dGdt,double** restrict dGdxp,double** restrict dGdxm,double** restrict dGdyp,double** restrict dGdym, double** restrict u, double** restrict v, double dx, double dy, double dt, int nx, int ny, int nghost);

	void levelset_advect_TVDRK3(double** restrict G,double** restrict G1,double** restrict G2,double** restrict dGdt,double** restrict dGdt1,double** restrict dGdt2,double** restrict dGdxp,double** restrict dGdxm,double** restrict dGdyp,double** restrict dGdym, double** restrict u, double** restrict v, double dx, double dy, double dt,double time, int nx, int ny, int nghost);

	void reinit_diff(double** restrict G, double** restrict G0, double** restrict dGdt,double** restrict dGdxp,double** restrict dGdxm,double** restrict dGdyp,double** restrict dGdym, double dx, double dy, int nx, int ny, int nghost);

	double reinit_advect_TVDRK3(double** restrict G, double** restrict G0, double** restrict G1,double** restrict G2,double** restrict dGdt,double** restrict dGdt1,double** restrict dGdt2,double** restrict dGdxp,double** restrict dGdxm,double** restrict dGdyp,double** restrict dGdym, double dx, double dy, double dt, int nx, int ny, int nghost);

	double get_vol(double** restrict G, double a, int nx, int ny, int nghost,double dx, double dy);

	double get_shape_err(double** restrict G, double** restrict Gt, double volGt, double a, int nx, int ny, int nghost,double dx, double dy);
#endif /* LVLSET_H_ */
