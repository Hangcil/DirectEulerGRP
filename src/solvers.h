#ifndef SOLVERS_H
#define SOLVERS_H

#include <vector>
#include <cmath>
#include <functional>
#include "Eigen/Eigen"

using std::vector;
using std::ofstream;
using std::min;
using std::sqrt;
using std::pow;
using scalarFun1d = std::function<double(double)>;
using vec3d = Eigen::Vector3d;
using mat2d = Eigen::Matrix2d;
using vec2d = Eigen::Vector2d;

enum waveType
{
	CRW,
	Shock,
	CD,
	Vacuum,
};

class RPSolver
{
public:
	RPSolver(const vec3d& lState, const vec3d& rState); //[rho,p,u]
	void setGamma(double gamma);
	void setTol(double tol);
	void solve();

	vec3d operator()(double x, double t);
	double rho_starL();
	double rho_starR();
	double p_star();
	double u_star();
	waveType wave1();
	waveType wave3();


protected:
	double NewtonMethod(const scalarFun1d& f, const scalarFun1d& df, double start);
	void initVar();
	double gamma = 1.403;
	double tol = 1.0e-10;
	double rho_L = 0.0, rho_R = 0.0, p_L = 0.0, p_R = 0.0, u_L = 0.0, u_R = 0.0, c_L = 0.0, c_R = 0.0;
	double mu_2 = 0.0, psi_L = 0.0, phi_R = 0.0, A_L = 0.0, A_R = 0.0, B_L = 0.0, B_R = 0.0;
	double rho_starL_ = 0.0, rho_starR_ = 0.0, p_star_ = 0.0, u_star_ = 0.0, c_starL = 0.0, c_starR = 0.0;
	waveType wave1_ = Shock, wave2 = CD, wave3_ = Shock;
};

class GRPSolver :public RPSolver
{
public:
	GRPSolver(const vec3d& lState, const vec3d& rState, const vec3d& lStateSlope, const vec3d& rStateSlope); //[rho,p,u]
	void solve();
	vec3d timeDerivative();

protected:
	double rho_t = 0.0, p_t = 0.0, u_t = 0.0;
	vec3d lStateSlope = vec3d::Zero();
	vec3d rStateSlope = vec3d::Zero();
};

#endif