#ifndef RPSOLVER_H
#define RPSOLVER_H

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

enum waveType
{
	CRW,
	Shock,
	CD,
	Vacuum,
};

class RPSol
{
public:
	vec3d lState = vec3d::Zero(), mlState = vec3d::Zero(), mrState = vec3d::Zero(), rState = vec3d::Zero();
	vec3d GodunovSol = vec3d::Zero();
	waveType wave1 = Shock, wave2 = CD, wave3 = Shock;
};

class RPSolver
{
public:
	RPSolver(const vec3d& lState, const vec3d& rState); //[rho,p,u]
	void setGamma(double gamma);
	RPSol solve(int level, int i);

protected:
	static double NewtonMethod(const scalarFun1d& f, const scalarFun1d& df, double start);
	vec3d lState = vec3d::Zero();
	vec3d rState = vec3d::Zero();
	double gamma = 1.403;
};

#endif