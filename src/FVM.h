#ifndef FINITEVOLUME_H
#define FINITEVOLUME_H

#include "solvers.h"

class FVM
{
public:
	FVM(const vector<vec3d>& U, double lBoundary, double rBoundary);//rho,p,u
	void setTimeAxis(double endTime, double timeStep);
	void setGamma(double gamma);
	vector<vector<vec3d>> solve();

protected:
	virtual void setBoundaryValuesAndSlopes() {}; // 4
	virtual void iterateOnce() {};
	vec3d rhoPUToConserVar(const vec3d& rhoPU);
	vec3d conserVarToRhoPU(const vec3d& conserVar);
	vec3d fluxFun(const vec3d& conserVar);
	vector<vec3d> U_n_;//rho,rho*u,rho*(e+u*u/2)
	vector<vec3d> U_n, U_Slope;//rho,p,u
	vector<vector<vec3d>> Us;
	double lBoundary = -1.0, rBoundary = 1.0;
	double endTime = 1.0, timeStep = 0.01;
	double h = 0.01;
	double gamma = 1.403;
	double currentTime = 0.0;
	int currentLevel = 0;
	int spatialSize = 0;
};

class FVM_Godunov :public FVM
{
public:
	FVM_Godunov(const vector<vec3d>& U, double lBoundary, double rBoundary);

protected:
	virtual void setBoundaryValuesAndSlopes() override; // 4
	virtual void iterateOnce() override;
};

class FVM_2ndRK :public FVM_Godunov
{
public:
	FVM_2ndRK(const vector<vec3d>& U, double lBoundary, double rBoundary);

protected:
	virtual void iterateOnce() override;
	void linearizeU0();
	vec3d minmod(const vec3d& v1, const vec3d& v2);
};

class FVM_GRP :public FVM_2ndRK
{
public:
	FVM_GRP(const vector<vec3d>& U, double lBoundary, double rBoundary);
	void setAlpha(double alpha);

protected:
	virtual void iterateOnce() override;
	vec3d minmod3(const vec3d& v1, const vec3d& v2, const vec3d& v3);
	double alpha = 1.9;
};



#endif