#ifndef FINITEVOLUME_H
#define FINITEVOLUME_H

#include "GRPSolver.h"
#include <iostream>

class FVMSolver
{
public:
    FVMSolver(const vector<vec3d> &U, double lBoundary, double rBoundary);//rho,p,u
    void setTimeAxis(double endTime, double timeStep);
    void setGamma(double gamma);
    void setAlpha(double alpha);
    vector<vector<vec3d>> solve();

protected:
    vector<vec3d> U_;//rho,rho*u,rho*(e+u*u/2)
    vector<vec3d> U, U_Slope;//rho,p,u
    vector<vector<vec3d>> Us;
    double lBoundary = -1.0, rBoundary = 1.0;
    double endTime = 1.0, timeStep = 0.01;
    double h = 0.01;
    double gamma = 1.403;
    double alpha = 1.9;
    double currentTime = 0.0;
    int currentLevel = 0;
    int spatialSize = 0;
    virtual void compleBoundary(); // 4
    vec3d minmod(const vec3d& v1, const vec3d& v2);
    vec3d minmod3(const vec3d& v1, const vec3d& v2, const vec3d& v3);
    void linearizeU0();
};

#endif