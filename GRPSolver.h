#ifndef GRPSOLVER_H
#define GRPSOLVER_H

#include "RPSolver.h"
using mat2d = Eigen::Matrix2d;
using vec2d = Eigen::Vector2d;

class GRPSol
{
public:
    RPSol ARPSol;
    vec3d timeDerivative = vec3d::Zero();
};

class GRPSolver
{
public:
    GRPSolver(const vec3d &lState, const vec3d &rState, const vec3d &lStateSlope, const vec3d &rStateSlope); //[rho,p,u]
    void setGamma(double gamma);
    double getGamma();
    GRPSol solve(int level, int i);

protected:
    vec3d lState = vec3d::Zero();
    vec3d rState = vec3d::Zero();
    vec3d lStateSlope = vec3d::Zero();
    vec3d rStateSlope = vec3d::Zero();
    double gamma = 1.403;
};

#endif