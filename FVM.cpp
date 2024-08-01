#include "FVM.h"

FVMSolver::FVMSolver(const vector<vec3d> &U, double lBoundary, double rBoundary) : U(U), lBoundary(lBoundary), rBoundary(rBoundary)
{
    spatialSize = int(U.size());
    U_Slope = vector<vec3d>(spatialSize);
    for (auto i = 0; i < spatialSize; i++)
    {
        vec3d U__;
        U__ << U[i](0), U[i](0)* U[i](2), U[i](1) / (gamma - 1.0) + U[i](0) * U[i](2) * U[i](2) / 2.0;
        this->U_.push_back(U__);
    }
    Us.push_back(U);
    h = (rBoundary - lBoundary) / double(spatialSize);
}

void FVMSolver::setTimeAxis(double endTime, double timeStep)
{
    this->endTime = endTime;
    this->timeStep = timeStep;
}

void FVMSolver::setGamma(double gamma)
{
    U_.clear();
    this->gamma = gamma;
    for (auto i = 0; i < spatialSize; i++)
    {
        vec3d U__;
        U__ << U[i](0), U[i](0)* U[i](2), U[i](1) / (gamma - 1.0) + U[i](0) * U[i](2) * U[i](2) / 2.0;
        this->U_.push_back(U__);
    }
}

void FVMSolver::setAlpha(double alpha)
{
    this->alpha = alpha >= 0.0 && alpha < 2.0 ? alpha : this->alpha;
}

void FVMSolver::compleBoundary()
{
    if (spatialSize <= 4)
    {
        return;
    }
    if (spatialSize <= 6)
    {
        U[0] = U[2];
        U[1] = U[2];
        U[spatialSize - 2] = U[spatialSize - 3];
        U[spatialSize - 1] = U[spatialSize - 3];
    }
    else
    {
        vec3d d_L = U[3] - U[2];
        U[0] = U[2] - 2.0 * d_L;
        U[1] = U[2] - d_L;
        vec3d d_R = U[spatialSize - 3] - U[spatialSize - 4];
        U[spatialSize - 1] = U[spatialSize - 3] + 2.0 * d_R;
        U[spatialSize - 2] = U[spatialSize - 3] + d_R;
        U_Slope[0] = d_L / h;
        U_Slope[1] = d_L / h;
        U_Slope[2] = d_L / h;
        U_Slope[spatialSize - 1] = d_R / h;
        U_Slope[spatialSize - 2] = d_R / h;
        U_Slope[spatialSize - 3] = d_R / h;


        vec3d d_L_ = U_[3] - U_[2];
        U_[0] = U_[2] - 2.0 * d_L_;
        U_[1] = U_[2] - d_L_;
        vec3d d_R_ = U_[U_.size() - 3] - U_[U_.size() - 4];
        U_[U_.size() - 1] = U_[U_.size() - 3] + 2.0 * d_R_;
        U_[U_.size() - 2] = U_[U_.size() - 3] + d_R_;
    }
}

vec3d FVMSolver::minmod(const vec3d &v1, const vec3d &v2)
{
    vec3d ret;
    for (auto i = 0; i <= 2; i++)
    {
        if (v1(i) * v2(i) <= 0.0)
        {
            ret(i) = 0.0;
        }
        else if (abs(v1(i)) <= abs(v2(i)))
        {
            ret(i) = v1(i);
        }
        else
        {
            ret(i) = v2(i);
        }
    }
    return ret;
}

vec3d FVMSolver::minmod3(const vec3d& v1, const vec3d& v2, const vec3d& v3)
{
    vec3d ret;
    for (auto i = 0; i <= 2; i++)
    {
        double sgn13 = v1(i) * v3(i), sgn23 = v2(i) * v3(i);
        if (sgn13 > 0.0 && sgn23 > 0.0)
        {
            double sgn = 1.0;
            if (v1(i) <= 0.0)
            {
                sgn = -1.0;
            }
            ret(i) = sgn * min(abs(v1(i)), min(abs(v2(i)), abs(v3(i))));
        }
        else
        {
            ret(i) = 0.0;
        }
    }
    return ret;
}

void FVMSolver::linearizeU0()
{
    for (auto j = 1; j < spatialSize - 1; j++)
    {
        U_Slope[j] = minmod(U[j + 1] - U[j], U[j] - U[j - 1]) / h;
    }
}

vector<vector<vec3d>> FVMSolver::solve()
{
    if (h <= 0 || timeStep <= 0 || endTime <= 0 || endTime < timeStep || gamma < 1.0)
    {
        return {};
    }

    currentTime = timeStep;
    linearizeU0();
    while (currentTime <= endTime)
    {
        currentLevel += 1;
        currentTime += timeStep;
        auto timeDeris = vector<vec3d>(spatialSize - 1);
        auto GodunovSols = vector<vec3d>(spatialSize - 1);
        for (auto i = 1; i < spatialSize - 2; i++) // i+1/2
        {
            vec3d U_Li = U[i] + U_Slope[i] * h / 2.0;
            vec3d U_Ri = U[i + 1] - U_Slope[i + 1] * h / 2.0;

            auto localGRPSolver = GRPSolver(U_Li, U_Ri, U_Slope[i], U_Slope[i + 1]);
            localGRPSolver.setGamma(gamma);
            auto localGRPSol = localGRPSolver.solve(currentLevel, i);
            timeDeris[i] = localGRPSol.timeDerivative;
            GodunovSols[i] = localGRPSol.ARPSol.GodunovSol;
        }
        auto U_nextLevel = vector<vec3d>(spatialSize);
        auto U_nextLevel_ = vector<vec3d>(spatialSize);
        auto Slopes = vector<vec3d>(spatialSize);
        for (auto i = 2; i < spatialSize - 2; i++)
        {
            vec3d U_L = GodunovSols[i - 1] + timeStep / 2.0 * timeDeris[i - 1];
            vec3d U_R = GodunovSols[i] + timeStep / 2.0 * timeDeris[i];

            vec3d U_L_ = GodunovSols[i - 1] + timeStep * timeDeris[i - 1];
            vec3d U_R_ = GodunovSols[i] + timeStep * timeDeris[i];
            Slopes[i] = (U_R_ - U_L_) / h;
     
            vec3d F_L;
            F_L << U_L(0) * U_L(2), U_L(0)* U_L(2)* U_L(2) + U_L(1), U_L(2)* (gamma * U_L(1) / (gamma - 1.0) + U_L(0) * U_L(2) * U_L(2) / 2.0);
            vec3d F_R;
            F_R << U_R(0) * U_R(2), U_R(0)* U_R(2)* U_R(2) + U_R(1), U_R(2)* (gamma * U_R(1) / (gamma - 1.0) + U_R(0) * U_R(2) * U_R(2) / 2.0);
            U_nextLevel_[i] = U_[i] - timeStep / h * (F_R - F_L);
            U_nextLevel[i](0) = U_nextLevel_[i](0);
            U_nextLevel[i](2) = U_nextLevel_[i](1) / U_nextLevel_[i](0);
            U_nextLevel[i](1) = (gamma - 1.0) * (U_nextLevel_[i](2) - U_nextLevel[i](2) * U_nextLevel_[i](1) / 2.0);
        }
        U = U_nextLevel;
        U_ = U_nextLevel_;
        compleBoundary();
        for (auto j = 3; j < spatialSize - 3; j++)
        {
            U_Slope[j] = minmod3(alpha * (U[j + 1] - U[j]) / h, alpha * (U[j] - U[j - 1]) / h, Slopes[j]);
            //U_Slope[j] = vec3d::Zero();
        }
        Us.push_back(U);
    }
    return Us;
}