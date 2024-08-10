#include "FVM.h"

FVM::FVM(const vector<vec3d>& U, double lBoundary, double rBoundary) : U_n(U), lBoundary(lBoundary), rBoundary(rBoundary)
{
    spatialSize = int(U.size());
    U_Slope = vector<vec3d>(spatialSize);
    for (auto i = 0; i < spatialSize; i++)
    {
        this->U_n_.push_back(rhoPUToConserVar(U[i]));
    }
    Us.push_back(U);
    h = (rBoundary - lBoundary) / double(spatialSize);
}

void FVM::setTimeAxis(double endTime, double timeStep)
{
    this->endTime = endTime;
    this->timeStep = timeStep;
}

void FVM::setGamma(double gamma)
{
    U_n_.clear();
    this->gamma = gamma;
    for (auto i = 0; i < spatialSize; i++)
    {
        this->U_n_.push_back(rhoPUToConserVar(U_n[i]));
    }
}

vec3d FVM::rhoPUToConserVar(const vec3d& rhoPU)
{
    return { rhoPU(0), 
        rhoPU(0) * rhoPU(2), 
        rhoPU(1) / (gamma - 1.0) + rhoPU(0) * rhoPU(2) * rhoPU(2) / 2.0 };
}

vec3d FVM::conserVarToRhoPU(const vec3d& conserVar)
{
    return { conserVar(0) ,
        (gamma - 1.0) * (conserVar(2) - conserVar(1) / conserVar(0) * conserVar(1) / 2.0) ,
        conserVar(1) / conserVar(0) };
}

vec3d FVM::fluxFun(const vec3d& conserVar)
{
    return { conserVar(1) ,
        (gamma - 1.0) * conserVar(2) + 0.5 * (3.0 - gamma) * conserVar(1) * conserVar(1) / conserVar(0) ,
        conserVar(1) / conserVar(0) * (gamma * conserVar(2) - 0.5 * (gamma - 1.0) * conserVar(1) * conserVar(1) / conserVar(0)) };
}

vector<vector<vec3d>> FVM::solve()
{
    if (h <= 0 || timeStep <= 0 || endTime <= 0 || endTime < timeStep || gamma < 1.0)
    {
        return {};
    }
    currentTime = timeStep;
    currentLevel = 1;
    while (currentTime <= endTime)
    {
        currentTime += timeStep;
        currentLevel++;
        iterateOnce();
        setBoundaryValuesAndSlopes();
        Us.push_back(U_n);
    }
    return Us;
}



FVM_Godunov::FVM_Godunov(const vector<vec3d>& U, double lBoundary, double rBoundary) :FVM(U, lBoundary, rBoundary)
{
}

void FVM_Godunov::setBoundaryValuesAndSlopes()
{
    if (spatialSize <= 4)
    {
        return;
    }
    if (spatialSize <= 6)
    {
        U_n[0] = U_n[2];
        U_n[1] = U_n[2];
        U_n[spatialSize - 2] = U_n[spatialSize - 3];
        U_n[spatialSize - 1] = U_n[spatialSize - 3];
        U_n_[0] = rhoPUToConserVar(U_n[0]);
        U_n_[1] = rhoPUToConserVar(U_n[1]);
        U_n_[spatialSize - 1] = rhoPUToConserVar(U_n[spatialSize - 1]);
        U_n_[spatialSize - 2] = rhoPUToConserVar(U_n[spatialSize - 2]);
    }
    else
    {
        vec3d d_L = U_n[3] - U_n[2];
        U_n[0] = U_n[2] - 2.0 * d_L;
        U_n[1] = U_n[2] - d_L;
        vec3d d_R = U_n[spatialSize - 3] - U_n[spatialSize - 4];
        U_n[spatialSize - 1] = U_n[spatialSize - 3] + 2.0 * d_R;
        U_n[spatialSize - 2] = U_n[spatialSize - 3] + d_R;
        U_n_[0] = rhoPUToConserVar(U_n[0]);
        U_n_[1] = rhoPUToConserVar(U_n[1]);
        U_n_[spatialSize - 1] = rhoPUToConserVar(U_n[spatialSize - 1]);
        U_n_[spatialSize - 2] = rhoPUToConserVar(U_n[spatialSize - 2]);

        U_Slope[0] = d_L / h;
        U_Slope[1] = d_L / h;
        U_Slope[2] = d_L / h;
        U_Slope[spatialSize - 1] = d_R / h;
        U_Slope[spatialSize - 2] = d_R / h;
        U_Slope[spatialSize - 3] = d_R / h;
    }
}

void FVM_Godunov::iterateOnce()
{
    auto GodunovSols = vector<vec3d>(spatialSize - 1);

    for (auto i = 1; i < spatialSize - 2; i++) // i+1/2
    {
        auto localRPSolver = RPSolver(U_n[i], U_n[i + 1]);
        localRPSolver.setGamma(gamma);
        auto localRPSol = localRPSolver.solve();
        GodunovSols[i] = localRPSol.GodunovSol;
    }

    for (auto i = 2; i < spatialSize - 2; i++)
    {
        vec3d F_L = fluxFun(rhoPUToConserVar(GodunovSols[i - 1]));
        vec3d F_R = fluxFun(rhoPUToConserVar(GodunovSols[i]));
        U_n_[i] = U_n_[i] - timeStep / h * (F_R - F_L);
        U_n[i] = conserVarToRhoPU(U_n_[i]);
    }
}



FVM_2ndRK::FVM_2ndRK(const vector<vec3d>& U, double lBoundary, double rBoundary) :FVM_Godunov(U, lBoundary, rBoundary)
{
    linearizeU0();
}

void FVM_2ndRK::linearizeU0()
{
    for (auto j = 1; j < spatialSize - 1; j++)
    {
        U_Slope[j] = minmod(U_n[j + 1] - U_n[j], U_n[j] - U_n[j - 1]) / h;
    }
}

vec3d FVM_2ndRK::minmod(const vec3d& v1, const vec3d& v2)
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

void FVM_2ndRK::iterateOnce()
{
    auto GodunovSols = vector<vec3d>(spatialSize - 1);

    for (auto i = 1; i < spatialSize - 2; i++) // i+1/2
    {
        vec3d U_n_i_L = U_n[i] + U_Slope[i] * h / 2.0;
        vec3d U_n_i_R = U_n[i + 1] - U_Slope[i + 1] * h / 2.0;

        auto localRPSolver = RPSolver(U_n_i_L, U_n_i_R);
        localRPSolver.setGamma(gamma);
        auto localRPSol = localRPSolver.solve();
        GodunovSols[i] = localRPSol.GodunovSol;
    }

    auto L = [&](const vec3d& u, int cell_ind) -> vec3d
        {
            vec3d F_L = fluxFun(rhoPUToConserVar(GodunovSols[cell_ind - 1]));
            vec3d F_R = fluxFun(rhoPUToConserVar(GodunovSols[cell_ind]));
            return (F_L - F_R) / h;
        };

    vector<vec3d> u1(spatialSize), u1_(spatialSize);
    for (auto i = 2; i < spatialSize - 2; i++)
    {
        u1_[i] = U_n_[i] + timeStep * L(U_n[i], i);
        u1[i](0) = u1_[i](0);
        u1[i](2) = u1_[i](1) / u1_[i](0);
        u1[i](1) = (gamma - 1.0) * (u1_[i](2) - u1[i](2) * u1_[i](1) / 2.0);
        U_n_[i] = 0.5 * U_n_[i] + 0.5 * u1_[i] + 0.5 * timeStep * L(u1[i], i);
        U_n[i] = conserVarToRhoPU(U_n_[i]);
    }

    for (auto j = 3; j < spatialSize - 3; j++)
    {
        U_Slope[j] = minmod((U_n[j + 1] - U_n[j]) / h, (U_n[j] - U_n[j - 1]) / h);
    }
}



FVM_GRP::FVM_GRP(const vector<vec3d>& U, double lBoundary, double rBoundary) :FVM_2ndRK(U, lBoundary, rBoundary)
{
}

void FVM_GRP::setAlpha(double alpha)
{
    this->alpha = alpha >= 0.0 && alpha < 2.0 ? alpha : this->alpha;
}

vec3d FVM_GRP::minmod3(const vec3d& v1, const vec3d& v2, const vec3d& v3)
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

void FVM_GRP::iterateOnce()
{
    auto timeDeris = vector<vec3d>(spatialSize - 1);
    auto GodunovSols = vector<vec3d>(spatialSize - 1);

    auto rhoPUTimeDeriToConserVarTimeDeri = [&](const vec3d& rhoPU, const vec3d& timeDeri, double k)->vec3d
        {
            vec3d ret;
            double rho = rhoPU(0), p = rhoPU(1), u = rhoPU(2);
            double rho_t = timeDeri(0), p_t = timeDeri(1), u_t = timeDeri(2);
            ret(0) = rho + k * rho_t;
            ret(1) = rho * u + k * (rho_t * u + rho * u_t);
            ret(2) = p / (gamma - 1.0) + rho * u * u / 2.0 + k * (p_t / (gamma - 1.0) + rho_t * u * u / 2.0 + rho * u * u_t);
            return ret;
        };

    for (auto i = 1; i < spatialSize - 2; i++) // i+1/2
    {
        vec3d U_n_i_L = U_n[i] + U_Slope[i] * h / 2.0;
        vec3d U_n_i_R = U_n[i + 1] - U_Slope[i + 1] * h / 2.0;

        auto localGRPSolver = GRPSolver(U_n_i_L, U_n_i_R, U_Slope[i], U_Slope[i + 1]);
        localGRPSolver.setGamma(gamma);
        auto localGRPSol = localGRPSolver.solve();
        timeDeris[i] = localGRPSol.timeDerivative;
        GodunovSols[i] = localGRPSol.GodunovSol;
    }

    auto Slopes = vector<vec3d>(spatialSize);
    for (auto i = 2; i < spatialSize - 2; i++)
    {
        vec3d U_L_n_p1 = GodunovSols[i - 1] + timeStep * timeDeris[i - 1];
        vec3d U_R_n_p1 = GodunovSols[i] + timeStep * timeDeris[i];
        Slopes[i] = (U_R_n_p1 - U_L_n_p1) / h;

        vec3d F_L = fluxFun(rhoPUTimeDeriToConserVarTimeDeri(GodunovSols[i - 1], timeDeris[i - 1], timeStep / 2.0));
        vec3d F_R = fluxFun(rhoPUTimeDeriToConserVarTimeDeri(GodunovSols[i], timeDeris[i], timeStep / 2.0));

        U_n_[i] = U_n_[i] - timeStep / h * (F_R - F_L);
        U_n[i] = conserVarToRhoPU(U_n_[i]);
    }

    for (auto j = 3; j < spatialSize - 3; j++)
    {
        U_Slope[j] = minmod3(alpha * (U_n[j + 1] - U_n[j]) / h, alpha * (U_n[j] - U_n[j - 1]) / h, Slopes[j]);
    }
}



