#include "RPSolver.h"

RPSolver::RPSolver(const vec3d& lState, const vec3d& rState) : lState(std::move(lState)), rState(std::move(rState))
{
}

void RPSolver::setGamma(double gamma)
{
    this->gamma = 1.0 < gamma && gamma <= 5.0 / 3.0 ? gamma : this->gamma;
}

double RPSolver::NewtonMethod(const scalarFun1d& f, const scalarFun1d& df, double start)
{
    double save = start < 0.0 ? 1.0e-6 : start, p = save - f(save) / df(save), cha = 2.0 * abs(p - save) / abs(p + save);
    int i = 0;
    while (cha >= 1.0e-6 && i <= 1000)
    {
        if (p < 0.0)
        {
            p = 1.0e-6;
        }
        i++;
        save = p;
        p = p - f(p) / df(p);
        cha = 2.0 * abs(p - save) / abs(p + save);
    }
    return p;
}

RPSol RPSolver::solve(int level, int i)
{
    RPSol ret;
    ret.lState = this->lState;
    ret.rState = this->rState;

    double rho_l = lState(0), rho_r = rState(0), p_l = lState(1), p_r = rState(1), u_l = lState(2), u_r = rState(2);
    double c_l = sqrt(gamma * p_l / rho_l), c_r = sqrt(gamma * p_r / rho_r);
    double mu = (gamma - 1.0) / (gamma + 1.0);
    double psi_l = u_l + 2.0 * c_l / (gamma - 1.0);
    double phi_r = u_r - 2.0 * c_r / (gamma - 1.0);
    double A_l = 2.0 / (gamma + 1.0) / rho_l, B_l = mu * p_l, A_r = 2.0 / (gamma + 1.0) / rho_r, B_r = mu * p_r;

    auto f_L = [&](double p) -> double
        {
            if (p > p_l)
            {
                return (p - p_l) * sqrt(A_l / (p + B_l));
            }
            else
            {
                return 2.0 * c_l / (gamma - 1.0) * (pow(p / p_l, (gamma - 1.0) / 2.0 / gamma) - 1.0);
            }
        };

    auto f_R = [&](double p) -> double
        {
            if (p > p_r)
            {
                return (p - p_r) * sqrt(A_r / (p + B_r));
            }
            else
            {
                return 2.0 * c_r / (gamma - 1.0) * (pow(p / p_r, (gamma - 1.0) / 2.0 / gamma) - 1.0);
            }
        };

    auto df_L = [&](double p) -> double
        {
            if (p > p_l)
            {
                return sqrt(A_l / (B_l + p)) * (1.0 - (p - p_l) / 2.0 / (B_l + p));
            }
            else
            {
                return 1.0 / rho_l / c_l * pow(p / p_l, -(gamma + 1.0) / 2.0 / gamma);
            }
        };

    auto df_R = [&](double p) -> double
        {
            if (p > p_r)
            {
                return sqrt(A_r / (B_r + p)) * (1.0 - (p - p_r) / 2.0 / (B_r + p));
            }
            else
            {
                return 1.0 / rho_r / c_r * pow(p / p_r, -(gamma + 1.0) / 2.0 / gamma);
            }
        };

    auto df = [&](double p) -> double
        {
            return df_L(p) + df_R(p);
        };

    auto f = [&](double p) -> double
        {
            return f_L(p) + f_R(p) + u_r - u_l;
        };

    double p_star = NewtonMethod(f, df, 0.5 * (p_l + p_r));
    double u_star = 0.5 * (u_l + u_r + f_R(p_star) - f_L(p_star));
    double rho_ml = 0.0, rho_mr = 0.0;

    if (p_star > p_l)
    {
        ret.wave1 = Shock;
        rho_ml = rho_l * (mu * p_l + p_star) / (mu * p_star + p_l);
    }
    else
    {
        ret.wave1 = CRW;
        rho_ml = rho_l * pow(p_star / p_l, 1.0 / gamma);
    }

    if (p_star > p_r)
    {
        ret.wave3 = Shock;
        rho_mr = rho_r * (mu * p_r + p_star) / (mu * p_star + p_r);
    }
    else
    {
        ret.wave3 = CRW;
        rho_mr = rho_r * pow(p_star / p_r, 1.0 / gamma);
    }
    ret.mlState = vec3d(rho_ml, p_star, u_star);
    ret.mrState = vec3d(rho_mr, p_star, u_star);

    if (ret.wave1 == CRW && ret.wave3 == CRW)
    {
        double c_ml = sqrt(gamma * p_star / rho_ml);
        double c_mr = sqrt(gamma * p_star / rho_mr);
        if (u_l - c_l >= 0.0)
        {
            ret.GodunovSol(0) = rho_l;
            ret.GodunovSol(1) = p_l;
            ret.GodunovSol(2) = u_l;
        }
        else if (u_r + c_r <= 0.0)
        {

            ret.GodunovSol(0) = rho_r;
            ret.GodunovSol(1) = p_r;
            ret.GodunovSol(2) = u_r;
        }
        else
        {
            if (u_star - c_ml >= 0)
            {
                double c_0_0 = mu * psi_l;
                double u_0_0 = c_0_0;
                double rho_0_0 = pow(c_0_0 * c_0_0 * pow(rho_l, gamma) / gamma / p_l, 1.0 / (gamma - 1.0));
                ret.GodunovSol(0) = rho_0_0;
                ret.GodunovSol(1) = rho_0_0 * c_0_0 * c_0_0 / gamma;
                ret.GodunovSol(2) = u_0_0;
            }
            else if (u_star + c_mr <= 0)
            {
                double c_0_0 = -mu * phi_r;
                double u_0_0 = -c_0_0;
                double rho_0_0 = pow(c_0_0 * c_0_0 * pow(rho_r, gamma) / gamma / p_r, 1.0 / (gamma - 1.0));
                ret.GodunovSol(0) = rho_0_0;
                ret.GodunovSol(1) = rho_0_0 * c_0_0 * c_0_0 / gamma;
                ret.GodunovSol(2) = u_0_0;
            }
            else
            {
                ret.GodunovSol(1) = p_star;
                ret.GodunovSol(2) = u_star;
                if (u_star >= 0.0)
                {
                    ret.GodunovSol(0) = rho_ml;
                }
                else
                {
                    ret.GodunovSol(0) = rho_mr;
                }
            }
        }
    }
    else if (ret.wave1 == Shock && ret.wave3 == Shock)
    {
        double ARPSol_sigmaR = 0.0;
        if (abs(rho_r - rho_mr) >= 10e-10)
        {
            ARPSol_sigmaR = (rho_r * u_r - rho_mr * u_star) / (rho_r - rho_mr);
        }
        else
        {
            ARPSol_sigmaR = u_r + c_r;
        }
        double ARPSol_sigmaL = 0.0;
        if (abs(rho_l - rho_ml) >= 10e-10)
        {
            ARPSol_sigmaL = (rho_l * u_l - rho_ml * u_star) / (rho_l - rho_ml);
        }
        else
        {
            ARPSol_sigmaL = u_l - c_l;
        }
        if (ARPSol_sigmaR <= 0.0)
        {
            ret.GodunovSol(0) = rho_r;
            ret.GodunovSol(1) = p_r;
            ret.GodunovSol(2) = u_r;
        }
        else if (ARPSol_sigmaL >= 0.0)
        {
            ret.GodunovSol(0) = rho_l;
            ret.GodunovSol(1) = p_l;
            ret.GodunovSol(2) = u_l;
        }
        else
        {
            ret.GodunovSol(1) = p_star;
            ret.GodunovSol(2) = u_star;
            if (u_star >= 0.0)
            {
                ret.GodunovSol(0) = rho_ml;
            }
            else
            {
                ret.GodunovSol(0) = rho_mr;
            }
        }
    }
    else if (ret.wave1 == Shock && ret.wave3 == CRW)
    {
        double ARPSol_sigma = 0.0;
        if (abs(rho_l - rho_ml) >= 10e-10)
        {
            ARPSol_sigma = (rho_l * u_l - rho_ml * u_star) / (rho_l - rho_ml);
        }
        else
        {
            ARPSol_sigma = u_l - c_l;
        }
        if (ARPSol_sigma >= 0)
        {
            ret.GodunovSol(0) = rho_l;
            ret.GodunovSol(1) = p_l;
            ret.GodunovSol(2) = u_l;
        }
        else if (u_r + c_r <= 0)
        {
            ret.GodunovSol(0) = rho_r;
            ret.GodunovSol(1) = p_r;
            ret.GodunovSol(2) = u_r;
        }
        else
        {
            double c_mr = sqrt(gamma * p_star / rho_mr);
            if (u_star + c_mr <= 0)
            {
                double c_0_0 = -mu * phi_r;
                double u_0_0 = -c_0_0;
                double rho_0_0 = pow(c_0_0 * c_0_0 * pow(rho_r, gamma) / gamma / p_r, 1.0 / (gamma - 1.0));
                ret.GodunovSol(0) = rho_0_0;
                ret.GodunovSol(1) = rho_0_0 * c_0_0 * c_0_0 / gamma;
                ret.GodunovSol(2) = u_0_0;
            }
            else
            {
                ret.GodunovSol(1) = p_star;
                ret.GodunovSol(2) = u_star;
                if (u_star >= 0.0)
                {
                    ret.GodunovSol(0) = rho_ml;
                }
                else
                {
                    ret.GodunovSol(0) = rho_mr;
                }
            }
        }
    }
    else
    {
        double ARPSol_sigma = 0.0;
        if (abs(rho_r - rho_mr) >= 10e-10)
        {
            ARPSol_sigma = (rho_r * u_r - rho_mr * u_star) / (rho_r - rho_mr);
        }
        else
        {
            ARPSol_sigma = u_r + c_r;
        }

        if (u_l - c_l >= 0.0)
        {
            ret.GodunovSol(0) = rho_l;
            ret.GodunovSol(1) = p_l;
            ret.GodunovSol(2) = u_l;
        }
        else if (ARPSol_sigma <= 0.0)
        {
            ret.GodunovSol(0) = rho_r;
            ret.GodunovSol(1) = p_r;
            ret.GodunovSol(2) = u_r;
        }
        else
        {
            double c_ml = sqrt(gamma * p_star / rho_ml);
            if (u_star - c_ml >= 0)
            {
                double c_0_0 = mu * psi_l;
                double u_0_0 = c_0_0;
                double rho_0_0 = pow(c_0_0 * c_0_0 * pow(rho_l, gamma) / gamma / p_l, 1.0 / (gamma - 1.0));
                ret.GodunovSol(0) = rho_0_0;
                ret.GodunovSol(1) = rho_0_0 * c_0_0 * c_0_0 / gamma;
                ret.GodunovSol(2) = u_0_0;
            }
            else
            {
                ret.GodunovSol(1) = p_star;
                ret.GodunovSol(2) = u_star;
                if (u_star >= 0.0)
                {
                    ret.GodunovSol(0) = rho_ml;
                }
                else
                {
                    ret.GodunovSol(0) = rho_mr;
                }
            }
        }
    }
    return ret;
}