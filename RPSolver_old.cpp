#include "RPSolver.h"

RPSolver::RPSolver(const vec3d& lState, const vec3d& rState) : lState(std::move(lState)), rState(std::move(rState))
{
}

void RPSolver::setGamma(double gamma)
{
    if (1.0 < gamma && gamma <= 5.0 / 3.0)
    {
        this->gamma = gamma;
    }
}

double RPSolver::bisection(const scalarFun1d& f, double start, double end)
{
    double a = start, b = end, c = 0.0;
    int i = 0;
    while (i <= 1000)
    {
        i++;
        c = (a + b) / 2.0;
        if (f(c) >= 1.0e-15)
        {
            b = c;
        }
        else if (f(c) <= -1.0e-15)
        {
            a = c;
        }
        else
        {
            return c;
        }
    }
    return c;
}

double RPSolver::NewtonIter(const scalarFun1d &f, const scalarFun1d &df, double start)
{
    return 0.0;
}

RPSol RPSolver::solve(int level, int i)
{
    RPSol ret;
    ret.lState = this->lState;
    ret.rState = this->rState;

    auto c_rho_p = [&](double rho, double p) -> double
        {
            return sqrt(gamma * p / rho);
        };
    double rho_l = lState(0), rho_r = rState(0), p_l = lState(1), p_r = rState(1), u_l = lState(2), u_r = rState(2);
    double c_l = c_rho_p(rho_l, p_l), c_r = c_rho_p(rho_r, p_r);
    double v_1 = 2.0 / (gamma - 1.0) * (c_l + c_r);
    double v_2l = 2.0 / (gamma - 1.0) * (c_l - c_r * sqrt(rho_r / rho_l) * pow(p_l / p_r, 0.5 / gamma));
    double v_2r = 2.0 / (gamma - 1.0) * (c_r - c_l * sqrt(rho_l / rho_r) * pow(p_r / p_l, 0.5 / gamma));
    double v_3l = -(p_l - p_r) * sqrt(2.0 / rho_r / ((gamma + 1.0) * p_l + (gamma - 1.0) * p_r));
    double v_3r = -(p_r - p_l) * sqrt(2.0 / rho_l / ((gamma + 1.0) * p_r + (gamma - 1.0) * p_l));
    double mu = (gamma - 1.0) / (gamma + 1.0);
    double psi_L = u_l + 2.0 * c_l / (gamma - 1.0);
    double phi_R = u_r - 2.0 * c_r / (gamma - 1.0);

    if ((lState - rState).norm() <= 1.0e-15)
    {
        ret.wave1 = CRW;
        ret.wave2 = Smooth;
        ret.wave3 = CRW;
        ret.mlState = lState;
        ret.mrState = rState;
        ret.GodunovSol(0) = rho_l;
        ret.GodunovSol(1) = p_l;
        ret.GodunovSol(2) = u_l;
    }
    else if (u_r - u_l >= v_1)
    {
        ret.wave1 = CRW;
        ret.wave2 = Vacuum;
        ret.wave3 = CRW;
        ret.mlState = vec3d(0.0, 0.0, u_l + 2.0 * c_l / (gamma - 1.0));
        ret.mrState = vec3d(0.0, 0.0, u_r - 2.0 * c_r / (gamma - 1.0));
    }
    else if ((p_l <= p_r && u_r - u_l >= v_2r) || (p_l >= p_r && u_r - u_l >= v_2l))
    {
        ret.wave1 = CRW;
        ret.wave3 = CRW;
        double p_star = pow(((gamma - 1.0) * (u_l - u_r) / 2.0 + c_l + c_r) / sqrt(gamma) / (sqrt(pow(p_l, 1.0 / gamma) / rho_l) + sqrt(pow(p_r, 1.0 / gamma) / rho_r)), 2.0 * gamma / (gamma - 1.0));
        double rho_mr = rho_r * pow(p_star / p_r, 1.0 / gamma);
        double rho_ml = rho_l * pow(p_star / p_l, 1.0 / gamma);
        double c_ml = sqrt(gamma * p_star / rho_ml);
        double c_mr = sqrt(gamma * p_star / rho_mr);
        double u_star = u_r + 2.0 * (c_mr - c_r) / (gamma - 1.0);
        ret.mlState = vec3d(rho_ml, p_star, u_star);
        ret.mrState = vec3d(rho_mr, p_star, u_star);
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
                double c_0_0 = mu * psi_L;
                double u_0_0 = c_0_0;
                double rho_0_0 = pow(c_0_0 * c_0_0 * pow(rho_l, gamma) / gamma / p_l, 1.0 / (gamma - 1.0));
                ret.GodunovSol(0) = rho_0_0;
                ret.GodunovSol(1) = rho_0_0 * c_0_0 * c_0_0 / gamma;
                ret.GodunovSol(2) = u_0_0;
            }
            else if (u_star + c_mr <= 0)
            {
                double c_0_0 = -mu * phi_R;
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
                if (u_star >= 0)
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
    else if ((p_l <= p_r && u_r - u_l <= v_3r) || (p_l >= p_r && u_r - u_l <= v_3l))
    {
        ret.wave1 = Shock;
        ret.wave3 = Shock;
        auto f = [&](double p_star) -> double
            {
                return u_r - u_l + (p_star - p_r) * sqrt(2.0 / rho_r / ((gamma + 1.0) * p_star + (gamma - 1.0) * p_r)) + (p_star - p_l) * sqrt(2.0 / rho_l / ((gamma + 1.0) * p_star + (gamma - 1.0) * p_l));
            };
        double end = max(p_l, p_r);
        for (auto i = 0; i <= 20; i++)
        {
            if (f(end) < 0.0)
            {
                end *= 2.0;
            }
            else
            {
                break;
            }
        }
        double p_star = bisection(f, max(p_l, p_r), end);
        double rho_mr = rho_r * (p_star + mu * p_r) / (p_r + mu * p_star);
        double rho_ml = rho_l * (p_star + mu * p_l) / (p_l + mu * p_star);
        double u_star = u_r + sqrt((p_star - p_r) * (1.0 / rho_r - 1.0 / rho_mr));

        ret.mlState = vec3d(rho_ml, p_star, u_star);
        ret.mrState = vec3d(rho_mr, p_star, u_star);

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
            if (u_star >= 0)
            {
                ret.GodunovSol(0) = rho_ml;
            }
            else
            {
                ret.GodunovSol(0) = rho_mr;
            }
        }
    }
    else if (p_l <= p_r && u_r - u_l >= v_3r && u_r - u_l <= v_2r)
    {
        ret.wave1 = Shock;
        ret.wave3 = CRW;
        auto f = [&](double p_star) -> double
            {
                return u_r - u_l - 2.0 / (gamma - 1.0) * c_r + (p_star - p_l) * sqrt(2.0 / rho_l / ((gamma + 1.0) * p_star + (gamma - 1.0) * p_l)) + 2.0 / (gamma - 1.0) * sqrt(gamma * pow(p_star, (gamma - 1.0) / gamma) * pow(p_r, 1.0 / gamma) / rho_r);
            };

        double p_star = bisection(f, p_l, p_r);
        double rho_mr = rho_r * pow(p_star / p_r, 1.0 / gamma);
        double rho_ml = rho_l * (p_star + mu * p_l) / (p_l + mu * p_star);
        double u_star = u_l - sqrt((p_star - p_l) * (1.0 / rho_l - 1.0 / rho_ml));
        ret.mlState = vec3d(rho_ml, p_star, u_star);
        ret.mrState = vec3d(rho_mr, p_star, u_star);

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
                double c_0_0 = -mu * phi_R;
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
                if (u_star >= 0)
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
        ret.wave1 = CRW;
        ret.wave3 = Shock;
        auto f = [&](double p_star) -> double
            {
                return u_r - u_l - 2.0 / (gamma - 1.0) * c_l + (p_star - p_r) * sqrt(2.0 / rho_r / ((gamma + 1.0) * p_star + (gamma - 1.0) * p_r)) + 2.0 / (gamma - 1.0) * sqrt(gamma * pow(p_star, (gamma - 1.0) / gamma) * pow(p_l, 1.0 / gamma) / rho_l);
            };

        double p_star = bisection(f, p_r, p_l);
        double rho_mr = rho_r * (p_star + mu * p_r) / (p_r + mu * p_star);
        double rho_ml = rho_l * pow(p_star / p_l, 1.0 / gamma);
        double debug1 = p_star - p_r, debug2 = 1.0 / rho_r - 1.0 / rho_mr,debug3= (p_star + mu * p_r) / (p_r + mu * p_star);
        double debug5 = p_star + mu * p_r - p_r - mu * p_star;
        bool debug4 = debug3 >= 1.0;
        double rho_mr_ = rho_r * debug3;
        double u_star = u_r + sqrt((p_star - p_r) * (1.0 / rho_r - 1.0 / rho_mr));

        ret.mlState = vec3d(rho_ml, p_star, u_star);
        ret.mrState = vec3d(rho_mr, p_star, u_star);

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
                double c_0_0 = mu * psi_L;
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
                if (u_star >= 0)
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