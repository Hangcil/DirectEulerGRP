#include "GRPSolver.h"

GRPSolver::GRPSolver(const vec3d &lState, const vec3d &rState, const vec3d &lStateSlope, const vec3d &rStateSlope) : lState(lState), rState(rState), lStateSlope(lStateSlope), rStateSlope(rStateSlope)
{
}

void GRPSolver::setGamma(double gamma)
{
    this->gamma = 1.0 < gamma && gamma <= 5.0 / 3.0 ? gamma : this->gamma;
}

double GRPSolver::getGamma()
{
    return gamma;
}

GRPSol GRPSolver::solve(int level, int i)
{
    RPSolver ARPSolver(lState, rState);
    ARPSolver.setGamma(gamma);
    auto ARPSol = ARPSolver.solve(level, i);
    GRPSol ret;

    double mu_2 = (gamma - 1.0) / (gamma + 1.0);
    double rho_L = lState(0), p_L = lState(1), u_L = lState(2), c_L = sqrt(gamma * p_L / rho_L);
    double rho_R = rState(0), p_R = rState(1), u_R = rState(2), c_R = sqrt(gamma * p_R / rho_R);
    double rho_starL = ARPSol.mlState(0), p_star = ARPSol.mlState(1), u_star = ARPSol.mlState(2), c_starL = sqrt(gamma * p_star / rho_starL);
    double rho_starR = ARPSol.mrState(0), c_starR = sqrt(gamma * p_star / rho_starR);
    double psi_L = u_L + 2.0 * c_L / (gamma - 1.0), phi_R = u_R - 2.0 * c_R / (gamma - 1.0);
    double rho_t = 0.0, u_t = 0.0, p_t = 0.0;
 
    if ((lState - rState).norm() <= 10e-10)
    {
        u_t = -0.5 * ((u_star + c_starL) * (lStateSlope(2) + lStateSlope(1) / rho_starL / c_starL) + (u_star - c_starL) * (rStateSlope(2) - rStateSlope(1) / rho_starL / c_starL));
        p_t = -0.5 * rho_starL * c_starL * ((u_star + c_starL) * (lStateSlope(2) + lStateSlope(1) / rho_starL / c_starL) - (u_star - c_starL) * (rStateSlope(2) - rStateSlope(1) / rho_starL / c_starL));
        if (u_star >= 0)
        {
            rho_t = (p_t + u_star * (lStateSlope(1) - c_starL * c_starL * lStateSlope(0))) / c_starL / c_starL;
        }
        else
        {
            rho_t = (p_t + u_star * (rStateSlope(1) - c_starL * c_starL * rStateSlope(0))) / c_starL / c_starL;
        }
        ret.ARPSol = ARPSol;
        ret.timeDerivative(0) = rho_t;
        ret.timeDerivative(1) = p_t;
        ret.timeDerivative(1) = u_t;
        return ret;
    }

    if (ARPSol.wave1 == CRW && ARPSol.wave3 == Shock)
    {
        double ARPSol_sigma = 0.0;
        if (abs(rho_R - rho_starR) >= 10e-10)
        {
            ARPSol_sigma = (rho_R * u_R - rho_starR * u_star) / (rho_R - rho_starR);
        }
        else
        {
            ARPSol_sigma = u_R + c_R;
        }

        if (u_L - c_L >= 0.0)
        {
            rho_t = -u_L * lStateSlope(0) - rho_L * lStateSlope(2);
            u_t = -u_L * lStateSlope(2) - lStateSlope(1) / rho_L;
            p_t = -u_L * lStateSlope(1) - rho_L * c_L * c_L * lStateSlope(2);
        }
        else if (ARPSol_sigma <= 0.0)
        {
            rho_t = -u_R * rStateSlope(0) - rho_R * rStateSlope(2);
            u_t = -u_R * rStateSlope(2) - rStateSlope(1) / rho_R;
            p_t = -u_R * rStateSlope(1) - rho_R * c_R * c_R * rStateSlope(2);
        }
        else
        {
            double T_LS_Slope_L = lStateSlope(1) / (gamma - 1.0) / rho_L - c_L * c_L * lStateSlope(0) / (gamma - 1.0) / rho_L;
            double temp1_d_L = (1.0 + mu_2) / (1.0 + 2.0 * mu_2) * pow(c_starL / c_L, 0.5 / mu_2) + mu_2 / (1.0 + 2.0 * mu_2) * pow(c_starL / c_L, (1.0 + mu_2) / mu_2);
            double temp2_d_L = c_L * pow(c_starL / c_L, 0.5 / mu_2) * (lStateSlope(2) + gamma * lStateSlope(1) / (gamma - 1.0) / rho_L / c_L - c_L * lStateSlope(0) / (gamma - 1.0) / rho_L);
            double a_L = 1.0, b_L = 1.0 / (rho_starL * c_starL);
            double d_L = temp1_d_L * T_LS_Slope_L - temp2_d_L;

            double H_1 = 0.5 * sqrt((1.0 - mu_2) / rho_R / (p_star + mu_2 * p_R)) * (p_star + (1.0 + 2.0 * mu_2) * p_R) / (p_star + mu_2 * p_R);
            double H_2 = -0.5 * sqrt((1.0 - mu_2) / rho_R / (p_star + mu_2 * p_R)) * ((2.0 + mu_2) * p_star + mu_2 * p_R) / (p_star + mu_2 * p_R);
            double H_3 = -(p_star - p_R) / (2.0 * rho_R) * sqrt((1.0 - mu_2) / rho_R / (p_star + mu_2 * p_R));
            double L_p_R = -1.0 / rho_R + (ARPSol_sigma - u_R) * H_2;
            double L_u_R = ARPSol_sigma - u_R - rho_R * c_R * c_R * H_2 - rho_R * H_3;
            double L_rho_R = (ARPSol_sigma - u_R) * H_3;
            double a_R = 1.0 + rho_starR * (ARPSol_sigma - u_star) * H_1;
            double b_R = (u_star - ARPSol_sigma) / (rho_starR * c_starR * c_starR) - H_1;
            double d_R = L_p_R * rStateSlope(1) + L_u_R * rStateSlope(2) + L_rho_R * rStateSlope(0);

            if (u_star - c_starL >= 0)
            {
                double c_0_0 = mu_2 * psi_L;
                double u_0_0 = c_0_0;
                double rho_0_0 = pow(c_0_0 * c_0_0 * pow(rho_L, gamma) / gamma / p_L, 1.0 / (gamma - 1.0));
                double temp1_d_L0 = (1.0 +  mu_2) / (1.0 + 2.0 * mu_2) * pow(c_0_0 / c_L, 0.5 / mu_2) + mu_2 / (1.0 + 2.0 * mu_2) * pow(c_0_0 / c_L, (1.0 + mu_2) / mu_2);
                double temp2_d_L0 = c_L * pow(c_0_0 / c_L, 0.5 / mu_2) * (lStateSlope(2) + gamma * lStateSlope(1) / (gamma - 1.0) / rho_L / c_L - c_L * lStateSlope(0) / (gamma - 1.0) / rho_L);
                u_t = temp1_d_L0 * T_LS_Slope_L - temp2_d_L0;
                p_t = rho_0_0 * u_0_0 * u_t;
                rho_t = 1.0 / c_0_0 / c_0_0 * (p_t + (gamma - 1.0) * rho_0_0 * u_0_0 * pow(c_0_0 / c_L, (1.0 + mu_2) / mu_2) * T_LS_Slope_L);
            }
            else
            {
                mat2d A;
                A << a_L, b_L, a_R, b_R;
                vec2d x;
                x << d_L, d_R;
                vec2d u_p_material_deri = A.fullPivLu().solve(x);
                double Du_Dt = u_p_material_deri(0);
                double Dp_Dt = u_p_material_deri(1);
                if (u_star >= 0)
                {
                    rho_t = 1.0 / c_starL / c_starL * (p_t + (gamma - 1.0) * rho_starL * u_star * pow(c_starL / c_L, (1.0 + mu_2) / mu_2) * T_LS_Slope_L);
                    u_t = Du_Dt + u_star * Dp_Dt / rho_starL / c_starL / c_starL;
                    p_t = Dp_Dt + rho_starL * u_star * Du_Dt;
                }
                else
                {
                    double H_1_ = rho_R * (1.0 - mu_2 * mu_2) * p_R / (p_R + mu_2 * p_star) / (p_R + mu_2 * p_star);
                    double H_2_ = rho_R * (mu_2 * mu_2 - 1.0) * p_star / (p_R + mu_2 * p_star) / (p_R + mu_2 * p_star);
                    double H_3_ = (p_star + mu_2 * p_R) / (p_R + mu_2 * p_star);
                    double g_rho_R = u_star - ARPSol_sigma;
                    double g_p_R = ARPSol_sigma / c_starR / c_starR - u_star * H_1_;
                    double g_u_R = u_star * rho_starR * (ARPSol_sigma - u_star) * H_1_;
                    double f_R = (ARPSol_sigma - u_R) * H_2_ * rStateSlope(1) + (ARPSol_sigma - u_R) * H_3_ * rStateSlope(0) - rho_R * (H_2_ * c_R * c_R + H_3_) * rStateSlope(2);
                    rho_t = (u_star * f_R - g_p_R * Dp_Dt - g_u_R * Du_Dt) / g_rho_R;
                    u_t = Du_Dt + u_star * Dp_Dt / rho_starR / c_starR / c_starR;
                    p_t = Dp_Dt + rho_starR * u_star * Du_Dt;
                }
            }
        }
    }

    else if (ARPSol.wave1 == Shock && ARPSol.wave3 == CRW)
    {
        double ARPSol_sigma = 0.0;
        if (abs(rho_L - rho_starL) >= 10e-10)
        {
            ARPSol_sigma = (rho_L * u_L - rho_starL * u_star) / (rho_L - rho_starL);
        }
        else
        {
            ARPSol_sigma = u_L - c_L;
        }

        if (ARPSol_sigma >= 0)
        {
            rho_t = -u_L * lStateSlope(0) - rho_L * lStateSlope(2);
            u_t = -u_L * lStateSlope(2) - 1.0 / rho_L * lStateSlope(1);
            p_t = -u_L * lStateSlope(1) - rho_L * c_L * c_L * lStateSlope(2);
        }
        else if (u_R + c_R <= 0)
        {
            rho_t = -u_R * rStateSlope(0) - rho_R * rStateSlope(2);
            u_t = -u_R * rStateSlope(2) - 1.0 / rho_R * rStateSlope(1);
            p_t = -u_R * rStateSlope(1) - rho_R * c_R * c_R * rStateSlope(2);
        }
        else
        {
            double T_RS_Slope_R = rStateSlope(1) / (gamma - 1.0) / rho_R - c_R * c_R * rStateSlope(0) / (gamma - 1.0) / rho_R;
            double temp1_d_R = (1.0 + mu_2) / (1.0 + 2.0 * mu_2) * pow(c_starR / c_R, 0.5 / mu_2) + mu_2 / (1.0 + 2.0 * mu_2) * pow(c_starR / c_R, (1.0 + mu_2) / mu_2);
            double temp2_d_R = c_R * pow(c_starR / c_R, 0.5 / mu_2) * (rStateSlope(2) - gamma * rStateSlope(1) / (gamma - 1.0) / rho_R / c_R + c_R * rStateSlope(0) / (gamma - 1.0) / rho_R);
            double a_R = 1.0, b_R = -1.0 / (rho_starR * c_starR);
            double d_R = temp1_d_R * T_RS_Slope_R + temp2_d_R;

            double H_1 = 0.5 * sqrt((1.0 - mu_2) / rho_L / (p_star + mu_2 * p_L)) * (p_star + (1.0 + 2.0 * mu_2) * p_L) / (p_star + mu_2 * p_L);
            double H_2 = -0.5 * sqrt((1.0 - mu_2) / rho_L / (p_star + mu_2 * p_L)) * ((2.0 + mu_2) * p_star + mu_2 * p_L) / (p_star + mu_2 * p_L);
            double H_3 = -(p_star - p_L) / (2.0 * rho_L) * sqrt((1.0 - mu_2) / rho_L / (p_star + mu_2 * p_L));
            double L_p_L = -1.0 / rho_L - (ARPSol_sigma - u_L) * H_2;
            double L_u_L = ARPSol_sigma - u_L + rho_L * c_L * c_L * H_2 + rho_L * H_3;
            double L_rho_L = -(ARPSol_sigma - u_L) * H_3;
            double a_L = 1.0 - rho_starL * (ARPSol_sigma - u_star) * H_1;
            double b_L = (u_star - ARPSol_sigma) / (rho_starL * c_starL * c_starL) + H_1;
            double d_L = L_p_L * lStateSlope(1) + L_u_L * lStateSlope(2) + L_rho_L * lStateSlope(0);

            if (u_star + c_starR <= 0)
            {
                double c_0_0 = -mu_2 * phi_R;
                double u_0_0 = -c_0_0;
                double rho_0_0 = pow(c_0_0 * c_0_0 * pow(rho_R, gamma) / gamma / p_R, 1.0 / (gamma - 1.0));
                double temp1_d_R0 = (1.0+ mu_2) / (1.0 + 2.0 * mu_2) * pow(c_0_0 / c_R, 0.5 / mu_2) + mu_2 / (1.0 + 2.0 * mu_2) * pow(c_0_0 / c_R, (1.0 + mu_2) / mu_2);
                double temp2_d_R0 = c_R * pow(c_0_0 / c_R, 0.5 / mu_2) * (rStateSlope(2) - gamma * rStateSlope(1) / (gamma - 1.0) / rho_R / c_R + c_R * rStateSlope(0) / (gamma - 1.0) / rho_R);
                u_t = temp1_d_R0 * T_RS_Slope_R + temp2_d_R0;
                p_t = -rho_0_0 * u_0_0 * u_t;
                rho_t = 1.0 / c_0_0 / c_0_0 * (p_t + (gamma - 1.0) * rho_0_0 * u_0_0 * pow(c_0_0 / c_R, (1.0 + mu_2) / mu_2) * T_RS_Slope_R);
            }
            else
            {
                mat2d A;
                A << a_L, b_L, a_R, b_R;
                vec2d x;
                x << d_L, d_R;
                vec2d u_p_material_deri = A.fullPivLu().solve(x);
                double Du_Dt = u_p_material_deri(0);
                double Dp_Dt = u_p_material_deri(1);
                if (u_star <= 0)
                {
                    rho_t = 1.0 / c_starR / c_starR * (p_t + (gamma - 1.0) * rho_starR * u_star * pow(c_starR / c_R, (1.0 + mu_2) / mu_2) * T_RS_Slope_R);
                    u_t = Du_Dt + u_star * Dp_Dt / rho_starR / c_starR / c_starR;
                    p_t = Dp_Dt + rho_starR * u_star * Du_Dt;
                }
                else
                {
                    double H_1_ = rho_L * (1.0 - mu_2 * mu_2) * p_L / (p_L + mu_2 * p_star) / (p_L + mu_2 * p_star);
                    double H_2_ = rho_L * (mu_2 * mu_2 - 1.0) * p_star / (p_L + mu_2 * p_star) / (p_L + mu_2 * p_star);
                    double H_3_ = (p_star + mu_2 * p_L) / (p_L + mu_2 * p_star);
                    double g_rho_L = u_star - ARPSol_sigma;
                    double g_p_L = ARPSol_sigma / c_starL / c_starL - u_star * H_1_;
                    double g_u_L = u_star * rho_starL * (ARPSol_sigma - u_star) * H_1_;
                    double f_L = (ARPSol_sigma - u_L) * H_2_ * lStateSlope(1) + (ARPSol_sigma - u_L) * H_3_ * lStateSlope(0) - rho_L * (H_2_ * c_L * c_L + H_3_) * lStateSlope(2);
                    rho_t = (u_star * f_L - g_p_L * Dp_Dt - g_u_L * Du_Dt) / g_rho_L;
                    u_t = Du_Dt + u_star * Dp_Dt / rho_starL / c_starL / c_starL;
                    p_t = Dp_Dt + rho_starL * u_star * Du_Dt;
                }
            }
        }
    }

    else if (ARPSol.wave1 == CRW && ARPSol.wave3 == CRW)
    {
        if (u_L - c_L >= 0.0)
        {
            rho_t = -u_L * lStateSlope(0) - rho_L * lStateSlope(2);
            u_t = -u_L * lStateSlope(2) - 1.0 / rho_L * lStateSlope(1);
            p_t = -u_L * lStateSlope(1) - rho_L * c_L * c_L * lStateSlope(2);
        }
        else if (u_R + c_R <= 0.0)
        {
            rho_t = -u_R * rStateSlope(0) - rho_R * rStateSlope(2);
            u_t = -u_R * rStateSlope(2) - 1.0 / rho_R * rStateSlope(1);
            p_t = -u_R * rStateSlope(1) - rho_R * c_R * c_R * rStateSlope(2);
        }
        else
        {
            double T_LS_Slope_L = lStateSlope(1) / (gamma - 1.0) / rho_L - c_L * c_L * lStateSlope(0) / (gamma - 1.0) / rho_L;
            double temp1_d_L = (1.0 + mu_2) / (1.0 + 2.0 * mu_2) * pow(c_starL / c_L, 0.5 / mu_2) + mu_2 / (1.0 + 2.0 * mu_2) * pow(c_starL / c_L, (1.0 + mu_2) / mu_2);
            double temp2_d_L = c_L * pow(c_starL / c_L, 0.5 / mu_2) * (lStateSlope(2) + gamma * lStateSlope(1) / (gamma - 1.0) / rho_L / c_L - c_L * lStateSlope(0) / (gamma - 1.0) / rho_L);
            double a_L = 1.0, b_L = 1.0 / (rho_starL * c_starL);
            double d_L = temp1_d_L * T_LS_Slope_L - temp2_d_L;

            double T_RS_Slope_R = rStateSlope(1) / (gamma - 1.0) / rho_R - c_R * c_R * rStateSlope(0) / (gamma - 1.0) / rho_R;
            double temp1_d_R = (1.0 + mu_2) / (1.0 + 2.0 * mu_2) * pow(c_starR / c_R, 0.5 / mu_2) + mu_2 / (1.0 + 2.0 * mu_2) * pow(c_starR / c_R, (1.0 + mu_2) / mu_2);
            double temp2_d_R = c_R * pow(c_starR / c_R, 0.5 / mu_2) * (rStateSlope(2) - gamma * rStateSlope(1) / (gamma - 1.0) / rho_R / c_R + c_R * rStateSlope(0) / (gamma - 1.0) / rho_R);
            double a_R = 1.0, b_R = -1.0 / (rho_starR * c_starR);
            double d_R = temp1_d_R * T_RS_Slope_R + temp2_d_R;

            if (u_star - c_starL >= 0)
            {
                double c_0_0 = mu_2 * psi_L;
                double u_0_0 = c_0_0;
                double rho_0_0 = pow(c_0_0 * c_0_0 * pow(rho_L, gamma) / gamma / p_L, 1.0 / (gamma - 1.0));
                double temp1_d_L0 = (1.0 + mu_2) / (1.0 + 2.0 * mu_2) * pow(c_0_0 / c_L, 0.5 / mu_2) + mu_2 / (1.0 + 2.0 * mu_2) * pow(c_0_0 / c_L, (1.0 + mu_2) / mu_2);
                double temp2_d_L0 = c_L * pow(c_0_0 / c_L, 0.5 / mu_2) * (lStateSlope(2) + gamma * lStateSlope(1) / (gamma - 1.0) / rho_L / c_L - c_L * lStateSlope(0) / (gamma - 1.0) / rho_L);
                u_t = temp1_d_L0 * T_LS_Slope_L - temp2_d_L0;
                p_t = rho_0_0 * u_0_0 * u_t;
                rho_t = 1.0 / c_0_0 / c_0_0 * (p_t + (gamma - 1.0) * rho_0_0 * u_0_0 * pow(c_0_0 / c_L, (1.0 + mu_2) / mu_2) * T_LS_Slope_L);
            }
            else if (u_star + c_starR <= 0)
            {
                double c_0_0 = -mu_2 * phi_R;
                double u_0_0 = -c_0_0;
                double rho_0_0 = pow(c_0_0 * c_0_0 * pow(rho_R, gamma) / gamma / p_R, 1.0 / (gamma - 1.0));
                double temp1_d_R0 = (1.0 + mu_2) / (1.0 + 2.0 * mu_2) * pow(c_0_0 / c_R, 0.5 / mu_2) + mu_2 / (1.0 + 2.0 * mu_2) * pow(c_0_0 / c_R, (1.0 + mu_2) / mu_2);
                double temp2_d_R0 = c_R * pow(c_0_0 / c_R, 0.5 / mu_2) * (rStateSlope(2) - gamma * rStateSlope(1) / (gamma - 1.0) / rho_R / c_R + c_R * rStateSlope(0) / (gamma - 1.0) / rho_R);
                u_t = temp1_d_R0 * T_RS_Slope_R + temp2_d_R0;
                p_t = -rho_0_0 * u_0_0 * u_t;
                rho_t = 1.0 / c_0_0 / c_0_0 * (p_t + (gamma - 1.0) * rho_0_0 * u_0_0 * pow(c_0_0 / c_R, (1.0 + mu_2) / mu_2) * T_RS_Slope_R);
            }
            else
            {
                mat2d A;
                A << a_L, b_L, a_R, b_R;
                vec2d x;
                x << d_L, d_R;
                vec2d u_p_material_deri = A.fullPivLu().solve(x);
                double Du_Dt = u_p_material_deri(0);
                double Dp_Dt = u_p_material_deri(1);
                if (u_star >= 0)
                {
                    rho_t = 1.0 / c_starL / c_starL * (p_t + (gamma - 1.0) * rho_starL * u_star * pow(c_starL / c_L, (1.0 + mu_2) / mu_2) * T_LS_Slope_L);
                    u_t = Du_Dt + u_star * Dp_Dt / rho_starL / c_starL / c_starL;
                    p_t = Dp_Dt + rho_starL * u_star * Du_Dt;
                }
                else
                {
                    rho_t = 1.0 / c_starR / c_starR * (p_t + (gamma - 1.0) * rho_starR * u_star * pow(c_starR / c_R, (1.0 + mu_2) / mu_2) * T_RS_Slope_R);
                    u_t = Du_Dt + u_star * Dp_Dt / rho_starR / c_starR / c_starR;
                    p_t = Dp_Dt + rho_starR * u_star * Du_Dt;
                }
            }
        }
    }

    else
    {
        double ARPSol_sigmaL = 0.0;
        if (abs(rho_L - rho_starL) >= 10e-10)
        {
            ARPSol_sigmaL = (rho_L * u_L - rho_starL * u_star) / (rho_L - rho_starL);
        }
        else
        {
            ARPSol_sigmaL = u_L - c_L;
        }
        double ARPSol_sigmaR = 0.0;
        if (abs(rho_R - rho_starR) >= 10e-10)
        {
            ARPSol_sigmaR = (rho_R * u_R - rho_starR * u_star) / (rho_R - rho_starR);
        }
        else
        {
            ARPSol_sigmaR = u_R + c_R;
        }

        if (ARPSol_sigmaR <= 0.0)
        {
            rho_t = -u_R * rStateSlope(0) - rho_R * rStateSlope(2);
            u_t = -u_R * rStateSlope(2) - 1.0 / rho_R * rStateSlope(1);
            p_t = -u_R * rStateSlope(1) - rho_R * c_R * c_R * rStateSlope(2);
        }
        else if (ARPSol_sigmaL >= 0)
        {
            rho_t = -u_L * lStateSlope(0) - rho_L * lStateSlope(2);
            u_t = -u_L * lStateSlope(2) - 1.0 / rho_L * lStateSlope(1);
            p_t = -u_L * lStateSlope(1) - rho_L * c_L * c_L * lStateSlope(2);
        }
        else
        {
            double H_1_R = 0.5 * sqrt((1.0 - mu_2) / rho_R / (p_star + mu_2 * p_R)) * (p_star + (1.0 + 2.0 * mu_2) * p_R) / (p_star + mu_2 * p_R);
            double H_2_R = -0.5 * sqrt((1.0 - mu_2) / rho_R / (p_star + mu_2 * p_R)) * ((2.0 + mu_2) * p_star + mu_2 * p_R) / (p_star + mu_2 * p_R);
            double H_3_R = -(p_star - p_R) / (2.0 * rho_R) * sqrt((1.0 - mu_2) / rho_R / (p_star + mu_2 * p_R));
            double L_p_R = -1.0 / rho_R + (ARPSol_sigmaR - u_R) * H_2_R;
            double L_u_R = ARPSol_sigmaR - u_R - rho_R * c_R * c_R * H_2_R - rho_R * H_3_R;
            double L_rho_R = (ARPSol_sigmaR - u_R) * H_3_R;
            double a_R = 1.0 + rho_starR * (ARPSol_sigmaR - u_star) * H_1_R;
            double b_R = (u_star - ARPSol_sigmaR) / (rho_starR * c_starR * c_starR) - H_1_R;
            double d_R = L_p_R * rStateSlope(1) + L_u_R * rStateSlope(2) + L_rho_R * rStateSlope(0);

            double H_1_L = 0.5 * sqrt((1.0 - mu_2) / rho_L / (p_star + mu_2 * p_L)) * (p_star + (1.0 + 2.0 * mu_2) * p_L) / (p_star + mu_2 * p_L);
            double H_2_L = -0.5 * sqrt((1.0 - mu_2) / rho_L / (p_star + mu_2 * p_L)) * ((2.0 + mu_2) * p_star + mu_2 * p_L) / (p_star + mu_2 * p_L);
            double H_3_L = -(p_star - p_L) / (2.0 * rho_L) * sqrt((1.0 - mu_2) / rho_L / (p_star + mu_2 * p_L));
            double L_p_L = -1.0 / rho_L - (ARPSol_sigmaL - u_L) * H_2_L;
            double L_u_L = ARPSol_sigmaL - u_L + rho_L * c_L * c_L * H_2_L + rho_L * H_3_L;
            double L_rho_L = -(ARPSol_sigmaL - u_L) * H_3_L;
            double a_L = 1.0 - rho_starL * (ARPSol_sigmaL - u_star) * H_1_L;
            double b_L = (u_star - ARPSol_sigmaL) / (rho_starL * c_starL * c_starL) + H_1_L;
            double d_L = L_p_L * lStateSlope(1) + L_u_L * lStateSlope(2) + L_rho_L * lStateSlope(0);

            mat2d A;
            A << a_L, b_L, a_R, b_R;
            vec2d x;
            x << d_L, d_R;
            vec2d u_p_material_deri = A.fullPivLu().solve(x);
            double Du_Dt = u_p_material_deri(0);
            double Dp_Dt = u_p_material_deri(1);
            if (u_star >= 0)
            {
                double H_1_ = rho_L * (1.0 - mu_2 * mu_2) * p_L / (p_L + mu_2 * p_star) / (p_L + mu_2 * p_star);
                double H_2_ = rho_L * (mu_2 * mu_2 - 1.0) * p_star / (p_L + mu_2 * p_star) / (p_L + mu_2 * p_star);
                double H_3_ = (p_star + mu_2 * p_L) / (p_L + mu_2 * p_star);
                double g_rho_L = u_star - ARPSol_sigmaL;
                double g_p_L = ARPSol_sigmaL / c_starL / c_starL - u_star * H_1_;
                double g_u_L = u_star * rho_starL * (ARPSol_sigmaL - u_star) * H_1_;
                double f_L = (ARPSol_sigmaL - u_L) * H_2_ * lStateSlope(1) + (ARPSol_sigmaL - u_L) * H_3_ * lStateSlope(0) - rho_L * (H_2_ * c_L * c_L + H_3_) * lStateSlope(2);
                rho_t = (u_star * f_L - g_p_L * Dp_Dt - g_u_L * Du_Dt) / g_rho_L;
                u_t = Du_Dt + u_star * Dp_Dt / rho_starL / c_starL / c_starL;
                p_t = Dp_Dt + rho_starL * u_star * Du_Dt;
            }
            else
            {
                double H_1_ = rho_R * (1.0 - mu_2 * mu_2) * p_R / (p_R + mu_2 * p_star) / (p_R + mu_2 * p_star);
                double H_2_ = rho_R * (mu_2 * mu_2 - 1.0) * p_star / (p_R + mu_2 * p_star) / (p_R + mu_2 * p_star);
                double H_3_ = (p_star + mu_2 * p_R) / (p_R + mu_2 * p_star);
                double g_rho_R = u_star - ARPSol_sigmaR;
                double g_p_R = ARPSol_sigmaR / c_starR / c_starR - u_star * H_1_;
                double g_u_R = u_star * rho_starR * (ARPSol_sigmaR - u_star) * H_1_;
                double f_R = (ARPSol_sigmaR - u_R) * H_2_ * rStateSlope(1) + (ARPSol_sigmaR - u_R) * H_3_ * rStateSlope(0) - rho_R * (H_2_ * c_R * c_R + H_3_) * rStateSlope(2);
                rho_t = (u_star * f_R - g_p_R * Dp_Dt - g_u_R * Du_Dt) / g_rho_R;
                u_t = Du_Dt + u_star * Dp_Dt / rho_starR / c_starR / c_starR;
                p_t = Dp_Dt + rho_starR * u_star * Du_Dt;
            }
        }
    }

    ret.ARPSol = ARPSol;
    ret.timeDerivative(0) = rho_t;
    ret.timeDerivative(1) = p_t;
    ret.timeDerivative(2) = u_t;
    return ret;
}