#include "solvers.h"

RPSolver::RPSolver(const vec3d& lState, const vec3d& rState)
{
	rho_L = lState(0), rho_R = rState(0), p_L = lState(1), p_R = rState(1), u_L = lState(2), u_R = rState(2);
	initVar();
}

void RPSolver::initVar()
{
	mu_2 = (gamma - 1.0) / (gamma + 1.0);
	c_L = sqrt(gamma * p_L / rho_L), c_R = sqrt(gamma * p_R / rho_R);
	psi_L = u_L + 2.0 * c_L / (gamma - 1.0);
	phi_R = u_R - 2.0 * c_R / (gamma - 1.0);
	A_L = 2.0 / (gamma + 1.0) / rho_L, B_L = mu_2 * p_L, A_R = 2.0 / (gamma + 1.0) / rho_R, B_R = mu_2 * p_R;
}

void RPSolver::setGamma(double gamma)
{
	this->gamma = 1.0 < gamma && gamma <= 5.0 / 3.0 ? gamma : this->gamma;
	initVar();
}

void RPSolver::setTol(double tol)
{
	this->tol = 0.0 < tol && tol <= 1.0e-4 ? tol : this->tol;
}

double RPSolver::NewtonMethod(const scalarFun1d& f, const scalarFun1d& df, double start)
{
	double save = start < 0.0 ? 1.0e-6 : start, p = save - f(save) / df(save), cha = 2.0 * abs(p - save) / abs(p + save);
	int i = 0;
	while (cha >= tol && i <= 1000)
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

void RPSolver::solve()
{
	auto f_L = [&](double p) -> double
		{
			if (p > p_L)
			{
				return (p - p_L) * sqrt(A_L / (p + B_L));
			}
			else
			{
				return 2.0 * c_L / (gamma - 1.0) * (pow(p / p_L, (gamma - 1.0) / 2.0 / gamma) - 1.0);
			}
		};

	auto f_R = [&](double p) -> double
		{
			if (p > p_R)
			{
				return (p - p_R) * sqrt(A_R / (p + B_R));
			}
			else
			{
				return 2.0 * c_R / (gamma - 1.0) * (pow(p / p_R, (gamma - 1.0) / 2.0 / gamma) - 1.0);
			}
		};

	auto df_L = [&](double p) -> double
		{
			if (p > p_L)
			{
				return sqrt(A_L / (B_L + p)) * (1.0 - (p - p_L) / 2.0 / (B_L + p));
			}
			else
			{
				return 1.0 / rho_L / c_L * pow(p / p_L, -(gamma + 1.0) / 2.0 / gamma);
			}
		};

	auto df_R = [&](double p) -> double
		{
			if (p > p_R)
			{
				return sqrt(A_R / (B_R + p)) * (1.0 - (p - p_R) / 2.0 / (B_R + p));
			}
			else
			{
				return 1.0 / rho_R / c_R * pow(p / p_R, -(gamma + 1.0) / 2.0 / gamma);
			}
		};

	auto df = [&](double p) -> double
		{
			return df_L(p) + df_R(p);
		};

	auto f = [&](double p) -> double
		{
			return f_L(p) + f_R(p) + u_R - u_L;
		};

	p_star_ = NewtonMethod(f, df, 0.5 * (p_L + p_R));
	u_star_ = 0.5 * (u_L + u_R + f_R(p_star_) - f_L(p_star_));

	if (p_star_ > p_L)
	{
		wave1_ = Shock;
		rho_starL_ = rho_L * (mu_2 * p_L + p_star_) / (mu_2 * p_star_ + p_L);
	}
	else
	{
		wave1_ = CRW;
		rho_starL_ = rho_L * pow(p_star_ / p_L, 1.0 / gamma);
	}

	if (p_star_ > p_R)
	{
		wave3_ = Shock;
		rho_starR_ = rho_R * (mu_2 * p_R + p_star_) / (mu_2 * p_star_ + p_R);
	}
	else
	{
		wave3_ = CRW;
		rho_starR_ = rho_R * pow(p_star_ / p_R, 1.0 / gamma);
	}
	c_starL = sqrt(gamma * p_star_ / rho_starL_);
	c_starR = sqrt(gamma * p_star_ / rho_starR_);
}

vec3d RPSolver::operator()(double x, double t)
{
	if (t <= 0.0)
	{
		return{ -1.0,-1.0,0.0 };
	}

	double beta = x / t;

	if (wave1_ == CRW)
	{
		if (beta <= u_L - c_L)
		{
			return { rho_L, p_L, u_L };
		}
		else if (beta <= u_star_ - c_starL)
		{
			double rho = pow(mu_2 * mu_2 * (psi_L - beta) * (psi_L - beta) * pow(rho_L, gamma) / gamma / p_L, 1.0 / (gamma - 1.0));
			return { rho, rho * mu_2 * mu_2 * (psi_L - beta) * (psi_L - beta) / gamma, mu_2 * psi_L + (1.0 - mu_2) * beta };
		}
		else if (beta <= u_star_)
		{
			return { rho_starL_, p_star_, u_star_ };
		}
	}
	else
	{
		double sigma_L = (rho_starL_ * u_star_ - rho_L * u_L) / (rho_starL_ - rho_L);
		if (beta <= sigma_L)
		{
			return { rho_L, p_L, u_L };
		}
		else if (beta <= u_star_)
		{
			return { rho_starL_, p_star_, u_star_ };
		}
	}

	if (wave3_ == CRW)
	{
		if (beta >= u_R + c_R)
		{
			return { rho_R, p_R, u_R };
		}
		else if (beta >= u_star_ + c_starR)
		{
			double rho = pow(mu_2 * mu_2 * (phi_R - beta) * (phi_R - beta) * pow(rho_R, gamma) / gamma / p_R, 1.0 / (gamma - 1.0));
			return { rho, rho * mu_2 * mu_2 * (phi_R - beta) * (phi_R - beta), mu_2 * phi_R + (1.0 - mu_2) * beta };
		}
		else if (beta >= u_star_)
		{
			return { rho_starR_, p_star_, u_star_ };
		}
	}
	else
	{
		double sigma_R = (rho_starR_ * u_star_ - rho_R * u_R) / (rho_starR_ - rho_R);
		if (beta >= sigma_R)
		{
			return { rho_R, p_R, u_R };
		}
		else if (beta >= u_star_)
		{
			return { rho_starR_, p_star_, u_star_ };
		}
	}

	return { -1.0, -1.0, 0.0 };
}

double RPSolver::rho_starL()
{
	return rho_starL_;
}

double RPSolver::rho_starR()
{
	return rho_starR_;
}

double RPSolver::p_star()
{
	return p_star_;
}

double RPSolver::u_star()
{
	return u_star_;
}

waveType RPSolver::wave1()
{
	return wave1_;
}

waveType RPSolver::wave3()
{
	return wave3_;
}

GRPSolver::GRPSolver(const vec3d& lState, const vec3d& rState, const vec3d& lStateSlope, const vec3d& rStateSlope) : RPSolver(lState, rState), lStateSlope(lStateSlope), rStateSlope(rStateSlope)
{
}

void GRPSolver::solve()
{
	this->RPSolver::solve();

	if (abs(rho_L - rho_R) + abs(p_L - p_R) + abs(u_L - u_R) <= tol)
	{
		u_t = -0.5 * ((u_star_ + c_starL) * (lStateSlope(2) + lStateSlope(1) / rho_starL_ / c_starL) + (u_star_ - c_starL) * (rStateSlope(2) - rStateSlope(1) / rho_starL_ / c_starL));
		p_t = -0.5 * rho_starL_ * c_starL * ((u_star_ + c_starL) * (lStateSlope(2) + lStateSlope(1) / rho_starL_ / c_starL) - (u_star_ - c_starL) * (rStateSlope(2) - rStateSlope(1) / rho_starL_ / c_starL));
		if (u_star_ >= 0)
		{
			rho_t = (p_t + u_star_ * (lStateSlope(1) - c_starL * c_starL * lStateSlope(0))) / c_starL / c_starL;
		}
		else
		{
			rho_t = (p_t + u_star_ * (rStateSlope(1) - c_starL * c_starL * rStateSlope(0))) / c_starL / c_starL;
		}
		return;
	}

	if (wave1_ == CRW && wave3_ == Shock)
	{
		double ARPSol_sigma = 0.0;
		if (abs(rho_R - rho_starR_) >= 10e-10)
		{
			ARPSol_sigma = (rho_R * u_R - rho_starR_ * u_star_) / (rho_R - rho_starR_);
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
			double a_L = 1.0, b_L = 1.0 / (rho_starL_ * c_starL);
			double d_L = temp1_d_L * T_LS_Slope_L - temp2_d_L;

			double H_1 = 0.5 * sqrt((1.0 - mu_2) / rho_R / (p_star_ + mu_2 * p_R)) * (p_star_ + (1.0 + 2.0 * mu_2) * p_R) / (p_star_ + mu_2 * p_R);
			double H_2 = -0.5 * sqrt((1.0 - mu_2) / rho_R / (p_star_ + mu_2 * p_R)) * ((2.0 + mu_2) * p_star_ + mu_2 * p_R) / (p_star_ + mu_2 * p_R);
			double H_3 = -(p_star_ - p_R) / (2.0 * rho_R) * sqrt((1.0 - mu_2) / rho_R / (p_star_ + mu_2 * p_R));
			double L_p_R = -1.0 / rho_R + (ARPSol_sigma - u_R) * H_2;
			double L_u_R = ARPSol_sigma - u_R - rho_R * c_R * c_R * H_2 - rho_R * H_3;
			double L_rho_R = (ARPSol_sigma - u_R) * H_3;
			double a_R = 1.0 + rho_starR_ * (ARPSol_sigma - u_star_) * H_1;
			double b_R = (u_star_ - ARPSol_sigma) / (rho_starR_ * c_starR * c_starR) - H_1;
			double d_R = L_p_R * rStateSlope(1) + L_u_R * rStateSlope(2) + L_rho_R * rStateSlope(0);

			if (u_star_ - c_starL >= 0)
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
			else
			{
				mat2d A;
				A << a_L, b_L, a_R, b_R;
				vec2d x;
				x << d_L, d_R;
				vec2d u_p_material_deri = A.fullPivLu().solve(x);
				double Du_Dt = u_p_material_deri(0);
				double Dp_Dt = u_p_material_deri(1);
				if (u_star_ >= 0)
				{
					rho_t = 1.0 / c_starL / c_starL * (p_t + (gamma - 1.0) * rho_starL_ * u_star_ * pow(c_starL / c_L, (1.0 + mu_2) / mu_2) * T_LS_Slope_L);
					u_t = Du_Dt + u_star_ * Dp_Dt / rho_starL_ / c_starL / c_starL;
					p_t = Dp_Dt + rho_starL_ * u_star_ * Du_Dt;
				}
				else
				{
					double H_1_ = rho_R * (1.0 - mu_2 * mu_2) * p_R / (p_R + mu_2 * p_star_) / (p_R + mu_2 * p_star_);
					double H_2_ = rho_R * (mu_2 * mu_2 - 1.0) * p_star_ / (p_R + mu_2 * p_star_) / (p_R + mu_2 * p_star_);
					double H_3_ = (p_star_ + mu_2 * p_R) / (p_R + mu_2 * p_star_);
					double g_rho_R = u_star_ - ARPSol_sigma;
					double g_p_R = ARPSol_sigma / c_starR / c_starR - u_star_ * H_1_;
					double g_u_R = u_star_ * rho_starR_ * (ARPSol_sigma - u_star_) * H_1_;
					double f_R = (ARPSol_sigma - u_R) * H_2_ * rStateSlope(1) + (ARPSol_sigma - u_R) * H_3_ * rStateSlope(0) - rho_R * (H_2_ * c_R * c_R + H_3_) * rStateSlope(2);
					rho_t = (u_star_ * f_R - g_p_R * Dp_Dt - g_u_R * Du_Dt) / g_rho_R;
					u_t = Du_Dt + u_star_ * Dp_Dt / rho_starR_ / c_starR / c_starR;
					p_t = Dp_Dt + rho_starR_ * u_star_ * Du_Dt;
				}
			}
		}
	}

	else if (wave1_ == Shock && wave3_ == CRW)
	{
		double ARPSol_sigma = 0.0;
		if (abs(rho_L - rho_starL_) >= 10e-10)
		{
			ARPSol_sigma = (rho_L * u_L - rho_starL_ * u_star_) / (rho_L - rho_starL_);
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
			double a_R = 1.0, b_R = -1.0 / (rho_starR_ * c_starR);
			double d_R = temp1_d_R * T_RS_Slope_R + temp2_d_R;

			double H_1 = 0.5 * sqrt((1.0 - mu_2) / rho_L / (p_star_ + mu_2 * p_L)) * (p_star_ + (1.0 + 2.0 * mu_2) * p_L) / (p_star_ + mu_2 * p_L);
			double H_2 = -0.5 * sqrt((1.0 - mu_2) / rho_L / (p_star_ + mu_2 * p_L)) * ((2.0 + mu_2) * p_star_ + mu_2 * p_L) / (p_star_ + mu_2 * p_L);
			double H_3 = -(p_star_ - p_L) / (2.0 * rho_L) * sqrt((1.0 - mu_2) / rho_L / (p_star_ + mu_2 * p_L));
			double L_p_L = -1.0 / rho_L - (ARPSol_sigma - u_L) * H_2;
			double L_u_L = ARPSol_sigma - u_L + rho_L * c_L * c_L * H_2 + rho_L * H_3;
			double L_rho_L = -(ARPSol_sigma - u_L) * H_3;
			double a_L = 1.0 - rho_starL_ * (ARPSol_sigma - u_star_) * H_1;
			double b_L = (u_star_ - ARPSol_sigma) / (rho_starL_ * c_starL * c_starL) + H_1;
			double d_L = L_p_L * lStateSlope(1) + L_u_L * lStateSlope(2) + L_rho_L * lStateSlope(0);

			if (u_star_ + c_starR <= 0)
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
				if (u_star_ <= 0)
				{
					rho_t = 1.0 / c_starR / c_starR * (p_t + (gamma - 1.0) * rho_starR_ * u_star_ * pow(c_starR / c_R, (1.0 + mu_2) / mu_2) * T_RS_Slope_R);
					u_t = Du_Dt + u_star_ * Dp_Dt / rho_starR_ / c_starR / c_starR;
					p_t = Dp_Dt + rho_starR_ * u_star_ * Du_Dt;
				}
				else
				{
					double H_1_ = rho_L * (1.0 - mu_2 * mu_2) * p_L / (p_L + mu_2 * p_star_) / (p_L + mu_2 * p_star_);
					double H_2_ = rho_L * (mu_2 * mu_2 - 1.0) * p_star_ / (p_L + mu_2 * p_star_) / (p_L + mu_2 * p_star_);
					double H_3_ = (p_star_ + mu_2 * p_L) / (p_L + mu_2 * p_star_);
					double g_rho_L = u_star_ - ARPSol_sigma;
					double g_p_L = ARPSol_sigma / c_starL / c_starL - u_star_ * H_1_;
					double g_u_L = u_star_ * rho_starL_ * (ARPSol_sigma - u_star_) * H_1_;
					double f_L = (ARPSol_sigma - u_L) * H_2_ * lStateSlope(1) + (ARPSol_sigma - u_L) * H_3_ * lStateSlope(0) - rho_L * (H_2_ * c_L * c_L + H_3_) * lStateSlope(2);
					rho_t = (u_star_ * f_L - g_p_L * Dp_Dt - g_u_L * Du_Dt) / g_rho_L;
					u_t = Du_Dt + u_star_ * Dp_Dt / rho_starL_ / c_starL / c_starL;
					p_t = Dp_Dt + rho_starL_ * u_star_ * Du_Dt;
				}
			}
		}
	}

	else if (wave1_ == CRW && wave3_ == CRW)
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
			double a_L = 1.0, b_L = 1.0 / (rho_starL_ * c_starL);
			double d_L = temp1_d_L * T_LS_Slope_L - temp2_d_L;

			double T_RS_Slope_R = rStateSlope(1) / (gamma - 1.0) / rho_R - c_R * c_R * rStateSlope(0) / (gamma - 1.0) / rho_R;
			double temp1_d_R = (1.0 + mu_2) / (1.0 + 2.0 * mu_2) * pow(c_starR / c_R, 0.5 / mu_2) + mu_2 / (1.0 + 2.0 * mu_2) * pow(c_starR / c_R, (1.0 + mu_2) / mu_2);
			double temp2_d_R = c_R * pow(c_starR / c_R, 0.5 / mu_2) * (rStateSlope(2) - gamma * rStateSlope(1) / (gamma - 1.0) / rho_R / c_R + c_R * rStateSlope(0) / (gamma - 1.0) / rho_R);
			double a_R = 1.0, b_R = -1.0 / (rho_starR_ * c_starR);
			double d_R = temp1_d_R * T_RS_Slope_R + temp2_d_R;

			if (u_star_ - c_starL >= 0)
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
			else if (u_star_ + c_starR <= 0)
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
				if (u_star_ >= 0)
				{
					rho_t = 1.0 / c_starL / c_starL * (p_t + (gamma - 1.0) * rho_starL_ * u_star_ * pow(c_starL / c_L, (1.0 + mu_2) / mu_2) * T_LS_Slope_L);
					u_t = Du_Dt + u_star_ * Dp_Dt / rho_starL_ / c_starL / c_starL;
					p_t = Dp_Dt + rho_starL_ * u_star_ * Du_Dt;
				}
				else
				{
					rho_t = 1.0 / c_starR / c_starR * (p_t + (gamma - 1.0) * rho_starR_ * u_star_ * pow(c_starR / c_R, (1.0 + mu_2) / mu_2) * T_RS_Slope_R);
					u_t = Du_Dt + u_star_ * Dp_Dt / rho_starR_ / c_starR / c_starR;
					p_t = Dp_Dt + rho_starR_ * u_star_ * Du_Dt;
				}
			}
		}
	}

	else
	{
		double ARPSol_sigmaL = 0.0;
		if (abs(rho_L - rho_starL_) >= 10e-10)
		{
			ARPSol_sigmaL = (rho_L * u_L - rho_starL_ * u_star_) / (rho_L - rho_starL_);
		}
		else
		{
			ARPSol_sigmaL = u_L - c_L;
		}
		double ARPSol_sigmaR = 0.0;
		if (abs(rho_R - rho_starR_) >= 10e-10)
		{
			ARPSol_sigmaR = (rho_R * u_R - rho_starR_ * u_star_) / (rho_R - rho_starR_);
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
			double H_1_R = 0.5 * sqrt((1.0 - mu_2) / rho_R / (p_star_ + mu_2 * p_R)) * (p_star_ + (1.0 + 2.0 * mu_2) * p_R) / (p_star_ + mu_2 * p_R);
			double H_2_R = -0.5 * sqrt((1.0 - mu_2) / rho_R / (p_star_ + mu_2 * p_R)) * ((2.0 + mu_2) * p_star_ + mu_2 * p_R) / (p_star_ + mu_2 * p_R);
			double H_3_R = -(p_star_ - p_R) / (2.0 * rho_R) * sqrt((1.0 - mu_2) / rho_R / (p_star_ + mu_2 * p_R));
			double L_p_R = -1.0 / rho_R + (ARPSol_sigmaR - u_R) * H_2_R;
			double L_u_R = ARPSol_sigmaR - u_R - rho_R * c_R * c_R * H_2_R - rho_R * H_3_R;
			double L_rho_R = (ARPSol_sigmaR - u_R) * H_3_R;
			double a_R = 1.0 + rho_starR_ * (ARPSol_sigmaR - u_star_) * H_1_R;
			double b_R = (u_star_ - ARPSol_sigmaR) / (rho_starR_ * c_starR * c_starR) - H_1_R;
			double d_R = L_p_R * rStateSlope(1) + L_u_R * rStateSlope(2) + L_rho_R * rStateSlope(0);

			double H_1_L = 0.5 * sqrt((1.0 - mu_2) / rho_L / (p_star_ + mu_2 * p_L)) * (p_star_ + (1.0 + 2.0 * mu_2) * p_L) / (p_star_ + mu_2 * p_L);
			double H_2_L = -0.5 * sqrt((1.0 - mu_2) / rho_L / (p_star_ + mu_2 * p_L)) * ((2.0 + mu_2) * p_star_ + mu_2 * p_L) / (p_star_ + mu_2 * p_L);
			double H_3_L = -(p_star_ - p_L) / (2.0 * rho_L) * sqrt((1.0 - mu_2) / rho_L / (p_star_ + mu_2 * p_L));
			double L_p_L = -1.0 / rho_L - (ARPSol_sigmaL - u_L) * H_2_L;
			double L_u_L = ARPSol_sigmaL - u_L + rho_L * c_L * c_L * H_2_L + rho_L * H_3_L;
			double L_rho_L = -(ARPSol_sigmaL - u_L) * H_3_L;
			double a_L = 1.0 - rho_starL_ * (ARPSol_sigmaL - u_star_) * H_1_L;
			double b_L = (u_star_ - ARPSol_sigmaL) / (rho_starL_ * c_starL * c_starL) + H_1_L;
			double d_L = L_p_L * lStateSlope(1) + L_u_L * lStateSlope(2) + L_rho_L * lStateSlope(0);

			mat2d A;
			A << a_L, b_L, a_R, b_R;
			vec2d x;
			x << d_L, d_R;
			vec2d u_p_material_deri = A.fullPivLu().solve(x);
			double Du_Dt = u_p_material_deri(0);
			double Dp_Dt = u_p_material_deri(1);
			if (u_star_ >= 0)
			{
				double H_1_ = rho_L * (1.0 - mu_2 * mu_2) * p_L / (p_L + mu_2 * p_star_) / (p_L + mu_2 * p_star_);
				double H_2_ = rho_L * (mu_2 * mu_2 - 1.0) * p_star_ / (p_L + mu_2 * p_star_) / (p_L + mu_2 * p_star_);
				double H_3_ = (p_star_ + mu_2 * p_L) / (p_L + mu_2 * p_star_);
				double g_rho_L = u_star_ - ARPSol_sigmaL;
				double g_p_L = ARPSol_sigmaL / c_starL / c_starL - u_star_ * H_1_;
				double g_u_L = u_star_ * rho_starL_ * (ARPSol_sigmaL - u_star_) * H_1_;
				double f_L = (ARPSol_sigmaL - u_L) * H_2_ * lStateSlope(1) + (ARPSol_sigmaL - u_L) * H_3_ * lStateSlope(0) - rho_L * (H_2_ * c_L * c_L + H_3_) * lStateSlope(2);
				rho_t = (u_star_ * f_L - g_p_L * Dp_Dt - g_u_L * Du_Dt) / g_rho_L;
				u_t = Du_Dt + u_star_ * Dp_Dt / rho_starL_ / c_starL / c_starL;
				p_t = Dp_Dt + rho_starL_ * u_star_ * Du_Dt;
			}
			else
			{
				double H_1_ = rho_R * (1.0 - mu_2 * mu_2) * p_R / (p_R + mu_2 * p_star_) / (p_R + mu_2 * p_star_);
				double H_2_ = rho_R * (mu_2 * mu_2 - 1.0) * p_star_ / (p_R + mu_2 * p_star_) / (p_R + mu_2 * p_star_);
				double H_3_ = (p_star_ + mu_2 * p_R) / (p_R + mu_2 * p_star_);
				double g_rho_R = u_star_ - ARPSol_sigmaR;
				double g_p_R = ARPSol_sigmaR / c_starR / c_starR - u_star_ * H_1_;
				double g_u_R = u_star_ * rho_starR_ * (ARPSol_sigmaR - u_star_) * H_1_;
				double f_R = (ARPSol_sigmaR - u_R) * H_2_ * rStateSlope(1) + (ARPSol_sigmaR - u_R) * H_3_ * rStateSlope(0) - rho_R * (H_2_ * c_R * c_R + H_3_) * rStateSlope(2);
				rho_t = (u_star_ * f_R - g_p_R * Dp_Dt - g_u_R * Du_Dt) / g_rho_R;
				u_t = Du_Dt + u_star_ * Dp_Dt / rho_starR_ / c_starR / c_starR;
				p_t = Dp_Dt + rho_starR_ * u_star_ * Du_Dt;
			}
		}
	}
}

vec3d GRPSolver::timeDerivative()
{
	return { rho_t,p_t,u_t };
}
