#include "basicFuns.h"

double maxPropagtingSpeed(const vector<vec3d>& U, double gamma)
{
	double cha_speed_save = 0.0, cha_speed = 0.0;
	for (auto& u : U)
	{
		double velo = u(2), rho = u(0), p = u(1), c = sqrt(gamma * p / rho);
		cha_speed = std::abs(velo) + c;
		if (cha_speed >= cha_speed_save)
		{
			cha_speed_save = cha_speed;
		}
	}
	return cha_speed_save;
}

vec3d rhoPUToConserVar(const vec3d& rhoPUV, double gamma)
{
	return { rhoPUV(0),
		     rhoPUV(0) * rhoPUV(2),
			 rhoPUV(1) / (gamma - 1.0) + rhoPUV(0) * rhoPUV(2) * rhoPUV(2) / 2.0};
}



void Lx_relGRP::operator()(vector<vec3d>& U, double k, double lB, double rB, double currentTime)
{
	int N = U.size();
	double h = (rB - lB) / double(N);
	U_ = vector<vec3d>(N + 4);

	for (auto i = 0; i < N; i++)
	{
		U_[i + 2] = rhoPUToConserVar(U[i], gamma);
	}

	addGhostCells(U, currentTime);
	U_[0] = rhoPUToConserVar(U[0], gamma);
	U_[1] = rhoPUToConserVar(U[1], gamma);
	U_[N + 2] = rhoPUToConserVar(U[N + 2], gamma);
	U_[N + 3] = rhoPUToConserVar(U[N + 3], gamma);


	auto minmod = [&](const vec3d& v1, const vec3d& v2) -> vec3d
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
					bool modified = false;
				bool modified = false;
				{
					ret(i) = v2(i);
				}
			}
			return ret;
		};
	auto minmod3 = [&](const vec3d& v1, const vec3d& v2, const vec3d& v3) -> vec3d
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
		};

	vector<vec3d> U_Slope(N + 4);
	for (auto j = 1; j < N + 3; j++)
	{
		U_Slope[j] = minmod3((U[j + 1] - U[j - 1]) / 2.0 / h, alpha * (U[j + 1] - U[j]) / h, alpha * (U[j] - U[j - 1]) / h);
	}

	auto timeDeris = vector<vec3d>(N + 3);
	auto GodunovSols_Rel = vector<vec3d>(N + 3);

	for (auto i = 1; i < N + 2; i++) // i+1/2
	{
		vec3d U_n_i_L = U[i] + U_Slope[i] * h / 2.0;
		vec3d U_n_i_R = U[i + 1] - U_Slope[i + 1] * h / 2.0;
		auto localGRPSolver_Rel = GRPSolver(U_n_i_L, U_n_i_R, U_Slope[i], U_Slope[i + 1]);
		localGRPSolver_Rel.setGamma(gamma);
		localGRPSolver_Rel.solve();
		timeDeris[i] = localGRPSolver_Rel.timeDerivative();
		GodunovSols_Rel[i] = localGRPSolver_Rel(0.0, 1.0);
	}

	auto fluxFun = [&](const vec3d& U, const vec3d& U_t) -> vec3d
		{
			double rho = U(0), rho_t = U_t(0);
			double p = U(1), p_t = U_t(1);
			double u = U(2), u_t = U_t(2);
			double e = p / rho / (gamma - 1.0), e_t = p_t / rho / (gamma - 1.0) - p * rho_t / rho / rho / (gamma - 1.0);
			double x = rho * u + k / 2.0 * (rho * u_t + u * rho_t);
			double y = rho * u * u + p + k / 2.0 * (rho_t * u * u + 2.0 * rho * u * u_t + p_t);
			double z = rho * u * (e + u * u / 2.0) + p * u + k / 2.0 * (rho_t * u * e + rho * u_t * e + rho * u * e_t + rho_t * u * u * u / 2.0 + 1.5 * rho * u * u * u_t + p_t * u + p * u_t);
			return { x, y, z };
		};
	for (auto i = 2; i < N + 2; i++)
	{
		vec3d F_L = fluxFun(GodunovSols_Rel[i - 1], timeDeris[i - 1]);
		vec3d F_R = fluxFun(GodunovSols_Rel[i], timeDeris[i]);

		U_[i] = U_[i] - k / h * (F_R - F_L);
		U[i] = { U_[i](0),
				(gamma - 1.0) * (U_[i](2) - U_[i](1) * U_[i](1) / U_[i](0) / 2.0),
				U_[i](1) / U_[i](0)};
	}
	U.erase(U.begin(), U.begin() + 2);
	U.erase(U.end() - 2, U.end());
}

void Lx_relGRP::setGhostCellStrategy(const string& strategy)
{
	this->ghostCellStrategy = strategy;
}

void Lx_relGRP::setCustomGhostCell(const std::function<void(vector<vec3d>&, double)>& customCellTreatment)
{
	addGhostCells = customCellTreatment;
}

void Lx_relGRP::setAlpha(double alpha)
{
	this->alpha = alpha;
}

void Lx_relGRP::setGamma(double gamma)
{
	this->gamma = gamma;
}

