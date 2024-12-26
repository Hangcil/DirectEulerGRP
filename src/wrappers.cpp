#include "wrappers.h"


FVM_1D::FVM_1D(const vector<vec3d>& U0, double lB, double rB) :U(U0), lB(lB), rB(rB)
{
	h = (rB - lB) / double(U0.size());
}

void FVM_1D::setEndingTime(double T)
{
	this->T = T;
}

void FVM_1D::setCFL(double CFL)
{
	this->CFL = CFL;
}

void FVM_1D::setAlpha(double alpha)
{
	this->alpha = alpha;
}

void FVM_1D::setGamma(double gamma)
{
	this->gamma = gamma;
}

void FVM_1D::setGhostCellStrategy(const string& strategy)
{
	this->ghostCellStrategy = strategy;
	useCustomGhostCell = false;
}

void FVM_1D::setCustomCellTreatment(const std::function<void(vector<vec3d>&, double)>& treatmentFun)
{
	customCellTreatment = treatmentFun;
	useCustomGhostCell = true;
}

vector<vec3d> FVM_1D::solve()
{
	double currentTime = 0.0;
	while (currentTime < T)
	{
		double c_max = maxPropagtingSpeed(U, gamma);
		double dt = CFL * h / c_max;
		dt = std::min(dt, T - currentTime);
		if (dt <= 1e-10 * T)
		{
			return U;
		}
		Lx_relGRP solver;
		if (!useCustomGhostCell)
		{
			solver.setGhostCellStrategy(ghostCellStrategy);
		}
		else
		{
			solver.setCustomGhostCell(customCellTreatment);
		}
		solver.setAlpha(alpha);
		solver.setGamma(gamma);
		solver(U, dt, lB, rB, currentTime);
		currentTime += dt;
	}
	return U;
}

