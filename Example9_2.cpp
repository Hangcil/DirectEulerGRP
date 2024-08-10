#include <iostream>
#include <fstream>
#include "FVM.h"

class solver9_2 :public FVM_GRP
{
public:
	solver9_2(const vector<vec3d>& U, double lBoundary, double rBoundary) :FVM_GRP(U, lBoundary, rBoundary) {}
protected:
	// setup the boundary conditions
	virtual void setBoundaryValuesAndSlopes() override
	{
		if (currentTime <= 50.0)
		{
			U_n[0] << 4.0, 4.0 / 3.0, -0.3;
			U_n[1] = U_n[0];
			U_n_[0] = rhoPUToConserVar(U_n[0]);
			U_n_[1] = rhoPUToConserVar(U_n[1]);
		}
		else if (currentTime <= 166)
		{
			U_n[0] << 4.0000021378484485, 1.3333345210271275, -0.30000039836446712;
			U_n[1] = U_n[0];
			U_n_[0] = rhoPUToConserVar(U_n[0]);
			U_n_[1] = rhoPUToConserVar(U_n[1]);
		}
		else
		{
			U_n[0] << 3.9999887500437716, 1.3333345210271275, -0.30000039836446712;
			U_n[1] = U_n[0];
			U_n_[0] = rhoPUToConserVar(U_n[0]);
			U_n_[1] = rhoPUToConserVar(U_n[1]);
		}
		U_n[spatialSize - 1] << 1.0, 1.0e-6, -1.3;
		U_n[spatialSize - 2] = U_n[U_n.size() - 1];
		U_n_[spatialSize - 1] = rhoPUToConserVar(U_n[spatialSize - 1]);
		U_n_[spatialSize - 2] = rhoPUToConserVar(U_n[spatialSize - 2]);

		U_Slope[0] = { 0.0,0.0,0.0 };
		U_Slope[1] = { 0.0,0.0,0.0 };
		U_Slope[2] = { 0.0,0.0,0.0 };
		U_Slope[spatialSize - 1] = { 0.0,0.0,0.0 };
		U_Slope[spatialSize - 2] = { 0.0,0.0,0.0 };
		U_Slope[spatialSize - 3] = { 0.0,0.0,0.0 };

	}
};

int main(int, char**)
{
	// analytic solution
	vec3d Ul(4.0, 4.0 / 3.0, -0.3);
	vec3d Ur(1.0, 1.0e-6, -1.3);
	auto solver = RPSolver(Ul, Ur);
	solver.setGamma(5.0 / 3.0);
	auto r1 = solver.solve();
	double sigmaL = (r1.mlState(0) * r1.mlState(2) - 4.0 * (-0.3)) / (r1.mlState(0) - 4.0);
	double sigmaR = (r1.mrState(0) * r1.mrState(2) - 1.0 * (-1.3)) / (r1.mrState(0) - 1.0);

	// numeric solution
	auto U = vector<vec3d>(100);
	for (auto i = 0; i < 100; i++)
	{
		if (i <= 19)
		{
			U[i] << 4.0, 4.0 / 3.0, -0.3;
		}
		else
		{
			U[i] << 1.0, 1.0e-6, -1.3;
		}
	}
	solver9_2 numericSolver(U, 0.0, 100.0);
	numericSolver.setGamma(5.0 / 3.0);
	numericSolver.setTimeAxis(2000.0, 0.2);
	numericSolver.setAlpha(1.0);
	auto r = numericSolver.solve();

	// write the numeric result to the files
	auto r_T = r[r.size() - 1];
	vector<double> rho, p, u;
	for (auto i = 0; i < r_T.size(); i++)
	{
		rho.push_back(r_T[i](0));
		p.push_back(r_T[i](1));
		u.push_back(r_T[i](2));
	}
	ofstream oFileRho("./test_rho.txt");
	if (oFileRho)
	{
		for (auto i = 0; i < 100; i++)
		{
			oFileRho << rho[i] << ' ';
		}
		oFileRho.close();
	}
	ofstream oFileP("./test_p.txt");
	if (oFileP)
	{
		for (auto i = 0; i < 100; i++)
		{
			oFileP << p[i] << ' ';
		}
		oFileP.close();
	}
	ofstream oFileU("./test_u.txt");
	if (oFileU)
	{
		for (auto i = 0; i < 100; i++)
		{
			oFileU << u[i] << ' ';
		}
		oFileU.close();
	}

	return 0;
}