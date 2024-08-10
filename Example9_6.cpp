#include <iostream>
#include <fstream>
#include "FVM.h"

int main(int, char**)
{
	// analytic solution
	vec3d Ul(1.0, 0.4, -2.0);
	vec3d Ur(1.0, 0.4, 2.0);
	auto solver = RPSolver(Ul, Ur);
	solver.setGamma(1.4);
	auto r1 = solver.solve();
	double c_ml = sqrt(1.4 * r1.mlState(1) / r1.mlState(0));
	double c_mr = sqrt(1.4 * r1.mrState(1) / r1.mrState(0));
	
	// numeric solution
	auto U = vector<vec3d>(100);
	for (auto i = 0; i < 100; i++)
	{
		if (i < 50)
		{
			U[i] << 1.0, 0.4, -2.0;
		}
		else
		{
			U[i] << 1.0, 0.4, 2.0;
		}
	}
	FVM_GRP numericSolver(U, 0.0, 100.0);
	numericSolver.setGamma(1.4);
	numericSolver.setTimeAxis(10.0, 0.2);
	numericSolver.setAlpha(1.9);
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