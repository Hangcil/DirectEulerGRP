#include <iostream>
#include <fstream>
#include "FVM.h"

int main(int, char**)
{
	// analytic solution
	vec3d Ul(2.0, 1.0e9, 0.0);
	vec3d Ur(0.001, 1.0, 0.0);
	auto solver = RPSolver(Ul, Ur);
	solver.setGamma(1.4);
	auto r1 = solver.solve(0, 0);

	// numeric solution
	auto U = vector<vec3d>(200);
	for (auto i = 0; i < 200; i++)
	{
		if (i < 100)
		{
			U[i] << 2.0, 1.0e9, 0.0;
		}
		else
		{
			U[i] << 0.001, 1.0, 0.0;
		}
	}
	FVMSolver numericSolver(U, -10.0, 10.0);
	numericSolver.setGamma(1.4);
	numericSolver.setTimeAxis(0.0001, 0.0000001);
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
		for (auto i = 0; i < 200; i++)
		{
			oFileRho << rho[i] << ' ';
		}
		oFileRho.close();
	}
	ofstream oFileP("./test_p.txt");
	if (oFileP)
	{
		for (auto i = 0; i < 200; i++)
		{
			oFileP << p[i] << ' ';
		}
		oFileP.close();
	}
	ofstream oFileU("./test_u.txt");
	if (oFileU)
	{
		for (auto i = 0; i < 200; i++)
		{
			oFileU << u[i] << ' ';
		}
		oFileU.close();
	}


	return 0;
}