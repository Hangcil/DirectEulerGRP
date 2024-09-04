#include <iostream>
#include <fstream>
#include "FVM.h"

int main(int, char**)
{
	// analytic solution
	vec3d Ul(1.0, 1.0, 0.0);
	vec3d Ur(0.125, 0.100, 0.0);
	auto solver = RPSolver(Ul, Ur);
	solver.setGamma(1.4);
	solver.solve();

	// numeric solution
	auto U = vector<vec3d>(100);
	for (auto i = 0; i < 100; i++)
	{
		if (i < 50)
		{
			U[i] << 1.0, 1.0, 0.0;
		}
		else
		{
			U[i] << 0.125, 0.100, 0.0;
		}
	}
	FVM_GRP numericSolver(U, 0.0, 100.0);
	numericSolver.setGamma(1.4);
	numericSolver.setTimeAxis(15.0, 0.2);
	auto r = numericSolver.solve();

	// write the numeric result to the files
	auto r_T = r[r.size() - 1];
	ofstream oFileRho("./test_rho.txt");
	ofstream oFileP("./test_p.txt");
	ofstream oFileU("./test_u.txt");
	if (oFileRho && oFileP && oFileU)
	{
		for (auto i = 0; i < r_T.size(); i++)
		{
			oFileRho << r_T[i](0) << ' ';
			oFileP << r_T[i](1) << ' ';
			oFileU << r_T[i](2) << ' ';
		}
		oFileRho.close();
		oFileP.close();
		oFileU.close();
	}

	return 0;
}
