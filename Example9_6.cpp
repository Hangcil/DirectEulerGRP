#include <iostream>
#include <fstream>
#include "wrappers.h"


int main(int, char**)
{
	// analytic solution
	vec3d Ul(1.0, 0.4, -2.0);
	vec3d Ur(1.0, 0.4, 2.0);
	auto solver = RPSolver(Ul, Ur);
	solver.solve();

	int n = 100;
	auto U0 = vector<vec3d>(n);
	for (auto i = 0; i < n; i++)
	{
		if (i < n / 2)
		{
			U0[i] << 1.0, 0.4, -2.0;
		}
		else
		{
			U0[i] << 1.0, 0.4, 2.0;
		}
	}
	FVM_1D numericSolver(U0, 0.0, 1.0);
	numericSolver.setEndingTime(0.1);
	auto r = numericSolver.solve();

	// write the numeric result to the files
	ofstream oFileRho("./test_rho.txt");
	ofstream oFileP("./test_p.txt");
	ofstream oFileU("./test_u.txt");
	if (oFileRho && oFileP && oFileU)
	{
		for (auto i = 0; i < r.size(); i++)
		{
			oFileRho << r[i](0) << ' ';
			oFileP << r[i](1) << ' ';
			oFileU << r[i](2) << ' ';
		}
		oFileRho.close();
		oFileP.close();
		oFileU.close();
	}

	ofstream oFileRho_("./test_rho_exa.txt");
	ofstream oFileP_("./test_p_exa.txt");
	ofstream oFileU_("./test_u_exa.txt");
	if (oFileRho_ && oFileP_ && oFileU_)
	{
		for (auto i = 0; i < 1000; i++)
		{
			oFileRho_ << solver(0.001 * double(i - 500), 0.1)(0) << ' ';
			oFileP_ << solver(0.001 * double(i - 500), 0.1)(1) << ' ';
			oFileU_ << solver(0.001 * double(i - 500), 0.1)(2) << ' ';
		}
		oFileRho_.close();
		oFileP_.close();
		oFileU_.close();
	}

	return 0;
}

