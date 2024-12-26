#include <iostream>
#include <fstream>
#include "wrappers.h"


int main(int, char**)
{
	// numeric solution
	int n = 1600;
	auto U = vector<vec3d>(n);
	for (auto i = 0; i < n; i++)
	{
		if (i < n / 10)
		{
			U[i] << 1.0, 1000.0, 0.0;
		}
		else if (i < n * 0.9)
		{
			U[i] << 1.0, 0.01, 0.0;
		}
		else
		{
			U[i] << 1.0, 100.0, 0.0;
		}
	}
	auto solver = FVM_1D(U, 0.0, 1.0);
	solver.setEndingTime(0.038);
	solver.setGhostCellStrategy("reflective");

	// write the numeric result to the files
	auto r_T = solver.solve();
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