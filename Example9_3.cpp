#include <iostream>
#include <fstream>
#include "FVM.h"

int main(int, char**)
{
	// numeric solution
	auto U = vector<vec3d>(70);
	for (auto i = 0; i < 70; i++)
	{
		if (i < 7)
		{
			U[i] << 2.8182, 5.0, 1.6064;
		}
		else if (i == 7)
		{
			vec3d U_local1, U_local2;
			U_local1 << 2.8182, 5.0, 1.6064;
			U_local2 << 1.0, 1.0, 0.0;
			U[i] = 0.5 * (1.1 * U_local1 + 0.9 * U_local2);
		}
		else if (i <= 19)
		{
			U[i] << 1.0, 1.0, 0.0;
		}
		else
		{
			U[i] << 0.3, 1.0, 0.0;
		}
	}
	FVM_GRP numericSolver(U, -40.0, 100.0);
	numericSolver.setGamma(1.4);
	numericSolver.setTimeAxis(30.0, 0.2);
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