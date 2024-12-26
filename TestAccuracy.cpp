#include <iostream>
#include <fstream>
#include <algorithm>
#include "wrappers.h"

double error_num(int N)
{
	double h = 2.0 / N;
	vector<vec3d> U0(N);
	for (auto i = 0; i < N; i++)
	{
		double x = (0.5 + double(i)) * h;
		U0[i] << 1.0 + 0.2 * sin(EIGEN_PI * x), 1.0, 1.0;
	}

	double T = 1.0;
	FVM_1D solver(U0, 0.0, 2.0);
	solver.setEndingTime(T);
	solver.setGhostCellStrategy("periodic");
	auto r = solver.solve();



	vector<vec3d> U_exa(N);
	for (auto i = 0; i < N; i++)
	{
		U_exa[i] << 1.0 - 0.2 * sin(EIGEN_PI * ((double(i) + 0.5) * h)), 1.0, 1.0;
	}

	double error_sum = 0.0;
	for (auto i = 0; i < N; i++)
	{
		error_sum += (r[i] - U_exa[i]).cwiseAbs()(0) * h;
	}

	return error_sum;
}

int main(int, char**)
{
	double e1 = error_num(25);
	double e2 = error_num(50);
	double a12 = log(e1 / e2) / log(2);
	double e3 = error_num(100);
	double a23 = log(e2 / e3) / log(2);
	double e4 = error_num(200);
	double a34 = log(e3 / e4) / log(2);
	double e5 = error_num(400);
	double a45 = log(e4 / e5) / log(2);
	/*double e6 = error_num(800);
	double a56 = log(e5 / e6) / log(2);
	double e7 = error_num(1600);
	double a67 = log(e6 / e7) / log(2);*/


	return 0;
}