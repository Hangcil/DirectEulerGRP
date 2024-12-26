#ifndef BASICFUNS_H
#define BASICFUNS_H

#include "solvers.h"
#include <utility>
#include <string>
#include <fstream>
using std::string;


double maxPropagtingSpeed(const vector<vec3d>& U, double gamma);
vec3d rhoPUToConserVar(const vec3d& rhoPUV, double gamma);

class Lx_relGRP
{
public:
	void operator()(vector<vec3d>& U, double k, double lB, double rB, double currentTime);
	void setGhostCellStrategy(const string& strategy);
	void setCustomGhostCell(const std::function<void(vector<vec3d>&, double)>& customCellTreatment);
	void setAlpha(double alpha);
	void setGamma(double gamma);
protected:
	string ghostCellStrategy = "flat";
	vector<vec3d> U_;
	double lb = 0.0, rb = 0.0;
	int N = 0;
	double h = 0;
	double CFL = 0.5;
	double alpha = 1.9;
	double gamma = 1.4;
	std::function<void(vector<vec3d>&, double)> addGhostCells = [&](vector<vec3d>& U, double currentTime)->void
		{
			if (ghostCellStrategy == "flat")
			{
				U.insert(U.begin(), U[0]);
				U.insert(U.begin(), U[0]);
				U.push_back(U.back());
				U.push_back(U.back());
			}
			else if (ghostCellStrategy == "reflective")
			{
				U.insert(U.begin(), U[0]);
				U.insert(U.begin(), U[2]);
				U[0](2) = -U[3](2);
				U[1](2) = -U[2](2);
				U.push_back(U.back());
				U.push_back(U[U.size() - 3]);
				U[U.size() - 1](2) = -U[U.size() - 4](2);
				U[U.size() - 2](2) = -U[U.size() - 3](2);
			}
			else if (ghostCellStrategy == "periodic")
			{
				auto temp0 = U[0], temp1 = U[1];
				U.insert(U.begin(), U.back());
				U.insert(U.begin(), U[U.size() - 2]);
				U.push_back(temp0);
				U.push_back(temp1);
			}
		};
};






#endif