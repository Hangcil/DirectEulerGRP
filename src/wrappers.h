#ifndef WRAPPERS_H
#define WRAPPERS_H

#include "basicFuns.h"

class FVM_1D
{
public:
	FVM_1D(const vector<vec3d>& U0, double lB, double rB);
	void setEndingTime(double T);
	void setCFL(double CFL);
	void setAlpha(double alpha);
	void setGamma(double gamma);
	void setGhostCellStrategy(const string& strategy);
	void setCustomCellTreatment(const std::function<void(vector<vec3d>&, double)>& treatmentFun);
	vector<vec3d> solve();
protected:
	double T = 1.0;
	double CFL = 0.5;
	double alpha = 1.9;
	double lB = 0.0, rB = 0.0, h = 0.0;
	double gamma = 1.4;
	vector<vec3d> U;
	string ghostCellStrategy = "flat";
	std::function<void(vector<vec3d>&, double)> customCellTreatment;
	bool useCustomGhostCell = false;
};



#endif