#include <iostream>
#include <fstream>
#include "FVM.h"

int main(int, char **)
{
    
    auto U = vector<vec3d>(101);
    for (auto i = 0; i < 101; i++)
    {
        vec3d U_local;
        if (i <= 50)
        {
            U_local << 1.0, 1.0, 0.0;
            U[i] = U_local;
        }
        else
        {
            U_local << 0.125, 0.100, 0.0;
            U[i] = U_local;
        }
    }
    FVMSolver numericSolver(U, 0.0, 100.0);
    numericSolver.setGamma(1.4);
    numericSolver.setTimeStep(15.0, 0.02);
    auto r = numericSolver.solve();

    vector<int> errorLevel;
    for (auto i = 0; i < r.size(); i++)
    {
        for (auto j = 0; j < r[0].size(); j++)
        {
            if (r[i][j](0) >= 1.001 || r[i][j](2) >= 1.001)
            {
                errorLevel.push_back(i);
                break;
            }
        }
    }
    
    auto r_20 = r[r.size() -1];
    vector<double> u_20;
    for (auto i = 0; i < r_20.size(); i++)
    {
        u_20.push_back(r_20[i](0));
    }
    

    /*
    vec3d lState(0.125, 0.99255583126550834, 4.4408920985006264e-17), rState(0.125, 0.099255583126550848, -0.0000000000000000);
    vec3d lStateSlope(0,0,0), rStateSlope(0,0,0);
    GRPSolver solver(lState, rState, lStateSlope, rStateSlope);
    solver.setGamma(1.4);
    auto GRPSol = solver.solve();
    */

    /*
    RPSolver solver(lState, rState);
    solver.setGamma(1.4);
    auto RiemannSol = solver.solve();
    double rho_ml = RiemannSol.mlState(0);
    double rho_mr = RiemannSol.mrState(0);
    double u_m = RiemannSol.mlState(2);
    double p_m = RiemannSol.mrState(1);
    cout << "[rho,p,u]_*l= " << rho_ml << ' ' << p_m << ' ' << u_m << endl;
    cout << "[rho,p,u]_*r= " << rho_mr << ' ' << p_m << ' ' << u_m << endl;
    */
    
    
    ofstream oFile("./test.txt");
    if (oFile)
    {
        for (auto i = 0; i < 101; i++)
        {
            oFile << u_20[i] << ' ';
        }
        oFile.close();
    }
    

    return 0;
}
