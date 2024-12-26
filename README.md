# The aim of this reposiry
This reposiry is only used for the academic purpose. It implements the direct Eulerian GRP scheme introduced in [A direct Eulerian GRP scheme for compressible fluid flows](https://www.sciencedirect.com/science/article/abs/pii/S0021999106000581) The algorithm only completes the 1D case, and the Strang splitting method should be further implemented based on this version.

# Functions that has been implemetned
1. The 1D Riemann solver for the Euler equations.
   
   Usage. Define the left and right states using ```vec3d ul; vec3d ur;```, where ```vec3d=Eigen::Vector3d``` is an alias of the Eigen vector. Define the Riemann solver by ```auto solver=RPSolver(ul,ur)```. To solve the Riemman problem, call ```solver.solve()```. The Riemann solution is obtained by calling an overloaded ```()```: ```vec3d solution=solver(x,t)```, where ```x```: spatial coordinate and ```t```: time.

   Options. Call ```solver.setGamma(gamma)``` to change the value of gamma (default value: 1.4) . Call ```solver.setTol(tol)``` to change the tolerance of the error of the Riemann solution, which requires a solution to the nonlinear equation.

2. The 1D direct Eulerian GRP solver  for the Euler equations.

   Usage. The class ```GRPSolver``` is directly inheritted from ```RPSolver```, which means it preserves the full capability of ```RPSolver```. Define the left and right slopes using ```vec3d ulSlope; vec3d urSlope;```, and define the solver using ```auto gSolver=GRPSolver(ul,ur,ulSlope,urSlope);```. To solve the generalized Riemman problem, call ```gSolver.solve()```. The Riemann solution is obtained by calling an overloaded ```()```: ```vec3d solution=gSolver(x,t)```, and the time derivatives is obtained by calling ```vec3d derivatives=gSolver.timeDerivatives()```.

   Options. Remains the same as Riemann solver.

3. The 1D finite volume algorithm.

   Usage. Define the initial data ```auto U0 = vector<vec3d>(n);```. Then define the FVM solver by ```FVM_1D solver(U0, lB, rB);```, where ```lB,rB``` represent the left/right boundaries of the spatial domain. Set the ending time by calling ```solver.setEndingTime(T);```, and solve with ```vector<vec3d> result=solver.solve()```.

   Options. Set the ghost cell strategy by ```solver.setGhostCellStrategy("reflective")```. Available strategies: "flat"(default), "reflective", "periodic". Set $\gamma$ by calling ```setGamma(gamma)```. Set CFL by calling ```setCFL(CFL)```. Set the parameter in the slope limiter by calling ```setAlpha(alpha)``` (1~2). You can set up your own ghost cell strategy by calling ```setCustomCellTreatment(fun);```, where ```fun``` is a function that adds values to the ghost cells.

# Compile & run
This repo only requires the Eigen library. Install it and the dependency should be OK. On Windows, revise the path to the Eigen in the ```CMakeLists.txt```. The ```Example__.cpp``` and ```TestAccuracy.cpp``` are both ready to execute.
