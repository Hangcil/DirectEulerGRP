# DirectEulerGRP
The GRP scheme is a second-order-accuracy finite volume method(FVM) that utilizes the time derivatives on the cell interfaces, and achieves stunning numeric results in a wide variety of CFD topics.

# What algorithm is this implementation for?
The GRP scheme is implemented by the guide of the article [A direct Eulerian GRP scheme for compressible fluid flows](https://www.sciencedirect.com/science/article/abs/pii/S0021999106000581?via%3Dihub). And we also provide Godunov and 2nd-order Runge-Kutta methods, which can be seen as the fundamentals for the GRP scheme.

# Test Results
We choose Leblanc's problem that consists of two constant initial states $(\rho,p,u)_l=(2.0, 10^9, 0.0)$ and $(\rho,p,u)_r=(0.001, 1.0, 0.0)$, and predict the solution at time $T=0.0001$. The spatial domain $[-10,10]$ is divided into 200 cells, but it's worth mentioning that, for a classical higher-order numerical scheme like WENO, Leblanc's problem costs thousands of cells to obtain a correct capture of the shock position. So, for comparation we also calculate with 2nd-order Runge-Kutta method, and the test results are represented as the following plot. The most impressive part is the high accuracy of the GRP scheme with merely 200 cells, which makes the GRP scheme possibly the most efficient numeric solver for certain topics in CFD.
![Leblanc's problem, GRP vs 2nd-order Runge-Kutta](https://github.com/Hangcil/DirectEulerGRP/blob/main/test.jpg)
