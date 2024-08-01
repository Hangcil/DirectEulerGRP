# DirectEulerGRP
The GRP scheme is a second-order-accuracy finite volume method(FVM) that utilizes the time derivatives on the cell interfaces, which produces stunning numeric results in a wide variety of CFD topics.

# What algorithm is this implementation for?
See the article [A direct Eulerian GRP scheme for compressible fluid flows.](https://www.sciencedirect.com/science/article/abs/pii/S0021999106000581?via%3Dihub)

# Test Results
We choose Leblanc's problem that consists two constant initial states $(\rho,p,u)_l=(2.0, 10^9, 0.0)$ and $(\rho,p,u)_r=(0.001, 1.0, 0.0)$, and predict the solution at time $T=0.0001$. The spatial domain $[-10,10]$ is divided into 200 cells, but it's worth mentioning that, in a typical higher-order numerical scheme, Leblanc's problem requires thousands of cells to obtain a correct capture of the shock. For comparation we also calculate with 2nd-order Runge-Kutta method, and the test results are represented as the following plot. The most impressive part is the accuracy of the GRP scheme with a considerably fewer dicretized mesh.
![Leblanc's problem, GRP vs 2nd-order Runge-Kutta](https://github.com/Hangcil/DirectEulerGRP/blob/main/test.jpg)
