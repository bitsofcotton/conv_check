# Konbu Check
Get feasible point from multiple linear constraints.

This program aims to check and gain the solvee from multiple set of linear constraints.
C++ and Eigen library needed, and to calculate more accurate, we may need gmp library, mpfr library and real.hpp library.
For older information, please refer http://sourceforge.net/projects/convcheck/

# Tips
In steps function, we cannot get inner point when ratio of diameter are large or inner points far enough from origin point.
This depends on huge intercepts on rank decreasing loop.

# Proof
Ax&lt;=b

[-b,P][t,x']&lt;=0

[-b+1&epsilon;,P][t,x']&lt;=1&epsilon;
