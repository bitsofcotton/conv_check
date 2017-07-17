# Konbu Check
Get feasible point from multiple linear constraints.

This program aims to check and gain the solvee from multiple set of linear constraints.
C++ and Eigen library needed, and to calculate more accurate, we may need a gmp library, a mpfr library and a real.hpp library.
For older information, please refer http://sourceforge.net/projects/convcheck/

# Tips
In steps function, we cannot get inner point when the ratio of diameter are large or inner points are far enough from origin point.
This depends on huge intercepts on rank decreasing loop.

# Bugs
If the accuracy or parameter configuration is not valid for the problem to be solved, the feasibility that this program checks will bugly, If original problem is good scaled, it's hard to be happen, but, scaling in konbu.hh will destroys such condition even if good, so if you don't need scaling, please comment out the corresponding code.

# Proof
Ax&lt;=b

[-b,P][t,x']&lt;=0
P is part of orthogonal matrix.

[-b'+1&epsilon;,P][t,x'+b'']&lt;=1&epsilon;
b' is orthogonal to P.

After loop we get :
P'[Q[t,x'']]&lt;=0

# Usage
    #include "konbu_init.h"
    ...
    int m; // number of rows;
    int n; // number of columns;
    Mat A(m, n);
    ...
    Vec result;
    bool fix_partial[A.rows()];
    LP<num_t> lp;
    bool feas = lp.inner(fix_partial
    result, A, b);

# Makefile
Please configure Makefile manually, there's -DACC_GMP=$bits option or -DDOUBLE option and so on, -DWITHOUT_EIGEN option, and compiler options such like -fopenmp, -pg, -I$where_eigen_lives, ...
