# Konbu Check
Get feasible point from multiple linear constraints.

This program aims to check and gain the solvee from multiple set of linear constraints.
C++ is needed, and to calculate faster, Eigen library is needed, and to calculate more accurate (in the most case we need), we may need a gmp library, a mpfr library and a real.hpp library.
For older information, please refer http://sourceforge.net/projects/convcheck/ .

# Tips
The shown string after 'err_error' or 'intercept' is the value depends on the problem and accuracy.
If the value >> 0 (especially >= 1), it is hard to solve in the accuracy.

# Bugs
If the accuracy or parameter configuration is not valid for the problem to be solved, the feasibility that this program checks will be bugly, If original problem is good scaled, it rarely happens.

# Parameters
We shall configure the parameters in LP<T>::LP() in konbu.hh. Description is in https://konbu.sakura.ne.jp/ .

# Proof
Ax&lt;=b

[-b,P][t,x']&lt;=0
P is part of orthogonal matrix.

[-b'+1&epsilon;,P][t,x'+b'']&lt;=1&epsilon;
b' is orthogonal to P.

After loop we get :
P'[Q[t,x'',0]]&lt;=0

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
    bool feas = lp.inner(fix_partial, result, A, b);

# Makefile
Please configure Makefile manually, there's -DACC_GMP=$bits option or -DACC_DOUBLE option, -DWITHOUT_EIGEN option, and compiler options such like -fopenmp, -pg, -I$where_eigen_lives, ...

