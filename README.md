# Konbu Check
Get feasible point from multiple linear constraints.

This program aims to check and gain the solvee from multiple set of linear constraints in O(mn^2) arithmetic and elementary function operations.
C++ is needed, and to calculate faster, Eigen library is needed, and to calculate more accurate (in the most case we need), we may need a gmp library, a mpfr library and a mpfr++ library is needed, or, with certain accuracy calculation, QD library is needed.
For older information, please refer http://sourceforge.net/projects/convcheck/ .

# Tips
The shown string after 'err_error' or 'intercept' is the value depends on the problem and accuracy.
If the value >> 0 (especially >= 1), it is hard to solve in the accuracy.  
If feasible region of the original problem is too tight, there's a possibility fails to get feasible point.
If then, please extend little more feasible region by adding A.row norm to b.

# Bugs
If the accuracy or parameter configuration is not valid for the problem to be solved, the feasibility that this program checks will be bugly, If original problem is good scaled and extended little, it rarely happens.

# Parameters
We shall configure the parameters in LP<T>::LP() in konbu.hh.
* threshold_feas   : QR decomposition errors.
* threshold_p0     : each loop P matrix errors.
* threshold_loop   : [-b+1&epsilon;,P][t,x]&leq;b, &epsilon;
* threshold_inner  : each loop checking inner or not.
* largest_intercept: box constraints that Ax&leq;b fix, not Ax&geq;b.
* largest_opt      : assuming optimal value max ratio.
* n_opt_steps      : number of loops to find optimal values.

# Context
Japanese Patent Draft : JP2014089683 . 

# Status
Freezed. But to make this program worthy, some user interface or some program corresponding to phenomenon needed. And waiting g++ for compatible with OpenACC with iris chips.
And deducting when the matrix multiplication and same order similar operations can be done in O(1), maybe this time, O(max(m, n)) arithmetic and matrix multiplication similar operations will be taken. 

# How to install
Please rewrite Makefile as the libraries enabled.

# Demos
http://services.limpid-intensity.info/konbu.php is a working sample.

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
    // if you use with mpfr, num_t::set_default_proc(BITS); is needed.
    // if you use with QD,   unsigned int old_cw; fpu_fix_start(&old_cw); is needed.
    ...
    int m; // number of rows;
    int n; // number of columns;
    Mat A(m, n);
    ...
    Vec result;
    bool fix_partial[A.rows()];
    LP<num_t> lp;
    bool feas = lp.inner(fix_partial, result, A, b);
    // if you use with QD,    fpu_fix_end(&old_cw); is needed.

# Makefile
Please configure Makefile manually, there's -DACC_GMP=$bits option or -DACC_QD_QDOUBLE option or -DACC_DOUBLE option, -DWITHOUT_EIGEN option, and compiler options such like -fopenmp, -pg, -I$where_eigen_lives, ...

# Little Tips
You can also use this program for pattern matching (finding smaller co-space ∞-norm) of Rn to Rn+m function, or, for finding possible solvees of PDEs (with 2013 memo last stage), (even if it is not only one solvee because of restrict equations # but if so, solvee we get in this shall not be useful, region or shape or rank will be useful).
If we are using 32 bit machine, it is limitted that we can solve the problem smaller than around 8k x 4k matrix (because of index type is integer, in fact, if we're using 8 bytes floating point number, and like most implementation pointer 1 bit for system, and matrix copy on the memory.).
And with -DWITHOUT_EIGEN option, we can use this without eigen library, (with simple pure C++ and OpenMP implementation) but it is dramatically slow. 

# Another download sites.
* https://ja.osdn.net/projects/conv-check/
* https://www.sourceforge.net/projects/convcheck/
* https://konbu.sakura.ne.jp/files/konbu_check-1.01.tar.gz
* http://files.limpid-intensity.info/files/konbu_check-1.01.tar.gz (preparing...)
