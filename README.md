# Konbu Check
Get feasible point from multiple linear constraints.

This program aims to check and gain the solvee from multiple set of linear constraints in O(mn^2) arithmetic and O(mn) sqrt function operations for m constraints n variable problem.
We need C++, and to calculate faster, needs Eigen library, and to calculate better accurate (N.B. in the most case we need), we may need a gmp library, a mpfr library and a mpfr++ library, or, with certain accuracy calculation, we need a QD library.

# Tips
The shown string after 'err_error' is the value depends on the problem and accuracy.
If the value >> 0 (especially >= 1), it is hard to solve in the accuracy.
If feasible region of the original problem is too tight, there's a possibility fails to get feasible point.
If then, please extend little more feasible region by changing threshold_loop parameter.  
And, parameter large can be tiny value for the problem with medium accuracy, so please configure this first.  
Include guard definition this uses seems high probability to conflict, please patch before to use.

# Bugs
If the accuracy or parameter configuration is not valid for the problem to be solved, the feasibility that
this program checks will be bugly, If a original problem is good scaled and extended little, it rarely happens.  
And, when we're gaining optimal value, in the parameter initialize function, we assume the largest ratio of
optimal value and variable upper (lower) bound. So if we gained seems to be not optimal, please configure the parameters.

# Parameters
We shall configure the parameters in Linner<T>::Linner() in konbu.hh.
* threshold_feas   : QR decomposition errors.
* threshold_p0     : each loop P matrix errors.
* threshold_loop   : [-b+1&epsilon;,P][t,x]&leq;b, &epsilon;
* threshold_inner  : each loop checking inner or not.
* large            : box constraints that Ax&leq;b fix, not Ax&geq;b.

# Context
Japanese Patent Draft : JP2014089683 . 

# How to install
Please rewrite Makefile as the libraries enabled.
There's -DACC_GMP=$bits option or -DACC_QD_QDOUBLE option or -DACC_DOUBLE option, -D_WITHOUT_EIGEN_ option, and compiler options such like -fopenmp, -pg, -I$where_eigen_lives, ...

# Proof
Ax&lt;=b

[-b,P][t,x']&lt;=0,
P is part of orthogonal matrix.

[-b'+1&epsilon;,P][t,x'+b'']&lt;=1&epsilon;,
b' is orthogonal to P.

After loop we get :
P'[Q[t,x'',0]]&lt;=0

And when loop, the intercept added is remains, but it is fixed in the last.
So we gain P^t z_0 with fixed intercept.

# Usage
    // if you need, please scope with namespace block, but include guard may harms.
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
    Linner<num_t> lp;
    bool feas = lp.inner(fix_partial, result, A, b);
    // if you use with QD,    fpu_fix_end(&old_cw); is needed.

# Little Tips
You can also use this program for pattern matching (finding smaller co-space ∞-norm) of Rn to Rn+m function, or, for finding possible solvees of PDEs (with 2013 memo last stage), (even if it is not only one solvee because of restrict equations # but if so, solvee we get in this shall not be useful, region or shape or rank will be useful).
If we are using 32 bit machine, it is limitted that we can solve the problem smaller than around 8k x 4k matrix (because of index type is integer, in fact, if we're using 8 bytes floating point number, and like most implementation pointer 1 bit for system, and matrix copy on the memory.).  
If, we run this program with over mn core MPUs, we can gain inner vector in O(n) time order.  
This program solves the feasibility better if the problem variable range is in certain ranges especially eg. [-&pi;, &pi;] and normalized.
And, if we can fix all in once the inner vector instead of fixing one by one, it's O(lg(n)*lg(mn)) time order but is seems not. 

# Another download sites.
* https://konbu.azurewebsites.net/ (Sample Site)
* https://sites.google.com/view/bitsofcotton/
* https://ja.osdn.net/projects/conv-check/
* https://www.sourceforge.net/projects/convcheck/

