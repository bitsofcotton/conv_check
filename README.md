# Konbu Check
Get feasible point from multiple linear constraints.

This program aims to check and gain the solvee from multiple set of linear constraints in O(mn^2) arithmetic and O(mn) sqrt function operations for m constraints n variable problem, but the accuracy is very severe.

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

When this condition, if P'^t -1 is inner point, it's ok, otherwise,
we add some salt on them: P''x+a&lt;=a, 0&lt;a, ||a||/||x||-&gt;0,
(if a is orthogonal to P'', this changes feasible region if
 we avoid left equation's a, simply a little extends them.).
then, scaling projection itself results some crossing point to fix
hyperplane. So extreme small a and any of the direction leads us
any of positive nonzero projection can be fixed if ||a||/||x||-&gt; 0.
we can choose P''^t -1 norm on ||x|| first, then, move ||a||&lt;||x||&epsilon;
amounts, and calculate them while that's illegal.
But in this, the accuracy is so severe.

And in loop orthogonality on P' is saved structure.

# Usage
    // if you need, please scope with namespace block, but include guard may harms.
    #include "konbu_init.h"
    ...
    // if you use with mpfr, num_t::set_default_proc(BITS); is needed.
    // if you use with QD,   unsigned int old_cw; fpu_fix_start(&old_cw); is needed.
    ...
    Mat A(/* some rows */, /* some cols */);
    Vec b(A.rows());
    ...
    const auto error(A * Linner<num_t>().inner(A, b) - b);
    ...
    // if you use with QD,    fpu_fix_end(&old_cw); is needed.

# Little Tips
This program solves the feasibility better if the problem variable range is in certain ranges especially eg. [-&pi;, &pi;] and normalized.
And, if we can fix all in once the inner vector instead of fixing one by one, it's O(lg(n)*lg(mn)) time order but is seems not. 
If we have the hyperplane that doesn't include origin point, the feasible point this program get is stable.

# Another download sites.
* https://konbu.azurewebsites.net/ (Sample Site)
* https://drive.google.com/drive/folders/1B71X1BMttL6yyi76REeOTNRrpopO8EAR?usp=sharing
* https://1drv.ms/u/s!AnqkwcwMjB_PaDIfXya_M3-aLXw?e=qzfKcU
* https://ja.osdn.net/users/bitsofcotton/
* https://www.sourceforge.net/projects/convcheck/

# Archive
This repository is archived, so without bugreport, will no change.

