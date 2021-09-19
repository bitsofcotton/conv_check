# Konbu Check
Get feasible point from multiple linear constraints with both side conditions.

This program aims to check and gain the solvee from multiple set of both side linear constraints in O(mn^2) arithmetic and O(mn) sqrt function operations for m constraints n variable problem.

# Context
Japanese Patent Draft : JP2014089683 . 

# Proof
b_l&lt;=Ax&lt;=b_u

-b&lt;=Ax+b'&lt;=b

|Ax+b'|&lt;=b

|A'x+b''|&lt;=1

so pass linear invariant on 0 &lt;= |A'' x' + 1| &lt;= |b'''| causes
calculation ok.

To calculate linear invariant, ||P' x''||_2 / ||x''||_2 -&gt; maximum condition in |P' x''|&lt;=1 epsilon, minimize sup_k |&lt;p', x''&gt;|\_k. But optimal condition fixes some on the index |P' x'|&lt;=1 each (optimal is on the some of a vertex != 0. i.e. some of line segment.), causes ||P'_partial x''||_2 / ||x''||_2 -&gt; maximum, so the first condition, find minimum of ||P'_partial x''||_2 / ||x''||_2 -&gt; maximim condition, This can be done by sort each |&lt;p',x''&gt;|, then fix them in ascendant order because fix one of the index causes orthogonalize original matrix and it's linear dependant in 2nd-norm condition. This method also finds minimum combination on |&lt;p'x,x''&gt;|.

# Usage
    // if you need, please scope with namespace block, but include guard may harms.
    #include "lieonn.hh"
    ...
    SimpleMatrix<num_t> A(/* some rows */, /* some cols */);
    SimpleVector<num_t> b(A.rows());
    ...
    const auto error(A * A.inner(b * num_t(0), b) - b);
    ...

# Another download sites.
* https://konbu.azurewebsites.net/ (Sample Site)
* https://drive.google.com/drive/folders/1B71X1BMttL6yyi76REeOTNRrpopO8EAR?usp=sharing
* https://1drv.ms/u/s!AnqkwcwMjB_PaDIfXya_M3-aLXw?e=qzfKcU
* https://ja.osdn.net/users/bitsofcotton/
* https://www.sourceforge.net/projects/convcheck/

# Archive
This repository is archived, so without bugreport, will no change.

