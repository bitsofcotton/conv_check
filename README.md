# Konbu Check
Get feasible point from multiple linear constraints with both side conditions.

This program aims to check and gain the solvee from multiple set of both side linear constraints in O(mn^2) arithmetic and O(mn) sqrt function operations for m constraints n variable problem.

# Context
Japanese Patent Draft : JP2014089683 . 

# Proof
b_l&lt;=Ax&lt;=b_u

|(2 / b_u)_each A x - 1 - (b_l / b_u)_each| &lt;= |1 - (b_l / b_u)_each|.

choose b' := 1 - (b_l / b_u)_each,
(we can choose |b_l| &lt; |b_u| if same sign, otherwise right hand side &gt;0)

if left hand &gt; b' or left hand &lt; - b', |A'x - 2| &lt; 0 but we don't meet this condition.

else if left hand side &gt; 0 case, it's also |A'x - 2| &lt; 0 but in reasonable condition, |A'x - 2| &lt; 2b' in opposite sign condition.

otherwise: - b' &lt; left hand &lt; 0 case, it's equivalent to
some external intercept: |A'x - 2| &lt; 2 b'.

Twice same method, A'x - 2 &lt; 2b' &lt;=&gt; A'x &lt; 2b' + 2 &lt;=&gt; |A'x| &lt; 2b' + 2, sign range is from feasible region.

So we avoid const. multiply, exchange |b_l| &lt; |b_u| condition, divide by b_u, then, divide by (2b' + 2 == 4 - 2 \* (b_l / b_u)) each the original matrix's linear invariant has to make a sense.

We ignore sign on x, so they're equivalent to |A'x|==2 we also ignore ratio on x, it's the problem linearInvariant A'.

To calculate linear invariant, ||P' x''||_2 / ||x''||_2 -&gt; maximum condition in |P' x''|&lt;=1 epsilon, minimize sup_k |&lt;p', x''&gt;|\_k.
But optimal condition, they fixes some on the index |P' x'|&lt;=1 each (optimal is on the some of a vertex != 0. i.e. some of line segment.), causes ||P'_partial x''||_2 / ||x''||_2 -&gt; maximum.
So the linear invariant condition is equivalent to the first condition, and find minimum of ||P'_partial x''||_2 / ||x''||_2 -&gt; maximim condition.
This can be done by sort each |&lt;p',x''&gt;|, then fix them in ascendant order because fix one of the index causes orthogonalize original matrix and it's linear dependant in 2nd-norm condition.

This method also finds minimum combination on |&lt;p'x,x''&gt;|_1.
(||Px||_2^2-||Px||_1^2 == ||x||_2^2-(Sum\<p\_k,x\>)^2==||x||_2^2-Sum\<p_k,x\>\<p_l,x\>==||x'||_2^2-Sum x'\_kx'\_l) == ||x'||_2^2-||x'||_1^2, d/dt (at(x)\*bt(x))^t==(at(x)\*bt(x))^(t-1)(at'(x)\*bt(x) + at(x)\*bt'(x)), t==4 and t==2 case, ratio is also t==2 case, limit condition).

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
2023/04/21 make/revert ProgramInvariant algorithm changed.

