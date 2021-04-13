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

|[-b''',P][t,x']|&lt;=1.

After loop we get :
| P'[Q[t',x'',0]] |&lt;=0

# Usage
    // if you need, please scope with namespace block, but include guard may harms.
    #include "konbu.hh"
    ...
    SimpleMatrix<num_t> A(/* some rows */, /* some cols */);
    SimpleVector<num_t> b(A.rows());
    ...
    const auto error(A * inner<num_t>(A, b * num_t(0), b) - b);
    ...

# Tips
This program optimize inner point in Ax &lt;= bt meaning, so t moves \[0, \infty\[ .
So if t &lt;= 1, it's unstable but this is specification of this.

# Another download sites.
* https://konbu.azurewebsites.net/ (Sample Site)
* https://drive.google.com/drive/folders/1B71X1BMttL6yyi76REeOTNRrpopO8EAR?usp=sharing
* https://1drv.ms/u/s!AnqkwcwMjB_PaDIfXya_M3-aLXw?e=qzfKcU
* https://ja.osdn.net/users/bitsofcotton/
* https://www.sourceforge.net/projects/convcheck/

# Archive
This repository is archived, so without bugreport, will no change.

