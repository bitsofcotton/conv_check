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

|\[-b''',P\]\[t,x'\]|&lt;=1, t == 1.

scaling b*_k := 1,
|P x' - 1|&lt;=1

To minimize: |&lt;p\_k,x&gt; - 1|, we choose orthogonality condition on (P P^t 1)'s minimum to fix, this also includes also -1 intercept condition with some ratio.

# Usage
    // if you need, please scope with namespace block, but include guard may harms.
    #include "simplelin.hh"
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

