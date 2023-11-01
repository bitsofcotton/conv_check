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


To calculate linear invariant, -1 &lt;= Px'' &lt;= 1, so we take:
-1 + epsilon 1 &lt; Px'' \pm epsilon 1 &lt; 1 - epsilon 1.

So each loop we fix least enough on effect to |epsilon| with maximizing ||x''||,
overall we get:
|epsilon x\* p \pm epsilon 1(orthogonal to p)| &lt;= 1 - epsilon 1.

We can configure epsilon after doing them.
Also, the effect to epsilon on this form is minimal in abstract value (not the ratio) because we fixing orthogonal ones.

So reverse path leads us to -1 &lt;= P x_0 &lt;= 1 on their x_0, this can satisfy the condition |epsilon| 1 (because first optimize gets minimal in abstract value).

This method can be applied to -1 &lt;= (&lt;p_k,x&gt;)^2 &lt;= 1 conditions (if all of them are squared ones).


If the first hypothesis fails, this is the case P is orthogonal to 1,
with flipping the sign of some each of the original A'' matrix row can
improve the condition, but this is not implemented.

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
2023/04/23 nand.cc fix.
2023/11/01 readme large fix causes a little improve on proof. nand.cc revertProgramInvariant change.

