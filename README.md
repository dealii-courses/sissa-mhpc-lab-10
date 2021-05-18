#  Lab 10 - Mixed problems: Stokes
## Theory and Practice of Finite Elements

**Luca Heltai** <luca.heltai@sissa.it>

Starter code documentation can be accessed here:

https://dealii-courses.github.io/sissa-mhpc-lab-10/

* * * * *
## Lab-10

1. Derive from the `BaseProblem` class a `BaseBlockProblem` class, replacing the
flat versions of the Linear alagebra classes with `Block` versions, initialized
with a number of blocks deduced at construction time, using a `blocking`
parameter, i.e., `blocking = {0, 0, 0, 1}` will refer to a system with two
blocks, where in the first block we group the first three components (`0,0,0`
-> `0` is the block index), and in the second block we group the fourth
component.

2. Make sure you renumber your degrees of freedom block-wise, and you
initialize the block matrices with the correct block structure. Take as an
example `step-32` of the library.

3. Create a `Stokes` problem, derived from the `BaseBlockProblem` class, and
assemble the mixed stokes system.

4. Use the `LinearOperator` class facilities to create an efficient Schur
complement solver for the Stokes system (follow the discussion on `step-60` for
the resolution of the block system).

5. Create a non-trivial manufactured solution starting from an arbitrary vector
valued smooth function, and take its curl (this is automatically divergence
free) as an expected velocity solution. Check the error convergence for
Taylor-Hood FESystem. What happens for the following choices of finite element
spaces, on a uniformly refined squared mesh?
   1. FE_Q(1)^d-FE_DGQ(0)
   2. FE_Q(2)^d-FE_DGP(1)
   3. FE_Q(2)^d-FE_Q(1)