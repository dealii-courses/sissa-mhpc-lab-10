#  Lab 07 - Adaptive FEM and shared memory parallelization
## Theory and Practice of Finite Elements

**Luca Heltai** <luca.heltai@sissa.it>

Starter code documentation can be accessed here:

https://dealii-courses.github.io/sissa-mhpc-lab-08/

* * * * *

## General Instructions

For each of the point below, extend the `Poisson` class with functions that
perform the indicated tasks, trying to minimize the amount of code you copy and
paste, possibly restructuring existing code by adding arguments to existing
functions, and generating wrappers similar to the `run` method (e.g.,
`run_exercise_3`).

Once you created a function that performs the given task, add it to the
`poisson-tester.cc` file, and make sure all the exercises are run through the
`gtest` executable, e.g., adding a test for each exercise, as in the following
snippet:

```C++
TEST_F(PoissonTester, Exercise3) {
   run_exercise_3();
}
```

By the end of this laboratory, you will have modified your Poisson code to run
in parallel using shared memory parallelization on multiple threads, and you
will have some knowledge of Task based parallelization
## Lab-07

1. Repeat the basic MPI commands from the file mpihello/main.cc and understand
how things work

2. Read carefully the documentation of `step-40`, and adapt your poisson solver
to become a parallel distributed solver

3. Make the solver run in "hypbrid" parallelization mode, i.e., according the
number of threads you specified in the parameter file, assemble things in
parallel using both MPI and multithread

2. Similar to what is shown in the lecture, visualize the view of the mesh from
each individual processor using ``GridOut::write_vtk`` and the "global" mesh.
Use 3 MPI tasks for this

3. Create a simple mesh (hyper_cube refined twice globally), run with two MPI
tasks and print locally owned, locally active, and locally relevant
``IndexSet`` for each task

4. Switch to release mode (``make release``), decide on a global refinement
level that takes in the order of 30-60 seconds to solve, and study assembly and
solve time with 1,2,4,8,12,16 MPI tasks. Which is the fastest, do the timings
make sense based on how many cores your machine has?

5. Play with the test problem by switching to 3d and changing the geometry to
something interesting. Your choice!