#  Lab 07 - Adaptive FEM and shared memory parallelization
## Theory and Practice of Finite Elements

**Luca Heltai** <luca.heltai@sissa.it>

Starter code documentation can be accessed here:

https://dealii-courses.github.io/sissa-mhpc-lab-07/

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

1. Instrument your `Poisson` solver with a `TimerOutput` class, and extract
timing information about each method of the `Poisson` class (see
https://www.dealii.org/current/doxygen/deal.II/classTimerOutput.html) using
scoped timers. Make a couple of test runs, and keep aside both the parameter
file you used to run the code, and the timing results

2. Read the documentation about shared memory parallelization here:
   https://www.dealii.org/current/doxygen/deal.II/group__threads.html

3. Add the parameter
   
    - `Maximum number of threads`
   
which is allowed to take an integer number. A value of `-1` means *choose
automatically*. For any other number, make sure you call `MultithreadInfo::set_thread_limit` with the given argument

4. Make a *coarse granined* parallelization using `Threads::Thread` or
`Threads::ThreadGroup` on each major task of your `Poisson` solver. Experiment
by changing the argument above, and report on the speed up of your code when
using multiple threads, compared with the results you obtained in step 1 above

5. Use `Threads::split_range` to assemble the system in parallel with as many
CPUS as you have available. Check that you do get an improvement

6. Using the default `ScratchData` and `CopyData` objects of deal.II (i.e.,
https://www.dealii.org/current/doxygen/deal.II/classMeshWorker_1_1ScratchData.html
and
https://www.dealii.org/current/doxygen/deal.II/structMeshWorker_1_1CopyData.html
replace the `Threads::split_range` assembly by one based on the use of the
`WorkStream::run` method, and compare the running times with your original
code

7. Document all methods and members of the `Poisson` class, using `Doxygen`
syntax. Make sure that the GitHub Action `documentation.yaml` is working, and
generates the documentation of your code automatically. If everything works,
you should be able to access it as
`https://username.github.io/sissa-mhpc-lab-07-username/` as soon as as the
action has completed (provided that you enable the Pages section of your
repository, and point to the root folder on the `gh-pages` branch, which is
generated automatically by the Doxygen action)