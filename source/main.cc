#include <deal.II/base/utilities.h>

#include "poisson.h"

/**
 * \mainpage Adaptive FEM for Poisson problem
 *
 * This is the starting code for laboratory number 7 of the course "Theory and
 * Practice of Finite Element methods".
 *
 * We solve the Poisson problem
 *
 * \[
 * \begin{split}
 * -\Delta u &= f \quad \text{ in } \Omega\\
 * n\cdot \nabla u &= g_N \quad \text{ on } \partial\Omega_N\\
 *  u &= g_D \quad \text{ on } \partial\Omega_D
 * \end{split}
 * \]
 *
 * on convex and Lipschitz domains $\Omega$, using the Finite Element method,
 * and the deal.II library (www.dealii.org).
 */

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize init(argc, argv);
  std::string                      par_name = "";
  if (argc > 1)
    par_name = argv[1];

  deallog.depth_console(2);
  Poisson<2> laplace_problem;
  laplace_problem.initialize(par_name);
  laplace_problem.run();
  return 0;
}