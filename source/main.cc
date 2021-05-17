#include <deal.II/base/utilities.h>

#include "base_problem.h"
#include "linear_elasticity.h"
#include "poisson.h"

/**
 * \mainpage Adaptive FEM
 *
 * This is the starting code for laboratory number 10 of the course "Theory and
 * Practice of Finite Element methods".
 */
int
main(int argc, char **argv)
{
  const std::string program_name(argv[0]);
  if (program_name.find("poisson") != std::string::npos)
    return run<Poisson<2>>(argc, argv);
  if (program_name.find("linear_elasticity") != std::string::npos)
    return run<LinearElasticity<2>>(argc, argv);
}