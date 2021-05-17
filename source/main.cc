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
  try
    {
      using namespace dealii;

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      std::string                      par_name = "";
      if (argc > 1)
        par_name = argv[1];

      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        deallog.depth_console(2);
      else
        deallog.depth_console(0);

      LinearElasticity<2> base_problem;
      base_problem.initialize(par_name);
      base_problem.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
