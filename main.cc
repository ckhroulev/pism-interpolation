#include <mpi.h>

#include "pism/util/petscwrappers/PetscInitializer.hh"
#include "pism/util/Context.hh"
#include "pism/util/ConfigInterface.hh"
#include "pism/util/Logger.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/pism_options.hh"

extern "C" {
#include "yac_interface.h"
}

int main(int argc, char **argv) {

  MPI_Comm com = MPI_COMM_WORLD;
  pism::petsc::Initializer petsc(argc, argv, "");

  int exit_code = 0;
  try {

    auto ctx = pism::context_from_options(com, "");

    auto config = ctx->config();

    pism::options::String input("-input", "name of the file describing the input grid");

    auto log = ctx->log();

    log->message(2, "Reading from %s...\n", input->c_str());

    {
      int instance_id = 0;
      yac_cinit_instance(&instance_id);

      yac_ccleanup_instance(instance_id);
    }

  } catch (...) {
    pism::handle_fatal_errors(com);
    exit_code = 1;
  }

  return exit_code;
}
