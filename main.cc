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



int create_field(int instance_id, const char* comp_name, int comp_id, const char *field_name) {

  int cyclic[] = {0, 0};
  int n_vertices[] = {2, 2};
  double grid_x[] = {0, 1, 0, 1};
  double grid_y[] = {-1, -1, 1, 1};
  int grid_id = 0;

  const int n_cells = 1;
  int n_points[] = {1, 1};
  double points_x[n_cells] = {0.5};
  double points_y[n_cells] = {0.0};

  yac_cdef_grid_curve2d(comp_name, n_vertices, cyclic, grid_x, grid_y, &grid_id);

  int point_id = 0;
  yac_cdef_points_curve2d(grid_id, n_points, YAC_LOCATION_CELL, points_x,
                          points_y, &point_id);

  int point_set_size = 1;
  int collection_size = 1;

  int target_field_id = 0;
  const char *time_step_length = "1";
  yac_cdef_field(field_name, comp_id, &point_id, point_set_size,
                 collection_size, time_step_length, YAC_TIME_UNIT_SECOND,
                 &target_field_id);

  return target_field_id;
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
      yac_cdef_calendar(YAC_PROLEPTIC_GREGORIAN);

      int interp_stack_id = 0;
      yac_cget_interp_stack_config(&interp_stack_id);
      // add average
      yac_cadd_interp_stack_config_average(interp_stack_id, YAC_AVG_DIST, 0);
      // constant if averaging failed
      yac_cadd_interp_stack_config_fixed(interp_stack_id, -9999);

      int src_lag = 0;
      int tgt_lag = 0;
      int mapping_side = 1;
      yac_cdef_couple_instance(
          instance_id,
          "input", "input", "source",
          "output", "output", "target",
          "1", YAC_TIME_UNIT_SECOND, YAC_REDUCTION_TIME_NONE,
          interp_stack_id, src_lag, tgt_lag, nullptr, mapping_side);

      /*
        Setup:

        1. Initialize an instance
        2. Define a component.
        3. Define a grid.
        4. Define points on the grid.
        5. Define a field using these points.

        Once this is done, use yac_cput() to send and yac_cget() to
        receive data.
        */

      // Create components: this has to be done using *one* call:
      const char *comp_names[2] = {"input", "output"};
      int comp_ids[2] = {0, 0};
      yac_cdef_comps_instance(instance_id, comp_names, 2, comp_ids);
      int input_id = comp_ids[0];
      int output_id = comp_ids[1];

      int input_field = create_field(instance_id, "input", input_id, "source");
      int target_field = create_field(instance_id, "output", output_id, "target");

      yac_cenddef_instance(instance_id);

      int collection_size = 1;
      int ierror = 0;
      {
        double send_field_data[2] = {3.0, 4.0};
        double *send_field_[1] = {&send_field_data[0]};
        double **send_field[1] = {&send_field_[0]};
        int info;
        yac_cput(input_field, collection_size, send_field, &info, &ierror);
      }

      {
        double recv_field_data[2] = {-1.0, -1.0};
        double *recv_field[1] = {&recv_field_data[0]};
        int info;
        yac_cget(target_field, collection_size, recv_field, &info, &ierror);
      }

      yac_ccleanup_instance(instance_id);
    }

  } catch (...) {
    pism::handle_fatal_errors(com);
    exit_code = 1;
  }

  return exit_code;
}
