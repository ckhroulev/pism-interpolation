#include <mpi.h>

#include "pism/util/ConfigInterface.hh"
#include "pism/util/Context.hh"
#include "pism/util/Logger.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/petscwrappers/PetscInitializer.hh"
#include "pism/util/pism_options.hh"

extern "C" {
#include "yac_interface.h"
}

int define_input_field(int instance_id, const char *comp_name, int comp_id,
                       const char *field_name) {

  int cyclic[] = {0, 0};
  int n_vertices[] = {3, 3};
  double grid_x[] = {-1, 0, 1, -1, 0, 1, -1, 0, 1};
  double grid_y[] = {-1, -1, -1, 0, 0, 0, 1, 1, 1};
  int grid_id = 0;

  const int n_cells = 4;
  int n_points[] = {2, 2};
  double points_x[n_cells] = {-0.5, 0.5, -0.5, 0.5};
  double points_y[n_cells] = {-0.5, -0.5, 0.5, 0.5};

  yac_cdef_grid_curve2d(comp_name, n_vertices, cyclic, grid_x, grid_y,
                        &grid_id);

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

int define_target_field(int instance_id, const char *comp_name, int comp_id,
                        const char *field_name) {

  int cyclic[] = {0, 0};
  int n_vertices[] = {2, 2};
  double grid_x[] = {-0.5, 0.5, -0.5, 0.5};
  double grid_y[] = {-0.5, -0.5, 0.5, 0.5};
  int grid_id = 0;

  const int n_cells = 1;
  int n_points[] = {1, 1};
  double points_x[n_cells] = {0};
  double points_y[n_cells] = {0};

  yac_cdef_grid_curve2d(comp_name, n_vertices, cyclic, grid_x, grid_y,
                        &grid_id);

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

    pism::options::String input("-input",
                                "name of the file describing the input grid");

    auto log = ctx->log();

    log->message(2, "Reading from %s...\n", input->c_str());

    {
      // Initialize an instance:
      int instance_id = 0;
      {
        yac_cinit_instance(&instance_id);
        yac_cdef_calendar(YAC_PROLEPTIC_GREGORIAN);
        // Note: zero-padding of months and days *is* required.
        yac_cdef_datetime_instance(instance_id, "-1-01-01", "+1-01-01");
      }

      // Define components: this has to be done using *one* call
      // (cannot call yac_cdef_comp?_instance() more than once)
      const char *comp_names[2] = {"input", "output"};
      int comp_ids[2] = {0, 0};
      yac_cdef_comps_instance(instance_id, comp_names, 2, comp_ids);
      int input_id = comp_ids[0];
      int output_id = comp_ids[1];

      // Define fields:
      int input_field =
          define_input_field(instance_id, "input", input_id, "source");
      int target_field =
          define_target_field(instance_id, "output", output_id, "target");

      // Define the interpolation stack:
      int interp_stack_id = 0;
      {
        yac_cget_interp_stack_config(&interp_stack_id);
        // add average
        yac_cadd_interp_stack_config_average(interp_stack_id, YAC_AVG_DIST,
                                             0 /* partial_coverage = false */);
        yac_cadd_interp_stack_config_nnn(interp_stack_id, YAC_NNN_AVG,
                                         1 /* one nearest neighbor */,
                                         1.0 /* scaling */);
        // constant if averaging failed
        yac_cadd_interp_stack_config_fixed(interp_stack_id, -9999);
      }

      // Define the coupling between fields:
      int src_lag = 0;
      int tgt_lag = 0;
      int mapping_side = 1;
      yac_cdef_couple_instance(
          instance_id,
          "input",              // input component name
          "input",              // input grid name
          "source",             // input field name
          "output",             // target component name
          "output",             // target grid name
          "target",             // target field name
          "1",                  // time step length in units below
          YAC_TIME_UNIT_SECOND, // time step length units
          YAC_REDUCTION_TIME_NONE, // reduction in time (for asynchronous coupling)
          interp_stack_id,
          src_lag,
          tgt_lag,
          nullptr,              // weight file name (disabled)
          mapping_side);

      yac_cenddef_instance(instance_id);

      int collection_size = 1;
      int ierror = 0;
      {
        double send_field_data[4] = {1.0, 2.0, 3.0, 4.0};
        double *send_field_[1] = {&send_field_data[0]};
        double **send_field[1] = {&send_field_[0]};
        int info;
        yac_cput(input_field, collection_size, send_field, &info, &ierror);
      }

      {
        double recv_field_data[1] = {-1.0};
        double *recv_field[1] = {&recv_field_data[0]};
        int info;
        yac_cget(target_field, collection_size, recv_field, &info, &ierror);

        log->message(2, "Got %f\n", recv_field_data[0]);
      }

      yac_ccleanup_instance(instance_id);
    }

  } catch (...) {
    pism::handle_fatal_errors(com);
    exit_code = 1;
  }

  return exit_code;
}
