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

/*!
 * Define a curvilinear Cartesian 2D grid on a sphere.
 *
 * `grid_x` and `grid_y` should define cell vertices (cell bounds)
 *
 * `point_x` and `point_y` correspond to cell centers (for grids used
 * for conservative interpolation). Set to `nullptr` to define a grid
 * using cell vertices as the point set.
 */
static int define_grid(const char *grid_name, int n_vertices[2], double *grid_x,
                       double *grid_y, int n_points[2], double *point_x,
                       double *point_y) {

  int cyclic[] = {0, 0};
  int grid_id = 0;

  yac_cdef_grid_curve2d(grid_name, n_vertices, cyclic, grid_x, grid_y,
                        &grid_id);

  int point_id = 0;
  if (n_points != nullptr and point_x != nullptr and point_y != nullptr) {
    yac_cdef_points_curve2d(grid_id, n_points, YAC_LOCATION_CELL, point_x,
                            point_y, &point_id);
  } else {
    yac_cdef_points_curve2d(grid_id, n_vertices, YAC_LOCATION_CORNER, grid_x,
                            grid_y, &point_id);
  }

  return point_id;
}

static int define_field(int comp_id, int point_id, const char *field_name) {

  const char *time_step_length = "1";
  const int point_set_size = 1;
  const int collection_size = 1;

  int field_id = 0;
  yac_cdef_field(field_name, comp_id, &point_id, point_set_size,
                 collection_size, time_step_length, YAC_TIME_UNIT_SECOND,
                 &field_id);

  return field_id;
}

int define_input_field(const char *comp_name, int comp_id,
                       const char *grid_name, const char *field_name) {

  int n_vertices[] = {3, 3};
  double grid_x[] = {-1, 0, 1, -1, 0, 1, -1, 0, 1};
  double grid_y[] = {-1, -1, -1, 0, 0, 0, 1, 1, 1};

  int n_points[] = {2, 2};
  double points_x[] = {-0.5, 0.5, -0.5, 0.5};
  double points_y[] = {-0.5, -0.5, 0.5, 0.5};

  int point_id = define_grid(grid_name, n_vertices, grid_x, grid_y, n_points,
                             points_x, points_y);

  return define_field(comp_id, point_id, field_name);
}

int define_target_field(const char *comp_name, int comp_id,
                        const char *grid_name, const char *field_name) {

  int n_vertices[] = {2, 2};
  double grid_x[] = {-0.25, 0.25, -0.25, 0.25};
  double grid_y[] = {-0.25, -0.25, 0.25, 0.25};

  int point_id = define_grid(grid_name, n_vertices, grid_x, grid_y, nullptr,
                             nullptr, nullptr);

  return define_field(comp_id, point_id, field_name);
}

static int define_average_interpolation(double missing_value) {
  int id = 0;
  yac_cget_interp_stack_config(&id);
  // add average
  yac_cadd_interp_stack_config_average(id, YAC_AVG_DIST,
                                       0 /* partial_coverage = false */);
  // add nearest neighbor
  yac_cadd_interp_stack_config_nnn(
      id, YAC_NNN_DIST, 1 /* one nearest neighbor */, 1.0 /* scaling */);
  // constant if the above failed
  yac_cadd_interp_stack_config_fixed(id, missing_value);

  return id;
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
      const char *comp_name = "interpolation";
      int comp_id = 0;
      yac_cdef_comp_instance(instance_id, comp_name, &comp_id);

      // Define fields:
      int input_field =
          define_input_field(comp_name, comp_id, "source", "source");
      int target_field =
          define_target_field(comp_name, comp_id, "target", "target");

      // Define the interpolation stack:
      int interp_stack_id = define_average_interpolation(-9999.0);

      // Define the coupling between fields:
      const int src_lag = 0;
      const int tgt_lag = 0;
      const int mapping_side = 1; // 1 means "mapping on source"
      yac_cdef_couple_instance(
          instance_id,
          "interpolation",      // input component name
          "source",             // input grid name
          "source",             // input field name
          "interpolation",      // target component name
          "target",             // target grid name
          "target",             // target field name
          "1",                  // time step length in units below
          YAC_TIME_UNIT_SECOND, // time step length units
          YAC_REDUCTION_TIME_NONE, // reduction in time (for asynchronous coupling)
          interp_stack_id,
          src_lag,
          tgt_lag,
          nullptr,              // weight file name (disabled)
          mapping_side);

      // free the interpolation stack config now that we defined the coupling
      yac_cfree_interp_stack_config(interp_stack_id);

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
        const int N = 4;
        double recv_field_data[N] = {-1.0, -1.0, -1.0, -1.0};
        double *recv_field[1] = {&recv_field_data[0]};
        int info;
        yac_cget(target_field, collection_size, recv_field, &info, &ierror);

        for (int k = 0; k < N; ++k) {
          log->message(2, "Got %f\n", recv_field_data[k]);
        }
      }

      yac_ccleanup_instance(instance_id);
    }

  } catch (...) {
    pism::handle_fatal_errors(com);
    exit_code = 1;
  }

  return exit_code;
}
