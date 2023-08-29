#include <cassert>
#include <memory>
#include <mpi.h>
#include <pism/util/array/Array.hh>
#include <vector>

#include "pism/util/ConfigInterface.hh"
#include "pism/util/Context.hh"
#include "pism/util/Grid.hh"
#include "pism/util/Logger.hh"
#include "pism/util/Units.hh"
#include "pism/util/array/Scalar.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/File.hh"
#include "pism/util/io/IO_Flags.hh"
#include "pism/util/petscwrappers/PetscInitializer.hh"
#include "pism/util/petscwrappers/Vec.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/projection.hh"

extern "C" {
#include "yac_interface.h"
}

/*!
 * Define a curvilinear Cartesian 2D grid on a sphere.
 *
 * `grid_lon` and `grid_lat` should define cell vertices (cell bounds).
 *
 * `cell_lon` and `cell_lat` correspond to cell centers (for grids used
 * for conservative interpolation). Set to `nullptr` to define a grid
 * using cell vertices as the point set.
 *
 * All coordinates are in radians.
 *
 * Returns the point ID that can be used to define a "field".
 */
static int define_grid(const char *grid_name, int n_vertices[2],
                       double *grid_lon, double *grid_lat, int n_points[2],
                       double *cell_lon, double *cell_lat, int *cell_global_index) {



  int cyclic[] = {0, 0};
  int grid_id = 0;

  yac_cdef_grid_curve2d(grid_name, n_vertices, cyclic, grid_lon, grid_lat,
                        &grid_id);

  yac_cset_global_index(cell_global_index, YAC_LOCATION_CELL, grid_id);

  assert(n_points != nullptr and cell_lon != nullptr and cell_lat != nullptr);

  int point_id = 0;
  yac_cdef_points_curve2d(grid_id, n_points, YAC_LOCATION_CELL, cell_lon,
                          cell_lat, &point_id);

  return point_id;
}

/*!
 * Grid definition using coordinates in radians.
 */
struct Grid {
  std::vector<double> lon;
  std::vector<double> lat;

  /*!
   *
   * Converts a Cartesian grid in a `projection` that uses coordinates
   * `x` and `y` in meters into the form that can be used to define a
   * curvilinear grid in YAC.
   *
   * The `projection` string has to use the format compatible with PROJ.
   */
  Grid(const std::vector<double> &x, const std::vector<double> &y,
       const std::string &projection) {

    int nrow = y.size();
    int ncol = x.size();
    int N = nrow * ncol;

    lon.resize(N);
    lat.resize(N);

    // convert from (row, col) to the linear index in "cell" arrays
    auto C = [ncol](int row, int col) { return row * ncol + col; };

    // convert from degrees to radians
    auto deg2rad = [](double degree) { return degree * M_PI / 180; };

    pism::LonLatCalculator mapping(projection);

    for (int row = 0; row < nrow; ++row) {
      for (int col = 0; col < ncol; ++col) {
        auto coords = mapping.lonlat(x[col], y[row]);

        lon[C(row, col)] = deg2rad(coords[0]);
        lat[C(row, col)] = deg2rad(coords[1]);
      }
    }
  }
};

/*!
 * Define the PISM grid. Each PE defines its own subdomain.
 */
static int define_grid(std::shared_ptr<const pism::Grid> grid, const std::string &grid_name,
                       const std::string &projection) {

  std::vector<double> x(grid->xm());
  std::vector<double> y(grid->ym());

  // Set x and y to coordinates of cell centers:
  {
    for (int k = 0; k < x.size(); ++k) {
      x[k] = grid->x(grid->xs() + k);
    }
    for (int k = 0; k < y.size(); ++k) {
      y[k] = grid->y(grid->ys() + k);
    }
  }

  // Compute lon,lat coordinates of cell centers:
  Grid cells(x, y, projection);
  int n_cells[2] = {(int)x.size(), (int)y.size()};

  // Shift x and y by half a grid spacing and add one more row and
  // column to get coordinates of cell corners:
  {
    double dx = grid->dx();
    double dy = grid->dy();

    double x_last = x.back() + 0.5 * dx;
    for (int k = 0; k < x.size(); ++k) {
      x[k] -= 0.5 * dx;
    }
    x.push_back(x_last);

    double y_last = y.back() + 0.5 * dy;
    for (int k = 0; k < y.size(); ++k) {
      y[k] -= 0.5 * dy;
    }
    y.push_back(y_last);
  }

  // Compute lon,lat coordinates of cell corners:
  Grid nodes(x, y, projection);
  int n_nodes[2] = {(int)x.size(), (int)y.size()};

  std::vector<int> cell_global_index(n_cells[0] * n_cells[1]);
  {
    int k = 0;
    for (auto p = grid->points(); p; p.next()) {
      int i = p.i(), j = p.j();
      cell_global_index[k] = j * grid->Mx() + i;
      ++k;
    }
  }

  return define_grid(grid_name.c_str(), n_nodes, nodes.lon.data(),
                     nodes.lat.data(), n_cells, cells.lon.data(),
                     cells.lat.data(), cell_global_index.data());
}

/*!
 * Define a "field" on a point set `point_id` in the component `comp_id`.
 */
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

static int define_interpolation(double missing_value) {
  int id = 0;
  yac_cget_interp_stack_config(&id);

  int order = 1;
  int enforce_conservation = 1;
  int partial_coverage = 0;

  // conservative
  yac_cadd_interp_stack_config_conservative(
      id, order, enforce_conservation, partial_coverage, YAC_CONSERV_DESTAREA);

  // average over source grid nodes containing a target point, weighted using barycentric
  yac_cadd_interp_stack_config_average(id, YAC_AVG_BARY, partial_coverage);

  // nearest neighbor
  int n_neighbors = 1;
  double scaling = 1.0;
  yac_cadd_interp_stack_config_nnn(id, YAC_NNN_DIST, n_neighbors, scaling);

  // constant if all of the above failed
  yac_cadd_interp_stack_config_fixed(id, missing_value);

  return id;
}

//! Get projection info from a NetCDF file. Returns a PROJ string.
std::string projection(MPI_Comm com, const std::string &filename,
                       pism::units::System::Ptr sys) {
  pism::File input_file(com, filename, pism::io::PISM_GUESS,
                        pism::io::PISM_READONLY);

  return pism::get_projection_info(input_file, "mapping", sys).proj;
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

    pism::options::String output_file("-o", "name of the output file");

    pism::options::String output("-output",
                                 "name of the file describing the output grid");

    auto log = ctx->log();

    auto input_grid = pism::Grid::FromFile(ctx, input, {"topg", "thk"},
                                           pism::grid::CELL_CENTER);

    log->message(2, "\nInput grid read from %s:\n", input->c_str());
    input_grid->report_parameters();

    auto output_grid = pism::Grid::FromFile(ctx, output, {"topg", "thk"},
                                            pism::grid::CELL_CENTER);
    log->message(2, "\nOutput grid read from %s:\n", output->c_str());
    output_grid->report_parameters();

    auto input_projection = projection(ctx->com(), input, ctx->unit_system());
    auto output_projection = projection(ctx->com(), output, ctx->unit_system());

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

      // Define grids:
      int input_grid_id = define_grid(input_grid, "source", input_projection);
      int output_grid_id =
          define_grid(output_grid, "target", output_projection);

      // Define fields:
      int source_field = define_field(comp_id, input_grid_id, "source");
      int target_field = define_field(comp_id, output_grid_id, "target");

      // Define the interpolation stack:
      double fill_value = -99999;
      int interp_stack_id = define_interpolation(fill_value);

      // Define the coupling between fields:
      const int src_lag = 0;
      const int tgt_lag = 0;
      const int mapping_side = 1; // 1 means "mapping on source"
      yac_cdef_couple_instance(
          instance_id,
          "interpolation",         // input component name
          "source",                // input grid name
          "source",                // input field name
          "interpolation",         // target component name
          "target",                // target grid name
          "target",                // target field name
          "1",                     // time step length in units below
          YAC_TIME_UNIT_SECOND,    // time step length units
          YAC_REDUCTION_TIME_NONE, // reduction in time (for asynchronous
                                   // coupling)
          interp_stack_id, src_lag, tgt_lag,
          nullptr, // weight file name (disabled)
          mapping_side);

      // free the interpolation stack config now that we defined the coupling
      yac_cfree_interp_stack_config(interp_stack_id);

      double start = 0;
      double end = 0;

      log->message(2, "Initializing interpolation... ");
      start = MPI_Wtime();
      yac_cenddef_instance(instance_id);
      end = MPI_Wtime();
      log->message(2, "done in %f seconds.\n", end - start);

      int collection_size = 1;
      int ierror = 0;
      {
        pism::array::Scalar input_field(input_grid, "topg");
        input_field.metadata().units("m").standard_name("bedrock_altitude");
        input_field.regrid(input, pism::io::CRITICAL);

        pism::petsc::VecArray array(input_field.vec());

        double *send_field_[1] = {array.get()};
        double **send_field[1] = {&send_field_[0]};
        int info;
        start = MPI_Wtime();
        yac_cput(source_field, collection_size, send_field, &info, &ierror);
      }
      {
        pism::array::Scalar output(output_grid, "bed_topography");
        output.metadata().units("m").standard_name("bedrock_altitude");
        output.metadata()["_FillValue"] = {fill_value};
        pism::petsc::VecArray array(output.vec());

        double *recv_field[1] = {array.get()};
        int info;
        yac_cget(target_field, collection_size, recv_field, &info, &ierror);
        end = MPI_Wtime();
        output.dump(output_file->c_str());
      }

      log->message(2, "Data transfer took %f seconds.\n", end - start);

      yac_ccleanup_instance(instance_id);
    }

  } catch (...) {
    pism::handle_fatal_errors(com);
    exit_code = 1;
  }

  return exit_code;
}
