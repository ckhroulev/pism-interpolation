#include <cassert>
#include <memory>
#include <mpi.h>
#include <pism/util/array/Array.hh>
#include <string>
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
 * Grid definition using coordinates in radians.
 */
struct LonLatGrid {
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
  LonLatGrid(const std::vector<double> &x, const std::vector<double> &y,
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

struct ProjectedGrid {

  // grid coordinates, in meters
  std::vector<double> x, y;

  // PROJ string defining the projection
  std::string projection;

  // Domain decomposition

  // first x index of the local subdomain
  int xs;
  // number of grid points in the x direction in the local subdomain
  int xm;
  // first y index of the local subdomain
  int ys;
  // number of grid points in the y direction in the local subdomain
  int ym;

  ProjectedGrid(const std::vector<double> &x_, const std::vector<double> &y_,
                const std::string &projection_) {
    x = x_;
    y = y_;
    projection = projection_;
    xs = 0;
    xm = x.size();
    ys = 0;
    ym = y.size();
  }

  void set_decomposition(int xs_, int xm_, int ys_, int ym_) {
    xs = xs_;
    xm = xm_;
    ys = ys_;
    ym = ym_;
  }
};

/*!
 * Define the PISM grid. Each PE defines its own subdomain.
 *
 * Returns the point ID that can be used to define a "field".
 */
static int define_grid(std::shared_ptr<const pism::Grid> grid,
                       const std::string &grid_name,
                       const std::string &projection) {

  ProjectedGrid info(grid->x(), grid->y(), projection);
  info.set_decomposition(grid->xs(), grid->xm(), grid->ys(), grid->ym());

  std::vector<double> x(info.xm);
  std::vector<double> y(info.ym);

  // Set x and y to coordinates of cell centers:
  {
    for (int k = 0; k < x.size(); ++k) {
      x[k] = info.x[info.xs + k];
    }
    for (int k = 0; k < y.size(); ++k) {
      y[k] = info.y[info.ys + k];
    }
  }

  // Compute lon,lat coordinates of cell centers:
  LonLatGrid cells(x, y, projection);
  int n_cells[2] = {(int)x.size(), (int)y.size()};

  // Shift x and y by half a grid spacing and add one more row and
  // column to get coordinates of cell corners:
  {
    double dx = info.x[1] - info.x[0];
    double dy = info.y[1] - info.y[0];

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
  LonLatGrid nodes(x, y, projection);
  int n_nodes[2] = {(int)x.size(), (int)y.size()};

  std::vector<int> cell_global_index(n_cells[0] * n_cells[1]);
  {
    int Mx = info.x.size();
    int k = 0;
    for (int j = info.ys; j < info.ys + info.ym; ++j) {
      for (int i = info.xs; i < info.xs + info.xm; ++i) {
        cell_global_index[k] = j * Mx + i;
        ++k;
      }
    }
  }

  int point_id = 0;
  {
    int cyclic[] = {0, 0};

    int grid_id = 0;

    yac_cdef_grid_curve2d(grid_name.c_str(), n_nodes, cyclic, nodes.lon.data(),
                          nodes.lat.data(), &grid_id);

    yac_cset_global_index(cell_global_index.data(), YAC_LOCATION_CELL, grid_id);

    yac_cdef_points_curve2d(grid_id, n_cells, YAC_LOCATION_CELL,
                            cells.lon.data(), cells.lat.data(), &point_id);
  }
  return point_id;
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

static int interpolation_fine_to_coarse(double missing_value) {
  int id = 0;
  yac_cget_interp_stack_config(&id);

  int order = 1;
  int enforce_conservation = 1;
  int partial_coverage = 0;

  // conservative
  yac_cadd_interp_stack_config_conservative(
      id, order, enforce_conservation, partial_coverage, YAC_CONSERV_DESTAREA);

  // average over source grid nodes containing a target point, weighted using
  // barycentric local coordinates
  yac_cadd_interp_stack_config_average(id, YAC_AVG_BARY, partial_coverage);

  // nearest neighbor
  int n_neighbors = 1;
  double scaling = 1.0;
  yac_cadd_interp_stack_config_nnn(id, YAC_NNN_DIST, n_neighbors, scaling);

  // constant if all of the above failed
  yac_cadd_interp_stack_config_fixed(id, missing_value);

  return id;
}

static int interpolation_coarse_to_fine(double missing_value) {
  int id = 0;
  yac_cget_interp_stack_config(&id);

  int partial_coverage = 0;

  // average over source grid nodes containing a target point, weighted using
  // barycentric local coordinates
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

  return input_file.read_text_attribute("PISM_GLOBAL", "proj");
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
      {
        yac_cinit();
        yac_cdef_calendar(YAC_PROLEPTIC_GREGORIAN);
        // Note: zero-padding of months and days *is* required.
        yac_cdef_datetime("-1-01-01", "+1-01-01");
      }

      // Define components: this has to be done using *one* call
      // (cannot call yac_cdef_comp?_instance() more than once)
      const int n_comps = 2;
      const char *comp_names[n_comps] = {"input", "output"};
      int comp_ids[n_comps] = {0, 0};
      yac_cdef_comps(comp_names, n_comps, comp_ids);

      // Define grids:
      int input_grid_id = define_grid(input_grid, "source", input_projection);
      int output_grid_id =
          define_grid(output_grid, "target", output_projection);

      // Define fields:
      int source_field = define_field(comp_ids[0], input_grid_id, "source");
      int target_field = define_field(comp_ids[1], output_grid_id, "target");

      // Define the interpolation stack:
      double fill_value = -99999;
      int interp_stack_id = 0;
      if (input_grid->dx() < output_grid->dx() or
          input_grid->dy() < output_grid->dy()) {
        interp_stack_id = interpolation_fine_to_coarse(fill_value);
      } else {
        interp_stack_id = interpolation_coarse_to_fine(fill_value);
      }

      // Define the coupling between fields:
      const int src_lag = 0;
      const int tgt_lag = 0;
      yac_cdef_couple("input",              // input component name
                      "source",             // input grid name
                      "source",             // input field name
                      "output",             // target component name
                      "target",             // target grid name
                      "target",             // target field name
                      "1",                  // time step length in units below
                      YAC_TIME_UNIT_SECOND, // time step length units
                      YAC_REDUCTION_TIME_NONE, // reduction in time (for
                                               // asynchronous coupling)
                      interp_stack_id, src_lag, tgt_lag);

      // free the interpolation stack config now that we defined the coupling
      yac_cfree_interp_stack_config(interp_stack_id);

      double start = 0;
      double end = 0;

      log->message(2, "Initializing interpolation... ");
      start = MPI_Wtime();
      yac_cenddef();
      end = MPI_Wtime();
      log->message(2, "done in %f seconds.\n", end - start);

      int collection_size = 1;
      int ierror = 0;

      pism::array::Scalar input_field(input_grid, "topg");
      input_field.metadata().units("m").standard_name("bedrock_altitude");
      input_field.regrid(input, pism::io::Default::Nil());

      pism::array::Scalar output(output_grid, "topg");
      output.metadata().units("m").standard_name("bedrock_altitude");
      output.metadata()["_FillValue"] = {fill_value};
      {
        int send_info = 0;
        int recv_info = 0;

        pism::petsc::VecArray input_array(input_field.vec());
        pism::petsc::VecArray output_array(output.vec());

        double *send_field_ = input_array.get();
        double **send_field[1] = {&send_field_};

        double *recv_field[1] = {output_array.get()};

        start = MPI_Wtime();
        yac_cexchange(source_field, target_field, collection_size, send_field,
                      recv_field, &send_info, &recv_info, &ierror);
        end = MPI_Wtime();
      }

      output.dump(output_file->c_str());

      log->message(2, "Data transfer took %f seconds.\n", end - start);

      yac_cfinalize();
    }

  } catch (...) {
    pism::handle_fatal_errors(com);
    exit_code = 1;
  }

  return exit_code;
}
