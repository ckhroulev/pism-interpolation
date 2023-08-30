#include <array> // std::array
#include <cassert>
#include <cmath> // M_PI
#include <stdexcept> // std::runtime_error
#include <string>    // std::string
#include <vector>    // std::vector

#include <mpi.h>

#include <proj.h>

extern "C" {
#include "yac_interface.h"
}

/*!
 * Utility class converting `x,y` coordinates in a projection to a `lon,lat`
 * pair.
 *
 * Requires the `PROJ` library.
 */
class LonLatCalculator {
public:
  LonLatCalculator(const std::string &spec) {
    m_pj = proj_create_crs_to_crs(PJ_DEFAULT_CTX, spec.c_str(), "EPSG:4326", 0);
    if (m_pj == 0) {
      throw std::runtime_error(
          "Failed to initialize projection transformation");
    }
  }

  ~LonLatCalculator() { proj_destroy(m_pj); }

  std::array<double, 2> lonlat(double x, double y) const {
    PJ_COORD in, out;

    in.xy = {x, y};
    out = proj_trans(m_pj, PJ_FWD, in);

    return {out.lp.phi, out.lp.lam};
  }

private:
  PJ *m_pj;
};

/*!
 * Grid definition using lon,lat coordinates in radians.
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

    // convert from (row, col) to the linear index lon and lat arrays below
    auto C = [ncol](int row, int col) { return row * ncol + col; };

    // convert from degrees to radians
    auto deg2rad = [](double degree) { return degree * M_PI / 180; };

    LonLatCalculator mapping(projection);

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
 * A regular Cartesian grid in a projected coordinate system.
 *
 * Includes domain decomposition info.
 */
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

// Creates an equally spaced 1D grid of size `N` starting at `start`
// and using grid spacing `step`.
std::vector<double> linspace(double start, double step, int N) {
  std::vector<double> result(N);
  for (int k = 0; k < N; ++k) {
    result[k] = start + k * step;
  }
  return result;
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
int define_grid(const char *grid_name, int n_vertices[2], double *grid_lon,
                double *grid_lat, int n_points[2], double *cell_lon,
                double *cell_lat, int *cell_global_index) {

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

ProjectedGrid source(int rank) {
  auto x = linspace(-678200, 900, 1760);
  auto y = linspace(-3371150, 900, 3040);
  const char *proj = "EPSG:3413";

  ProjectedGrid result(x, y, proj);

  switch (rank) {
  case 0:
    result.set_decomposition(0, 880, 0, 1520);
    break;
  case 1:
    result.set_decomposition(880, 880, 0, 1520);
    break;
  case 2:
    result.set_decomposition(0, 880, 1520, 1520);
    break;
  case 3:
    result.set_decomposition(880, 880, 1520, 1520);
    break;
  }

  return result;
}

ProjectedGrid target(int rank) {

  auto x = linspace(-800000, 5000, 301);
  auto y = linspace(-3400000, 5000, 561);
  const char *proj =
      "+proj=stere +lat_0=90 +lat_ts=71 +lon_0=-39 +k=1 +x_0=0 +y_0=0 "
      "+ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs";

  ProjectedGrid result(x, y, proj);

  switch (rank) {
  case 0: result.set_decomposition(0, 301, 0, 141); break;
  case 1: result.set_decomposition(0, 301, 141, 140); break;
  case 2: result.set_decomposition(0, 301, 281, 140); break;
  case 3: result.set_decomposition(0, 301, 421, 140); break;
  }

  return result;
}

/*!
 * Define the a regular grid. Each PE defines its own subdomain.
 *
 * Returns the point ID that can be used to define a "field".
 */
int define_grid(const ProjectedGrid &info, const std::string &grid_name) {

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
  LonLatGrid cells(x, y, info.projection);
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
  LonLatGrid nodes(x, y, info.projection);
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

  return define_grid(grid_name.c_str(), n_nodes, nodes.lon.data(),
                     nodes.lat.data(), n_cells, cells.lon.data(),
                     cells.lat.data(), cell_global_index.data());
}

/*!
 * Define a "field" on a point set `point_id` in the component `comp_id`.
 *
 * Returns the field ID that can be used to define a "couple".
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

/*!
 * Define an interpolation stack.
 *
 * Returns the interpolation stack ID.
 */
int define_interpolation(double missing_value) {
  int id = 0;
  yac_cget_interp_stack_config(&id);

  int order = 1;
  int enforce_conservation = 1;
  int partial_coverage = 0;

  // conservative
  yac_cadd_interp_stack_config_conservative(
      id, order, enforce_conservation, partial_coverage, YAC_CONSERV_DESTAREA);

  // average over source grid nodes containing a target point, weighted using
  // barycentric
  yac_cadd_interp_stack_config_average(id, YAC_AVG_BARY, partial_coverage);

  // nearest neighbor
  int n_neighbors = 1;
  double scaling = 1.0;
  yac_cadd_interp_stack_config_nnn(id, YAC_NNN_DIST, n_neighbors, scaling);

  // constant if all of the above failed
  yac_cadd_interp_stack_config_fixed(id, missing_value);

  return id;
}

int main(int argc, char **argv) {

  MPI_Init(&argc, &argv);
  {
    MPI_Comm com = MPI_COMM_WORLD;

    int rank = 0;
    MPI_Comm_rank(com, &rank);
    int size = 0;
    MPI_Comm_size(com, &size);

    assert(size == 4);

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
      const int n_comps = 2;
      const char *comp_names[n_comps] = {"input", "output"};
      int comp_ids[n_comps] = {0, 0};
      yac_cdef_comps_instance(instance_id, comp_names, n_comps, comp_ids);

      // Define grids:
      ProjectedGrid input = source(rank);
      int input_grid_id = define_grid(input, "source");

      ProjectedGrid output = target(rank);
      int output_grid_id = define_grid(output, "target");

      // Define fields:
      int source_field = define_field(comp_ids[0], input_grid_id, "source");
      int target_field = define_field(comp_ids[1], output_grid_id, "target");

      // Define the interpolation stack:
      double fill_value = -99999;
      int interp_stack_id = define_interpolation(fill_value);

      // Define the coupling between fields:
      const int src_lag = 0;
      const int tgt_lag = 0;
      const int mapping_side = 1; // 1 means "mapping on source"
      const char *weight_file_name = nullptr;
      yac_cdef_couple_instance(
          instance_id,
          "input",                 // input component name
          "source",                // input grid name
          "source",                // input field name
          "output",                // target component name
          "target",                // target grid name
          "target",                // target field name
          "1",                     // time step length in units below
          YAC_TIME_UNIT_SECOND,    // time step length units
          YAC_REDUCTION_TIME_NONE, // reduction in time (for asynchronous
                                   // coupling)
          interp_stack_id, src_lag, tgt_lag, weight_file_name, mapping_side);

      // free the interpolation stack config now that we defined the coupling
      yac_cfree_interp_stack_config(interp_stack_id);

      yac_cenddef_instance(instance_id);

      int collection_size = 1;
      int ierror = 0;
      {
        std::vector<double> input_array(input.xm * input.ym);
        // the contents of input_array don't matter

        double *send_field_[1] = {input_array.data()};
        double **send_field[1] = {&send_field_[0]};
        int info;
        yac_cput(source_field, collection_size, send_field, &info, &ierror);
      }

      {
        std::vector<double> output_array(output.xm * output.ym);

        double *recv_field[1] = {output_array.data()};
        int info;
        yac_cget(target_field, collection_size, recv_field, &info, &ierror);
      }

      yac_ccleanup_instance(instance_id);
    }
  }
  MPI_Finalize();

  return 0;
}
