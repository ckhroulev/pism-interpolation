#include "pism/util/Context.hh"
#include "pism/util/Grid.hh"
#include "pism/util/Logger.hh"
#include "pism/util/array/Scalar.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/io/File.hh"
#include "pism/util/io/IO_Flags.hh"
#include "pism/util/petscwrappers/PetscInitializer.hh"
#include "pism/util/pism_options.hh"

#include "YACInterpolation.hh"

int main(int argc, char **argv) {

  MPI_Comm com = MPI_COMM_WORLD;
  pism::petsc::Initializer petsc(argc, argv, "");

  double fill_value = -99999;

  int exit_code = 0;
  try {

    auto ctx = pism::context_from_options(com, "");

    pism::options::String input_filename(
        "-input", "name of the file describing the input grid");

    pism::options::String output_filename("-o", "name of the output file");

    pism::options::String output_grid_filename(
        "-output", "name of the file describing the output grid");

    auto output_grid = pism::Grid::FromFile(ctx, output_grid_filename, {"topg"},
                                            pism::grid::CELL_CENTER);

    {
      pism::File output_file(ctx->com(), output_grid_filename, pism::io::PISM_GUESS,
                            pism::io::PISM_READONLY);

      output_grid->set_mapping_info(mapping(output_file, ctx->unit_system()));
    }

    auto log = ctx->log();

    pism::File input_file(ctx->com(), input_filename, pism::io::PISM_GUESS, pism::io::PISM_READONLY);

    YACInterpolation interp(*output_grid, input_file, "topg");

    log->message(2, "Output:\n");
    output_grid->report_parameters();

    pism::array::Scalar target(output_grid, "topg");
    target.metadata().units("m").standard_name("bedrock_altitude");
    target.metadata()["_FillValue"] = {fill_value};

    interp.regrid(input_file, pism::io::Default::Nil(), target);

    target.dump(output_filename->c_str());

  } catch (...) {
    pism::handle_fatal_errors(com);
    exit_code = 1;
  }

  return exit_code;
}
