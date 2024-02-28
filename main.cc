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
#include <memory>
#include <pism/util/pism_utilities.hh>

int main(int argc, char **argv) {

  MPI_Comm com = MPI_COMM_WORLD;
  pism::petsc::Initializer petsc(argc, argv, "");

  int exit_code = 0;
  try {

    auto ctx = pism::context_from_options(com, "");

    auto log = ctx->log();

    pism::options::String input_filename(
        "-input", "name of the file describing the input grid");

    pism::options::String variables("-v", "comma separated list of variables");

    pism::options::String output_grid_filename(
        "-output", "name of the file describing the output grid");

    auto output_grid = pism::Grid::FromFile(ctx, output_grid_filename, {"topg"},
                                            pism::grid::CELL_CENTER);

    {
      pism::File F(ctx->com(), output_grid_filename, pism::io::PISM_GUESS,
                   pism::io::PISM_READONLY);

      output_grid->set_mapping_info(mapping(F, ctx->unit_system()));
    }

    log->message(2, "Output grid:\n");
    output_grid->report_parameters();

    pism::File input_file(ctx->com(), input_filename, pism::io::PISM_GUESS, pism::io::PISM_READONLY);

    pism::array::Scalar target(output_grid, "temporary_storage");

    std::map<std::string, std::shared_ptr<YACInterpolation>> maps;

    for (const auto &v : pism::split(variables, ',')) {
      auto grid_name = YACInterpolation::grid_name(input_file, v, ctx->unit_system());

      auto map = maps[grid_name];
      if (map == nullptr) {
        map = std::make_shared<YACInterpolation>(*output_grid, input_file, v);
        maps[grid_name] = map;
      }
      target.metadata()
        .set_name(v)
        .units(input_file.read_text_attribute(v, "units"));

      map->regrid(input_file, pism::io::Default::Nil(), target);

      auto filename = pism::printf("o_%s.nc", v.c_str());
      target.dump(filename.c_str());
    }
  } catch (...) {
    pism::handle_fatal_errors(com);
    exit_code = 1;
  }

  return exit_code;
}
