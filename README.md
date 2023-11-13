### Flexible interpolation support in PISM

This repository contains the code I wrote to see if we can use
[YAC](https://dkrz-sw.gitlab-pages.dkrz.de/yac/index.html) to add the
ability to perform conservative interpolation "on the fly". This would
reduce the number of pre-processing steps required to set up a simulation.

Requires [PISM](https://github.com/pism/pism), YAC,
[YAXT](https://gitlab.dkrz.de/dkrz-sw/yaxt), and NetCDF. (YAXT is not
used directly but is a dependency of YAC. Similarly, NetCDF is
required by PISM.)
