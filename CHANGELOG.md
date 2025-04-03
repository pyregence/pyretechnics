# Version 2025.4.3

## Changes for Developers
- Removed deprecated `distutils` dependency in `setup.py`.
- Be explicit about all Numpy buffer types created in `random.py`.
- Only apply custom compiler flags when performing builds under `make.sh`.

## Changes for Users
- This release should work on a Windows system after being compiled
  with the Microsoft C++ compiler (cl.exe) and should respect the
  `CFLAGS` settings on computers that have set `CC="gcc"`.

# Version 2025.4.2

## Changes for Developers
- Restrict custom compiler flags to only be used when distutils.sysconfig.get_config_var("CC") == "gcc".
- Removed project.license entry from pyproject.toml to eliminate conflict between pypa-build and setuptools.

## Changes for Users
- This release should be C compiler agnostic when building wheels.

# Version 2025.3.27

## Changes for Developers
- Updated the `license` field definition in `pyproject.toml` to conform to PEP 639.
- Updated the `upload-testpypi` and `upload-pypi` commands in `make.sh` to only upload source builds.
- Updated object types in *.py and *.pxd files to more specific Python types where possible.
- Fixed optional typing mistake in burn_cells.pxd.
- Corrected Numpy ndarray types in all *.py and *.pxd files.
- Replaced all `x == None` expressions with `x is None`.
- Removed unnecessary `@cy.exceptval` decorators.
- Addressed index performance hints in `space_time_cube.py` by type-hinting more `SpaceTimeCube` methods.

## Changes for Users
- Added a `SpreadState` class to `pyretechnics.eulerian_level_set`.
- Updated `spread_fire_phi_field` to take and return a `SpreadState` object instead of an `output_matrices` dict.

# Version 2025.3.8

## Changes for Developers
- Major round of efficiency optimizations (many of them non-trivial and relying heavily on low-level programming and C semantics), improving throughput by over 1400x.
- Overhauled the high-level structure of the fire spread code. This was done for optimization purposes, but it also paves the way to using other fire behavior equations than Rothermel’s, as this change largely decouples the fire spread logic from the fire behavior calculations.
- Redefined frontier cells as burnable cells with a burnable neighbor of opposite phi sign. This optimization dramatically improves performance on landscapes with a lot of non-burnable cells, restoring near linear scaling of runtime to the number of cells burned.
- Fixed numerical stability issues (related to flux-limiting).
- The phi matrix is now padded by 2 cells on each side. It was already planned to do breaking changes to this part of the API (because we don’t want to let users create output matrices themselves), so it’s not a huge regression.
- Added a new `num_tracked_cells` metric to the return results from `spread_fire_with_phi_field`. This is mostly only useful for testing purposes.
- Added `build-cython`, `build-cython-profiled`, `benchmark`, `profile`, and `snakeviz` commands to `make.sh`.

## Changes for Users
- The `fuel_model_table` dictionary in `pyretechnics.fuel_models` is no longer public. Fuel models should now be looked up with `get_fuel_model(fm_number)`, `list_fuel_model_numbers()`, and `list_fuel_models()` instead.
- More fuel model properties (in `pyretechnics.fuel_models`) are now pre-calculated at module load time.
- The default value for `surface_lw_ratio_model` is now "behave" rather than "rothermel" to better agree with both ELMFIRE and Behave7.
- The `calc_combined_fire_behavior` function in `pyretechnics.crown_fire` was updated to allow for passive crown fires with a crown spread rate of 0 since this is mathematically possible and matches FlamMap's behavior.
- Added `burn_all_cells_as_head_fire` and `burn_all_cells_toward_azimuth` functions to `burn_cells.py`.
- The conversion coefficient between 10m and 20ft wind speeds has been updated to match Behave7.
- Floating point precision has been changed from 64 bit doubles to 32 bit floats for memory efficiency.
- Functions that return fire behavior metrics as dictionaries may contain a `dphi_dt` attribute that can be ignored. Also, the `fire_type` attribute in these dictionaries has changed type from string to integer as follows:
  - "unburned" -> 0
  - "surface" -> 1
  - "passive_crown" -> 2
  - "active_crown" -> 3
- Added an optional `cube_refresh_rates` parameter to `spread_fire_with_phi_field` in `pyretechnics.eulerian_level_set`, which makes Pyretechnics interpolate between temporal band readings for each layer at its specified refresh rate.
- The following default values have been set for unspecified `spot_config` parameters in `pyretechnics.eulerian_level_set.spread_fire_with_phi_field`:
  - `firebrands_per_unit_heat`: 5e-11
  - `downwind_distance_mean`: 11.7
  - `crosswind_distance_stdev`: 10.0
  - `decay_distance`: 500.0

# Version 2024.11.7

## Changes for Developers
- Added `org-eval` command to `make.sh` to complete Emacs-free developer tooling.
- Added separate `upload-pypi` and `upload-testpypi` commands to `make.sh`.
- Updated woven LP documentation for Github Pages compatibility.

## Changes for Users
- Added LP section "How to Spread a Fire, Pause, Fork, and Continue".
- Renamed `heat_per_firebrand` to `firebrands_per_unit_heat` in `spotting_config`.
- Allow `SpaceTimeCube` and `LazySpaceTimeCube` look-up functions to accept `None` for open index ranges.
- Return `random_generator` from `pyretechnics.eulerian_level_set.spread_fire_with_phi_field`, so it can be used to re-seed subsequent runs.
- Exposed `surface_lw_ratio_model` parameter in `burn_cell_*` and `spread_fire_*` functions.
- Renamed `max_length_to_width_ratio` parameter to `crown_max_lw_ratio` in all functions that use it for clarity.

# Version 2024.10.31

Initial beta release, containing functions for:
- dataset unification (via `SpaceTimeCube` and `LazySpaceTimeCube` classes)
- fuel model and moisture equations
- surface fire equations
- crown fire equations
- spot fire equations
- burning cells on a grid (surface fire + crown fire)
- fire spread algorithms
- 3D vector operations
- units conversion functions
- LANDFIRE data loading functions
