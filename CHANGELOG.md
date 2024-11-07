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

# Version 2024.11.07

Changes for Developers:
- Added `org-eval` command to `make.sh` to complete Emacs-free developer tooling.
- Added separate `upload-pypi` and `upload-testpypi` commands to `make.sh`.
- Updated woven LP documentation for Github Pages compatibility.

Changes for Users:
- Added LP section "How to Spread a Fire, Pause, Fork, and Continue".
- Renamed `heat_per_firebrand` to `firebrands_per_unit_heat` in `spotting_config`.
- Allow `SpaceTimeCube` and `LazySpaceTimeCube` look-up functions to accept `None` for open index ranges.
- Return `random_generator` from `pyretechnics.eulerian_level_set.spread_fire_with_phi_field`, so it can be used to re-seed subsequent runs.
- Exposed `surface_lw_ratio_model` parameter in `burn_cell_*` and `spread_fire_*` functions.
- Renamed `max_length_to_width_ratio` parameter to `crown_max_lw_ratio` in all functions that use it for clarity.
