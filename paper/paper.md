---
title: 'Pyretechnics: A High-Performance Python Library for Wildland Fire Simulation'
tags:
  - Python
  - Cython
  - wildland fire
  - fire modeling
  - fire spread
  - fire behavior
  - fire suppression
  - Rothermel equations
  - level-set method
authors:
  - name: Gary W. Johnson
    orcid: 0000-0001-8609-9256
    affiliation: 1
  - name: Valentin Waeselynck
    orcid: 0009-0007-4577-7912
    affiliation: 1
affiliations:
 - name: Spatial Informatics Group, LLC
   index: 1
date: 23 November 2026
bibliography: paper.bib
---

# Summary

`Pyretechnics` is a Python library (implemented in Cython for
performance) for simulating wildland fire spread and behavior. It can
be used by itself or in combination with other Python libraries and
tools as part of larger data science workflows. Its public API enables
users to perform a wide variety of scientific analyses from simple,
targeted questions like "how does the spread rate vary when the wind
speed changes" all the way up to large scale Monte Carlo simulations
of fire spread and behavior over delimited regions of space-time,
producing outputs such as annual burn probabilities and fire type and
intensity distributions.

Because individual fires are events with clearly defined start and end
times and which only burn a finite region of space, Pyretechnics
models their evolution from ignition point or perimeter to final
burned area as the rendering of a three-dimensional shape in
space-time, wherein the X-Y plane represents space and the Z index
represents time.

![Example Fireline Intensity in a Space-Time Cube](pics/els_fireline_intensity_main_cube.png){ width=40% }

Users wishing to simulate fire spread and behavior begin by loading
rasters (2D or 3D arrays) of landscape measurements, such as slope,
aspect, wind speed and direction, surface and canopy fuels (e.g.,
grasses, shrubs, downed timber, and trees), and fuel moisture, into
Numpy arrays and then wrapping these in `SpaceTimeCube` objects (from
`pyretechnics.space_time_cube`) to initialize the simulation
environment. They then call functions to set the ignition point(s) or
perimeter(s) and launch the spread algorithm. The results of these
function calls are new Numpy arrays, containing per-cell values for
the fire type, spread rate, spread direction, flame length, fireline
intensity, and time of arrival over the final burned area.

Importantly, the outputs of the fire spread functions can be fed back
into themselves to resume a "paused" simulation from where it left off
and allow the fire to continue burning until the user's next stop
condition is reached. This design enables user-extensible workflows
which alternate between Pyretechnics' fire spread algorithms and
user-supplied modifications to the landscape or model parameters, such
as for modeling suppression or feedback loops between the fire and the
surrounding environment.

# Statement of need

`Pyretechnics` implements a wide variety of physical and empirical
equations, which together can be used to predict the evolution of a
fire front and compute related fire behavior metrics (e.g., fire type,
spread rate and direction, fireline intensity, and flame length).

These computations involve 2 sorts of components:

1. A theoretical *fire behavior model*, which specifies the equations
   that model the fire. This includes among other things the Rothermel
   model with elliptical Huygens wavelets.

2. *Numerical algorithms*, which approximate the above by discrete
   data structures that computers can process. An example of this is
   the Eulerian level set method, which Pyretechnics uses for fire
   front tracking in spread simulations.

At the core of Pyretechnics is an algorithm for propagating fire
fronts. This involves:

1. The Rothermel model [@Rothermel1972], which computes the
   spread rate and behavior metrics of surface fire in the heading
   direction only.

2. Equations for *elliptical dimensions*, enabling the model to
   calculate the spread rate and behavior metrics of surface fire in
   the direction *normal* to the fire front, by applying the Huygens
   principle [@Anderson1982ModellingSpread;
   @KnightColeman1993PerimeterExpansionHuygens;
   @Richards1995fire2dSpread].

3. Crown fire equations, which extend the above to model the
   combustion of the canopy, which can make the fire front spread
   faster and more intensely. Currently, Pyretechnics implements the
   equations of Cruz and Van Wagner [@Cruz2005; @VanWagner1977].

4. A spotting module, which models the ignition of long-range spot
   fires using equations of Sardoy
   [@Sardoy2008FirebrandDistribution].

5. A numerical solver for tracking the fire front. Currently, that's
   the Eulerian level set method [@Giga2002SurfaceEvolution;
   @Osher2004LevelSetMethods], similarly to ELMFIRE
   ([https://elmfire.io/](https://elmfire.io/)).

The inputs of Pyretechnics are rasters of fuels, topography, and
weather with (at most) 2 spatial dimensions and 1 temporal dimension.
Fuels inputs include a categorical raster of fuels (following the
classification of Scott and Burgan [@ScottBurgan2005]), and 4
canopy-related metrics, which are used both to resolve the midflame
wind speed of surface fire and to compute crown fire behavior.
Topography inputs include DEM-derived slope and aspect rasters.
Weather inputs include wind speed and direction as well as fuel
moistures for 6 surface fuel classes and 1 canopy fuel class.

Pyretechnics puts a strong emphasis on computing efficiency at the
expense of some physical accuracy. It is suitable for simulating
hundreds of thousands of fire seasons on the scale of a country, with
a pixel resolution on the order of 50m. Pyretechnics does not
currently compute any feedback of the fire on the atmosphere, nor does
it automatically simulate suppression by firefighters or extinction.
However, users can inject firelines onto the landscape via Python
scripting (see `pyretechnics.suppression`), and a weather-based spread
rate adjustment can be used to slow or stop fire spread in response to
changing conditions.

# State of the field

Pyretechnics' modeling capabilities are similar in scope to those of
FARSITE [@Finney2004FARSITE], Prometheus [@Tymstra2010Prometheus], or
ELMFIRE [@Lautenberger2013elmfire]. However, its ergonomics are quite
different because Pyretechnics is meant to be used as a Python library
rather than as a standalone program with a declarative config file
format.

This distinction provides a number of important benefits, including:

1. **Flexibility:** Rather than being constrained by the fixed set of
   options in these models' config files (each of which essentially
   specifies a unique and non-trivial DSL for fire modeling),
   Pyretechnics lets its users craft their experimental designs using
   the full power of the Python programming language. In particular,
   its users can create analytical workflows that its authors did not
   anticipate because they can combine its public API functions in
   novel ways. Notably, its design enables users to interleave their
   own computations into its fire spread simulations by passing
   control back and forth between Pyretechnics and the user.

2. **Ease of Installation and Use:** Pyretechnics can be easily
   imported into any existing Python environment and combined with
   other Python libraries, tools, and frameworks as part of larger
   data science workflows or compute pipelines. There is no software
   to manually compile or install, no config file syntax to learn, and
   no custom input or output file formats to work around. Just load
   input data into Numpy arrays, call Python functions on them, and
   get Python collections, scalar values, and Numpy arrays back.

3. **Efficiency:** Because fire modeling is computationally intensive,
   Pyretechnics is implemented almost entirely in Cython for
   performance, and its API design makes it easy to parallelize
   simulations using popular Python libraries, such as
   `multiprocessing`, `dask`, and `ray`. As part of larger
   Python-based compute pipelines, Pyretechnics also benefits from
   being able to keep intermediate calculations in memory and pass
   them directly between other Python library functions in the runtime
   rather than needing to serialize intermediate results out to disk
   and then deserialize them back into the next step of a computation.

4. **File Format Support:** Unlike most existing fire models, which
   only support a few input and output file formats (that are
   frequently incompatible with one another), Pyretechnics only
   expects its functions' input data to be provided as Numpy arrays.
   Since Python already provides libraries that can be used to read
   all prominent GIS formats into Numpy arrays, modify them after
   reading, or even generate new Numpy arrays programmatically,
   Pyretechnics can easily be used with any file formats supported by
   Python. For the same reason, Pyretechnics outputs can be written
   into any output file format, whether spatial or aspatial, that is
   supported by Python.

# Software design

The Pyretechnics library provides Python modules (written in the
Python-compatible syntax of Cython for performance) that implement
many of the fundamental equations of fire spread and behavior as
published in the wildland fire science literature. These fall into the
following categories:

- **pyretechnics.fuel_models:** Fuel Model and Moisture Definitions
- **pyretechnics.surface_fire:** Surface Fire Equations
- **pyretechnics.crown_fire:** Crown Fire Equations
- **pyretechnics.spot_fire:** Spot Fire Equations
- **pyretechnics.burn_cells:** Burning Cells on a Grid
- **pyretechnics.eulerian_level_set:** Fire Spread Algorithm (ELMFIRE): Eulerian Level Set Method

Each of these modules builds upon the functions and equations
implemented in those listed above it, while exposing a public API for
users who want to focus their analyses on these low-level elements of
fire modeling, perform verification or validation experiments with
them, or even build their own higher level spread algorithms on top of
them.

This offers a great deal of flexibility over existing fire modeling
applications, which typically read their inputs, simulate fire spread
and behavior in its entirety, and write their outputs back to disk.
Parameterizing these models is usually done by setting parameters in
their configuration files, and experimenting with or reusing their
internal functions may involve significantly more development effort.

In addition to these core wildland fire science modules, Pyretechnics
provides a module for unifying 0D (constant), 1D (temporal), 2D
(spatial), and 3D (spatiotemporal) input datasets (of potentially
different resolutions) for fuels, topography, wind, and moisture
called `pyretechnics.space_time_cube`.

The user of this library simply loads in their input datasets as Numpy
arrays, wraps them in `SpaceTimeCube` objects to make all of their
dimensions and resolutions match, and feeds them into the functions
for burning cells (`pyretechnics.burn_cells`) or spreading a fire from
a point or existing burn scar (`pyretechnics.eulerian_level_set`). The
cell-burning functions will return the fire behavior metrics for the
burned cell (i.e., fire type, spread rate, spread direction, fireline
intensity, flame length) as a Python dictionary. The fire spreading
functions return a dictionary containing these same metrics (+ time of
arrival) in Numpy arrays, covering the simulation area.

If a user wishes to model active fire suppression in a spread
simulation, they can either use the functions in the
`pyretechnics.suppression` module or write their own Python functions
to burn firelines into the input fuel rasters.

Since the library's main functions take and return Numpy arrays, users
are free to create them as needed for their projects. They can load
arrays from numerous geospatial file formats using the Python
`rasterio` library. They could grab data from a PostGIS raster using
`psycopg2` or load in NetCDF files with the `netCDF4` library. They
could also create arrays using Numpy or Scipy functions directly
and/or use them to apply random noise, value perturbations, or other
preprocessing steps on the arrays before initiating fire spread
simulations.

Similarly, any Python library or algorithm can be used to post-process
the output arrays into aggregate arrays (e.g., annual/conditional burn
frequency), visualize them (e.g., with `matplotlib`), and/or write
them out as raster files, CSVs, or any other format required by the
user.

Notably, Pyretechnics is written as a Literate Program, which can be
read online here: https://pyregence.github.io/pyretechnics

Each Python module has its own dedicated chapter, with two subsections:
- For Developers
- For Users

The `For Developers` section contains the **literate programming**
implementation of the module's equations and algorithms for those
interested in understanding the science and engineering decisions used
within it.

The `For Users` section contains **notebook programming** examples of
how to use the public functions in each module along with their
computed outputs. Users can read these sections and cut-and-paste
their example scripts into their own Python files or Jupyter notebooks
to quickly start using the library.

# Research impact statement

## TODO: Include evidence of realized impact (publications, external use, integrations) or credible near-term significance (benchmarks, reproducible materials, community-readiness signals). The evidence should be compelling and specific, not aspirational.

# AI usage disclosure

No generative AI tools were used in the development of this software,
the writing of this manuscript, or the preparation of supporting
materials.

# Acknowledgements

We acknowledge support from David Saah and Spatial Informatics Group,
LLC in the creation and ongoing development of this software project
and technical inspiration from Chris Lautenberger of CloudFire Inc.

# References
