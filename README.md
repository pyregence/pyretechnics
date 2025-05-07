# Pyretechnics

A Python library for simulating fire behavior in a variety of ways.

Latest API Documentation: https://pyregence.github.io/pyretechnics/

Conda installation instructions: https://www.anaconda.com/docs/tools/working-with-conda/packages/pip-install

## Overview

The Pyretechnics library provides modules that implement the fundamental equations used in most operational wildland fire behavior models like GridFire, ELMFIRE, FlamMap, FARSITE, FSIM, and BehavePlus. These fall into the following categories:

- pyretechnics.fuel_models: Fuel Model and Moisture Definitions
- pyretechnics.surface_fire: Surface Fire Equations
- pyretechnics.crown_fire: Crown Fire Equations
- pyretechnics.spot_fire: Spot Fire Equations
- pyretechnics.burn_cells: Burning Cells on a Grid
- pyretechnics.eulerian_level_set: Fire Spread Algorithm (ELMFIRE)

In addition, it provides a module for unifying 0D (constant), 1D (temporal), 2D (spatial), and 3D (spatiotemporal) input datasets (of potentially different resolutions) for fuels, topography, wind, and moisture called `pyretechnics.space_time_cube`.

The user of this library simply loads in their input datasets as Numpy arrays, wraps them in `SpaceTimeCube` objects to make all of their dimensions and resolutions match, and feeds them into the functions for burning cells (`pyretechnics.burn_cells`) or spreading a fire from a point or existing burn scar (`pyretechnics.eulerian_level_set`).
The cell-burning functions will return the fire behavior metrics for the burned cell (i.e., fire type, spread rate, spread direction, fireline intensity, flame length) as a Python dictionary. The fire spreading functions return a dictionary containing these same metrics (+ time of arrival) in Numpy arrays, covering the simulation area.

Since the library's main functions take and return Numpy arrays, you are free to create them as needed for your own projects. You can load arrays from numerous geospatial file formats using the Python `rasterio` library. You could grab data from a PostGIS raster using `psycopg2` or load in NetCDF files with the `netCDF4` library. You could also create your own arrays using Numpy or Scipy functions based on your needs. And, of course, you could use these same libraries to apply random noise, value perturbations, or other preprocessing steps on the arrays before initiating fire spread simulations.

Similary, you can use any Python library or algorithm to post-process the output arrays into aggregate arrays (e.g., annual/conditional burn frequency), visualize them (e.g., with `matplotlib`), and/or write them out as raster files, CSVs, or any other format that you can think of.

Notably, Pyretechnics is written as a Literate Program, which you can read online here: https://pyregence.github.io/pyretechnics

Each module has its own dedicated chapter, with two subsections:
- For Developers
- For Users

The `For Developers` sections are the literate programming implementation of the equations and algorithms included in that module. Read these sections if you want to understand the science and engineering behind that module.

The `For Users` sections are notebook programming examples of *how to use* the public functions in each module along with their computed outputs. If you just want to use the Pyretechnics library, you can jump to these sections and cut-and-paste the example scripts into your own Python file or Jupyter notebook to get yourself going quickly.

## Design Principles

### Free and Open Source Software

To promote open review and encourage collaborative development of the
various algorithms implemented in Pyretechnics, this library is
released under a free and open source license. See the
[License](#license) section in this document for more information.

### Reproducible Research

One of the most persistently challenging aspects of software
development is the fact that the environments in which software is
built vary from one person's computer to another, including but not
limited to different types and versions of operating systems,
applications, libraries, and services installed. Since most software,
including Pyretechnics, is written so as to rely on external libraries
and applications at both build and run time, it is necessary to be
able to easily reproduce the computer environments needed for these
steps on new users' machines.

To do so, Pyretechnics uses [GNU Guix](https://guix.gnu.org) to
automatically install software dependencies and create ephemeral
development environments without requiring root privileges or
interfering with the main package manager of the underlying operating
system.

If you are running a GNU/Linux distribution on your computer, the
easiest way to install Guix is to simply follow the "Binary
Installation" instructions in the manual here:

https://guix.gnu.org/manual/en/html_node/Binary-Installation.html

This will add the `guix` command as an auxiliary package, environment,
and container manager on your machine.

If you are not running GNU/Linux, you will need to run Guix System in
a virtual machine. This is a complete GNU/Linux distribution that uses
Guix as its package manager and Shepherd as its service manager. You
can find instructions on getting this installed and running through
QEMU in the "Running Guix in a Virtual Machine" section of the Guix
manual here:

https://guix.gnu.org/manual/en/html_node/Running-Guix-in-a-VM.html

1.  Creating a Reproducible Development Environment

    Once installed, you can tell Guix to download all the necessary
    dependencies for Pyretechnics development and enter a shell in
    which they are available by running this command from the root
    directory of this repository:

    ```sh
    ./make.sh shell
    ```

    On its first invocation, this command will download the necessary
    software packages and install them under `/gnu/store`. When this
    is done, you will be dropped into a shell environment in which the
    environment variables have been automatically configured to point
    to these newly downloaded packages.
    
    On subsequent calls to `./make.sh shell`, you should be dropped
    directly into the shell environment without the need to install
    any new software unless the [guix.scm](guix.scm) or
    [channels.scm](channels.scm) files have been updated.

2.  Authorizing Guix to Automatically Read guix.scm

    The first time that you run `./make.sh shell`, you will be
    prompted to authorize Guix to read the [guix.scm](guix.scm) file
    in this repository for instructions on what to download and how to
    set up your ephemeral shell environment. Assuming you have set the
    `PYRETECHNICS_REPO` shell variable to the directory containing
    this repository, you can do so with this command:

    ```sh
    echo $PYRETECHNICS_REPO >> $HOME/.config/guix/shell-authorized-directories
    ```

3.  Exiting the Reproducible Development Environment

    You can always exit from the shell with this command:

    ```sh
    exit
    ```

4.  Running the Test Suite

    You can run the Pyretechnics library's test suite by invoking
    `pytest` through the Guix shell like so:

    ```sh
    ./make.sh test
    ```

5.  Building the Pyretechnics Library with Guix

    To build the Pyretechnics library, including running its tests,
    constructing a Python wheel, and unpacking it into the Guix
    /gnu/store directory, simply run this command:

    ```sh
    ./make.sh build-guix
    ```

6.  Building the Pyretechnics Library as a Distribution

    To create a *dist* folder containing source (*.tar.gz*) and built (*.whl*)
    distributions of the Pyretechnics library, you can run this command:

    ```sh
    ./make.sh build-dist
    ```

7.  Uploading the Built Distribution to TestPyPI

    To upload the built distribution to
    [TestPyPI](https://test.pypi.org/), you can use this command:

    ```sh
    ./make.sh upload-testpypi
    ```

    You will be prompted for a username and password. For the
    username, use `__token__`. For the password, use the TestPyPI API
    token value that you created
    [here](https://test.pypi.org/manage/account/#api-tokens), including the
    `pypi-` prefix.

8.  Uploading the Built Distribution to PyPI

    To upload the built distribution to [PyPI](https://pypi.org/), you
    can use this command:

    ```sh
    ./make.sh upload-pypi
    ```

    You will be prompted for a username and password. For the
    username, use `__token__`. For the password, use the PyPI API
    token value that you created
    [here](https://pypi.org/manage/account/#api-tokens), including the
    `pypi-` prefix.

9.  Installing the Pyretechnics Library with Guix

    You have two options for installing the Pyretechnics library locally:

    First, you can simply install it into a temporary shell
    environment like so:

    ```sh
    ./make.sh install-shell
    ```

    You can leave this shell by typing `exit`.

    Your second option is to install the Pyretechnics library into
    your Guix profile with this command:

    ```sh
    ./make.sh install-guix
    ```

    Next, you will need to invoke the following Bash commands in your
    shell to make the newly installed library available via
    `$GUIX_PYTHONPATH`. This environment variable is referenced
    automatically by the Guix-installed Python package.

    ```sh
    GUIX_PROFILE="$HOME/.guix-profile"
    . "$GUIX_PROFILE/etc/profile"
    ```

    It is recommended that you add these two lines to your
    `$HOME/.bash_profile`, so that they are run automatically each
    time you login.

10. Using the Pyretechnics Library

    Once you have installed the library into a temporary shell
    environment, installed it into your Guix profile, or downloaded it
    from [PyPI](https://pypi.org/), you should be able to launch
    `python` and load the library as follows:

    ```python
    import pyretechnics
    ```

### Literate Programming

Pyretechnics has been coded as a [literate
program](https://en.wikipedia.org/wiki/Literate_programming) using
Emacs' [org-mode](http://orgmode.org/worg/org-contrib/babel/). This
means that its source code is embedded within a larger document, which
explains the rationale behind the code using text, equations,
citations, tables, and figures. The reason for this approach is to
make the science and engineering decisions within Pyretechnics fully
transparent to our users, whether or not they feel confident reading
source code directly. We believe that this better promotes the goals
of open science than open source software alone.

What this means practically is that the [org](org) directory in this
software repository contains a single literate program document called
[org/pyretechnics.org](org/pyretechnics.org), which is used to
automatically generate all of the other source code and documentation
files within this repository.

1.  Regenerating HTML Documentation

    The latest HTML documentation can always be found in
    [docs/index.html](docs/index.html).

    After editing [org/pyretechnics.org](org/pyretechnics.org), you can
    regenerate the HTML documentation by running this command:

    ```sh
    ./make.sh weave
    ```

2.  Regenerating the Source Code Tree

    After editing [org/pyretechnics.org](org/pyretechnics.org), you can
    regenerate all the source code files in this repository by running
    this command:

    ```sh
    ./make.sh tangle
    ```

3.  Bringing Source Code File Edits Back into the Literate Program

    If you edit a source code file directly, its changes can be
    automatically incorporated back into the literate program by
    running this command:

    ```sh
    ./make.sh detangle
    ```

4.  Running All the Source Code Blocks in the Literate Program

    If your changes would impact the results of the example code
    blocks in the literate program, then you can run them again to
    update their results in
    [org/pyretechnics.org](org/pyretechnics.org) with this command:

    ```sh
    ./make.sh org-eval
    ```

## Contact

### Authors

- Gary W. Johnson, PhD
  - Email: gjohnson@sig-gis.com
  - Web: https://sig-gis.com

- Valentin Waeselynck
  - Email: vwaeselynck@sig-gis.com
  - Web: https://sig-gis.com

- Chris Lautenberger, PhD, PE
  - Email: chris@cloudfire.com
  - Web: https://cloudfire.com

- David Saah, PhD
  - Email: dsaah@sig-gis.com
  - Web: https://sig-gis.com

<a id="license"></a>
## License

Copyright © 2023-2025 Spatial Informatics Group, LLC.

Pyretechnics is distributed by Spatial Informatics Group, LLC. under
the terms of the Eclipse Public License version 2.0 (EPLv2). See
[LICENSE](LICENSE) in this directory for more information.
