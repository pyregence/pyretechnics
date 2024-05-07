# Pyretechnics

A Python library for simulating fire behavior in a variety of ways.

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

Visit <https://guix.gnu.org/en/download/> to find instructions for
downloading and installing Guix onto your computer or into a VM.

1.  Creating a Reproducible Development Environment

    Once installed, you can tell Guix to download all the necessary
    dependencies for Pyretechnics development and enter a shell in
    which they are available by running this command from the root
    directory of this repository:

    ```sh
    guix shell
    ```

    If you would like this shell environment to be isolated from your
    host machine within a container, you can use this command instead:

    ```sh
    guix shell --container --network --link-profile -S /usr/bin/env=bin/env --share=$HOME/.ssh
    ```

    On their first invocations, these commands will download the
    necessary software packages and install them under `/gnu/store`.
    When this is done, they will drop you into a shell environment in
    which the environment variables have been automatically configured
    to point to these newly downloaded packages.
    
    On subsequent calls to `guix shell`, you should be dropped
    directly into the shell environment without the need to install
    any new software unless the [guix.scm](guix.scm) file has been
    updated or you have updated your Guix version by running `guix
    pull`.

2.  Authorizing Guix to Automatically Read guix.scm

    The first time that you run `guix shell`, you will be prompted to
    authorize Guix to read the [guix.scm](guix.scm) file in this
    repository for instructions on what to download and how to set up
    your ephemeral shell environment. Assuming you have set the
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
    guix shell -D -f guix.scm -- pytest -vv
    ```

5.  Building and Installing the Pyretechnics Library with Guix

    To build the Pyretechnics library, including running its tests and
    constructing Python wheels, simply run this command:

    ```sh
    guix build -f guix.scm
    ```

    Should you wish to install the package into your Guix profile, you
    can run this:

    ```sh
    guix install -f guix.scm
    ```

    Log out and back in to activate your profile's default environment
    variables. You should now be able to import the library's modules
    from a Python interpreter with the usual `import` command.

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
    [doc/pyretechnics.html](doc/pyretechnics.html).
    
    After editing [org/pyretechnics.org](org/pyretechnics.org), you can
    regenerate the HTML documentation by running this command:

    ```sh
    guix shell -D -f guix.scm -- ./org/weave.el
    ```

2.  Regenerating the Source Code Tree

    After editing [org/pyretechnics.org](org/pyretechnics.org), you can
    regenerate all the source code files in this repository by running
    this command:

    ```sh
    guix shell -D -f guix.scm -- ./org/tangle.el
    ```

3.  Bringing Source Code File Edits Back into the Literate Program

    If you edit a source code file directly, its changes can be
    automatically incorporated back into the literate program by
    running this command:

    ```sh
    guix shell -D -f guix.scm -- ./org/detangle.el
    ```

## Contact

### Authors

- Gary W. Johnson, PhD
  - Email: gjohnson@sig-gis.com
  - Web: <https://sig-gis.com>

- Valentin Waeselynck
  - Email: vwaeselynck@sig-gis.com
  - Web: <https://sig-gis.com>

- Chris Lautenberger, PhD, PE
  - Email: chris@cloudfire.ai
  - Web: <https://cloudfire.ai>

<a id="license"></a>
## License

Copyright © 2023-2024 Spatial Informatics Group, LLC.

Pyretechnics is distributed by Spatial Informatics Group, LLC. under
the terms of the Eclipse Public License version 2.0 (EPLv2). See
[LICENSE.txt](LICENSE.txt) in this directory for more information.
