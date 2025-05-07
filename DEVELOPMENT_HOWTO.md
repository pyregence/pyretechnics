# Setting up Your Development Environment

The `pyretechnics` repo brings along its own opinionated development environment through the magic of Guix. This is documented in the repo's `README.md` file.

What this means for you is that you need to install Guix on your machine in order to get going. Otherwise, you are going to have to recreate the entire development experience from scratch, which I would not recommend. Note that Guix is a Linux application, so this becomes quite easy to set up on a Linux machine but quite difficult on a Mac or Windows machine. In these cases, your best bet is to spin up a Linux VM to work in, which could either be running the Guix System distribution or another Linux distribution like Ubuntu that you can install Guix into as a package.

Once this is done, your main development tool is the `make.sh` shell script in the top level of the repository:

```
user@host pyretechnics $ ./make.sh
Usage: ./make.sh <command>

  where <command> is one of:

  - shell
  - test
  - benchmark
  - profile
  - snakeviz
  - build-cython
  - build-cython-profiled
  - build-guix
  - build-dist
  - upload-testpypi
  - upload-pypi
  - install-shell
  - install-guix
  - weave
  - tangle
  - detangle
  - org-eval

  See README.md for more information.
```

Each of these commands will launch a Guix environment (similar to a Docker container) containing the set of packages specified in `guix.scm`'s `native-inputs` and `propagated-inputs` fields and will then invoke some command before exiting the environment. Note, of course, that this is not the `make` command but a simple shell script with a similar name. The reason for this is to avoid a bootstrapping dependency for yet another development tool. You just need Guix and a shell. The rest should be handled by this setup.

# File Locations

- org/: Contains `pyretechnics.org` (the Literate Program source) and several Elisp scripts to enable the standard LP workflow. (Don't worry. Guix brings along its own minimal Emacs to run them.)

- docs/: Contains `index.html` (the LP woven documentation as shown on https://pyregence.github.io/pyretechnics/). Regenerate this from `pyretechnics.org` with `./make.sh weave`.

- src/: Contains source code for the Python modules. Regenerate this from `pyretechnics.org` with `./make.sh tangle` or apply changes made to these files back to `pyretechnics.org` with `./make.sh detangle`.

- test/: Contains source code for the Python module tests. Regenerate this from `pyretechnics.org` with `./make.sh tangle` or apply changes made to these files back to `pyretechnics.org` with `./make.sh detangle`.

- prof/: Contains `spread_fire.py`, which runs the simple 8-hour surface fire simulation example from [here](https://pyregence.github.io/pyretechnics/#how-to-spread-a-fire-from-a-point-ignition). Run it with `./make.sh benchmark` or profile it with `./make.sh profile`. The second command will write `prof/spread_fire.prof`, which can be viewed with `./make.sh snakeviz`.

# Development Workflow

If you want a persistent shell in the Guix environment, you can run `./make.sh shell`. You can think of this like running the `./activate.sh` script for a Python virtual environment or running `docker run -it some_container bash`. You can enter the same environment in a kernel isolated container with `./make.sh container-shell`. In either case, typing `exit` will leave this environment and return you to your original calling shell. Please note that you don't need to be in this Guix shell in order to do development work. It is simply provided as a convenience for developers who want to ensure reproducibility between their environment and the one in which the project's various build commands are run.

## Without Emacs

If you typically work outside of Emacs, your workflow is going to look like this:

1. Open some `*.py` or `*.pxd` files in `src/` or `test/` in your editor, eval them, and make some edits. (Either connect to a Python process running under `./make.sh shell` or spin up a new Python virtual environment using `requirements.txt`.)

2. Run `./make.sh build-cython` to regenerate the `*.c`, `*.html`, and `*.so` files from your changed `*.py` and `*.pxd` files. These new files will appear under `src/pyretechnics/`, `src/build/`, `test/pyretechnics/`, and `test/build/`.

3. [Optional] Run `./make.sh benchmark` to see how fast your new code can run the demo 8-hour simulation.

4. Run `./make.sh profile` to regenerate `prof/spread_fire.prof` by running your new code. Remember that the reported runtime is going to be about 2x what you see from `./make.sh benchmark` because of the slowdown from `cProfile`.

5. Run `./make.sh snakeviz` (in a separate terminal) to spin up a server for browsing `prof/spread_fire.prof` in your web browser. You only have to do this once and can then refresh your browser window on each successive `./make.sh profile` call.

6. Explore the profiling results in your web browser to determine what you want to work on next.

7. Run `./make.sh detangle` to copy your changes to the `*.py` and `*.pxd` files back into `org/pyretechnics.org`.

8. Run `./make.sh org-eval` to rerun all the example scripts inside `org/pyretechnics.org`. Then check your `git status` to see if anything changed in that file or the binary image files that it generates. If so, you broke something (unless it's just a tiny change in floating point values due to changing floating point functions or variable precisions).

9. Run `./make.sh weave` to regenerate the woven documentation in `org/index.html`.

10. Commit and push your changes!

## With Emacs

If you typically work in Emacs, you can either use the workflow above or make your edits directly to the literate program in `org/pyretechnics.org` like so:

1. Open `org/pyretechnics.org` in Emacs. The new buffer should open in `Org` major mode.

2. Move your cursor inside any code block (bounded by `#+begin_src ... #+end_src`) and press `C-c '` (`org-edit-special`) to open a temporary buffer in Python major mode.

3. Use `C-c C-p` (`run-python`) to open a new Python buffer bound to the literate program's outline section session.

4. Use `C-c C-c` (`org-ctrl-c-ctrl-c`) to evaluate the Python code from your buffer in the connected Python interpreter.

5. Make some edits and evaluate them with `C-c C-c` (`org-ctrl-c-ctrl-c`) as needed.

6. Press `C-c '` (`org-edit-special`) to return to `org/pyretechnics.org` with your new changes applied to the code block.

7. Press `C-c C-c` (`org-ctrl-c-ctrl-c`) to re-run the code block in the connected Python interpreter and embed its results (if any) in the literate program below the source code block.

8. Save the file and run `C-c C-v C-t` (`org-babel-tangle`) to regenerate the `*.py` and `*.pxd` files in `src/`, `test/`, and `prof/`.

9. Use the same steps 2-10 from the "Without Emacs" section above (skipping step 7 [detangle], of course).

Happy hacking!
