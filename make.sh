#!/usr/bin/env bash

USAGE=$(cat <<EOM
Usage: $0 <command>

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
EOM
)

if [[ $# == 0 ]]
then
    echo "$USAGE"
    exit 0
fi

if ! [[ -x "$(command -v guix)" ]]
then
    echo "Error: guix is not installed." >&2
    exit 1
fi

RUN_DIR=$(pwd)
SCRIPT_DIR=$(dirname "$0")
CMD=$1
ARGS=${@:2:$#}

export PYTHONPATH="src:test"

TIME_MACHINE="guix time-machine -C channels.scm"
DEV_CONTAINER_SHELL="shell -CN -D -f guix.scm --link-profile -S /usr/bin/env=bin/env --preserve=PYTHONPATH"
CONTAINER_SHELL="shell -CN -f guix.scm --link-profile"

# Ensure current directory contains channels.scm and guix.scm
cd $SCRIPT_DIR

case $CMD in

    "shell")
        $TIME_MACHINE -- $DEV_CONTAINER_SHELL $ARGS
        ;;

    "test")
        $TIME_MACHINE -- $DEV_CONTAINER_SHELL -- pytest -vv $ARGS
        ;;

    "benchmark")
        $TIME_MACHINE -- $DEV_CONTAINER_SHELL -- python prof/spread_fire.py
        ;;

    "profile")
        $TIME_MACHINE -- $DEV_CONTAINER_SHELL -- python -m cProfile -o prof/spread_fire.prof prof/spread_fire.py
        ;;

    "snakeviz")
        $TIME_MACHINE -- $DEV_CONTAINER_SHELL -- snakeviz -s prof/spread_fire.prof
        ;;

    "build-cython")
        $TIME_MACHINE -- $DEV_CONTAINER_SHELL -- env PYR_SET_GCC_ARGS=1 python setup.py build_ext --inplace
        ;;

    "build-cython-profiled")
        $TIME_MACHINE -- $DEV_CONTAINER_SHELL -- env PYR_SET_GCC_ARGS=1 PYR_PROFILE_CYTHON=1 python setup.py build_ext --inplace
        ;;

    "build-guix")
        $TIME_MACHINE -- build -f guix.scm $ARGS
        ;;

    "build-dist")
        rm -rf dist && $TIME_MACHINE -- $DEV_CONTAINER_SHELL -- env PYR_SET_GCC_ARGS=1 python -m build $ARGS
        ;;

    "upload-testpypi")
        $TIME_MACHINE -- $DEV_CONTAINER_SHELL -- python -m twine upload --repository testpypi dist/*.tar.gz $ARGS
        ;;

    "upload-pypi")
        $TIME_MACHINE -- $DEV_CONTAINER_SHELL -- python -m twine upload dist/*.tar.gz $ARGS
        ;;

    "install-shell")
        $TIME_MACHINE -- $CONTAINER_SHELL coreutils python-wrapper $ARGS
        ;;

    "install-guix")
        $TIME_MACHINE -- install -L .guix/modules python-wrapper python-pyretechnics $ARGS
        ;;

    "weave")
        $TIME_MACHINE -- $DEV_CONTAINER_SHELL -- ./org/weave.el $ARGS
        ;;

    "tangle")
        $TIME_MACHINE -- $DEV_CONTAINER_SHELL -- ./org/tangle.el $ARGS
        ;;

    "detangle")
        $TIME_MACHINE -- $DEV_CONTAINER_SHELL -- ./org/detangle.el $ARGS
        ;;

    "org-eval")
        $TIME_MACHINE -- $DEV_CONTAINER_SHELL -- ./org/eval.el $ARGS
        ;;

    *)
        echo "$USAGE"
        ;;

esac

# Return to the directory where make.sh was run
cd $RUN_DIR
