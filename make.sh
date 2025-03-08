#!/usr/bin/env bash

USAGE=$(cat <<EOM
Usage: $0 <command>

  where <command> is one of:

  - shell
  - container-shell
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
  - install-user
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

# Ensure current directory contains channels.scm and guix.scm
cd $SCRIPT_DIR

case $CMD in

    "shell")
        guix time-machine -C channels.scm -- shell -D -f guix.scm $ARGS
        ;;

    "container-shell")
        guix time-machine -C channels.scm -- shell -D -f guix.scm --container --network --link-profile -S /usr/bin/env=bin/env --preserve=PYTHONPATH --share=$HOME/.ssh $ARGS
        ;;

    "test")
        guix time-machine -C channels.scm -- shell -D -f guix.scm -- pytest -vv $ARGS
        ;;

    "benchmark")
        guix time-machine -C channels.scm -- shell -D -f guix.scm -- python prof/spread_fire.py
        ;;

    "profile")
        guix time-machine -C channels.scm -- shell -D -f guix.scm -- python -m cProfile -o prof/spread_fire.prof prof/spread_fire.py
        ;;

    "snakeviz")
        guix time-machine -C channels.scm -- shell -D -f guix.scm -- snakeviz prof/spread_fire.prof
        ;;

    "build-cython")
        guix time-machine -C channels.scm -- shell -D -f guix.scm -- python setup.py build_ext --inplace
        ;;

    "build-cython-profiled")
        guix time-machine -C channels.scm -- shell -D -f guix.scm -- env PROFILE_CYTHON=1 python setup.py build_ext --inplace
        ;;

    "build-guix")
        guix time-machine -C channels.scm -- build -f guix.scm $ARGS
        ;;

    "build-dist")
        rm -rf dist && guix time-machine -C channels.scm -- shell -D -f guix.scm -- python setup.py sdist $ARGS
        ;;

    "upload-testpypi")
        guix time-machine -C channels.scm -- shell -D -f guix.scm -- python -m twine upload --repository testpypi dist/*.tar.gz $ARGS
        ;;

    "upload-pypi")
        guix time-machine -C channels.scm -- shell -D -f guix.scm -- python -m twine upload dist/*.tar.gz $ARGS
        ;;

    "install-shell")
        guix time-machine -C channels.scm -- shell -f guix.scm $ARGS
        ;;

    "install-user")
        guix time-machine -C channels.scm -- install -L .guix/modules python-wrapper python-pyretechnics $ARGS
        ;;

    "weave")
        guix time-machine -C channels.scm -- shell -D -f guix.scm -- ./org/weave.el $ARGS
        ;;

    "tangle")
        guix time-machine -C channels.scm -- shell -D -f guix.scm -- ./org/tangle.el $ARGS
        ;;

    "detangle")
        guix time-machine -C channels.scm -- shell -D -f guix.scm -- ./org/detangle.el $ARGS
        ;;

    "org-eval")
        guix time-machine -C channels.scm -- shell -D -f guix.scm -- ./org/eval.el $ARGS
        ;;

    *)
        echo "$USAGE"
        ;;

esac

# Return to the directory where make.sh was run
cd $RUN_DIR
