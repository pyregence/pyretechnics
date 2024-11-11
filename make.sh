#!/usr/bin/env bash

USAGE=$(cat <<EOM
Usage: $0 <command>

  where <command> is one of:

  - shell
  - container-shell
  - test
  - build-cython
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

CMD=$1
ARGS=${@:2:$#}

case $CMD in

    "shell")
        guix time-machine -C channels.scm -- shell -D -f guix.scm $ARGS
        ;;

    "container-shell")
        guix time-machine -C channels.scm -- shell -D -f guix.scm --container --network --link-profile -S /usr/bin/env=bin/env --share=$HOME/.ssh $ARGS
        ;;

    "test")
        guix time-machine -C channels.scm -- shell -D -f guix.scm -- pytest -vv $ARGS
        ;;

    "build-cython")
        guix time-machine -c channels.scm -- shell -D -f guix.scm -- cythonize -X profile=True -i -a src/pyretechnics/*.py
        ;;

    "build-guix")
        guix time-machine -C channels.scm -- build -f guix.scm $ARGS
        ;;

    "build-dist")
        guix time-machine -C channels.scm -- shell -D -f guix.scm -- python -m build $ARGS
        ;;

    "upload-testpypi")
        guix time-machine -C channels.scm -- shell -D -f guix.scm -- python -m twine upload --repository testpypi dist/* $ARGS
        ;;

    "upload-pypi")
        guix time-machine -C channels.scm -- shell -D -f guix.scm -- python -m twine upload dist/* $ARGS
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
