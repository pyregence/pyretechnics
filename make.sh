#!/usr/bin/env bash

USAGE=$(cat <<EOM
Usage: $0 <command>

  where <command> is one of:

  - shell
  - container-shell
  - test
  - build
  - install
  - weave
  - tangle
  - detangle

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

case $1 in

    "shell")
        guix time-machine -C channels.scm -- shell -D -f guix.scm
        ;;

    "container-shell")
        guix time-machine -C channels.scm -- shell -D -f guix.scm --container --network --link-profile -S /usr/bin/env=bin/env --share=$HOME/.ssh
        ;;

    "test")
        guix time-machine -C channels.scm -- shell -D -f guix.scm -- pytest -vv
        ;;

    "build")
        guix time-machine -C channels.scm -- build -f guix.scm
        ;;

    "install")
        guix time-machine -C channels.scm -- install -L .guix/modules python-wrapper python-pyretechnics
        ;;

    "weave")
        guix time-machine -C channels.scm -- shell -D -f guix.scm -- ./org/weave.el
        ;;

    "tangle")
        guix time-machine -C channels.scm -- shell -D -f guix.scm -- ./org/tangle.el
        ;;

    "detangle")
        guix time-machine -C channels.scm -- shell -D -f guix.scm -- ./org/detangle.el
        ;;

    *)
        echo "$USAGE"
        ;;

esac
