#!/bin/sh

SCRIPT_DIR=`dirname "$(readlink -f "$0")"`

docker run --rm --volume $SCRIPT_DIR/../paper:/data --user $(id -u):$(id -g) --env JOURNAL=joss openjournals/inara
