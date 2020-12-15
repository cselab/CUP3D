#!/bin/bash

SCRIPTPATH="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

PYTHONPATH="$SCRIPTPATH/../../build:$PYTHONPATH" python -m unittest "$@"
