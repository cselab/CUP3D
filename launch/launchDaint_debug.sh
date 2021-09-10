#!/bin/bash

set -e

WCLOCK=${WCLOCK:-00:30:00}
PARTITION=${PARTITION:-debug}
NNODE=${NNODE:-4}
_SCRIPTPATH=`pwd`"/launchDaint_debug.sh"

source launchDaint.sh

cp $_SCRIPTPATH $FOLDER
