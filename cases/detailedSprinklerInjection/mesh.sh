#!/bin/sh
#set -x 

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

./clean.sh

runApplication blockMesh

echo "setSet -batch createSamplePlane.setSet"
runApplication setSet -batch system/createSamplePlane.setSet -time 0 -noZero -constant

runApplication topoSet -dict system/topoSetDict.thrust

cp 0/ph_rgh.orig 0/ph_rgh

decomposePar -force

# -----------------------------------------------------------------------------
