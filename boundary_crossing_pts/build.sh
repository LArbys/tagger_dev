#!/bin/bash

homedir=$PWD
echo "Starting from $homedir"

source setup.sh
source configure.sh

echo "LARLITE: ${LARLITE_BASEDIR}"
echo "LARLITE: ${GEO2D_BASEDIR}"
echo "LARCV: ${LARCV_BASEDIR}"
echo "LAROPENCV: ${LAROPENCV_BASEDIR}"
echo "LARLITECV: ${LARLITECV_BASEDIR}"

cd $LARLITE_BASEDIR
make -j4 || return 1

cd $LARLITE_BASEDIR/UserDev/BasicTool
make -j4 || return 1

cd $LARLITE_BASEDIR/UserDev/SelectionTool/OpT0Finder
make -j4 || return 1

cd $LARLITE_BASEDIR/UserDev/RecoTool
make -j4 || return 1

cd $LARLITE_BASEDIR/UserDev/RecoTool/ClusterRecoUtil
make -j4 || return 1

cd $GEO2D_BASEDIR
make -j4 || return 1

cd $LAROPENCV_BASEDIR
make -j4 || return 1

cd $LARCV_BASEDIR
make -j4 || return 1

cd $LARLITECV_BASEDIR
make -j4 || return 1

cd $homedir

echo "DONE"
