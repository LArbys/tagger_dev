#!/bin/bash

if [ -z ${FQ_TAGGERDEV_BASEDIR+x} ]; then
    export FQ_TAGGERDEV_BASEDIR=$PWD
fi

# setup environment variables
source $FQ_TAGGERDEV_BASEDIR/setup.sh

# setup larlite
source $FQ_TAGGERDEV_BASEDIR/larlite/config/setup.sh

# setup laropencv
source $FQ_TAGGERDEV_BASEDIR/LArOpenCV/setup_laropencv.sh

# setup Geo2D
source $FQ_TAGGERDEV_BASEDIR/Geo2D/config/setup.sh

# setup LArCV
source $FQ_TAGGERDEV_BASEDIR/LArCV/configure.sh

# setup larlitecv
source $FQ_TAGGERDEV_BASEDIR/larlitecv/configure.sh

export LD_LIBRARY_PATH=.:${LD_LIBRARY_PATH}
