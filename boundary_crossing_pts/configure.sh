#!/bin/bash

if [ -z ${BMT_TAGGERDEV_BASEDIR+x} ]; then
    export BMT_TAGGERDEV_BASEDIR=$PWD
fi

# setup environment variables
source $BMT_TAGGERDEV_BASEDIR/setup.sh

# setup larlite
source $BMT_TAGGERDEV_BASEDIR/larlite/config/setup.sh

# setup laropencv
#source $BMT_TAGGERDEV_BASEDIR/LArOpenCV/setup_laropencv.sh

# setup Geo2D
#source $BMT_TAGGERDEV_BASEDIR/Geo2D/config/setup.sh

# setup LArCV
source $BMT_TAGGERDEV_BASEDIR/LArCV/configure.sh

# setup larlitecv
source $BMT_TAGGERDEV_BASEDIR/larlitecv/configure.sh

export LD_LIBRARY_PATH=.:${LD_LIBRARY_PATH}
