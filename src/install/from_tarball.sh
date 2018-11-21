#!/usr/bin/env bash

mkdir -p src/

pushd src/ || exit

# popoolation
wget \
    --continue \
    http://downloads.sourceforge.net/project/popoolation/popoolation_1.2.2.zip
unzip popoolation_1.2.2.zip

# popoolation2
wget \
    --continue \
    http://downloads.sourceforge.net/project/popoolation2/popoolation2_1201.zip
unzip popoolation2_1201.zip

popd || exit
