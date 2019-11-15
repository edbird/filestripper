#!/bin/bash
mkdir -p build
cd build
#export PKG_CONFIG_PATH="/usr/lib/x86_64-linux-gnu/pkgconfig:/usr/share/pkgconfig"
cmake -DCMAKE_BUILD_TYPE=Debug ..
RET="$?"
if [ "$RET" -eq 0 ]
then
    make
    RET="$?"
else
    echo "cmake failed"
fi
cd ..
#return "$RET"
