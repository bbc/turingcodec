#!/bin/bash -xe

cleanup() {
    cat $TRAVIS_BUILD_DIR/turingcodec/build/CMakeFiles/CMakeOutput.log
}
trap cleanup SIGHUP SIGINT SIGTERM

mkdir -p build

cd build
cmake -DCMAKE_CXX_COMPILER=$(which ${CXX:-g++}) -DCMAKE_C_COMPILER=$(which ${CC:-gcc}) ..
cat $TRAVIS_BUILD_DIR/build/CMakeFiles/CMakeOutput.log || true
cat $TRAVIS_BUILD_DIR/build/CMakeFiles/CMakeError.log || true
make VERBOSE=1
cd ..

./build/turing/turing signature ./test/
./build/turing/turing testdecode
