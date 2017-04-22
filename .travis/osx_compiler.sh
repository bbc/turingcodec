#!/bin/bash -x

case "$CC" in
    gcc-4.9)
        brew install homebrew/versions/gcc49
        ;;
    gcc-5)
        brew install homebrew/versions/gcc5
        ;;
    gcc-6)
        brew install homebrew/versions/gcc6
        ;;
esac
