# HAVOC - Handcoded Assembly Optimisations for VideO Codecs

## Introduction

Havoc is a collection of permissively-licensed optimised DSP functions for video coding. Instead of using intrinsics or a build-time assembler, Havoc makes use of a just-in-time (JIT) assembler. 

Havoc has a C-language API and makes use of C++11 internally.


## Clone and build

Because the project includes Xbyak as a submodule, be sure to `git clone --recursive` or to use `git submodule` when making a working copy.

CMake is used to build the project for Linux and Windows. 


### 64-bit Linux 

The following is tested on Ubuntu 14.04 with g++ 4.8.4. First configure and create a makefile.

    mkdir buildLinux
    cd buildLinux
    cmake -DCMAKE_C_COMPILER=/usr/bin/gcc -DCMAKE_CXX_COMPILER=/usr/bin/g++ ..

Then use the generated makefile to build the project. A static library and a test executable are emitted - this performs various self tests.

```
make
src/bin/havoc
```


### 64-bit Windows (Visual Studio 2015)

Download the CMake GUI for Windows and start it.  Enter your Havoc working copy's root folder in "where is the source code." Choose any other folder, for example a `buildWindows` subdir  for "where to build the binaries."  

Click "Configure" and select your version of Visual Studio - be sure to select a "Win64" generator - and click Finish. Click "Generate" to create Visual Studio project and solution files that can be opened in your IDE.


## Just-in-Time assembler

Havoc uses just-in-time assembler [Xbyak](https://github.com/herumi/xbyak), . This means that the executable code is assembled at runtime, just before it runs. JIT assembly has the following advantages:

* All sorts of interesting optimisations and adaptations are possible. For example, video-stream-specific parameters could be embedded in the code and code could be tuned for a specific machine's measured characteristics.
* No special tool (assembler) is needed to build the project. That means one less dependency to find, version and install.
* Developers can write assembly code. Xbyak offers a syntactically sweet domain-specific language having a syntax very similar to assembly code.
* Preprocesing can be done with idiomatic C++ constructs such as functions and loops - there is no need for complex assembler-specific macros. Parameters are regular C++ variables and regular debugging can be used to debug the assembly.

There are a few drawbacks to the current JIT implementation for which workaround are possible:

* Debug symbols are unavailable - debuggers cannot know or show the name of the JIT-assembled function that is running. This also affects profiler - only the code address is known, not the function name.
* Assembly takes a small but non-zero amount of time during program initialisation.


## Continuous integration

The Linux build of this project is continuously integrated by Travis-CI. Because Travis-CI is hosted on Amazon EC2's older C3 instances, it does not support the AVX2 instruction set.  To get around this, an emulator (Intel SDE) is used to emulate a more recent processor. This allows all optimised functions to be tested, slowly.  The project CI has the following limitations:

* Only 64-bit Linux is tested currently
* Performance results are inaccurate (and completely wrong for AVX2 functions as these are emulated)

The results and output of the Travis-CI process can be seen here: [https://travis-ci.org/kupix/havoc](https://travis-ci.org/kupix/havoc)


## Credits

Some functions are derived from the [f265](http://vantrix.com/f-265/) project, others from Google's [WebM](http://www.webmproject.org/) project. Further functions were written from scratch. See COPYING and the source code for further details.

