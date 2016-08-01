# Turing codec

**The turing codec (see http://turingcodec.org) is a software HEVC codec**

The code is published under the GPL version 2 licence.  Commercial support and intellectual property rights for the Turing codec are also available under a proprietary license. 
For more information, contact us at info @ turingcodec.org.

## Building for the first time

If any problems are discovered during this process, please edit the project's `README.md` file to change or improve these instructions.  This way, others will not have to suffer the same issues.

First, clone the code from GitHub:

```
git clone https://github.com/bbc/turingcodec.git
cd turingcodec
git checkout {stable|develop}
```

### Microsoft Windows (with Visual Studio 2015)
#### Prequisites
  * CMake for Windows (version 3.4.1 used when writing these instructions).

#### Building the project executable
* Open CMake for Windows
 * Set "Where is the source code" as the folder containing this README.md file
 * Set "Where to build the binaries" to any new folder location you prefer
 * Click 'Configure' and select the Win64 configuration for Visual Studio 2015
 * Click 'Generate'
 * Project file "Turing.sln" should be emitted into the new build folder location - open this in Visual Studio and build either as `Debug` or `Release

### Linux
#### Prerequisites
CMake

Run the following commands to install necessary tools and libraries:


```
sudo apt-get update
sudo apt-get install git-core
sudo apt-get install g++
sudo apt-get install make
```

Make a local clone of this repository and `cd` to its root folder.

#### Build and test

```
mkdir -p build/debug build/release
cd build/release
cmake ../../
```

Build the Turing codec executable:

```
make
```
