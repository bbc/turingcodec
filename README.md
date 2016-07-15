# Turing codec

**The turing codec (see http://turingcodec.org) is a software HEVC codec**

The code is published under the GPL version 2 licence.  Commercial support and intellectual property rights for the Turing codec are also available under a proprietary license. 
For more information, contact us at info @ turingcodec.org.

## Code line policy
General points:

  * read about ["GitFlow"](http://nvie.com/posts/a-successful-git-branching-model/).
  * create a branch for each task - most likely the most suitable branch point will be the current tip of the `develop` branch.
  * commit work-in-progress to that branch using a meaningful but succinct commit message. Push to github for backup.
  * keep changes minimal and on-topic.

### Policy for merging a branch to `develop`
  * ensure new changes are appropriately:
    * mature,
    * complete,
    * tested with tests built into the code,
    * documented and
    * usable via command options/API.
  * initiate a GitHub pull request to inform others  https://help.github.com/articles/using-pull-requests.
  * before pushing merged changes to `develop`
    * ensure a new clone builds on Windows x64 Release, Windows x64 Debug and gcc
    * ensure `turing-exe smoke` runs without error in Windows x64 Release, Windows x64 Debug and gcc

Note that `turing-exe smoke` and the `smoke()` function are part of the test plan. This test is intended to run in 5-10 minutes or so and should be adapted to cover as many features as possible within that time.


## Building for the first time

If any problems are discovered during this process, please edit the project's `README.md` file to change or improve these instructions.  This way, others will not have to suffer the same issues.

First, clone the code from GitHub:

```
git clone https://github.com/bbc/turingcodec.git
cd turingcodec
git checkout {stable|develop}
```

### Microsoft Windows (with Visual Studio 2013 or 2015)
#### Prequisites
  * CMake for Windows (version 3.4.1 used when writing these instructions).

#### Building the project executable
* Open CMake for Windows
 * Set "Where is the source code" as the folder containing this README.md file
 * Set "Where to build the binaries" to any new folder location you prefer
 * Click 'Configure' and select the Win64 configuration corresponding to your Visual Studio installation (either VS 2013 or VS 2015)
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

Build the TUring codec executable:

```
make
```
