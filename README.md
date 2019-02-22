# OpenSMOKEppPlastics
C++ library for numerical modeling of thermal degradation of plastics

### Compulsory external dependencies ###
- Eigen (http://eigen.tuxfamily.org/index.php?title=Main_Page)
- CMake (https://cmake.org/)

### Optional external dependencies ###
- Intel MKL (https://software.intel.com/en-us/intel-mkl)
- CppCheck (http://cppcheck.sourceforge.net/)
- GoogleTest (https://github.com/google/googletest)
- RapidXML (http://rapidxml.sourceforge.net/)
- Boost C++ (http://www.boost.org/)

Source code and compilation
---------------------------
This section provides the basic instructions for compiling OpenSMOKEppPlastics under cmake. In the wiki and discussion forum you may find more detailed instructions for specific platforms. 

### Obtaining the source code ###
The source code of OpenSMOKEppPlastics project is available at:  
https://github.com/acuoci/OpenSMOKEppPlastics

### 1. Create directory where to download thing ###
`mkdir opensmokeppplastics`  
`cd opensmokeppplastics`  

### 2. Clone the git repository ###
`git clone git://www.github.com/acuoci/OpenSMOKEppPlastics`

### 3. Create build directory ###
`mkdir build`  
`cd build`  

### 4. Let cmake generate makefiles ###
We use the CMake build system, but only to build the documentation and unit-tests, and to automate installation. You can invoke the `cmake` command, together with proper options. Section 8 (below) reports the complete list of options which are available.  
The minimum set of options to be provided consists in specyfing the path where to install the library, the build type (Release vs Debug), and the path to the Eigen3 library:  
`cmake ../opensmokeppplastics/ -DCMAKE_INSTALL_PREFIX:PATH=/path/where/to/install -DCMAKE_BUILD_TYPE=Release -DEIGEN3=/path/to/eigen-3.3`

It is recommended to compile also the tests, to check the correct compilation of the library. In this case, you need to relace the above instruction with the following one, in which the path to the GoogleTest library is provided:   
`cmake ../opensmokeppplastics/ -DCMAKE_INSTALL_PREFIX:PATH=/path/where/to/install -DCMAKE_BUILD_TYPE=Release -DEIGEN3=/path/to/eigen-3.3 -DENABLE_GTEST=ON -DGTEST_ROOT=/path/to/googletest`

Alternatively, you can tune the compilation parameters graphically with:  
`ccmake .`  
or  
`cmake-gui .`  

### 5. Compile ### 
If no errors are produced, you can proceed to compilation:  
`make install`  
or alternatively, you can compile in parallel (2 procs, for example)  
`make -j2 install`

### 6. Testing the compilation ###
If the compilation of tests was enabled (see Section 4), you can proceed to testing (2 procs, for example):  
`ctest -j2`

### 7. Installation ### 
You can now proceed to the final installation to the folder specified in Section 4:  
`make install`  

### Alternative/Additional compilation options ### 
#### Intel MKL ####
It is highly recommanded, for improved speed, to enable the support to Intel MKL libraries:  
`-DENABLE_MKL=ON -DMKL_DIR=/path/to/intel/mkl`  

#### Compilers ####
You can specify the compilers to be used through the following options (for example, Intel compilers):  
`-DCMAKE_C_COMPILER=icc`  
`-DCMAKE_CXX_COMPILER=icpc`  
`-DCMAKE_Fortran_COMPILER=ifort`  

#### Doxygen documentation ####
Doxygen documentation can be generated through the following option:  
`-DENABLE_DOXYGEN=ON`

#### CppCheck ####
Code can be checked through the CppCheck utility by enabling the following option at compilation time:  
`-DENABLE_CPPCHECK=ON`  

After compilation, `cppcheck` utility can be automatically invoked through:  
`make analysis`

#### Code coverage ####
Code coveraged can be enabled through the following option at compilation time:  
`-DENABLE_CODECOV=ON`  

After compilation, `gcov` can be automatically invoked through:  
`make test_codecov`

The code coverage is available only if gcc compilers are adopted.  
