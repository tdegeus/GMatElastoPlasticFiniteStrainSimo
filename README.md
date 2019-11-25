
# GMatElastoPlasticFiniteStrainSimo

[![Travis](https://travis-ci.com/tdegeus/GMatElastoPlasticFiniteStrainSimo.svg?branch=master)](https://travis-ci.com/tdegeus/GMatElastoPlasticFiniteStrainSimo)

Simo elasto-plastic model. An overview of the theory can be found in `docs/` in particular in this [PDF](docs/readme.pdf).

# Contents

<!-- MarkdownTOC levels="1,2" -->

- [Implementation](#implementation)
- [Installation](#installation)
    - [C++ headers](#c-headers)
    - [Python module](#python-module)
- [Compiling](#compiling)
    - [By hand](#by-hand)
    - [Using pkg-config](#using-pkg-config)
    - [Using `CMakeLists.txt`](#using-cmakeliststxt)
- [References / Credits](#references--credits)

<!-- /MarkdownTOC -->

# Implementation

The headers are meant to be self-explanatory, please check them out:

* [Cartesian3d.h](include/GMatElastoPlasticFiniteStrainSimo/Cartesian3d.h)

Only a tiny example is presented here, that is meant to understand the code's structure:

```cpp
#include <GMatElastoPlasticFiniteStrainSimo/Cartesian3d.h>

int main()
{
    // a single material point
    // - create class
    GMatElastoPlasticFiniteStrainSimo::Cartesian3d::Elastic elastic(K, G);
    GMatElastoPlasticFiniteStrainSimo::Cartesian3d::LinearHardening plastic(K, G, tauy0, H);
    // - compute stress [allocate result]
    //   (N.B. returns the Cauchy stress)
    Sig = elastic.Stress(F);
    ...
    // - compute stress [no allocation]
    //   (N.B. returns the Cauchy stress)
    elastic.stress(F, Sig); 
    ...

    // a "matrix" of material points
    // - create class
    GMatElastoPlasticFiniteStrainSimo::Cartesian3d::Elastic matrix(nelem, nip);
    // - set material
    matrix.setElastic(I, K, G);
    matrix.setLinearHardening(I, K, G, tauy0, H);
    // - compute stress [allocate result]
    //   (N.B. returns the Cauchy stress)
    Sig = matrix.Stress(F);
    ...
    // - compute stress [no allocation]
    //   (N.B. returns the Cauchy stress)
    matrix.stress(F, Sig); 
    ...
}
```

# Installation

## C++ headers

### Using conda

```bash
conda install -c conda-forge gmatelastoplasticfinitestrainsimo
```

### From source

```bash
# Download GMatElastoPlasticFiniteStrainSimo
git checkout https://github.com/tdegeus/GMatElastoPlasticFiniteStrainSimo.git
cd GMatElastoPlasticFiniteStrainSimo

# Install headers, CMake and pkg-config support
cmake .
make install
```

## Python module

### Using conda

> Warning: this has the disadvantage of xsimd optimisation being switched off

```bash
conda install -c conda-forge python-gmatelastoplasticfinitestrainsimo
```

### From source

> To get the prerequisites you *can* use conda
> 
> ```bash
> conda install -c conda-forge pyxtensor
> conda install -c conda-forge xsimd
> ```

```bash
# Download GMatElastoPlasticFiniteStrainSimo
git checkout https://github.com/tdegeus/GMatElastoPlasticFiniteStrainSimo.git
cd GMatElastoPlasticFiniteStrainSimo

# Compile and install the Python module
python setup.py build
python setup.py install
```

# Compiling

## By hand

Presuming that the compiler is `c++`, compile using:

```
c++ -I/path/to/GMatElastoPlasticFiniteStrainSimo/include ...
```

## Using pkg-config

Presuming that the compiler is `c++`, compile using:

```
c++ `pkg-config --cflags GMatElastoPlasticFiniteStrainSimo` ...
```

## Using `CMakeLists.txt`

Using *GMatElastoPlasticFiniteStrainSimo* the `CMakeLists.txt` can be as follows

```cmake
cmake_minimum_required(VERSION 3.1)

project(example)

find_package(xtensor REQUIRED)
find_package(GMatElastoPlasticFiniteStrainSimo REQUIRED)

add_executable(example example.cpp)

target_link_libraries(example
    PRIVATE
    xtensor
    GMatElastoPlasticFiniteStrainSimo)
```

Compilation can then proceed using 

```bash
cmake .
make
```

# References / Credits

*   The model is described in *M.G.D. Geers (2004). Finite strain logarithmic hyperelasto-plasticity with softening: a strongly non-local implicit gradient framework. Computer Methods in Applied Mechanics and Engineering, 193(30–32), 3377–3401, [doi: 10.1016/j.cma.2003.07.014](https://doi.org/10.1016/j.cma.2003.07.014)*.

*   [xtensor](https://github.com/QuantStack/xtensor) is used under the hood.

*   The eigenvalue and eigenvector decomposition is done by an algorithm due to [Joachim Kopp](https://www.mpi-hd.mpg.de/personalhomes/globes/3x3/). See also: *J. Kopp (2008). Efficient numerical diagonalization of hermitian 3x3 matrices. Int. J. Mod. Phys. C. 19:523-548, arXiv: [physics/0610206](https://arxiv.org/abs/physics/0610206)*.

