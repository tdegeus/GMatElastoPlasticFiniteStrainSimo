
# GMatElastoPlasticFiniteStrainSimo

[![Travis](https://travis-ci.com/tdegeus/GMatElastoPlasticFiniteStrainSimo.svg?branch=master)](https://travis-ci.com/tdegeus/GMatElastoPlasticFiniteStrainSimo)

Simo elasto-plastic model. An overview of the theory can be found in `docs/` in particular in this [PDF](docs/readme.pdf).

# Contents

<!-- MarkdownTOC levels="1,2" -->

- [Disclaimer](#disclaimer)
- [Implementation](#implementation)
    - [Overview](#overview)
    - [Conventions](#conventions)
    - [Example](#example)
    - [Debugging](#debugging)
- [Installation](#installation)
    - [C++ headers](#c-headers)
    - [Python module](#python-module)
- [Compiling](#compiling)
    - [Using CMake](#using-cmake)
    - [By hand](#by-hand)
    - [Using pkg-config](#using-pkg-config)
- [References / Credits](#references--credits)

<!-- /MarkdownTOC -->

# Disclaimer

This library is free to use under the [MIT license](https://github.com/tdegeus/GMatElastoPlasticFiniteStrainSimo/blob/master/LICENSE). Any additions are very much appreciated, in terms of suggested functionality, code, documentation, testimonials, word-of-mouth advertisement, etc. Bug reports or feature requests can be filed on [GitHub](https://github.com/tdegeus/GMatElastoPlasticFiniteStrainSimo). As always, the code comes with no guarantee. None of the developers can be held responsible for possible mistakes.

Download: [.zip file](https://github.com/tdegeus/GMatElastoPlasticFiniteStrainSimo/zipball/master) | [.tar.gz file](https://github.com/tdegeus/GMatElastoPlasticFiniteStrainSimo/tarball/master).

(c - [MIT](https://github.com/tdegeus/GMatElastoPlasticFiniteStrainSimo/blob/master/LICENSE)) T.W.J. de Geus (Tom) | tom@geus.me | www.geus.me | [github.com/tdegeus/GMatElastoPlasticFiniteStrainSimo](https://github.com/tdegeus/GMatElastoPlasticFiniteStrainSimo)

# Implementation

## Overview

The headers are meant to be self-explanatory, please inspect them:

* [Cartesian3d.h](include/GMatElastoPlasticFiniteStrainSimo/Cartesian3d.h)

## Conventions

Naming conventions

+   Functions whose name starts with a capital (e.g. `Stress`) return their result (allocating it internally).
+   Functions whose name starts with a small (e.g. `stress`) write to the, fully allocated, (last) input argument (avoiding re-allocation, but making the user responsible to do it properly).

Storage conventions

+   Scalar
    ```cpp
    double
    ```

+   2nd-order tensor
    ```cpp
    Tensor2 = xt::xtensor_fixed<double, xt::xshape<3,3>>
    ```

+   List *(a)* of second order tensors *(i,j)* : *A(a,i,j)*
    ```cpp
    xt::xtensor<double,3>
    ```

+   Matrix *(a,b)* of second order tensors *(i,j)* : *A(a,b,i,j)*
    ```cpp
    xt::xtensor<double,4>
    ```

## Example

Only a partial example is presented here, that is meant to understand the code's structure:

```cpp
#include <GMatElastoPlasticFiniteStrainSimo/Cartesian3d.h>

int main()
{
    // a single material point
    // - construct
    GMatElastoPlasticFiniteStrainSimo::Cartesian3d::Elastic elastic(K, G);
    GMatElastoPlasticFiniteStrainSimo::Cartesian3d::LinearHardening plastic(K, G, tauy0, H);
    // - set deformation gradient tensor (follows e.g. from FEM discretisation)
    GMatElastoPlasticFiniteStrainSimo::Tensor2 F;
    // - compute stress [allocate result]
    //   (N.B. returns the Cauchy stress)
    GMatElastoPlasticFiniteStrainSimo::Tensor2 Sig = elastic.Stress(F);
    ...
    // - compute stress [no allocation]
    //   (N.B. returns the Cauchy stress)
    elastic.stress(F, Sig); 
    ...

    // a matrix, of shape [nelem, nip], of material points
    // - construct
    GMatElastoPlasticFiniteStrainSimo::Cartesian3d::Elastic matrix(nelem, nip);
    // - set material
    matrix.setElastic(I, K, G);
    matrix.setLinearHardening(I, K, G, tauy0, H);
    // - set deformation gradient tensor (follows e.g. from FEM discretisation)
    xt::xtensor<double,4> F = xt::empty<double>({nelem, nip, 3ul, 3ul});
    ... 
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
>   See [Cartesian2d.h](include/GMatElastoPlasticFiniteStrainSimo/Cartesian3d.h) for more details.

## Debugging

To enable assertions define `GMATELASTOPLASTICFINITESTRAINSIMO_ENABLE_ASSERT` **before** including *GMatElastoPlasticFiniteStrainSimo* for the first time. 

Using *CMake* this can be done using the `GMatElastoPlasticFiniteStrainSimo::assert` target (see [below](#using-cmake)).

>   To also enable assertions of *xtensor* also define `XTENSOR_ENABLE_ASSERT` **before** including *xtensor* (and *GMatElastoPlasticFiniteStrainSimo*) for the first time. 
>   
>   Using *CMake* all assertions are enabled using the `GMatElastoPlasticFiniteStrainSimo::debug` target (see [below](#using-cmake)).

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

```bash
conda install -c conda-forge python-gmatelastoplasticfinitestrainsimo
```

Note that *xsimd* and hardware optimisations are not enabled. To enable them you have to compile on your system, as is discussed next.

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

>   You can also use a single call:
>   
>   ```
>   python -m pip install .
>   ```
>   
>   However potentially leading to less readable output.

# Compiling

## Using CMake

### Example

Using *GMatElastoPlasticFiniteStrainSimo* your `CMakeLists.txt` can be as follows

```cmake
cmake_minimum_required(VERSION 3.1)
project(example)
find_package(GMatElastoPlasticFiniteStrainSimo REQUIRED)
add_executable(example example.cpp)
target_link_libraries(example PRIVATE GMatElastoPlasticFiniteStrainSimo)
```

### Targets

The following targets are available:

*   `GMatElastoPlasticFiniteStrainSimo`
    Includes *GMatElastoPlasticFiniteStrainSimo* and the *xtensor* dependency.

*   `GMatElastoPlasticFiniteStrainSimo::assert`
    Enables assertions by defining `GMATELASTOPLASTICFINITESTRAINSIMO_ENABLE_ASSERT`.

*   `GMatElastoPlasticFiniteStrainSimo::debug`
    Enables all assertions by defining `GMATELASTOPLASTICFINITESTRAINSIMO_ENABLE_ASSERT` and `XTENSOR_ENABLE_ASSERT`.

*   `GMatElastoPlasticFiniteStrainSimo::compiler_warings`
    Enables compiler warnings (generic).

### Optimization

It is advised to think about compiler optimization and about enabling *xsimd. In *CMake* this can be done using the `xtensor::optimize` and `xtensor::use_xsimd` targets. The above example then becomes:

```cmake
cmake_minimum_required(VERSION 3.1)
project(example)
find_package(GMatElastoPlasticFiniteStrainSimo REQUIRED)
add_executable(example example.cpp)
target_link_libraries(example PRIVATE 
    GMatElastoPlasticFiniteStrainSimo 
    xtensor::optimize 
    xtensor::use_xsimd)
```

See the [documentation of xtensor](https://xtensor.readthedocs.io/en/latest/) concerning optimization.

## By hand

Presuming that the compiler is `c++`, compile using:

```
c++ -I/path/to/GMatElastoPlasticFiniteStrainSimo/include ...
```

Note that you have to take care of the *xtensor* dependency, the C++ version, optimization, enabling *xsimd*, ...

## Using pkg-config

Presuming that the compiler is `c++`, compile using:

```
c++ `pkg-config --cflags GMatElastoPlasticFiniteStrainSimo` ...
```

Note that you have to take care of the *xtensor* dependency, the C++ version, optimization, enabling *xsimd*, ...

# References / Credits

*   The model is described in *M.G.D. Geers (2004). Finite strain logarithmic hyperelasto-plasticity with softening: a strongly non-local implicit gradient framework. Computer Methods in Applied Mechanics and Engineering, 193(30–32), 3377–3401, [doi: 10.1016/j.cma.2003.07.014](https://doi.org/10.1016/j.cma.2003.07.014)*.

*   [xtensor](https://github.com/QuantStack/xtensor) is used under the hood.

*   The eigenvalue and eigenvector decomposition is done by an algorithm due to [Joachim Kopp](https://www.mpi-hd.mpg.de/personalhomes/globes/3x3/). See also: *J. Kopp (2008). Efficient numerical diagonalization of hermitian 3x3 matrices. Int. J. Mod. Phys. C. 19:523-548, arXiv: [physics/0610206](https://arxiv.org/abs/physics/0610206)*.

