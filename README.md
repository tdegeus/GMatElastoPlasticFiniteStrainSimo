# GMatElastoPlasticFiniteStrainSimo

Simo elasto-plastic model. An overview of the theory can be found in `docs/theory` in particular in this (PDF)[docs/theory/main.pdf].

# Contents

<!-- MarkdownTOC -->

- [Implementation](#implementation)
- [Installation](#installation)
    - [Linux / macOS](#linux--macos)
        - [Install system-wide \(depends on your privileges\)](#install-system-wide-depends-on-your-privileges)
        - [Install in custom location \(user\)](#install-in-custom-location-user)
- [References / Credits](#references--credits)

<!-- /MarkdownTOC -->

# Implementation

```cpp
#include <GMatElastoPlasticFiniteStrainSimo/Cartesian3d.h>

int main()
{
    // a single material point
    // - create class
    GMatElastoPlasticFiniteStrainSimo::Cartesian3d::Elastic elastic(K, G);
    GMatElastoPlasticFiniteStrainSimo::Cartesian3d::LinearHardening plastic(K, G, tauy0, H);
    // - compute stress [allocate result]
    Tau = elastic.Stress(F);
    ...
    // - compute stress [no allocation]
    elastic.stress(F, Tau); 
    ...

    // a "matrix" of material points
    // - create class
    GMatElastoPlasticFiniteStrainSimo::Cartesian3d::Elastic matrix(nelem, nip);
    // - set material
    matrix.setElastic(I, K, G);
    matrix.setLinearHardening(I, K, G, tauy0, H);
    // - compute stress [allocate result]
    Tau = matrix.Stress(F);
    ...
    // - compute stress [no allocation]
    matrix.stress(F, Tau); 
    ...
}
```

# Installation

## Linux / macOS

### Install system-wide (depends on your privileges)

1.  Proceed to a (temporary) build directory. For example

    ```bash
    cd /path/to/GMatElastoPlasticFiniteStrainSimo
    mkdir build
    cd build
    ```

2.  'Install' `GMatElastoPlasticFiniteStrainSimo`:

    ```bash
    cmake .. 
    make install
    ```

> One usually does not need any compiler arguments after following this protocol.

### Install in custom location (user)

1.  Proceed to a (temporary) build directory. For example

    ```bash
    cd /path/to/GMatElastoPlasticFiniteStrainSimo
    mkdir build
    cd build
    ```

2.  'Install' `GMatElastoPlasticFiniteStrainSimo`, to install it in a custom location

    ```bash
    mkdir /custom/install/path
    cmake .. -DCMAKE_INSTALL_PREFIX:PATH=/custom/install/path
    make install
    ```

3.  Add the following path to your ``~/.bashrc`` (or ``~/.zshrc``):

    ```bash
    export PKG_CONFIG_PATH=/custom/install/path/share/pkgconfig:$PKG_CONFIG_PATH
    export CPLUS_INCLUDE_PATH=$HOME/custom/install/path/include:$CPLUS_INCLUDE_PATH
    ```

> One usually has to inform the CMake or the compiler about `${CPLUS_INCLUDE_PATH}`.

# References / Credits

*   [xtensor](https://github.com/QuantStack/xtensor) is used under the hood.

*   The model is described in *M.G.D. Geers (2004). Finite strain logarithmic hyperelasto-plasticity with softening: a strongly non-local implicit gradient framework. Computer Methods in Applied Mechanics and Engineering, 193(30–32), 3377–3401.* https://doi.org/10.1016/j.cma.2003.07.014

*   The eigenvalue and eigenvector decomposition is done by an algorithm due to [Joachim Kopp](https://www.mpi-hd.mpg.de/personalhomes/globes/3x3/). See also *J. Kopp (2008). Efficient numerical diagonalization of hermitian 3x3 matrices. Int. J. Mod. Phys. C. 19:523-548. arXiv: [physics/0610206](https://arxiv.org/abs/physics/0610206)*.

