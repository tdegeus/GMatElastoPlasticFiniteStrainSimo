/* =================================================================================================

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastoPlasticQPot

================================================================================================= */

#ifndef GMATELASTOPLASTICFINITESTRAINSIMO_CONFIG_H
#define GMATELASTOPLASTICFINITESTRAINSIMO_CONFIG_H

// -------------------------------------------------------------------------------------------------

#include <tuple>
#include <stdexcept>
#include <limits>
#include <math.h>
#include <iostream>
#include <vector>
#include <tuple>
#include <xtensor/xarray.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xnoalias.hpp>
#include <xtensor/xfixed.hpp>
#include <xtensor/xadapt.hpp>
#include <xtensor/xview.hpp>
#include <xtensor/xio.hpp>
#include <xtensor/xoperation.hpp>
#include <xtensor/xsort.hpp>
#include <xtensor/xmath.hpp>

// -------------------------------------------------------------------------------------------------

#ifndef NDEBUG
#define GMATELASTOPLASTICFINITESTRAINSIMO_ENABLE_ASSERT
#endif

#ifdef GMATELASTOPLASTICFINITESTRAINSIMO_ENABLE_ASSERT
#define GMATELASTOPLASTICFINITESTRAINSIMO_ASSERT(expr) GMATELASTOPLASTICFINITESTRAINSIMO_ASSERT_IMPL(expr, __FILE__, __LINE__)
#define GMATELASTOPLASTICFINITESTRAINSIMO_ASSERT_IMPL(expr, file, line)                                                   \
    if (!(expr))                                                                                                          \
    {                                                                                                                     \
        throw std::runtime_error(std::string(file) + ':' + std::to_string(line) + ": assertion failed (" #expr ") \n\t"); \
    }
#else
#define GMATELASTOPLASTICFINITESTRAINSIMO_ASSERT(expr)
#endif

// -------------------------------------------------------------------------------------------------

#define GMATELASTOPLASTICFINITESTRAINSIMO_WORLD_VERSION 0
#define GMATELASTOPLASTICFINITESTRAINSIMO_MAJOR_VERSION 0
#define GMATELASTOPLASTICFINITESTRAINSIMO_MINOR_VERSION 2

#define GMATELASTOPLASTICFINITESTRAINSIMO_VERSION_AT_LEAST(x,y,z) \
  (GMATELASTOPLASTICFINITESTRAINSIMO_WORLD_VERSION>x || (GMATELASTOPLASTICFINITESTRAINSIMO_WORLD_VERSION>=x && \
  (GMATELASTOPLASTICFINITESTRAINSIMO_MAJOR_VERSION>y || (GMATELASTOPLASTICFINITESTRAINSIMO_MAJOR_VERSION>=y && \
                                         GMATELASTOPLASTICFINITESTRAINSIMO_MINOR_VERSION>=z))))

#define GMATELASTOPLASTICFINITESTRAINSIMO_VERSION(x,y,z) \
  (GMATELASTOPLASTICFINITESTRAINSIMO_WORLD_VERSION==x && \
   GMATELASTOPLASTICFINITESTRAINSIMO_MAJOR_VERSION==y && \
   GMATELASTOPLASTICFINITESTRAINSIMO_MINOR_VERSION==z)

// -------------------------------------------------------------------------------------------------

#include <xtensor/xtensor_simd.hpp>

#ifndef GMATELASTOPLASTICFINITESTRAINSIMO_ALIGNMENT
    #ifdef XTENSOR_USE_XSIMD
        #define GMATELASTOPLASTICFINITESTRAINSIMO_ALIGNMENT XSIMD_DEFAULT_ALIGNMENT
    #else
        #define GMATELASTOPLASTICFINITESTRAINSIMO_ALIGNMENT 0
    #endif
#endif

// -------------------------------------------------------------------------------------------------

#endif
