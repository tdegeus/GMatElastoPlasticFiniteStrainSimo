/*

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastoPlasticFiniteStrainSimo

*/

#ifndef GMATELASTOPLASTICFINITESTRAINSIMO_CARTESIAN3D_HPP
#define GMATELASTOPLASTICFINITESTRAINSIMO_CARTESIAN3D_HPP

#include "Cartesian3d.h"

namespace GMatElastoPlasticFiniteStrainSimo {
namespace Cartesian3d {

template <class T, class U>
inline void epseq(const T& A, U& ret)
{
    GMatTensor::Cartesian3d::norm_deviatoric(A, ret);
    ret *= std::sqrt(2.0 / 3.0);
}

template <class T>
inline auto Epseq(const T& A)
{
    return xt::eval(std::sqrt(2.0 / 3.0) * GMatTensor::Cartesian3d::Norm_deviatoric(A));
}

template <class T, class U>
inline void sigeq(const T& A, U& ret)
{
    GMatTensor::Cartesian3d::norm_deviatoric(A, ret);
    ret *= std::sqrt(1.5);
}

template <class T>
inline auto Sigeq(const T& A)
{
    return xt::eval(std::sqrt(1.5) * GMatTensor::Cartesian3d::Norm_deviatoric(A));
}

template <class T, class U>
inline void strain(const T& A, U& ret)
{
    GMatTensor::Cartesian3d::A2_dot_A2T(A, ret);
    GMatTensor::Cartesian3d::logs(ret, ret);
    ret *= 0.5;
}

template <class T>
inline auto Strain(const T& A)
{
    auto ret = GMatTensor::Cartesian3d::A2_dot_A2T(A);
    GMatTensor::Cartesian3d::logs(ret, ret);
    ret *= 0.5;
    return ret;
}

} // namespace Cartesian3d
} // namespace GMatElastoPlasticFiniteStrainSimo

#endif
