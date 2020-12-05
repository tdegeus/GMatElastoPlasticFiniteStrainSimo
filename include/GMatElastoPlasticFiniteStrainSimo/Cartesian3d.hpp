/*

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastoPlasticFiniteStrainSimo

*/

#ifndef GMATELASTOPLASTICFINITESTRAINSIMO_CARTESIAN3D_HPP
#define GMATELASTOPLASTICFINITESTRAINSIMO_CARTESIAN3D_HPP

#include "Cartesian3d.h"

namespace GMatElastoPlasticFiniteStrainSimo {
namespace Cartesian3d {

namespace detail {

    template <class T>
    struct equiv_impl : GMatTensor::Cartesian3d::detail::equiv_impl<T>
    {
        using value_type = typename T::value_type;
        using GMatTensor::Cartesian3d::detail::equiv_impl<T>::rank;
        using GMatTensor::Cartesian3d::detail::equiv_impl<T>::toMatrixSize;
        using GMatTensor::Cartesian3d::detail::equiv_impl<T>::toMatrixShape;
        using GMatTensor::Cartesian3d::detail::equiv_impl<T>::toShape;

        static void strain_no_alloc(const T& A, xt::xtensor<value_type, rank>& ret)
        {
            GMATTENSOR_ASSERT(xt::has_shape(A, toShape(A.shape())));
            GMATTENSOR_ASSERT(xt::has_shape(A, B.shape()));
            #pragma omp parallel for
            for (size_t i = 0; i < toMatrixSize(A.shape()); ++i) {
                std::array<double, 9> B;
                std::array<double, 3> val;
                std::array<double, 9> vec;
                GMatTensor::Cartesian3d::pointer::A2_dot_A2T(&A.data()[i * 9], &B[0]);
                GMatTensor::Cartesian3d::pointer::eigs(&B[0], &vec[0], &val[0]);
                for (size_t j = 0; j < 3; ++j) {
                    val[j] = 0.5 * std::log(val[j]);
                }
                GMatTensor::Cartesian3d::pointer::from_eigs(&vec[0], &val[0], &ret.data()[i * 9]);
            }
        }

        static auto strain_alloc(const T& A)
        {
            xt::xtensor<value_type, rank> B = xt::empty<value_type>(A.shape());
            strain_no_alloc(A, B);
            return B;
        }
    };

} // namespace detail

template <class T, class U>
inline void epseq(const T& A, U& B)
{
    GMatTensor::Cartesian3d::equivalent_deviatoric(A, B);
    B *= std::sqrt(2.0 / 3.0);
}

template <class T>
inline auto Epseq(const T& A)
{
    return xt::eval(std::sqrt(2.0 / 3.0) * GMatTensor::Cartesian3d::Equivalent_deviatoric(A));
}

template <class T, class U>
inline void sigeq(const T& A, U& B)
{
    GMatTensor::Cartesian3d::equivalent_deviatoric(A, B);
    B *= std::sqrt(1.5);
}

template <class T>
inline auto Sigeq(const T& A)
{
    return xt::eval(std::sqrt(1.5) * GMatTensor::Cartesian3d::Equivalent_deviatoric(A));
}

template <class T, class U>
inline void strain(const T& A, U& B)
{
    return detail::equiv_impl<T>::strain_no_alloc(A, B);
}

template <class T>
inline auto Strain(const T& A)
{
    return detail::equiv_impl<T>::strain_alloc(A);
}

} // namespace Cartesian3d
} // namespace GMatElastoPlasticFiniteStrainSimo

#endif
