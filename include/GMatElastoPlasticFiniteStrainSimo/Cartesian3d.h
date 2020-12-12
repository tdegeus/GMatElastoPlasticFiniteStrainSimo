/*

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastoPlasticFiniteStrainSimo

*/

#ifndef GMATELASTOPLASTICFINITESTRAINSIMO_CARTESIAN3D_H
#define GMATELASTOPLASTICFINITESTRAINSIMO_CARTESIAN3D_H

#include <GMatTensor/Cartesian3d.h>

#include "config.h"

namespace GMatElastoPlasticFiniteStrainSimo {
namespace Cartesian3d {

// Unit tensors

using GMatTensor::Cartesian3d::O2;
using GMatTensor::Cartesian3d::O4;
using GMatTensor::Cartesian3d::I2;
using GMatTensor::Cartesian3d::II;
using GMatTensor::Cartesian3d::I4;
using GMatTensor::Cartesian3d::I4rt;
using GMatTensor::Cartesian3d::I4s;
using GMatTensor::Cartesian3d::I4d;

// Tensor decomposition

using GMatTensor::Cartesian3d::hydrostatic;
using GMatTensor::Cartesian3d::Hydrostatic;
using GMatTensor::Cartesian3d::deviatoric;
using GMatTensor::Cartesian3d::Deviatoric;

// Equivalent strain

template <class T, class U>
inline void epseq(const T& A, U& ret);

template <class T>
inline auto Epseq(const T& A);

// Equivalent stress

template <class T, class U>
inline void sigeq(const T& A, U& ret);

template <class T>
inline auto Sigeq(const T& A);

// Deformation gradient tensor -> strain = 0.5 * ln(B), B = F . F^T

template <class T, class U>
inline void strain(const T& A, U& ret);

template <class T>
inline auto Strain(const T& A);

// Material point

class Elastic
{
public:
    Elastic() = default;
    Elastic(double K, double G);

    double K() const;
    double G() const;

    template <class T> void setDefGrad(const T& arg, bool compute_tangent = true);
    template <class T> void defGrad(T& ret) const;
    template <class T> void stress(T& ret) const;
    template <class T> void tangent(T& ret) const;

    template <class T> void setDefGradPtr(const T* arg, bool compute_tangent = true);
    template <class T> void defGradPtr(T* ret) const;
    template <class T> void stressPtr(T* ret) const;
    template <class T> void tangentPtr(T* ret) const;

    xt::xtensor<double, 2> DefGrad() const;
    xt::xtensor<double, 2> Stress() const;
    xt::xtensor<double, 4> Tangent() const;

private:
    double m_K; // bulk modulus
    double m_G; // shear modulus
    std::array<double, 9> m_F;   // deformation gradient tensor [xx, xy, xz, yx, yy, yz, zx, zy, zz]
    std::array<double, 9> m_Sig; // stress tensor ,,
    xt::xtensor<double, 4> m_C;  // tangent
};

// Material point

class LinearHardening
{
public:
    LinearHardening() = default;
    LinearHardening(double K, double G, double tauy0, double H);

    double K() const;
    double G() const;
    double tauy0() const;
    double H() const;
    double epsp() const;
    void increment(); // update history

    template <class T> void setDefGrad(const T& arg, bool compute_tangent = true);
    template <class T> void defGrad(T& ret) const;
    template <class T> void stress(T& ret) const;
    template <class T> void tangent(T& ret) const;

    template <class T> void setDefGradPtr(const T* arg, bool compute_tangent = true);
    template <class T> void defGradPtr(T* ret) const;
    template <class T> void stressPtr(T* ret) const;
    template <class T> void tangentPtr(T* ret) const;

    xt::xtensor<double, 2> DefGrad() const;
    xt::xtensor<double, 2> Stress() const;
    xt::xtensor<double, 4> Tangent() const;

private:
    double m_K;      // bulk modulus
    double m_G;      // shear modulus
    double m_tauy0;  // initial yield stress
    double m_H;      // hardening modulus
    double m_epsp;   // plastic strain
    double m_epsp_t; // ,, previous time step
    std::array<double, 9> m_F;    // deformation gradient tensor [xx, xy, xz, yx, yy, yz, zx, zy, zz]
    std::array<double, 9> m_F_t;  // ,, previous time step
    std::array<double, 9> m_Be;   // elastic Finger tensor ,,
    std::array<double, 9> m_Be_t; // ,, previous time step
    std::array<double, 9> m_Sig;  // stress tensor ,,
    xt::xtensor<double, 4> m_C;   // tangent
};

// Material identifier

struct Type {
    enum Value {
        Unset,
        Elastic,
        LinearHardening,
    };
};

// Array of material points

template <size_t N>
class Array : public GMatTensor::Cartesian3d::Array<N>
{
public:
    using GMatTensor::Cartesian3d::Array<N>::rank;

    // Constructors

    Array() = default;
    Array(const std::array<size_t, N>& shape);

    // Overloaded methods:
    // - "shape"
    // - unit tensors: "I2", "II", "I4", "I4rt", "I4s", "I4d"

    // Type

    xt::xtensor<size_t, N> type() const;
    xt::xtensor<size_t, N> isElastic() const;
    xt::xtensor<size_t, N> isPlastic() const;
    xt::xtensor<size_t, N> isLinearHardening() const;

    // Parameters

    xt::xtensor<double, N> K() const;
    xt::xtensor<double, N> G() const;
    xt::xtensor<double, N> epsp() const;

    // Set parameters for a batch of points
    // (uniform for all points specified: that have "I(i, j) == 1")

    void setElastic(
        const xt::xtensor<size_t, N>& I,
        double K,
        double G);

    void setLinearHardening(
        const xt::xtensor<size_t, N>& I,
        double K,
        double G,
        double tauy0,
        double H);

    // Set parameters for a batch of points:
    // each to the same material, but with different parameters:
    // the matrix "idx" refers to a which entry to use: "K(idx)", "G(idx)", or "epsy(idx,:)"

    void setElastic(
        const xt::xtensor<size_t, N>& I,
        const xt::xtensor<size_t, N>& idx,
        const xt::xtensor<double, 1>& K,
        const xt::xtensor<double, 1>& G);

    void setLinearHardening(
        const xt::xtensor<size_t, N>& I,
        const xt::xtensor<size_t, N>& idx,
        const xt::xtensor<double, 1>& K,
        const xt::xtensor<double, 1>& G,
        const xt::xtensor<double, 1>& tauy0,
        const xt::xtensor<double, 1>& H);

    // Set strain tensor, get the response

    void setDefGrad(const xt::xtensor<double, N + 2>& arg, bool compute_tangent = true);
    void defGrad(xt::xtensor<double, N + 2>& ret) const;
    void stress(xt::xtensor<double, N + 2>& ret) const;
    void tangent(xt::xtensor<double, N + 4>& ret) const;
    void epsp(xt::xtensor<double, N>& ret) const;
    void increment(); // update history variables

    // Auto-allocation of the functions above

    xt::xtensor<double, N + 2> DefGrad() const;
    xt::xtensor<double, N + 2> Stress() const;
    xt::xtensor<double, N + 4> Tangent() const;
    xt::xtensor<double, N> Epsp() const;

    // Get copy or reference to the underlying model at on point

    auto getElastic(const std::array<size_t, N>& index) const;
    auto getLinearHardening(const std::array<size_t, N>& index) const;
    auto* refElastic(const std::array<size_t, N>& index);
    auto* refLinearHardening(const std::array<size_t, N>& index);

private:
    // Material vectors
    std::vector<Elastic> m_Elastic;
    std::vector<LinearHardening> m_LinearHardening;

    // Identifiers for each matrix entry
    xt::xtensor<size_t, N> m_type;  // type (e.g. "Type::Elastic")
    xt::xtensor<size_t, N> m_index; // index from the relevant material vector (e.g. "m_Elastic")

    // Shape
    using GMatTensor::Cartesian3d::Array<N>::m_ndim;
    using GMatTensor::Cartesian3d::Array<N>::m_stride_tensor2;
    using GMatTensor::Cartesian3d::Array<N>::m_stride_tensor4;
    using GMatTensor::Cartesian3d::Array<N>::m_size;
    using GMatTensor::Cartesian3d::Array<N>::m_shape;
    using GMatTensor::Cartesian3d::Array<N>::m_shape_tensor2;
    using GMatTensor::Cartesian3d::Array<N>::m_shape_tensor4;
};

} // namespace Cartesian3d
} // namespace GMatElastoPlasticFiniteStrainSimo


#include "Cartesian3d.hpp"
#include "Cartesian3d_Array.hpp"
#include "Cartesian3d_Elastic.hpp"
#include "Cartesian3d_LinearHardening.hpp"

#endif
