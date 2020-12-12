/*

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastoPlasticFiniteStrainSimo

*/

#ifndef GMATELASTOPLASTICFINITESTRAINSIMO_CARTESIAN3D_ELASTIC_HPP
#define GMATELASTOPLASTICFINITESTRAINSIMO_CARTESIAN3D_ELASTIC_HPP

#include "Cartesian3d.h"

namespace GMatElastoPlasticFiniteStrainSimo {
namespace Cartesian3d {

inline Elastic::Elastic(double K, double G) : m_K(K), m_G(G)
{
    m_C = xt::empty<double>({3, 3, 3, 3});
}

inline double Elastic::K() const
{
    return m_K;
}

inline double Elastic::G() const
{
    return m_G;
}

template <class T>
inline void Elastic::setDefGradPtr(const T* arg, bool tangent)
{
    namespace GT = GMatTensor::Cartesian3d::pointer;

    using eigvals_type = std::array<double, 3>;
    using tensor2_type = decltype(m_F);
    using tensor4_type = decltype(m_C);

    std::copy(arg, arg + 9, &m_F[0]);

    // volume change ratio
    double J = GT::Det(&m_F[0]);

    // Finger tensor
    tensor2_type Be;
    GT::A2_dot_A2T(&m_F[0], &Be[0]);

    // eigenvalue decomposition of the trial "Be"
    eigvals_type Be_val;
    tensor2_type vec;
    GT::eigs(&Be[0], &vec[0], &Be_val[0]);

    // logarithmic strain "Eps := 0.5 ln(Be)" (in diagonalised form)
    eigvals_type Eps_val;
    for (size_t j = 0; j < 3; ++j) {
        Eps_val[j] = 0.5 * std::log(Be_val[j]);
    }

    // decompose strain (in diagonalised form)
    eigvals_type Epsd_val;
    double epsm = (Eps_val[0] + Eps_val[1] + Eps_val[2]) / 3.0;
    for (size_t j = 0; j < 3; ++j) {
        Epsd_val[j] = Eps_val[j] - epsm;
    }

    // Cauchy stress (in diagonalised form)
    eigvals_type Sig_val;
    for (size_t j = 0; j < 3; ++j) {
        Sig_val[j] = (3.0 * m_K * epsm + 2.0 * m_G * Epsd_val[j]) / J;
    }

    // compute Cauchy stress, in original coordinate frame
    GT::from_eigs(&vec[0], &Sig_val[0], m_Sig.data());

    if (!tangent) {
        return;
    }

    // unit tensors
    auto II = Cartesian3d::II();
    auto I4d = Cartesian3d::I4d();
    auto I4s = Cartesian3d::I4s();

    // 'linearisation' of the constitutive response
    // Use that "Tau := Ce : Eps = 0.5 * Ce : ln(Be)"
    auto dTau_dlnBe = xt::eval(0.5 * m_K * II + m_G * I4d);
    auto dlnBe_dBe = GMatTensor::Cartesian3d::O4();
    for (size_t m = 0; m < 3; ++m) {
        for (size_t n = 0; n < 3; ++n) {

            double gc = 2.0 * (Eps_val[n] - Eps_val[m]) / (Be_val[n] - Be_val[m]);

            if (Be_val[m] == Be_val[n]) {
                gc = 1.0 / Be_val[m];
            }

            for (size_t i = 0; i < 3; ++i) {
                for (size_t j = 0; j < 3; ++j) {
                    for (size_t k = 0; k < 3; ++k) {
                        for (size_t l = 0; l < 3; ++l) {
                            dlnBe_dBe(i, j, k, l) += gc *
                                vec[i * 3 + m] * vec[j * 3 + n] * vec[k * 3 + m] * vec[l * 3 + n];
                        }
                    }
                }
            }
        }
    }

    // linearization of "Be"
    // Use that "dBe = 2 * (I4s . Be) : LT" (where "LT" refers to "L_\delta^T")
    // Hence: "dBe_dLT = 2 * (I4s * Be)"
    auto dBe_dLT = tensor4_type::from_shape(m_C.shape());
    GT::A4_dot_B2(I4s.data(), &Be[0], dBe_dLT.data());
    dBe_dLT *= 2.0;

    // material tangent stiffness
    // Kmat = dTau_dlnBe : dlnBe_dBe : dBe_dLT
    auto Kmat = tensor4_type::from_shape(m_C.shape());
    GT::A4_ddot_B4_ddot_C4(dTau_dlnBe.data(), dlnBe_dBe.data(), dBe_dLT.data(), Kmat.data());

    // geometrically non-linear tangent
    // Kgeo = -I4rt . Tau
    auto nI4rt = xt::eval(-1.0 * Cartesian3d::I4rt());
    auto Kgeo = tensor4_type::from_shape(m_C.shape());
    GT::A4_dot_B2(nI4rt.data(), &m_Sig[0], Kgeo.data());

    // combine tangents:
    xt::noalias(m_C) = Kgeo + Kmat / J;
}

template <class T>
inline void Elastic::defGradPtr(T* ret) const
{
    std::copy(m_F.cbegin(), m_F.cend(), ret);
}

template <class T>
inline void Elastic::stressPtr(T* ret) const
{
    std::copy(m_Sig.cbegin(), m_Sig.cend(), ret);
}

template <class T>
inline void Elastic::tangentPtr(T* ret) const
{
    std::copy(m_C.cbegin(), m_C.cend(), ret);
}

template <class T>
inline void Elastic::setDefGrad(const T& arg, bool tangent)
{
    GMATELASTOPLASTICFINITESTRAINSIMO_ASSERT(xt::has_shape(arg, {3, 3}));
    return this->setDefGradPtr(arg.data(), tangent);
}

template <class T>
inline void Elastic::defGrad(T& ret) const
{
    GMATELASTOPLASTICFINITESTRAINSIMO_ASSERT(xt::has_shape(ret, {3, 3}));
    return this->defGradPtr(ret.data());
}

template <class T>
inline void Elastic::stress(T& ret) const
{
    GMATELASTOPLASTICFINITESTRAINSIMO_ASSERT(xt::has_shape(ret, {3, 3}));
    return this->stressPtr(ret.data());
}

template <class T>
inline void Elastic::tangent(T& ret) const
{
    GMATELASTOPLASTICFINITESTRAINSIMO_ASSERT(xt::has_shape(ret, {3, 3, 3, 3}));
    return this->tangentPtr(ret.data());
}

inline xt::xtensor<double, 2> Elastic::DefGrad() const
{
    xt::xtensor<double, 2> ret = xt::empty<double>({3, 3});
    this->defGradPtr(ret.data());
    return ret;
}

inline xt::xtensor<double, 2> Elastic::Stress() const
{
    xt::xtensor<double, 2> ret = xt::empty<double>({3, 3});
    this->stressPtr(ret.data());
    return ret;
}

inline xt::xtensor<double, 4> Elastic::Tangent() const
{
    return m_C;
}

} // namespace Cartesian3d
} // namespace GMatElastoPlasticFiniteStrainSimo

#endif
