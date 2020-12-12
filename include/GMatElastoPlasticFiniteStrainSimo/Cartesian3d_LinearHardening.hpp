/*

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastoPlasticFiniteStrainSimo

*/

#ifndef GMATELASTOPLASTICFINITESTRAINSIMO_CARTESIAN3D_LINEARHARDENING_HPP
#define GMATELASTOPLASTICFINITESTRAINSIMO_CARTESIAN3D_LINEARHARDENING_HPP

#include "Cartesian3d.h"

namespace GMatElastoPlasticFiniteStrainSimo {
namespace Cartesian3d {

inline LinearHardening::LinearHardening(double K, double G, double tauy0, double H)
    : m_K(K), m_G(G), m_tauy0(tauy0), m_H(H)
{
    namespace GT = GMatTensor::Cartesian3d::pointer;
    m_epsp = 0.0;
    m_epsp_t = 0.0;
    GT::I2(&m_F[0]);
    GT::I2(&m_F_t[0]);
    GT::I2(&m_Be[0]);
    GT::I2(&m_Be_t[0]);
    m_C = xt::empty<double>({3, 3, 3, 3});
}

inline double LinearHardening::K() const
{
    return m_K;
}

inline double LinearHardening::G() const
{
    return m_G;
}

inline double LinearHardening::tauy0() const
{
    return m_tauy0;
}

inline double LinearHardening::H() const
{
    return m_H;
}

inline double LinearHardening::epsp() const
{
    return m_epsp;
}

inline void LinearHardening::increment()
{
    m_epsp_t = m_epsp;
    std::copy(&m_F[0], &m_F[0] + 9, &m_F_t[0]);
    std::copy(&m_Be[0], &m_Be[0] + 9, &m_Be_t[0]);
}

template <class T>
inline void LinearHardening::setDefGradPtr(const T* arg, bool tangent)
{
    namespace GT = GMatTensor::Cartesian3d::pointer;

    using eigvals_type = std::array<double, 3>;
    using tensor2_type = decltype(m_F);
    using tensor4_type = decltype(m_C);

    std::copy(arg, arg + 9, &m_F[0]);

    // volume change ratio
    double J = GT::Det(&m_F[0]);

    // inverse of "F_t"
    tensor2_type Finv_t;
    GT::Inv(m_F_t.data(), &Finv_t[0]);

    // incremental deformation gradient tensor
    tensor2_type f;
    GT::A2_dot_B2(&m_F[0], &Finv_t[0], &f[0]);

    // trial elastic Finger tensor (symmetric)
    // this assumes "f" to lead to only elastic deformation, which is corrected below if needed
    GT::A2_dot_B2_dot_C2T(&f[0], &m_Be_t[0], &f[0], &m_Be[0]);

    // copy trial elastic Finger tensor (not updated by the return map)
    auto Be = m_Be;

    // eigenvalue decomposition of the trial "Be"
    eigvals_type Be_val;
    tensor2_type vec;
    GT::eigs(&Be[0], &vec[0], &Be_val[0]);

    // logarithmic strain "Eps := 0.5 ln(Be)" (in diagonalised form)
    eigvals_type Epse_val;
    for (size_t j = 0; j < 3; ++j) {
        Epse_val[j] = 0.5 * std::log(Be_val[j]);
    }

    // decompose strain (in diagonalised form)
    eigvals_type Epsed_val;
    double epsem = (Epse_val[0] + Epse_val[1] + Epse_val[2]) / 3.0;
    for (size_t j = 0; j < 3; ++j) {
        Epsed_val[j] = Epse_val[j] - epsem;
    }

    // decomposed trial Kirchhoff stress (in diagonalised form), and trial equivalent stress
    eigvals_type Taud_val;
    double taum = 3.0 * m_K * epsem;
    for (size_t j = 0; j < 3; ++j) {
        Taud_val[j] = 2.0 * m_G * Epsed_val[j];
    }
    double taueq = std::sqrt(1.5 *
        (std::pow(Taud_val[0], 2.0) + std::pow(Taud_val[1], 2.0) + std::pow(Taud_val[2], 2.0)));

    // evaluate the yield surface
    double phi = taueq - (m_tauy0 + m_H * m_epsp_t);

    // (direction of) plastic flow
    double dgamma = 0.0;
    tensor2_type N;

    // return map
    if (phi > 0) {
        // - plastic flow
        dgamma = phi / (3.0 * m_G + m_H);
        // - direction of plastic flow
        eigvals_type N_val;
        for (size_t j = 0; j < 3; ++j) {
            N_val[j] = 1.5 * Taud_val[j] / taueq;
        }
        GT::from_eigs(&vec[0], &N_val[0], &N[0]);
        // - update trial stress and elastic strain (only the deviatoric part)
        eigvals_type e;
        for (size_t j = 0; j < 3; ++j) {
            Taud_val[j] *= (1.0 - 3.0 * m_G * dgamma / taueq);
            Epsed_val[j] = Taud_val[j] / (2.0 * m_G);
            e[j] = std::exp(2.0 * (epsem +  Epsed_val[j]));
        }
        // - update elastic Finger tensor, in original coordinate frame
        GT::from_eigs(&vec[0], &e[0], &m_Be[0]);
        // - update equivalent plastic strain
        m_epsp = m_epsp_t + dgamma;
    }

    // compute Cauchy stress, in original coordinate frame
    eigvals_type Sig_val;
    for (size_t j = 0; j < 3; ++j) {
        Sig_val[j] = (taum + Taud_val[j]) / J;
    }
    GT::from_eigs(&vec[0], &Sig_val[0], &m_Sig[0]);

    if (!tangent) {
        return;
    }

    // unit tensors
    auto II = Cartesian3d::II();
    auto I4d = Cartesian3d::I4d();
    auto I4s = Cartesian3d::I4s();

    // linearisation of the constitutive response
    auto dTau_dlnBe = tensor4_type::from_shape(m_C.shape());

    if (phi <= 0) {
        // - Use that "Tau := Ce : Eps = 0.5 * Ce : ln(Be)"
        xt::noalias(dTau_dlnBe) = 0.5 * m_K * II + m_G * I4d;
    }
    else {
        // - Directions of plastic flow
        auto NN = tensor4_type::from_shape(m_C.shape());
        GT::A2_dyadic_B2(&N[0], &N[0], NN.data());
        // - Temporary constants
        double a0;
        double a1 = m_G / (m_H + 3.0 * m_G);
        if (dgamma != 0.0) {
            a0 = dgamma * m_G / taueq;
        }
        else {
            a0 = 0.0;
        }
        // - Elasto-plastic tangent
        xt::noalias(dTau_dlnBe) = (0.5 * (m_K - 2.0 / 3.0 * m_G) + a0 * m_G) * II +
            (1.0 - 3.0 * a0) * m_G * I4s + 2.0 * m_G * (a0 - a1) * NN;
    }

    auto dlnBe_dBe = GMatTensor::Cartesian3d::O4();

    for (size_t m = 0; m < 3; ++m) {
        for (size_t n = 0; n < 3; ++n) {

            double gc = (std::log(Be_val[n]) - std::log(Be_val[m])) / (Be_val[n] - Be_val[m]);

            if (Be_val[m] == Be_val[n]) {
                gc = 1.0 / Be_val[m];
            }

            for (size_t l = 0; l < 3; ++l) {
                for (size_t k = 0; k < 3; ++k) {
                    for (size_t j = 0; j < 3; ++j) {
                        for (size_t i = 0; i < 3; ++i) {
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
inline void LinearHardening::defGradPtr(T* ret) const
{
    std::copy(m_F.cbegin(), m_F.cend(), ret);
}

template <class T>
inline void LinearHardening::stressPtr(T* ret) const
{
    std::copy(m_Sig.cbegin(), m_Sig.cend(), ret);
}

template <class T>
inline void LinearHardening::tangentPtr(T* ret) const
{
    std::copy(m_C.cbegin(), m_C.cend(), ret);
}

template <class T>
inline void LinearHardening::setDefGrad(const T& arg, bool tangent)
{
    GMATELASTOPLASTICFINITESTRAINSIMO_ASSERT(xt::has_shape(arg, {3, 3}));
    return this->setDefGradPtr(arg.data(), tangent);
}

template <class T>
inline void LinearHardening::defGrad(T& ret) const
{
    GMATELASTOPLASTICFINITESTRAINSIMO_ASSERT(xt::has_shape(ret, {3, 3}));
    return this->defGradPtr(ret.data());
}

template <class T>
inline void LinearHardening::stress(T& ret) const
{
    GMATELASTOPLASTICFINITESTRAINSIMO_ASSERT(xt::has_shape(ret, {3, 3}));
    return this->stressPtr(ret.data());
}

template <class T>
inline void LinearHardening::tangent(T& ret) const
{
    GMATELASTOPLASTICFINITESTRAINSIMO_ASSERT(xt::has_shape(ret, {3, 3, 3, 3}));
    return this->tangentPtr(ret.data());
}

inline xt::xtensor<double, 2> LinearHardening::DefGrad() const
{
    xt::xtensor<double, 2> ret = xt::empty<double>({3, 3});
    this->defGradPtr(ret.data());
    return ret;
}

inline xt::xtensor<double, 2> LinearHardening::Stress() const
{
    xt::xtensor<double, 2> ret = xt::empty<double>({3, 3});
    this->stressPtr(ret.data());
    return ret;
}

inline xt::xtensor<double, 4> LinearHardening::Tangent() const
{
    return m_C;
}

} // namespace Cartesian3d
} // namespace GMatElastoPlasticFiniteStrainSimo

#endif
