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

template <class T>
inline void LinearHardening::stress(const Tensor2& F, T&& Sig)
{
    // volume change ratio
    double J = det(F);

    // inverse of "F_t"
    Tensor2 Finv_t;
    inv(m_F_t, Finv_t);

    // incremental deformation gradient tensor
    Tensor2 f;
    A2_dot_B2(F, Finv_t, f);

    // trial elastic Finger tensor (symmetric)
    // this assumes "f" to lead to only elastic deformation, which is corrected below if needed
    A2_dot_B2_dot_C2T(f, m_Be_t, f, m_Be);

    // eigenvalue decomposition of the trial "Be"
    Tensor2 vec;
    Vector Be_val;
    eig(m_Be, vec, Be_val);

    // trial logarithmic strain "Eps := 0.5 * ln(Be)" (in diagonalised form)
    Vector Epse_val = 0.5 * xt::log(Be_val);

    // decompose trail strain (in diagonalised form)
    double epsem = xt::sum(Epse_val)[0] / 3.0;
    Vector Epsed_val = Epse_val - epsem;

    // decomposed trial Kirchhoff stress (in diagonalised form), and trial equivalent stress
    double taum = 3.0 * m_K * epsem;
    Vector Taud_val = 2.0 * m_G * Epsed_val;
    double taueq = std::sqrt(1.5 * xt::sum(xt::pow(Taud_val, 2.0))[0]);

    // evaluate the yield surface
    double phi = taueq - (m_tauy0 + m_H * m_epsp_t);

    // return map
    if (phi > 0) {
        // - plastic flow
        double dgamma = phi / (3.0 * m_G + m_H);
        // - update trial stress (only the deviatoric part)
        Taud_val *= (1.0 - 3.0 * m_G * dgamma / taueq);
        // - update elastic strain (only the deviatoric part)
        xt::noalias(Epsed_val) = Taud_val / (2.0 * m_G);
        // - update elastic Finger tensor, in original coordinate frame
        inv_eig(vec, xt::exp(2.0 * (epsem +  Epsed_val)), m_Be);
        // - update equivalent plastic strain
        m_epsp = m_epsp_t + dgamma;
    }

    // compute Cauchy stress, in original coordinate frame
    inv_eig(vec, (taum + Taud_val) / J, Sig);

    // store history
    xt::noalias(m_F) = F;
}

inline Tensor2 LinearHardening::Stress(const Tensor2& F)
{
    Tensor2 Sig;
    this->stress(F, Sig);
    return Sig;
}

template <class T, class S>
inline void LinearHardening::tangent(const Tensor2& F, T&& Sig, S&& C)
{
    // volume change ratio
    double J = det(F);

    // inverse of "F_t"
    Tensor2 Finv_t;
    inv(m_F_t, Finv_t);

    // incremental deformation gradient tensor
    Tensor2 f;
    A2_dot_B2(F, Finv_t, f);

    // trial elastic Finger tensor (symmetric)
    // this assumes "f" to lead to only elastic deformation, which is corrected below if needed
    A2_dot_B2_dot_C2T(f, m_Be_t, f, m_Be);
    // copy trial elastic Finger tensor (not updated by the return map)
    Tensor2 Be = m_Be;

    // eigenvalue decomposition of the trial "Be"
    Tensor2 vec;
    Vector Be_val;
    eig(m_Be, vec, Be_val);

    // trial logarithmic strain "Eps := 0.5 * ln(Be)" (in diagonalised form)
    Vector Epse_val = 0.5 * xt::log(Be_val);

    // decompose trail strain (in diagonalised form)
    double epsem = xt::sum(Epse_val)[0] / 3.0;
    Vector Epsed_val = Epse_val - epsem;

    // decomposed trial Kirchhoff stress (in diagonalised form), and trial equivalent stress
    double taum = 3.0 * m_K * epsem;
    Vector Taud_val = 2.0 * m_G * Epsed_val;
    double taueq = std::sqrt(1.5 * xt::sum(xt::pow(Taud_val, 2.0))[0]);

    // evaluate the yield surface
    double phi = taueq - (m_tauy0 + m_H * m_epsp_t);

    // (direction of) plastic flow
    double dgamma = 0.0;
    Tensor2 N;

    // return map
    if (phi > 0) {
        // - plastic flow
        dgamma = phi / (3.0 * m_G + m_H);
        // - direction of plastic flow
        inv_eig(vec, 1.5 * Taud_val / taueq, N);
        // - update trial stress (only the deviatoric part)
        Taud_val *= (1.0 - 3.0 * m_G * dgamma / taueq);
        // - update elastic strain (only the deviatoric part)
        xt::noalias(Epsed_val) = Taud_val / (2.0 * m_G);
        // - update elastic Finger tensor, in original coordinate frame
        inv_eig(vec, xt::exp(2.0 * (epsem +  Epsed_val)), m_Be);
        // - update equivalent plastic strain
        m_epsp = m_epsp_t + dgamma;
    }

    // compute Kirchhoff/Cauchy stress, in original coordinate frame
    Tensor2 Tau;
    inv_eig(vec, taum + Taud_val, Tau);
    xt::noalias(Sig) = Tau / J;

    // store history
    xt::noalias(m_F) = F;

    // unit tensors
    auto II = Cartesian3d::II();
    auto I4d = Cartesian3d::I4d();
    auto I4s = Cartesian3d::I4s();
    auto I4rt = Cartesian3d::I4rt();

    // local variables
    Tensor4 dTau_dlnBe;
    Tensor4 dlnBe_dBe;
    Tensor4 dBe_dLT;
    Tensor4 Kmat;
    Tensor4 Kgeo;
    double gc;

    // linearisation of the constitutive response
    if (phi <= 0) {
        // - Use that "Tau := Ce : Eps = 0.5 * Ce : ln(Be)"
        dTau_dlnBe = 0.5 * m_K * II + m_G * I4d;
    }
    else {
        // - Directions of plastic flow
        Tensor4 NN;
        A2_dyadic_B2(N, N, NN);
        // - Temporary constants
        double a0;
        double a1 = m_G / (m_H + 3.0 * m_G);
        if (dgamma != 0.0)
          a0 = dgamma * m_G / taueq;
        else
          a0 = 0.0;
        // - Elasto-plastic tangent
        dTau_dlnBe = (0.5 * (m_K - 2.0 / 3.0 * m_G) + a0 * m_G) * II +
                     (1.0 - 3.0 * a0) * m_G * I4s + 2.0 * m_G * (a0 - a1) * NN;
    }

    // linearization of the logarithmic strain
    dlnBe_dBe.fill(0.0);

    for (size_t m = 0; m < 3; ++m) {
        for (size_t n = 0; n < 3; ++n) {

            if (Be_val(m) == Be_val(n)) {
                gc = 1.0 / Be_val(m);
            }
            else {
                gc = (std::log(Be_val(n)) - std::log(Be_val(m))) / (Be_val(n) - Be_val(m));
            }

            for (size_t l = 0; l < 3; ++l) {
                for (size_t k = 0; k < 3; ++k) {
                    for (size_t j = 0; j < 3; ++j) {
                        for (size_t i = 0; i < 3; ++i) {
                            dlnBe_dBe(i, j, k, l) +=
                                gc * vec(i, m) * vec(j, n) * vec(k, m) * vec(l, n);
                        }
                    }
                }
            }
        }
    }

    // linearization of "Be"
    // Use that "dBe = 2 * (I4s . Be) : LT" (where "LT" refers to "L_\delta^T")
    // Hence: "dBe_dLT = 2 * (I4s * Be)"
    A4_dot_B2(I4s, Be, dBe_dLT);
    dBe_dLT *= 2.0;

    // material tangent stiffness
    // Kmat = dTau_dlnBe : dlnBe_dBe : dBe_dLT
    A4_ddot_B4_ddot_C4(dTau_dlnBe, dlnBe_dBe, dBe_dLT, Kmat);

    // geometrically non-linear tangent
    // Kgeo = -I4rt . Tau
    A4_dot_B2(-I4rt, Tau, Kgeo);

    // combine tangents:
    xt::noalias(C) = (Kgeo + Kmat) / J;
}

inline std::tuple<Tensor2, Tensor4> LinearHardening::Tangent(const Tensor2& F)
{
    Tensor2 Sig;
    Tensor4 C;
    this->tangent(F, Sig, C);
    return std::make_tuple(Sig, C);
}

inline void LinearHardening::increment()
{
    m_epsp_t = m_epsp;
    xt::noalias(m_F_t) = m_F;
    xt::noalias(m_Be_t) = m_Be;
}

inline double LinearHardening::epsp() const
{
    return m_epsp;
}

} // namespace Cartesian3d
} // namespace GMatElastoPlasticFiniteStrainSimo

#endif
