/*

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastoPlasticFiniteStrainSimo

*/

#ifndef GMATELASTOPLASTICFINITESTRAINSIMO_CARTESIAN3D_ELASTIC_HPP
#define GMATELASTOPLASTICFINITESTRAINSIMO_CARTESIAN3D_ELASTIC_HPP

#include "Cartesian3d.h"

namespace GMatElastoPlasticFiniteStrainSimo {
namespace Cartesian3d {

// -------------------------------------------------------------------------------------------------

inline Elastic::Elastic(double K, double G) : m_K(K), m_G(G)
{
}

// -------------------------------------------------------------------------------------------------

inline double Elastic::K() const
{
  return m_K;
}

// -------------------------------------------------------------------------------------------------

inline double Elastic::G() const
{
  return m_G;
}

// -------------------------------------------------------------------------------------------------

template <class T>
inline void Elastic::stress(const Tensor2& F, T&& Sig) const
{
  // volume change ratio
  double J = det(F);

  // Finger tensor
  Tensor2 Be;
  finger(F, Be);

  // eigenvalue decomposition of "Be"
  Tensor2 vec;
  Vector  Be_val;
  eig(Be, vec, Be_val);

  // logarithmic strain "Eps := 0.5 ln(Be)" (in diagonalised form)
  Vector Eps_val = 0.5 * xt::log(Be_val);

  // decompose strain (in diagonalised form)
  double epsm     = xt::sum(Eps_val)[0] / 3.0;
  Vector Epsd_val = Eps_val - epsm;

  // Kirchhoff stress (in diagonalised form)
  Vector Tau_val = 3.0 * m_K * epsm + 2.0 * m_G * Epsd_val;

  // compute Cauchy stress, in original coordinate frame
  inv_eig(vec, Tau_val / J, Sig);
}

// -------------------------------------------------------------------------------------------------

inline Tensor2 Elastic::Stress(const Tensor2& F) const
{
  Tensor2 Sig;
  this->stress(F, Sig);
  return Sig;
}

// -------------------------------------------------------------------------------------------------

template <class T, class S>
inline void Elastic::tangent(const Tensor2& F, T&& Sig, S&& C) const
{
  // volume change ratio
  double J = det(F);

  // Finger tensor
  Tensor2 Be;
  finger(F, Be);

  // eigenvalue decomposition of "Be"
  Tensor2 vec;
  Vector  Be_val;
  eig(Be, vec, Be_val);

  // logarithmic strain "Eps := 0.5 ln(Be)" (in diagonalised form)
  Vector Eps_val = 0.5 * xt::log(Be_val);

  // decompose strain (in diagonalised form)
  double epsm     = xt::sum(Eps_val)[0] / 3.0;
  Vector Epsd_val = Eps_val - epsm;

  // Kirchhoff stress (in diagonalised form)
  Vector Tau_val = 3.0 * m_K * epsm + 2.0 * m_G * Epsd_val;

  // compute Kirchhoff/Cauchy stress, in original coordinate frame
  Tensor2 Tau;
  inv_eig(vec, Tau_val, Tau);
  xt::noalias(Sig) = Tau / J;

  // unit tensors
  auto II   = Cartesian3d::II();
  auto I4d  = Cartesian3d::I4d();
  auto I4s  = Cartesian3d::I4s();
  auto I4rt = Cartesian3d::I4rt();

  // local variables
  Tensor4 dTau_dlnBe, dlnBe_dBe, dBe_dLT, Kmat, Kgeo;

  // 'linearisation' of the constitutive response
  // Use that "Tau := Ce : Eps = 0.5 * Ce : ln(Be)"
  dTau_dlnBe = 0.5 * m_K * II + m_G * I4d;

  // linearization of the logarithmic strain
  dlnBe_dBe.fill(0.0);

  for (size_t m = 0; m < 3; ++m) {
    for (size_t n = 0; n < 3; ++n) {

      double gc = 2.0 * (Eps_val(n) - Eps_val(m)) / (Be_val(n) - Be_val(m));

      if (Be_val(m) == Be_val(n))
        gc = 1.0 / Be_val(m);

      for (size_t i = 0; i < 3; ++i)
        for (size_t j = 0; j < 3; ++j)
          for (size_t k = 0; k < 3; ++k)
            for (size_t l = 0; l < 3; ++l)
              dlnBe_dBe(i,j,k,l) += gc * vec(i,m) * vec(j,n) * vec(k,m) * vec(l,n);
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

// -------------------------------------------------------------------------------------------------

inline std::tuple<Tensor2,Tensor4> Elastic::Tangent(const Tensor2& F) const
{
  Tensor2 Sig;
  Tensor4 C;
  this->tangent(F, Sig, C);
  return std::make_tuple(Sig, C);
}

// -------------------------------------------------------------------------------------------------

inline void Elastic::increment()
{
}

// -------------------------------------------------------------------------------------------------

inline double Elastic::epsp() const
{
  return 0.0;
}

// -------------------------------------------------------------------------------------------------

}} // namespace ...

#endif
