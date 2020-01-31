/*

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastoPlasticFiniteStrainSimo

*/

#ifndef GMATELASTOPLASTICFINITESTRAINSIMO_CARTESIAN3D_HPP
#define GMATELASTOPLASTICFINITESTRAINSIMO_CARTESIAN3D_HPP

#include "Cartesian3d.h"

namespace GMatElastoPlasticFiniteStrainSimo {
namespace Cartesian3d {

// -------------------------------------------------------------------------------------------------

inline Tensor2 I2()
{
  return Tensor2({{1.0, 0.0, 0.0},
                  {0.0, 1.0, 0.0},
                  {0.0, 0.0, 1.0}});
}

// -------------------------------------------------------------------------------------------------

inline Tensor4 II()
{
  Tensor4 out;
  out.fill(0.0);

  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      for (size_t k = 0; k < 3; ++k)
        for (size_t l = 0; l < 3; ++l)
          if (i == j and k == l)
            out(i,j,k,l) = 1.0;

  return out;
}

// -------------------------------------------------------------------------------------------------

inline Tensor4 I4()
{
  Tensor4 out;
  out.fill(0.0);

  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      for (size_t k = 0; k < 3; ++k)
        for (size_t l = 0; l < 3; ++l)
          if (i == l and j == k)
            out(i,j,k,l) = 1.0;

  return out;
}

// -------------------------------------------------------------------------------------------------

inline Tensor4 I4rt()
{
  Tensor4 out;
  out.fill(0.0);

  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      for (size_t k = 0; k < 3; ++k)
        for (size_t l = 0; l < 3; ++l)
          if (i == k and j == l)
            out(i,j,k,l) = 1.0;

  return out;
}

// -------------------------------------------------------------------------------------------------

inline Tensor4 I4s()
{
  return 0.5 * ( I4() + I4rt() );
}

// -------------------------------------------------------------------------------------------------

inline Tensor4 I4d()
{
  return I4s() - II() / 3.0;
}

// -------------------------------------------------------------------------------------------------

inline double Hydrostatic(const Tensor2& A)
{
  return trace(A) / 3.0;
}

// -------------------------------------------------------------------------------------------------

inline Tensor2 Deviatoric(const Tensor2& A)
{
  return A - trace(A) / 3.0 * I2();
}

// -------------------------------------------------------------------------------------------------

inline Tensor2 Strain(const Tensor2 &F)
{
  // Finger tensor: B = F . F^T
  Tensor2 B;
  finger(F, B);

  // eigenvalue decomposition of "B"
  Tensor2 vec;
  Vector  B_val;
  eig(B, vec, B_val);

  // logarithmic strain "Eps == 0.5 * ln(Be)" (in diagonalised form)
  Vector Eps_val = 0.5 * xt::log(B_val);

  // reconstruct strain in original coordinate frame
  Tensor2 Eps;
  inv_eig(vec, Eps_val, Eps);

  return Eps;
}

// -------------------------------------------------------------------------------------------------

inline double Epseq(const Tensor2& Eps)
{
  Tensor2 Epsd = Eps - trace(Eps) / 3.0 * I2();
  return std::sqrt(2.0/3.0 * A2_ddot_B2(Epsd,Epsd));
}

// -------------------------------------------------------------------------------------------------

inline double Sigeq(const Tensor2& Sig)
{
  Tensor2 Sigd = Sig - trace(Sig) / 3.0 * I2();
  return std::sqrt(1.5 * A2_ddot_B2(Sigd,Sigd));
}

// -------------------------------------------------------------------------------------------------

inline void hydrostatic(const xt::xtensor<double,3>& A, xt::xtensor<double,1>& Am)
{
  GMATELASTOPLASTICFINITESTRAINSIMO_ASSERT(A.shape() ==\
    std::decay_t<decltype(A)>::shape_type({Am.shape(0), 3, 3}));

  #pragma omp parallel
  {
    #pragma omp for
    for (size_t e = 0; e < A.shape(0); ++e) {
      auto Ai = xt::adapt(&A(e,0,0), xt::xshape<3,3>());
      Am(e) = trace(Ai) / 3.0;
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline void deviatoric(const xt::xtensor<double,3>& A, xt::xtensor<double,3>& Ad)
{
  GMATELASTOPLASTICFINITESTRAINSIMO_ASSERT(A.shape() ==\
    std::decay_t<decltype(A)>::shape_type({Ad.shape(0), 3, 3}));

  #pragma omp parallel
  {
    Tensor2 unit = I2();
    #pragma omp for
    for (size_t e = 0; e < A.shape(0); ++e) {
      auto Ai  = xt::adapt(&A (e,0,0), xt::xshape<3,3>());
      auto Aid = xt::adapt(&Ad(e,0,0), xt::xshape<3,3>());
      xt::noalias(Aid) = Ai - trace(Ai) / 3.0 * unit;
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline void strain(const xt::xtensor<double,3>& F, xt::xtensor<double,3>& Eps)
{
  GMATELASTOPLASTICFINITESTRAINSIMO_ASSERT(F.shape() ==\
    std::decay_t<decltype(F)>::shape_type({Eps.shape(0), 3, 3}));

  #pragma omp parallel
  {
    Tensor2 B, vec;
    Vector  B_val, Eps_val;
    #pragma omp for
    for (size_t e = 0; e < F.shape(0); ++e) {
      auto Fi   = xt::adapt(&F  (e,0,0), xt::xshape<3,3>());
      auto Epsi = xt::adapt(&Eps(e,0,0), xt::xshape<3,3>());
      finger(Fi, B);
      eig(B, vec, B_val);
      Eps_val = 0.5 * xt::log(B_val);
      inv_eig(vec, Eps_val, Epsi);
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline void epseq(const xt::xtensor<double,3>& A, xt::xtensor<double,1>& Aeq)
{
  GMATELASTOPLASTICFINITESTRAINSIMO_ASSERT(A.shape() ==\
    std::decay_t<decltype(A)>::shape_type({Aeq.shape(0), 3, 3}));

  #pragma omp parallel
  {
    Tensor2 unit = I2();
    #pragma omp for
    for (size_t e = 0; e < A.shape(0); ++e) {
      auto Ai  = xt::adapt(&A(e,0,0), xt::xshape<3,3>());
      auto Aid = Ai - trace(Ai) / 3.0 * unit;
      Aeq(e) = std::sqrt(2.0/3.0 * A2_ddot_B2(Aid,Aid));
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline void sigeq(const xt::xtensor<double,3>& A, xt::xtensor<double,1>& Aeq)
{
  GMATELASTOPLASTICFINITESTRAINSIMO_ASSERT(A.shape() ==\
    std::decay_t<decltype(A)>::shape_type({Aeq.shape(0), 3, 3}));

  #pragma omp parallel
  {
    Tensor2 unit = I2();
    #pragma omp for
    for (size_t e = 0; e < A.shape(0); ++e) {
      auto Ai  = xt::adapt(&A(e,0,0), xt::xshape<3,3>());
      auto Aid = Ai - trace(Ai) / 3.0 * unit;
      Aeq(e) = std::sqrt(1.5 * A2_ddot_B2(Aid,Aid));
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,1> Hydrostatic(const xt::xtensor<double,3>& A)
{
  xt::xtensor<double,1> Am = xt::empty<double>({A.shape(0)});
  Cartesian3d::hydrostatic(A, Am);
  return Am;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,3> Deviatoric(const xt::xtensor<double,3>& A)
{
  xt::xtensor<double,3> Ad = xt::empty<double>(A.shape());
  Cartesian3d::deviatoric(A, Ad);
  return Ad;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,3> Strain(const xt::xtensor<double,3>& F)
{
  xt::xtensor<double,3> Eps = xt::empty<double>(F.shape());
  Cartesian3d::strain(F, Eps);
  return Eps;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,1> Epseq(const xt::xtensor<double,3>& A)
{
  xt::xtensor<double,1> Aeq = xt::empty<double>({A.shape(0)});
  Cartesian3d::epseq(A, Aeq);
  return Aeq;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,1> Sigeq(const xt::xtensor<double,3>& A)
{
  xt::xtensor<double,1> Aeq = xt::empty<double>({A.shape(0)});
  Cartesian3d::sigeq(A, Aeq);
  return Aeq;
}

// -------------------------------------------------------------------------------------------------

inline void hydrostatic(const xt::xtensor<double,4>& A, xt::xtensor<double,2>& Am)
{
  GMATELASTOPLASTICFINITESTRAINSIMO_ASSERT(A.shape() ==\
    std::decay_t<decltype(A)>::shape_type({Am.shape(0), Am.shape(1), 3, 3}));

  #pragma omp parallel
  {
    #pragma omp for
    for (size_t e = 0; e < A.shape(0); ++e) {
      for (size_t q = 0; q < A.shape(1); ++q) {
        auto Ai = xt::adapt(&A(e,q,0,0), xt::xshape<3,3>());
        Am(e,q) = trace(Ai) / 3.0;
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline void deviatoric(const xt::xtensor<double,4>& A, xt::xtensor<double,4>& Ad)
{
  GMATELASTOPLASTICFINITESTRAINSIMO_ASSERT(A.shape() ==\
    std::decay_t<decltype(A)>::shape_type({Ad.shape(0), Ad.shape(1), 3, 3}));

  #pragma omp parallel
  {
    Tensor2 unit = I2();
    #pragma omp for
    for (size_t e = 0; e < A.shape(0); ++e) {
      for (size_t q = 0; q < A.shape(1); ++q) {
        auto Ai  = xt::adapt(&A (e,q,0,0), xt::xshape<3,3>());
        auto Aid = xt::adapt(&Ad(e,q,0,0), xt::xshape<3,3>());
        xt::noalias(Aid) = Ai - trace(Ai) / 3.0 * unit;
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline void strain(const xt::xtensor<double,4>& F, xt::xtensor<double,4>& Eps)
{
  GMATELASTOPLASTICFINITESTRAINSIMO_ASSERT(F.shape() ==\
    std::decay_t<decltype(F)>::shape_type({Eps.shape(0), Eps.shape(1), 3, 3}));

  #pragma omp parallel
  {
    Tensor2 B, vec;
    Vector  B_val, Eps_val;
    #pragma omp for
    for (size_t e = 0; e < F.shape(0); ++e) {
      for (size_t q = 0; q < F.shape(1); ++q) {
        auto Fi   = xt::adapt(&F  (e,q,0,0), xt::xshape<3,3>());
        auto Epsi = xt::adapt(&Eps(e,q,0,0), xt::xshape<3,3>());
        finger(Fi, B);
        eig(B, vec, B_val);
        Eps_val = 0.5 * xt::log(B_val);
        inv_eig(vec, Eps_val, Epsi);
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline void epseq(const xt::xtensor<double,4>& A, xt::xtensor<double,2>& Aeq)
{
  GMATELASTOPLASTICFINITESTRAINSIMO_ASSERT(A.shape() ==\
    std::decay_t<decltype(A)>::shape_type({Aeq.shape(0), Aeq.shape(1), 3, 3}));

  #pragma omp parallel
  {
    Tensor2 unit = I2();
    #pragma omp for
    for (size_t e = 0; e < A.shape(0); ++e) {
      for (size_t q = 0; q < A.shape(1); ++q) {
        auto Ai  = xt::adapt(&A(e,q,0,0), xt::xshape<3,3>());
        auto Aid = Ai - trace(Ai) / 3.0 * unit;
        Aeq(e,q) = std::sqrt(2.0/3.0 * A2_ddot_B2(Aid,Aid));
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline void sigeq(const xt::xtensor<double,4>& A, xt::xtensor<double,2>& Aeq)
{
  GMATELASTOPLASTICFINITESTRAINSIMO_ASSERT(A.shape() ==\
    std::decay_t<decltype(A)>::shape_type({Aeq.shape(0), Aeq.shape(1), 3, 3}));

  #pragma omp parallel
  {
    Tensor2 unit = I2();
    #pragma omp for
    for (size_t e = 0; e < A.shape(0); ++e) {
      for (size_t q = 0; q < A.shape(1); ++q) {
        auto Ai  = xt::adapt(&A(e,q,0,0), xt::xshape<3,3>());
        auto Aid = Ai - trace(Ai) / 3.0 * unit;
        Aeq(e,q) = std::sqrt(1.5 * A2_ddot_B2(Aid,Aid));
      }
    }
  }
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,2> Hydrostatic(const xt::xtensor<double,4>& A)
{
  xt::xtensor<double,2> Am = xt::empty<double>({A.shape(0), A.shape(1)});
  Cartesian3d::hydrostatic(A, Am);
  return Am;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,4> Deviatoric(const xt::xtensor<double,4>& A)
{
  xt::xtensor<double,4> Ad = xt::empty<double>(A.shape());
  Cartesian3d::deviatoric(A, Ad);
  return Ad;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,4> Strain(const xt::xtensor<double,4>& F)
{
  xt::xtensor<double,4> Eps = xt::empty<double>(F.shape());
  Cartesian3d::strain(F, Eps);
  return Eps;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,2> Epseq(const xt::xtensor<double,4>& A)
{
  xt::xtensor<double,2> Aeq = xt::empty<double>({A.shape(0), A.shape(1)});
  Cartesian3d::epseq(A, Aeq);
  return Aeq;
}

// -------------------------------------------------------------------------------------------------

inline xt::xtensor<double,2> Sigeq(const xt::xtensor<double,4>& A)
{
  xt::xtensor<double,2> Aeq = xt::empty<double>({A.shape(0), A.shape(1)});
  Cartesian3d::sigeq(A, Aeq);
  return Aeq;
}

// -------------------------------------------------------------------------------------------------

template <class U>
inline double trace(const U& A)
{
  return A(0,0) + A(1,1) + A(2,2);
}

// -------------------------------------------------------------------------------------------------

template <class U, class V, class W>
inline void A2_dyadic_B2(const U& A, const V& B, W& C)
{
  C.fill(0.0);

  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      for (size_t k = 0; k < 3; ++k)
        for (size_t l = 0; l < 3; ++l)
          C(i,j,k,l) += A(i,j) * B(k,l);
}

// -------------------------------------------------------------------------------------------------

template <class U, class V>
inline double A2_ddot_B2(const U& A, const V& B)
{
  return A(0,0) * B(0,0) + 2.0 * A(0,1) * B(0,1) + 2.0 * A(0,2) * B(0,2) +
         A(1,1) * B(1,1) + 2.0 * A(1,2) * B(1,2) +
         A(2,2) * B(2,2);
}

// -------------------------------------------------------------------------------------------------

template <class U, class V, class W>
inline void A2_dot_B2(const U& A, const V& B, W& C)
{
  C(0,0) = A(0,0) * B(0,0) + A(0,1) * B(1,0) + A(0,2) * B(2,0);
  C(0,1) = A(0,0) * B(0,1) + A(0,1) * B(1,1) + A(0,2) * B(2,1);
  C(0,2) = A(0,0) * B(0,2) + A(0,1) * B(1,2) + A(0,2) * B(2,2);
  C(1,0) = A(1,0) * B(0,0) + A(1,1) * B(1,0) + A(1,2) * B(2,0);
  C(1,1) = A(1,0) * B(0,1) + A(1,1) * B(1,1) + A(1,2) * B(2,1);
  C(1,2) = A(1,0) * B(0,2) + A(1,1) * B(1,2) + A(1,2) * B(2,2);
  C(2,0) = A(2,0) * B(0,0) + A(2,1) * B(1,0) + A(2,2) * B(2,0);
  C(2,1) = A(2,0) * B(0,1) + A(2,1) * B(1,1) + A(2,2) * B(2,1);
  C(2,2) = A(2,0) * B(0,2) + A(2,1) * B(1,2) + A(2,2) * B(2,2);
}

// -------------------------------------------------------------------------------------------------

template <class U, class V, class W>
inline void A4_dot_B2(const U& A, const V& B, W& C)
{
  C.fill(0.0);

  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      for (size_t k = 0; k < 3; ++k)
        for (size_t l = 0; l < 3; ++l)
          for (size_t m = 0; m < 3; ++m)
            C(i,j,k,m) += A(i,j,k,l) * B(l,m);
}

// -------------------------------------------------------------------------------------------------

template <class U, class V, class W, class X>
inline void A4_ddot_B4_ddot_C4(const U& A, const V& B, const W& C, X& D)
{
  D.fill(0.0);

  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      for (size_t k = 0; k < 3; ++k)
        for (size_t l = 0; l < 3; ++l)
          for (size_t m = 0; m < 3; ++m)
            for (size_t n = 0; n < 3; ++n)
              for (size_t o = 0; o < 3; ++o)
                for (size_t p = 0; p < 3; ++p)
                  D(i,j,o,p) += A(i,j,k,l) * B(l,k,m,n) * C(n,m,o,p);
}

// -------------------------------------------------------------------------------------------------

template <class U, class V, class W, class X>
inline void A2_dot_B2_dot_C2T(const U& A, const V& B, const W& C, X& D)
{
  D.fill(0.0);

  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      for (size_t h = 0; h < 3; ++h)
        for (size_t l = 0; l < 3; ++l)
          D(i,l) += A(i,j) * B(j,h) * C(l,h);
}

// -------------------------------------------------------------------------------------------------

template <class U>
inline double det(const U& A)
{
  return ( A(0,0) * A(1,1) * A(2,2) +
           A(0,1) * A(1,2) * A(2,0) +
           A(0,2) * A(1,0) * A(2,1) ) -
         ( A(0,2) * A(1,1) * A(2,0) +
           A(0,1) * A(1,0) * A(2,2) +
           A(0,0) * A(1,2) * A(2,1) );
}

// -------------------------------------------------------------------------------------------------

template <class U, class V>
inline void inv(const U& A, V& Ainv)
{
  auto D = det(A);

  Ainv(0,0) = (A(1,1) * A(2,2) - A(1,2) * A(2,1)) / D;
  Ainv(0,1) = (A(0,2) * A(2,1) - A(0,1) * A(2,2)) / D;
  Ainv(0,2) = (A(0,1) * A(1,2) - A(0,2) * A(1,1)) / D;
  Ainv(1,0) = (A(1,2) * A(2,0) - A(1,0) * A(2,2)) / D;
  Ainv(1,1) = (A(0,0) * A(2,2) - A(0,2) * A(2,0)) / D;
  Ainv(1,2) = (A(0,2) * A(1,0) - A(0,0) * A(1,2)) / D;
  Ainv(2,0) = (A(1,0) * A(2,1) - A(1,1) * A(2,0)) / D;
  Ainv(2,1) = (A(0,1) * A(2,0) - A(0,0) * A(2,1)) / D;
  Ainv(2,2) = (A(0,0) * A(1,1) - A(0,1) * A(1,0)) / D;
}

// -------------------------------------------------------------------------------------------------

template <class U, class V>
inline void finger(const U& F, V& B)
{
  B(0,0) = F(0,0) * F(0,0) + F(0,1) * F(0,1) + F(0,2) * F(0,2);
  B(0,1) = F(0,0) * F(1,0) + F(0,1) * F(1,1) + F(0,2) * F(1,2);
  B(0,2) = F(0,0) * F(2,0) + F(0,1) * F(2,1) + F(0,2) * F(2,2);
  B(1,1) = F(1,0) * F(1,0) + F(1,1) * F(1,1) + F(1,2) * F(1,2);
  B(1,2) = F(1,0) * F(2,0) + F(1,1) * F(2,1) + F(1,2) * F(2,2);
  B(2,2) = F(2,0) * F(2,0) + F(2,1) * F(2,1) + F(2,2) * F(2,2);
  B(1,0) = B(0,1);
  B(2,0) = B(0,2);
  B(2,1) = B(1,2);
}

// -------------------------------------------------------------------------------------------------

}} // namespace ...

#endif
