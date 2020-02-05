
#include <catch2/catch.hpp>

#define EQ(a,b) REQUIRE_THAT((a), Catch::WithinAbs((b), 1.e-12));

#include <GMatElastoPlasticFiniteStrainSimo/Cartesian3d.h>

namespace GM = GMatElastoPlasticFiniteStrainSimo::Cartesian3d;


TEST_CASE("GMatElastoPlasticFiniteStrainSimo::Cartesian3d", "Cartesian3d.h")
{

double K = 12.3;
double G = 45.6;

double gamma = 0.02;

GM::Tensor2 F;
F.fill(0.0);
F(0,0) = 1.0 + gamma;
F(1,1) = 1.0 / (1.0 + gamma);
F(2,2) = 1.0;


SECTION("Elastic")
{
    GM::Elastic mat(K, G);

    auto Sig = mat.Stress(F);

    EQ(Sig(0,0), G * +2.0 * std::log(1.0 + gamma));
    EQ(Sig(1,1), G * -2.0 * std::log(1.0 + gamma));
    EQ(Sig(2,2), 0.0);
    EQ(Sig(0,1), 0.0);
    EQ(Sig(0,2), 0.0);
    EQ(Sig(1,0), 0.0);
    EQ(Sig(1,2), 0.0);
    EQ(Sig(2,0), 0.0);
    EQ(Sig(2,1), 0.0);
}


SECTION("Matrix")
{
    size_t nelem = 3;
    size_t nip = 2;

    GM::Matrix mat(nelem, nip);

    // all rows elastic
    {
        xt::xtensor<size_t,2> I = xt::ones<size_t>({nelem, nip});
        mat.setElastic(I, K, G);
    }

    xt::xtensor<double,4> f = xt::empty<double>({nelem, nip, 3ul, 3ul});

    for (size_t e = 0; e < nelem; ++e) {
        for (size_t q = 0; q < nip; ++q) {
            xt::view(f, e, q) = F;
        }
    }

    auto sig = mat.Stress(f);

    for (size_t e = 0; e < nelem; ++e) {
        for (size_t q = 0; q < nip; ++q) {
            EQ(sig(e,q,0,0), G * +2.0 * std::log(1.0 + gamma));
            EQ(sig(e,q,1,1), G * -2.0 * std::log(1.0 + gamma));
            EQ(sig(e,q,0,1), 0.0);
            EQ(sig(e,q,1,0), 0.0);
        }
    }

    REQUIRE(xt::allclose(xt::view(sig, xt::all(), xt::all(), 2, xt::all()), 0.0));
    REQUIRE(xt::allclose(xt::view(sig, xt::all(), xt::all(), xt::all(), 2), 0.0));
}


}
