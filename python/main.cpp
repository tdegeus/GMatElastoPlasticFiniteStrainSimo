/*

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastoPlasticFiniteStrainSimo

*/

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#define FORCE_IMPORT_ARRAY
#include <xtensor-python/pytensor.hpp>

#include <GMatElastoPlasticFiniteStrainSimo/version.h>
#include <GMatElastoPlasticFiniteStrainSimo/Cartesian3d.h>

namespace py = pybind11;

template <class S, class T>
auto construct_Array(T& self)
{
    namespace SM = GMatElastoPlasticFiniteStrainSimo::Cartesian3d;

    self.def(py::init<std::array<size_t, S::rank>>(), "Array of material points.", py::arg("shape"))

        .def("shape", &S::shape, "Shape of array.")
        .def("I2", &S::I2, "Array with 2nd-order unit tensors.")
        .def("II", &S::II, "Array with 4th-order tensors = dyadic(I2, I2).")
        .def("I4", &S::I4, "Array with 4th-order unit tensors.")
        .def("I4rt", &S::I4rt, "Array with 4th-order right-transposed unit tensors.")
        .def("I4s", &S::I4s, "Array with 4th-order symmetric projection tensors.")
        .def("I4d", &S::I4d, "Array with 4th-order deviatoric projection tensors.")
        .def("K", &S::K, "Array with K.")
        .def("G", &S::G, "Array with G.")
        .def("increment", &S::increment, "Increment history variables.")
        .def("type", &S::type, "Array with material types.")

        .def(
            "setElastic",
            py::overload_cast<
                const xt::xtensor<size_t, S::rank>&,
                const xt::xtensor<size_t, S::rank>&,
                const xt::xtensor<double, 1>&,
                const xt::xtensor<double, 1>&>(&S::setElastic),
            "Set specific entries 'Elastic'.",
            py::arg("I"),
            py::arg("idx"),
            py::arg("K"),
            py::arg("G"))

        .def(
            "setLinearHardening",
            py::overload_cast<
                const xt::xtensor<size_t, S::rank>&,
                const xt::xtensor<size_t, S::rank>&,
                const xt::xtensor<double, 1>&,
                const xt::xtensor<double, 1>&,
                const xt::xtensor<double, 1>&,
                const xt::xtensor<double, 1>&>(&S::setLinearHardening),
            "Set specific entries 'LinearHardening'.",
            py::arg("I"),
            py::arg("idx"),
            py::arg("K"),
            py::arg("G"),
            py::arg("tauy0"),
            py::arg("H"))

        .def(
            "setElastic",
            py::overload_cast<const xt::xtensor<size_t, S::rank>&, double, double>(
                &S::setElastic),
            "Set specific entries 'Elastic'.",
            py::arg("I"),
            py::arg("K"),
            py::arg("G"))

        .def(
            "setLinearHardening",
            py::overload_cast<
                const xt::xtensor<size_t, S::rank>&,
                double,
                double,
                double,
                double>(&S::setLinearHardening),
            "Set specific entries 'LinearHardening'.",
            py::arg("I"),
            py::arg("K"),
            py::arg("G"),
            py::arg("tauy0"),
            py::arg("H"))

        .def(
            "setDefGrad",
            &S::setDefGrad,
            "Set deformation gradient tensors (computes stress and optionally tangent).",
            py::arg("F"),
            py::arg("compute_tangent") = true)

        .def("DefGrad", &S::DefGrad, "Get deformation gradient tensors.")
        .def("Stress", &S::Stress, "Get stress tensors.")
        .def("Tangent", &S::Tangent, "Get stiffness tensors.")
        .def("Epsp", &S::Epsp, "Array with plastic strains.")
        .def("getElastic", &S::getElastic, "Returns underlying Elastic model.")
        .def("getLinearHardening", &S::getLinearHardening, "Returns underlying LinearHardening model.")

        .def("__repr__", [](const S&) { return "<GMatElastoPlasticFiniteStrainSimo.Cartesian3d.Array>"; });
}

template <class S, class T>
void add_strain_overloads(T& module)
{
    module.def(
        "Strain",
        static_cast<S (*)(const S&)>(&GMatElastoPlasticFiniteStrainSimo::Cartesian3d::Strain<S>),
        "Strain from of a(n) (array of) deformation gradient tensor(s).",
        py::arg("A"));
}

template <class S, class T>
void add_deviatoric_overloads(T& module)
{
    module.def(
        "Deviatoric",
        static_cast<S (*)(const S&)>(&GMatElastoPlasticFiniteStrainSimo::Cartesian3d::Deviatoric<S>),
        "Deviatoric part of a(n) (array of) tensor(s).",
        py::arg("A"));
}

template <class R, class S, class T>
void add_hydrostatic_overloads(T& module)
{
    module.def(
        "Hydrostatic",
        static_cast<R (*)(const S&)>(&GMatElastoPlasticFiniteStrainSimo::Cartesian3d::Hydrostatic<S>),
        "Hydrostatic part of a(n) (array of) tensor(s).",
        py::arg("A"));
}

template <class R, class S, class T>
void add_epseq_overloads(T& module)
{
    module.def(
        "Epseq",
        static_cast<R (*)(const S&)>(
            &GMatElastoPlasticFiniteStrainSimo::Cartesian3d::Epseq<S>),
        "Equivalent strain of a(n) (array of) tensor(s).",
        py::arg("A"));
}

template <class R, class S, class T>
void add_sigeq_overloads(T& module)
{
    module.def(
        "Sigeq",
        static_cast<R (*)(const S&)>(
            &GMatElastoPlasticFiniteStrainSimo::Cartesian3d::Sigeq<S>),
        "Equivalent stress of a(n) (array of) tensor(s).",
        py::arg("A"));
}

PYBIND11_MODULE(_GMatElastoPlasticFiniteStrainSimo, m)
{
    xt::import_numpy();

    m.doc() = "Elasto-plastic material model";

    m.def("version",
          &GMatElastoPlasticFiniteStrainSimo::version,
          "Return version string.");

    m.def("version_dependencies",
          &GMatElastoPlasticFiniteStrainSimo::version_dependencies,
          "Return list of strings.");

    // -----------------------------
    // GMatElastoPlasticFiniteStrainSimo.Cartesian3d
    // -----------------------------

    py::module sm = m.def_submodule("Cartesian3d", "3d Cartesian coordinates");

    namespace SM = GMatElastoPlasticFiniteStrainSimo::Cartesian3d;

    // Unit tensors

    sm.def("I2", &SM::I2, "Second order unit tensor.");
    sm.def("II", &SM::II, "Fourth order tensor with the result of the dyadic product II.");
    sm.def("I4", &SM::I4, "Fourth order unit tensor.");
    sm.def("I4rt", &SM::I4rt, "Fourth right-transposed order unit tensor.");
    sm.def("I4s", &SM::I4s, "Fourth order symmetric projection tensor.");
    sm.def("I4d", &SM::I4d, "Fourth order deviatoric projection tensor.");

    // Tensor algebra

    add_strain_overloads<xt::xtensor<double, 4>>(sm);
    add_strain_overloads<xt::xtensor<double, 3>>(sm);
    add_strain_overloads<xt::xtensor<double, 2>>(sm);
    add_deviatoric_overloads<xt::xtensor<double, 4>>(sm);
    add_deviatoric_overloads<xt::xtensor<double, 3>>(sm);
    add_deviatoric_overloads<xt::xtensor<double, 2>>(sm);
    add_hydrostatic_overloads<xt::xtensor<double, 2>, xt::xtensor<double, 4>>(sm);
    add_hydrostatic_overloads<xt::xtensor<double, 1>, xt::xtensor<double, 3>>(sm);
    add_hydrostatic_overloads<xt::xtensor<double, 0>, xt::xtensor<double, 2>>(sm);
    add_epseq_overloads<xt::xtensor<double, 2>, xt::xtensor<double, 4>>(sm);
    add_epseq_overloads<xt::xtensor<double, 1>, xt::xtensor<double, 3>>(sm);
    add_epseq_overloads<xt::xtensor<double, 0>, xt::xtensor<double, 2>>(sm);
    add_sigeq_overloads<xt::xtensor<double, 2>, xt::xtensor<double, 4>>(sm);
    add_sigeq_overloads<xt::xtensor<double, 1>, xt::xtensor<double, 3>>(sm);
    add_sigeq_overloads<xt::xtensor<double, 0>, xt::xtensor<double, 2>>(sm);

    // Material point: Elastic

    py::class_<SM::Elastic>(sm, "Elastic")

        .def(py::init<double, double>(), "Linear elastic material point.", py::arg("K"), py::arg("G"))

        .def("K", &SM::Elastic::K, "Returns the bulk modulus.")
        .def("G", &SM::Elastic::G, "Returns the shear modulus.")

        .def(
            "setDefGrad",
            &SM::Elastic::setDefGrad<xt::xtensor<double, 2>>,
            "Set current deformation gradient tensor (computes stress and optionally tangent).",
            py::arg("F"),
            py::arg("compute_tangent") = true)

        .def("DefGrad", &SM::Elastic::DefGrad, "Returns deformation gradient tensor.")
        .def("Stress", &SM::Elastic::Stress, "Returns stress tensor.")
        .def("Tangent", &SM::Elastic::Tangent, "Returns tangent stiffness.")

        .def("__repr__", [](const SM::Elastic&) { return "<GMatElastic.Cartesian3d.Elastic>"; });

    // Material point: LinearHardening

    py::class_<SM::LinearHardening>(sm, "LinearHardening")

        .def(py::init<double, double, double, double>(),
            "Elasto-plastic with linear hardening",
            py::arg("K"),
            py::arg("G"),
            py::arg("tauy0"),
            py::arg("H"))

        .def("K", &SM::LinearHardening::K, "Returns the bulk modulus.")
        .def("G", &SM::LinearHardening::G, "Returns the shear modulus.")
        .def("tauy0", &SM::LinearHardening::tauy0, "Returns the initial yield stress.")
        .def("H", &SM::LinearHardening::H, "Returns the hardening modulus.")
        .def("epsp", &SM::LinearHardening::epsp, "Returns the plastic strain.")
        .def("increment", &SM::LinearHardening::increment, "Updates history.")

        .def(
            "setDefGrad",
            &SM::LinearHardening::setDefGrad<xt::xtensor<double, 2>>,
            "Set current deformation gradient tensor (computes stress and optionally tangent).",
            py::arg("F"),
            py::arg("compute_tangent") = true)

        .def("DefGrad", &SM::LinearHardening::DefGrad, "Returns deformation gradient tensor.")
        .def("Stress", &SM::LinearHardening::Stress, "Returns stress tensor.")
        .def("Tangent", &SM::LinearHardening::Tangent, "Returns tangent stiffness.")

        .def("__repr__", [](const SM::LinearHardening&) {
            return "<GMatElastoPlasticFiniteStrainSimo.Cartesian3d.LinearHardening>";
        });

    py::module smm = sm.def_submodule("Type", "Type enumerator");

    py::enum_<SM::Type::Value>(smm, "Type")
        .value("Unset", SM::Type::Unset)
        .value("Elastic", SM::Type::Elastic)
        .value("LinearHardening", SM::Type::LinearHardening)
        .export_values();

    // Array

    py::class_<SM::Array<1>> array1d(sm, "Array1d");
    py::class_<SM::Array<2>> array2d(sm, "Array2d");
    py::class_<SM::Array<3>> array3d(sm, "Array3d");

    construct_Array<SM::Array<1>>(array1d);
    construct_Array<SM::Array<2>>(array2d);
    construct_Array<SM::Array<3>>(array3d);
}
