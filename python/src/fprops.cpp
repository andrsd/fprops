#include <pybind11/pybind11.h>

#include "SinglePhaseFluidProperties.h"
#include "Air.h"
#include "IdealGas.h"
#include "Nitrogen.h"
#include "Helium.h"
#include "CarbonDioxide.h"

using namespace fprops;

namespace py = pybind11;

PYBIND11_MODULE(fprops, m)
{
    m.doc() = "pybind11 plugin for fprops";

    py::class_<SinglePhaseFluidProperties> spfp(m, "SinglePhaseFluidProperties");

    py::class_<State>(spfp, "State")
        .def(py::init<>())
        .def_readwrite("u", &State::u)
        .def_readwrite("v", &State::v)
        .def_readwrite("rho", &State::rho)
        .def_readwrite("p", &State::p)
        .def_readwrite("T", &State::T)
        .def_readwrite("mu", &State::mu)
        .def_readwrite("cp", &State::cp)
        .def_readwrite("cv", &State::cv)
        .def_readwrite("rho", &State::rho)
        .def_readwrite("s", &State::s)
        .def_readwrite("k", &State::k)
        .def_readwrite("h", &State::h)
        .def_readwrite("w", &State::w);

    py::class_<Air>(m, "Air")
        .def(py::init())
        .def("rho_T", &Air::rho_T)
        .def("rho_p", &Air::rho_p)
        .def("p_T", &Air::p_T)
        .def("v_u", &Air::v_u);

    py::class_<IdealGas>(m, "IdealGas")
        .def(py::init<double, double>())
        .def("get_gamma", &IdealGas::get_gamma)
        .def("get_specific_gas_constant", &IdealGas::get_specific_gas_constant)
        .def("set_mu", &IdealGas::set_mu)
        .def("set_k", &IdealGas::set_k)
        .def("rho_T", &IdealGas::rho_T)
        .def("rho_p", &IdealGas::rho_p)
        .def("p_T", &IdealGas::p_T)
        .def("v_u", &IdealGas::v_u)
        .def("h_s", &IdealGas::h_s);

    py::class_<Nitrogen>(m, "Nitrogen")
        .def(py::init())
        .def("rho_T", &Nitrogen::rho_T)
        .def("rho_p", &Nitrogen::rho_p)
        .def("p_T", &Nitrogen::p_T)
        .def("v_u", &Nitrogen::v_u);

    py::class_<Helium>(m, "Helium")
        .def(py::init())
        .def("rho_T", &Helium::rho_T)
        .def("rho_p", &Helium::rho_p)
        .def("p_T", &Helium::p_T)
        .def("v_u", &Helium::v_u);

    py::class_<CarbonDioxide>(m, "CarbonDioxide")
        .def(py::init())
        .def("rho_T", &CarbonDioxide::rho_T)
        .def("rho_p", &CarbonDioxide::rho_p)
        .def("p_T", &CarbonDioxide::p_T)
        .def("v_u", &CarbonDioxide::v_u);
}
