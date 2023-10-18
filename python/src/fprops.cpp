#include <pybind11/pybind11.h>

#include "SinglePhaseFluidProperties.h"
#include "Air.h"
#include "IdealGas.h"
#include "Nitrogen.h"
#include "Helium.h"
#include "CarbonDioxide.h"

using namespace fprops;

namespace py = pybind11;

PYBIND11_MODULE(pyfprops, m)
{
    m.doc() = "pybind11 plugin for fprops";

    py::class_<SinglePhaseFluidProperties> spfp(m, "SinglePhaseFluidProperties");

    py::class_<Props>(spfp, "Props")
        .def(py::init<>())
        .def_readwrite("u", &Props::u)
        .def_readwrite("v", &Props::v)
        .def_readwrite("rho", &Props::rho)
        .def_readwrite("p", &Props::p)
        .def_readwrite("T", &Props::T)
        .def_readwrite("mu", &Props::mu)
        .def_readwrite("cp", &Props::cp)
        .def_readwrite("cv", &Props::cv)
        .def_readwrite("rho", &Props::rho)
        .def_readwrite("s", &Props::s)
        .def_readwrite("k", &Props::k)
        .def_readwrite("h", &Props::h)
        .def_readwrite("w", &Props::w);

    py::class_<Air>(m, "Air")
        .def(py::init())
        .def("rho_T", &Air::rho_T)
        .def("rho_p", &Air::rho_p)
        .def("p_T", &Air::p_T)
        .def("v_u", &Air::v_u);

    py::class_<IdealGas>(m, "IdealGas")
        .def(py::init<double, double>())
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
