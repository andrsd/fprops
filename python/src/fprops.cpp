#include <pybind11/pybind11.h>

#include "SinglePhaseFluidProperties.h"
#include "Air.h"
#include "IdealGas.h"
#include "Nitrogen.h"
#include "Helium.h"

using namespace fprops;

namespace py = pybind11;

PYBIND11_MODULE(pyfprops, m)
{
    m.doc() = "pybind11 plugin for fprops";

    py::class_<SinglePhaseFluidProperties> spfp(m, "SinglePhaseFluidProperties");

    py::class_<SinglePhaseFluidProperties::Props>(spfp, "Props")
        .def(py::init<>())
        .def_readwrite("u", &SinglePhaseFluidProperties::Props::u)
        .def_readwrite("v", &SinglePhaseFluidProperties::Props::v)
        .def_readwrite("rho", &SinglePhaseFluidProperties::Props::rho)
        .def_readwrite("p", &SinglePhaseFluidProperties::Props::p)
        .def_readwrite("T", &SinglePhaseFluidProperties::Props::T)
        .def_readwrite("mu", &SinglePhaseFluidProperties::Props::mu)
        .def_readwrite("cp", &SinglePhaseFluidProperties::Props::cp)
        .def_readwrite("cv", &SinglePhaseFluidProperties::Props::cv)
        .def_readwrite("rho", &SinglePhaseFluidProperties::Props::rho)
        .def_readwrite("s", &SinglePhaseFluidProperties::Props::s)
        .def_readwrite("k", &SinglePhaseFluidProperties::Props::k)
        .def_readwrite("h", &SinglePhaseFluidProperties::Props::h)
        .def_readwrite("w", &SinglePhaseFluidProperties::Props::w);

    py::class_<Air>(m, "Air").def(py::init()).def("p_T", &Air::p_T).def("v_u", &Air::v_u);

    py::class_<IdealGas>(m, "IdealGas")
        .def(py::init<double, double>())
        .def("set_mu", &IdealGas::set_mu)
        .def("set_k", &IdealGas::set_k)
        .def("p_T", &IdealGas::p_T)
        .def("v_u", &IdealGas::v_u);

    py::class_<Nitrogen>(m, "Nitrogen")
        .def(py::init())
        .def("p_T", &Nitrogen::p_T)
        .def("v_u", &Nitrogen::v_u);

    py::class_<Helium>(m, "Helium")
        .def(py::init())
        .def("p_T", &Nitrogen::p_T)
        .def("v_u", &Nitrogen::v_u);
}
