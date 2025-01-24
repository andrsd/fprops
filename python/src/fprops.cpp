// SPDX-FileCopyrightText: 2022 David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

#include <pybind11/pybind11.h>

#include "fprops/single_phase_fluid_properties.h"
#include "fprops/ideal_gas.h"
#include "fprops/helmholtz.h"
#include "fprops/air.h"
#include "fprops/ammonia.h"
#include "fprops/nitrogen.h"
#include "fprops/helium.h"
#include "fprops/carbon_dioxide.h"
#include "fprops/oxygen.h"
#include "fprops/methane.h"
#include "version.h"

using namespace fprops;

namespace py = pybind11;

PYBIND11_MODULE(fprops, m)
{
    m.doc() = "pybind11 plugin for fprops";
    py::setattr(m, "version", py::str(FPROPS_VERSION));

    py::class_<SinglePhaseFluidProperties>(m, "SinglePhaseFluidProperties")
        .def("rho_T", &SinglePhaseFluidProperties::rho_T)
        .def("rho_p", &SinglePhaseFluidProperties::rho_p)
        .def("p_T", &SinglePhaseFluidProperties::p_T)
        .def("v_u", &SinglePhaseFluidProperties::v_u)
        .def("h_s", &SinglePhaseFluidProperties::h_s);

    py::class_<State>(m, "State")
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
        .def_readwrite("w", &State::w)
        .def("__repr__", &State::to_string);

    py::class_<Helmholtz, SinglePhaseFluidProperties>(m, "Helmholtz");

    py::class_<Air, Helmholtz>(m, "Air").def(py::init());

    py::class_<Ammonia, Helmholtz>(m, "Ammonia").def(py::init());

    py::class_<IdealGas, SinglePhaseFluidProperties>(m, "IdealGas")
        .def(py::init<double, double>())
        .def("gamma", &IdealGas::gamma)
        .def("molar_mass", &IdealGas::molar_mass)
        .def("specific_gas_constant", &IdealGas::R_specific)
        .def("cp", &IdealGas::cp)
        .def("cv", &IdealGas::cv)
        .def("mu", &IdealGas::mu)
        .def("k", &IdealGas::k)
        .def("set_mu", &IdealGas::set_mu)
        .def("set_k", &IdealGas::set_k);

    py::class_<Nitrogen, Helmholtz>(m, "Nitrogen").def(py::init());

    py::class_<Helium, Helmholtz>(m, "Helium").def(py::init());

    py::class_<CarbonDioxide, Helmholtz>(m, "CarbonDioxide").def(py::init());

    py::class_<Oxygen, Helmholtz>(m, "Oxygen").def(py::init());

    py::class_<Methane, Helmholtz>(m, "Methane").def(py::init());
}
