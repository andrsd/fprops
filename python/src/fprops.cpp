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

    py::class_<Air>(m, "Air")
        .def(py::init())
        .def("rho_T", &Air::rho_T)
        .def("rho_p", &Air::rho_p)
        .def("p_T", &Air::p_T)
        .def("v_u", &Air::v_u)
        .def("h_s", &Air::h_s);

    py::class_<Ammonia>(m, "Ammonia")
        .def(py::init())
        .def("rho_T", &Ammonia::rho_T)
        .def("rho_p", &Ammonia::rho_p)
        .def("p_T", &Ammonia::p_T)
        .def("v_u", &Ammonia::v_u)
        .def("h_s", &Ammonia::h_s);

    py::class_<IdealGas>(m, "IdealGas")
        .def(py::init<double, double>())
        .def("rho_T", &IdealGas::rho_T)
        .def("rho_p", &IdealGas::rho_p)
        .def("p_T", &IdealGas::p_T)
        .def("v_u", &IdealGas::v_u)
        .def("h_s", &IdealGas::h_s)
        .def("gamma", &IdealGas::gamma)
        .def("molar_mass", &IdealGas::molar_mass)
        .def("specific_gas_constant", &IdealGas::R_specific)
        .def("cp", &IdealGas::cp)
        .def("cv", &IdealGas::cv)
        .def("mu", &IdealGas::mu)
        .def("k", &IdealGas::k)
        .def("set_mu", &IdealGas::set_mu)
        .def("set_k", &IdealGas::set_k);

    py::class_<Nitrogen>(m, "Nitrogen")
        .def(py::init())
        .def("rho_T", &Nitrogen::rho_T)
        .def("rho_p", &Nitrogen::rho_p)
        .def("p_T", &Nitrogen::p_T)
        .def("v_u", &Nitrogen::v_u)
        .def("h_s", &Nitrogen::h_s);

    py::class_<Helium>(m, "Helium")
        .def(py::init())
        .def("rho_T", &Helium::rho_T)
        .def("rho_p", &Helium::rho_p)
        .def("p_T", &Helium::p_T)
        .def("v_u", &Helium::v_u)
        .def("h_s", &Helium::h_s);

    py::class_<CarbonDioxide>(m, "CarbonDioxide")
        .def(py::init())
        .def("rho_T", &CarbonDioxide::rho_T)
        .def("rho_p", &CarbonDioxide::rho_p)
        .def("p_T", &CarbonDioxide::p_T)
        .def("v_u", &CarbonDioxide::v_u)
        .def("h_s", &CarbonDioxide::h_s);

    py::class_<Oxygen>(m, "Oxygen")
        .def(py::init())
        .def("rho_T", &Oxygen::rho_T)
        .def("rho_p", &Oxygen::rho_p)
        .def("p_T", &Oxygen::p_T)
        .def("v_u", &Oxygen::v_u)
        .def("h_s", &Oxygen::h_s);

    py::class_<Methane>(m, "Methane")
        .def(py::init())
        .def("rho_T", &Methane::rho_T)
        .def("rho_p", &Methane::rho_p)
        .def("p_T", &Methane::p_T)
        .def("v_u", &Methane::v_u)
        .def("h_s", &Methane::h_s);
}
