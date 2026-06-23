/*
 * Author:
 * 2025/11/23 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include "Constants.h"
#include "UV.h"

namespace py = pybind11;

void cstrUV(py::module_&m)
{
    py::class_<LNLib::UV>(m, "UV")
        .def(py::init<>())
        .def(py::init<double, double>())
        .def("SetU", &LNLib::UV::SetU)
        .def("GetU", &LNLib::UV::GetU)
        .def("SetV", &LNLib::UV::SetV)
        .def("GetV", &LNLib::UV::GetV)
        .def_property("U",
            [](const LNLib::UV& self) { return self.U(); },
            [](LNLib::UV& self, double val) { self.U() = val; })
        .def_property("V",
            [](const LNLib::UV& self) { return self.V(); },
            [](LNLib::UV& self, double val) { self.V() = val; })
        .def("IsZero", &LNLib::UV::IsZero, py::arg("epsilon") = LNLib::Constants::DoubleEpsilon)
        .def("IsUnit", &LNLib::UV::IsUnit, py::arg("epsilon") = LNLib::Constants::DoubleEpsilon)
        .def("IsAlmostEqualTo", &LNLib::UV::IsAlmostEqualTo)
        .def("Length", &LNLib::UV::Length)
        .def("SqrLength", &LNLib::UV::SqrLength)
        .def("Normalize", &LNLib::UV::Normalize)
        .def("Add", &LNLib::UV::Add)
        .def("Substract", &LNLib::UV::Substract)
        .def("Negative", &LNLib::UV::Negative)
        .def("DotProduct", &LNLib::UV::DotProduct)
        .def("CrossProduct", &LNLib::UV::CrossProduct)
        .def("Distance", &LNLib::UV::Distance)
        .def("__getitem__", [](const LNLib::UV& self, int index) { return self[index]; })
        .def("__setitem__", [](LNLib::UV& self, int index, double val) { self[index] = val; })
        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(py::self * double())
        .def(double() * py::self)
        .def(py::self / double())
        .def("__xor__", [](const LNLib::UV& a, const LNLib::UV& b) { return a ^ b; })
        .def(-py::self);
}



