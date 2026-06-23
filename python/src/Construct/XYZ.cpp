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
#include "XYZ.h"

namespace py = pybind11;

void cstrXYZ(py::module_&m)
{
    py::class_<LNLib::XYZ>(m, "XYZ")
        .def(py::init<>())
        .def(py::init<double, double, double>())
        .def("SetX", &LNLib::XYZ::SetX)
        .def("GetX", &LNLib::XYZ::GetX)
        .def("SetY", &LNLib::XYZ::SetY)
        .def("GetY", &LNLib::XYZ::GetY)
        .def("SetZ", &LNLib::XYZ::SetZ)
        .def("GetZ", &LNLib::XYZ::GetZ)
        .def_property("X",
            [](const LNLib::XYZ& self) { return self.X(); },
            [](LNLib::XYZ& self, double val) { self.X() = val; })
        .def_property("Y",
            [](const LNLib::XYZ& self) { return self.Y(); },
            [](LNLib::XYZ& self, double val) { self.Y() = val; })
        .def_property("Z",
            [](const LNLib::XYZ& self) { return self.Z(); },
            [](LNLib::XYZ& self, double val) { self.Z() = val; })
        .def("IsZero", &LNLib::XYZ::IsZero, py::arg("epsilon") = LNLib::Constants::DoubleEpsilon)
        .def("IsUnit", &LNLib::XYZ::IsUnit, py::arg("epsilon") = LNLib::Constants::DoubleEpsilon)
        .def("IsAlmostEqualTo", &LNLib::XYZ::IsAlmostEqualTo)
        .def("Length", &LNLib::XYZ::Length)
        .def("SqrLength", &LNLib::XYZ::SqrLength)
        .def("AngleTo", &LNLib::XYZ::AngleTo)
        .def("Normalize", &LNLib::XYZ::Normalize)
        .def("Add", &LNLib::XYZ::Add)
        .def("Substract", &LNLib::XYZ::Substract)
        .def("Negative", &LNLib::XYZ::Negative)
        .def("DotProduct", &LNLib::XYZ::DotProduct)
        .def("CrossProduct", &LNLib::XYZ::CrossProduct)
        .def("Distance", &LNLib::XYZ::Distance)
        .def_static("CreateRandomOrthogonal", &LNLib::XYZ::CreateRandomOrthogonal)
        .def("__getitem__", [](const LNLib::XYZ& self, int index) { return self[index]; })
        .def("__setitem__", [](LNLib::XYZ& self, int index, double val) { self[index] = val; })
        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(py::self * double())
        .def(double() * py::self)
        .def(py::self / double())
        .def("__xor__", [](const LNLib::XYZ& a, const LNLib::XYZ& b) { return a ^ b; })
        .def(-py::self);
}