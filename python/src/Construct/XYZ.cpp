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
        .def("X", py::overload_cast<>(&LNLib::XYZ::X, py::const_))
        .def("X", py::overload_cast<>(&LNLib::XYZ::X))
        .def("Y", py::overload_cast<>(&LNLib::XYZ::Y, py::const_))
        .def("Y", py::overload_cast<>(&LNLib::XYZ::Y))
        .def("Z", py::overload_cast<>(&LNLib::XYZ::Z, py::const_))
        .def("Z", py::overload_cast<>(&LNLib::XYZ::Z))
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
        .def("__getitem__", [](const LNLib::XYZ& xyz, int index) {
        if (index < 0 || index >= 3) throw py::index_error();
        return xyz[index];
            })
        .def("__setitem__", [](LNLib::XYZ& xyz, int index, double value) {
        if (index < 0 || index >= 3) throw py::index_error();
        xyz[index] = value;
            })
        .def("__add__", [](const LNLib::XYZ& xyz1, const LNLib::XYZ& xyz2) {
        return xyz1 + xyz2;
            })
        .def("__sub__", [](const LNLib::XYZ& xyz1, const LNLib::XYZ& xyz2) {
        return xyz1 - xyz2;
            })
        .def("__mul__", [](const LNLib::XYZ& xyz1, const LNLib::XYZ& xyz2) {
        return xyz1 * xyz2;
            })
        .def("__mul__", [](const LNLib::XYZ& xyz, double d) {
        return xyz * d;
            })
        .def("__rmul__", [](const LNLib::XYZ& xyz, double d) {
        return xyz * d;
            })
        .def("__truediv__", [](const LNLib::XYZ& xyz, double d) {
        return xyz / d;
            })
        .def("__xor__", [](const LNLib::XYZ& xyz1, const LNLib::XYZ& xyz2) {
        return xyz1 ^ xyz2;
            })
        .def("__neg__", [](const LNLib::XYZ& xyz) {
        return -xyz;
            })
        .def("__repr__", [](const LNLib::XYZ& xyz) {
        return "XYZ(" + std::to_string(xyz.GetX()) + ", " + std::to_string(xyz.GetY()) + ", " + std::to_string(xyz.GetZ()) + ")";
            });
}